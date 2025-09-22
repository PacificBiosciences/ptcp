version 1.0

workflow paraphase {
  meta {
    title: "Paraphase workflow"
    summary: "Haplotype phasing, genotyping, and F8 inversion calling workflow"
    description: "A workflow for phasing haplotypes using Paraphase followed by genotyping and F8 inversion calling. The workflow performs targeted phasing of specific genomic regions and generates phased BAMs, VCFs, and JSON reports."
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
    }
    sex: {
      help: "Biological sex of the sample",
      label: "Sample sex",
      choices: ["M", "F"]
    }
    mapped_bam: {
      help: "Input BAM file containing aligned reads",
      label: "Mapped BAM"
    }
    mapped_bam_bai: {
      help: "Index file for the input BAM",
      label: "Mapped BAM index"
    }
    ref_fasta: {
      help: "Reference genome in FASTA format",
      label: "Reference FASTA"
    }
    ref_index: {
      help: "Index file for the reference FASTA",
      label: "Reference FASTA index"
    }
    config_file: {
      help: "Paraphase configuration file",
      label: "Paraphase config"
    }
    genome_version: {
      help: "Reference genome version",
      label: "Genome version"
    }
    annotation_vcf: {
      help: "Optional VCF file for havanno annotation of paraphase VCF outputs",
      label: "Paraphase annotation VCF"
    }
    # Outputs
    paraphase_bam: {
      help: "Output phased BAM from Paraphase",
      label: "Phased BAM"
    }
    paraphase_bam_bai: {
      help: "Index for the phased BAM",
      label: "Phased BAM index"
    }
    paraphase_json: {
      help: "JSON report from Paraphase",
      label: "Paraphase report"
    }
    paraphase_vcfs: {
      help: "VCF files containing phased variants",
      label: "Phased VCFs"
    }
    f8_vcf: {
      help: "VCF containing F8 inversion calls",
      label: "F8 inversion VCF"
    }
    f8_json: {
      help: "JSON report for F8 inversion calls",
      label: "F8 inversion report"
    }
    havanno_json: {
      help: "Array of JSON files containing havanno annotations (optional)",
      label: "Paraphase havanno annotations"
    }
  }

  input {
    String sample_name
    String sex
    File mapped_bam
    File mapped_bam_bai
    File ref_fasta
    File ref_index
    File config_file
    String genome_version
    File? annotation_vcf
    
    String docker_smrttools = "quay.io/pacbio/smrttools@sha256:01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b"
  }

  call run_paraphase {
    input:
      sample_name               = sample_name,
      mapped_bam                = mapped_bam,
      mapped_bam_bai            = mapped_bam_bai,
      ref_fasta                 = ref_fasta,
      ref_index                 = ref_index,
      config_file               = config_file,
      genome_version            = genome_version,
      annotation_vcf            = annotation_vcf,
      docker_smrttools          = docker_smrttools
  }

  call call_f8 {
    input:
      sex                       = sex,
      paraphase_bam             = run_paraphase.paraphase_bam,
      paraphase_bam_bai         = run_paraphase.paraphase_bam_bai,
      out_prefix                = "~{sample_name}",
      genome_version            = genome_version,
      docker_smrttools          = docker_smrttools
  }

  output {
    File paraphase_bam          = run_paraphase.paraphase_bam
    File paraphase_bam_bai      = run_paraphase.paraphase_bam_bai
    File paraphase_json         = run_paraphase.paraphase_json
    Array[File] paraphase_vcfs  = run_paraphase.paraphase_vcfs
    File? havanno_json          = run_paraphase.havanno_json
    File f8_vcf                 = call_f8.vcf
    File f8_json                = call_f8.json
  }
}

task run_paraphase {
  meta {
    title: "Paraphase"
    summary: "Performs haplotype phasing and genotyping using Paraphase"
    description: "Uses Paraphase to phase haplotypes in targeted genomic regions. Generates phased BAM files, VCFs, and JSON reports."
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
    }
    mapped_bam: {
      help: "Input BAM file containing aligned reads",
      label: "Input BAM"
    }
    mapped_bam_bai: {
      help: "Index file for the input BAM",
      label: "Input BAM index"
    }
    ref_fasta: {
      help: "Reference genome in FASTA format",
      label: "Reference FASTA"
    }
    ref_index: {
      help: "Index file for the reference FASTA",
      label: "Reference FASTA index"
    }
    config_file: {
      help: "Paraphase configuration file",
      label: "Paraphase config"
    }
    genome_version: {
      help: "Reference genome version",
      label: "Genome version"
    }
    annotation_vcf: {
      help: "Optional VCF file for havanno annotation",
      label: "Annotation VCF"
    }
    threads: {
      help: "Number of CPU threads to use (default: 8)",
      label: "CPU threads"
    }
    mem_gb: {
      help: "Memory allocation in gigabytes (default: 16)",
      label: "Memory (GB)"
    }
  }

  input {
    String sample_name
    File mapped_bam
    File mapped_bam_bai
    File ref_fasta
    File ref_index
    File config_file
    String genome_version
    File? annotation_vcf

    Int threads = 4
    Int mem_gb = 8

    String docker_smrttools
  }

  String out_dir = "~{sample_name}_paraphase"

  Int disk_size = ceil((size(mapped_bam, 'GB') + size(ref_fasta, 'GB')) * 2 + 20)

  command <<<
    set -e

    paraphase \
      --bam ~{mapped_bam} \
      --reference ~{ref_fasta} \
      --out ~{out_dir} \
      --genome ~{genome_version} \
      --threads ~{threads} \
      --config ~{config_file} \
      --write-nocalls-in-vcf \
      --targeted

    ~{if defined(annotation_vcf) then
      "havanno --variant-vcf " + annotation_vcf + " --paraphase-dir " + out_dir + " > " + sample_name + ".havanno.json"
      else ""}
  >>>

  output {
    File paraphase_bam         = "~{out_dir}/~{sample_name}.paraphase.bam"
    File paraphase_bam_bai     = "~{out_dir}/~{sample_name}.paraphase.bam.bai"
    File paraphase_json        = "~{out_dir}/~{sample_name}.paraphase.json"
    Array[File] paraphase_vcfs = glob("~{out_dir}/~{sample_name}_paraphase_vcfs/*")
    File? havanno_json         = "~{sample_name}.havanno.json"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}

task call_f8 {
  meta {
    title: "F8 inversion caller"
    summary: "Calls F8 inversions from phased BAM"
    description: "Analyzes phased BAM files to detect and genotype F8 gene inversions, generates both VCF and JSON results."
  }

  parameter_meta {
    sex: {
      help: "Biological sex of the sample",
      label: "Sample sex",
      choices: ["M", "F"]
    }
    paraphase_bam: {
      help: "Input phased BAM file",
      label: "Phased BAM"
    }
    paraphase_bam_bai: {
      help: "Index for the phased BAM",
      label: "Phased BAM index"
    }
    out_prefix: {
      help: "Prefix for output files",
      label: "Output prefix"
    }
    genome_version: {
      help: "Reference genome version",
      label: "Genome version"
    }
    threads: {
      help: "Number of CPU threads to use (default: 1)",
      label: "CPU threads"
    }
    mem_gb: {
      help: "Memory allocation in gigabytes (default: 4)",
      label: "Memory (GB)"
    }
  }

  input {
    String sex
    File paraphase_bam
    File paraphase_bam_bai
    String genome_version

    String out_prefix

    String docker_smrttools

    Int threads = 1
    Int mem_gb  = 4
  }

  String sample_sex = if select_first([sex, "F"]) == "M" then "M" else "F"

  Int disk_size = ceil(size(paraphase_bam, 'GB') * 2 + 10)

  command <<<
    set -e

    f8_inversion.py \
      --bam ~{paraphase_bam} \
      --prefix ~{out_prefix} \
      --genome ~{genome_version} \
      --json \
      --sex ~{sample_sex} \
      --out out/

    mv out/~{out_prefix}.f8inversion.vcf ~{out_prefix}.f8inversion.vcf
    mv out/~{out_prefix}.f8inversion.json ~{out_prefix}.f8inversion.json
  >>>

  output {
    File vcf  = "~{out_prefix}.f8inversion.vcf"
    File json = "~{out_prefix}.f8inversion.json"
  }


  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}
