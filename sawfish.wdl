version 1.0

workflow sawfish {
  meta {
    title: "Sawfish workflow"
    summary: "Structural variant discovery and calling using Sawfish"
    description: "Workflow for discovering and genotyping structural variants using Sawfish. The workflow performs discovery of candidate structural variants followed by calling outputting an VCF with SV calls."
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
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
    genome_version: {
      help: "Reference genome version",
      label: "Genome version"
    }
    vcf: {
      help: "Output VCF containing structural variants",
      label: "Structural variants VCF"
    }
    vcf_tbi: {
      help: "Index for the structural variants VCF",
      label: "Structural variants VCF index"
    }
    minimap2_bam: {
      help: "BAM file with reads remapped using minimap2",
      label: "Minimap2 BAM"
    }
    minimap2_bam_bai: {
      help: "Index for the minimap2 remapped BAM file",
      label: "Minimap2 BAM index"
    }
  }

  input {
    String sample_name
    File mapped_bam
    File mapped_bam_bai
    File ref_fasta
    File ref_index
    String genome_version
    String docker_smrttools
  }

  call run_sawfish {
    input:
      sample_name = sample_name,
      mapped_bam = mapped_bam,
      mapped_bam_bai = mapped_bam_bai,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      genome_version = genome_version,
      docker_smrttools = docker_smrttools
  }

  output {
    File vcf = run_sawfish.vcf
    File vcf_tbi = run_sawfish.vcf_tbi
    File minimap2_bam = run_sawfish.minimap2_bam
    File minimap2_bam_bai = run_sawfish.minimap2_bam_bai
  }
}

task run_sawfish {
  meta {
    title: "Sawfish SV Caller"
    summary: "Discovers and calls structural variants using Sawfish this is restricted to HBA"
    description: "Performs structural variant discovery and joint calling using Sawfish. The task runs in two steps: discovery of candidate SVs followed by joint calling to produce a final VCF with genotypes. This is restricted to HBA"
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
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
    threads: {
      help: "Number of CPU threads to use (default: 16)",
      label: "CPU threads"
    }
    mem_gb: {
      help: "Memory allocation in gigabytes (default: 32)",
      label: "Memory (GB)"
    }
    vcf: {
      help: "Output VCF containing structural variants",
      label: "Structural variants VCF"
    }
    vcf_tbi: {
      help: "Index for the structural variants VCF",
      label: "Structural variants VCF index"
    }
    minimap2_bam: {
      help: "BAM file with reads remapped using minimap2",
      label: "Minimap2 BAM"
    }
    minimap2_bam_bai: {
      help: "Index for the minimap2 remapped BAM file",
      label: "Minimap2 BAM index"
    }
  }

  input {
    String sample_name
    File mapped_bam
    File mapped_bam_bai
    File ref_fasta
    File ref_index
    String genome_version

    Int threads = 16
    Int mem_gb = 32

    String docker_smrttools
  }

  String discover_dirname = "~{sample_name}_sawfish_discover"
  String call_dirname = "~{sample_name}_sawfish"
  String minimap2_bam_filename = "~{sample_name}.minimap2.bam"

  Int disk_size = ceil((size(mapped_bam, 'GB') + size(ref_fasta, 'GB')) * 3 + 20)

  command <<<
    set -e

    if [[ "~{genome_version}" == "38" ]]; then
      HBA_REGION="chr16:126123-209491"
    elif [[ "~{genome_version}" == "37" ]] || [[ "~{genome_version}" == "19" ]]; then
      HBA_REGION="16:176122-259490"
    else
      echo "Warning: Unrecognized genome version ~{genome_version}. Defaulting to 38"
      HBA_REGION="chr16:126123-209491"
    fi

    samtools view -h ~{mapped_bam} -e "[rq]>=0.99" ${HBA_REGION} -T ~{ref_fasta} | \
    samtools fastq -T MM,ML - | \
    minimap2 -a -y -x map-hifi -Y -E 1,0 ~{ref_fasta} - | \
    samtools view -b -h | \
    samtools sort > ~{minimap2_bam_filename}

    samtools index ~{minimap2_bam_filename}

    sawfish discover \
      --disable-cnv \
      --threads ~{threads} \
      --ref ~{ref_fasta} \
      --bam ~{minimap2_bam_filename} \
      --output-dir ~{discover_dirname}

    sawfish joint-call \
      --threads ~{threads} \
      --sample ~{discover_dirname} \
      --output-dir ~{call_dirname}

    mv ~{call_dirname}/genotyped.sv.vcf.gz ~{call_dirname}/~{sample_name}.sv.vcf.gz
    mv ~{call_dirname}/genotyped.sv.vcf.gz.tbi ~{call_dirname}/~{sample_name}.sv.vcf.gz.tbi
  >>>

  output {
    File vcf = "~{call_dirname}/~{sample_name}.sv.vcf.gz"
    File vcf_tbi = "~{call_dirname}/~{sample_name}.sv.vcf.gz.tbi"
    File minimap2_bam = "~{minimap2_bam_filename}"
    File minimap2_bam_bai = "~{minimap2_bam_filename}.bai"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}
