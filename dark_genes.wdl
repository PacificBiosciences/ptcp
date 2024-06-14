version 1.0

import "./common.wdl" as Common

workflow dark_genes {
  meta {
    description: "Cluster HiFi reads with fixed start/end positions, (e.g. from amplicon sequencing data), align to a reference, and call variants."
  }

  input {
    String sample_name
    String? sex
    File unmapped_bam

    Array[File] guides
    Array[File] guide_indices

    File ref_fasta
    File ref_index
    File ref_mmi

    String log_level        = "INFO"
    String docker_smrttools = "quay.io/pacbio/smrttools@sha256:01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b"
  }

  call bam2hififastq {
    # extract hifi reads and convert to fastq
    input:
      sample_name      = sample_name,
      bam              = unmapped_bam,
      docker_smrttools = docker_smrttools
  }

  call pbaa_cluster as cluster_by_groups {
    # cluster reads by groups of guides, e.g. all F8_inv22 reads grouped together
    # we use very relaxed settings to include as many reads as possible
    input:
      fastq                  = bam2hififastq.fastq,
      fastq_index            = bam2hififastq.fastq_index,
      guides                 = guides,
      guide_indices          = guide_indices,
      prefix                 = sample_name,
      pile_size              = 100,
      max_reads_per_guide    = 5000,
      max_amplicon_size      = 15000,
      min_cluster_frequency  = 0.001,
      min_cluster_read_count = 2,
      skip_chimera_detection = true,
      skip_consensus         = true,
      log_level              = "FATAL",
      docker_smrttools       = docker_smrttools
  }

  # for each group of guides
  scatter (index in range(length(guides))) {
    String guide_name = basename(guides[index], ".fasta")
    call extract_reads_from_fastq as extract_guide_groups {
      # extract all reads that clustered for this guide
      input:
        info             = cluster_by_groups.info,
        guide_name       = guide_name,
        fastq            = bam2hififastq.fastq,
        fastq_index      = bam2hififastq.fastq_index,
        docker_smrttools = docker_smrttools
    }

    call pbaa_cluster {
      # use pbaa to cluster reads for this guide and generate consensuses
      input:
        fastq                  = extract_guide_groups.clustered_fastq,
        fastq_index            = extract_guide_groups.clustered_fastq_index,
        guides                 = [guides[index]],
        guide_indices          = [guide_indices[index]],
        prefix                 = "~{sample_name}.~{guide_name}",
        pile_size              = 100,
        max_reads_per_guide    = 1000,
        max_consensus_reads    = 10,
        max_amplicon_size      = 15000,
        min_cluster_frequency  = 0.005,
        min_cluster_read_count = 2,
        skip_chimera_detection = true,
        log_level              = log_level,
        docker_smrttools       = docker_smrttools
    }

    call extract_reads_from_bam as extract_clustered_reads {
      # extract all reads that clustered for this guide
      input:
        info             = pbaa_cluster.info,
        bam              = unmapped_bam,
        docker_smrttools = docker_smrttools
    }

    call Common.pbmm2_align_hifi as align_clusters {
      # align consensus sequences to the masked reference
      input:
        sample_name      = sample_name,
        unmapped         = pbaa_cluster.passed_cluster_sequences,
        ref_fasta        = ref_mmi,
        out_prefix       = "~{sample_name}.~{guide_name}.consensus",
        log_level        = log_level,
        docker_smrttools = docker_smrttools
    }

    call Common.pbmm2_align_hifi as align_reads {
      # align hifi reads to the masked reference
      input:
        sample_name      = sample_name,
        unmapped         = extract_clustered_reads.clustered_bam,
        ref_fasta        = ref_mmi,
        out_prefix       = "~{sample_name}.~{guide_name}.hifi_reads",
        log_level        = log_level,
        docker_smrttools = docker_smrttools
    }

    call paint_bam {
      # add HP and YC tags to the aligned hifi reads
      input:
        info             = pbaa_cluster.info,
        bam              = align_reads.mapped_bam,
        bam_index        = align_reads.mapped_bam_bai,
        out_prefix       = "~{sample_name}.~{guide_name}.painted_hifi_reads",
        docker_smrttools = docker_smrttools
    }
  }

  call merge_bams as merge_consensus_bams {
    # merge all consensus bams for this sample
    input:
      bams             = align_clusters.mapped_bam,
      bam_indices      = align_clusters.mapped_bam_bai,
      out_prefix       = "~{sample_name}.consensus",
      docker_smrttools = docker_smrttools
  }

  call minipileup {
    # call variants from each consensus individually, and merge the vcfs
    input:
      sample_name      = sample_name,
      bam              = merge_consensus_bams.bam,
      bam_index        = merge_consensus_bams.bam_bai,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      out_prefix       = "~{sample_name}",
      docker_smrttools = docker_smrttools
  }

  call merge_bams as merge_painted_bams {
    # merge all painted bams for this sample
    input:
      bams             = paint_bam.painted_bam,
      bam_indices      = paint_bam.painted_bam_bai,
      out_prefix       = "~{sample_name}.painted_hifi_reads",
      docker_smrttools = docker_smrttools
  }

  call call_f8 as call_f8_pbaa {
    # call F8 inversion genotypes from the painted bams
    input:
      sex              = sex,
      bam              = merge_painted_bams.bam,
      bam_index        = merge_painted_bams.bam_bai,
      out_prefix       = "~{sample_name}",
      docker_smrttools = docker_smrttools
  }

  output {
    File cluster_info      = cluster_by_groups.info
    File painted_bam       = merge_painted_bams.bam
    File painted_bam_bai   = merge_painted_bams.bam_bai
    File consensus_bam     = merge_consensus_bams.bam
    File consensus_bam_bai = merge_consensus_bams.bam_bai
    File minipileup_vcf    = minipileup.vcf
    File f8_pbaa_vcf       = call_f8_pbaa.vcf
    File f8_pbaa_json      = call_f8_pbaa.json
  }
}

task make_pbaa_guides {
  input {
    File ref_fasta
    File ref_index

    File guide_bed

    String docker_smrttools
  }

  Int threads = 1
  Int mem_gb  = 4

  command <<<
    make_guides.py ~{guide_bed} ~{ref_fasta}
    for i in ./*.fasta; do samtools faidx "$i"; done
  >>>

  output {
    Array[File] guides        = glob("*.fasta")
    Array[File] guide_indices = glob("*.fasta.fai")
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task mask_reference {
  input {
    File ref_fasta
    File ref_index

    File mask_bed

    String docker_smrttools
  }

  Int threads = 1
  Int mem_gb  = 4

  String out_prefix = basename(ref_fasta, ".fasta")

  command <<<
    bedtools maskfasta \
      -fi ~{ref_fasta} \
      -bed ~{mask_bed} \
      -fo ~{out_prefix}.masked.fasta

    samtools faidx ~{out_prefix}.masked.fasta
  >>>

  output {
    File masked_ref   = "~{out_prefix}.masked.fasta"
    File masked_index = "~{out_prefix}.masked.fasta.fai"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task bam2hififastq {
  input {
    String sample_name
    File bam

    String docker_smrttools
  }

  Int threads = 8
  Int mem_gb  = 4

  command <<<
    extracthifi --num-threads ~{threads} \
      ~{bam} ~{sample_name}.hifi_only.bam

    samtools fastq --threads ~{threads - 1} \
      ~{sample_name}.hifi_only.bam \
      > ~{sample_name}.fastq
    samtools fqidx ~{sample_name}.fastq
  >>>

  output {
    File fastq       = "~{sample_name}.fastq"
    File fastq_index = "~{sample_name}.fastq.fai"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task pbaa_cluster {
  input {
    File fastq
    File fastq_index

    Array[File] guides
    Array[File] guide_indices

    String prefix

    Int   filter                  = 3
    Int   trim_ends               = 5
    Int   pile_size               = 30
    Float min_var_frequency       = 0.3
    Int   max_alignments_per_read = 1000
  
    Int max_reads_per_guide = 500
    Int iterations          = 9
    Int seed                = 1984

    Int max_consensus_reads = 100

    Int   max_amplicon_size = 15000
    Float min_read_qv       = 20
    String? off_target_groups
    Float min_cluster_frequency  = 0.1
    Int   min_cluster_read_count = 5
    Float max_uchime_score       = 1

    Boolean skip_chimera_detection = false
    Boolean skip_consensus         = false

    String log_level = "INFO"
    String docker_smrttools
  }

  Int threads = 16
  Int mem_gb  = 32

  command <<<
    if [ ~{length(guides)} -gt 1 ];
    then
      cat ~{sep=" " guides} > guide.fasta
      samtools faidx guide.fasta
    else
      cp ~{guides[0]} guide.fasta
      cp ~{guide_indices[0]} guide.fasta.fai
    fi

    pbaa cluster \
      --num-threads ~{threads} \
      --log-level ~{log_level} \
      --filter ~{filter} \
      --trim-ends ~{trim_ends} \
      --pile-size ~{pile_size} \
      --min-var-frequency ~{min_var_frequency} \
      --max-alignments-per-read ~{max_alignments_per_read} \
      --max-reads-per-guide ~{max_reads_per_guide} \
      --iterations ~{iterations} \
      --seed ~{seed} \
      --max-consensus-reads ~{max_consensus_reads} \
      --max-amplicon-size ~{max_amplicon_size} \
      --min-read-qv ~{min_read_qv} \
      ~{if defined(off_target_groups) then "--off-target-groups " + off_target_groups else ""} \
      --min-cluster-frequency ~{min_cluster_frequency} \
      --min-cluster-read-count ~{min_cluster_read_count} \
      --max-uchime-score ~{max_uchime_score} \
      ~{if skip_chimera_detection then "--skip-chimera-detection" else ""} \
      ~{if skip_consensus then "--skip-consensus" else ""} \
      guide.fasta ~{fastq} "~{prefix}_pbaa"
  >>>

  output {
    File passed_cluster_sequences = "~{prefix}_pbaa_passed_cluster_sequences.fasta"
    File failed_cluster_sequences = "~{prefix}_pbaa_failed_cluster_sequences.fasta"
    File info                     = "~{prefix}_pbaa_read_info.txt"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task extract_reads_from_fastq {
  input {
    File info

    String? guide_name

    File fastq
    File fastq_index

    String docker_smrttools
  }

  Int threads = 1
  Int mem_gb  = 4

  command <<<
    if ~{defined(guide_name)}; then
      awk -F ' ' -v "guide=~{guide_name}" '$2 == guide {print $1}' ~{info} > clustered_holes.txt
    else
      cut -d' ' -f1 ~{info} > clustered_holes.txt
    fi

    samtools fqidx \
      --fai-idx ~{fastq_index} \
      --region-file clustered_holes.txt \
      ~{fastq} > clustered_hifi.fastq
    samtools fqidx clustered_hifi.fastq
  >>>

  output {
    File clustered_fastq       = "clustered_hifi.fastq"
    File clustered_fastq_index = "clustered_hifi.fastq.fai"
    File clustered_list        = "clustered_holes.txt"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task extract_reads_from_bam {
  input {
    File info

    File bam

    String docker_smrttools
  }

  Int threads = 1
  Int mem_gb  = 4

  command <<<
    cut -d' ' -f1 ~{info} > clustered_holes.txt

    samtools view --bam --with-header \
      --qname-file clustered_holes.txt \
      --output clustered_hifi.bam \
      ~{bam}
  >>>

  output {
    File clustered_bam     = "clustered_hifi.bam"
    File clustered_list    = "clustered_holes.txt"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task paint_bam {
  input {
    File info

    File bam
    File bam_index

    String out_prefix

    String docker_smrttools
  }

  Int threads = 2
  Int mem_gb  = 4

  command <<<
    pbaa bampaint ~{info} ~{bam} ~{out_prefix}.painted.bam
    samtools index -@ ~{threads - 1} ~{out_prefix}.painted.bam
  >>>

  output {
    File painted_bam     = "~{out_prefix}.painted.bam"
    File painted_bam_bai = "~{out_prefix}.painted.bam.bai"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task minipileup {
  input {
    String sample_name
    File bam
    File bam_index

    File ref_fasta
    File ref_index

    String out_prefix

    String docker_smrttools
  }

  Int threads = 1
  Int mem_gb  = 4

  command <<<
    samtools view --no-header ~{bam} \
      | cut -f1 | sort -u \
      > ~{out_prefix}.consensus_ids.txt

    while IFS= read -r consensus_id || [[ -n "$consensus_id" ]]; do
      echo "$consensus_id" > read_id.txt
      CLUSTER="${consensus_id//sample-clustered_hifi_guide-/}"
      echo "~{out_prefix}.${consensus_id}.bam	~{sample_name}.${CLUSTER}" > samples.txt
      samtools view --bam \
        --qname-file read_id.txt \
        --output "~{out_prefix}.${consensus_id}.bam" \
        ~{bam}
      minipileup -c \
        -f ~{ref_fasta} \
        "~{out_prefix}.${consensus_id}.bam" \
        | bcftools reheader --samples samples.txt - \
        > "~{out_prefix}.${CLUSTER}.vcf"
      rm read_id.txt samples.txt "~{out_prefix}.${consensus_id}.bam"
    done < ~{out_prefix}.consensus_ids.txt

    bcftools merge --no-index \
      --output ~{out_prefix}.consensus.vcf \
      ./*.vcf
  >>>

  output {
    File vcf = "~{out_prefix}.consensus.vcf"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task merge_bams {
  input {
    Array[File] bams
    Array[File] bam_indices

    String out_prefix

    String docker_smrttools
  }

  Int threads = 2
  Int mem_gb  = 4

  command <<<
    samtools merge --threads ~{threads - 1} ~{out_prefix}.bam ~{sep=' ' bams}
    samtools index -@ ~{threads - 1} ~{out_prefix}.bam
  >>>

  output {
    File bam     = "~{out_prefix}.bam"
    File bam_bai = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: mem_gb + " GB"
  }
}

task call_f8 {
  input {
    String? sex
    File bam
    File bam_index

    String out_prefix

    String docker_smrttools
  }

  Int threads = 1
  Int mem_gb  = 4

  String sample_sex = if select_first([sex, "F"]) == "M" then "M" else "F"

  command <<<
    f8_inversion.py \
      --bam ~{bam} \
      --prefix ~{out_prefix} \
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
    memory: mem_gb + " GB"
  }
}