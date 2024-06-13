version 1.0

workflow trgt {
  input {
    String sample_name
    String sex
    File mapped_bam
    File mapped_bam_bai

    File ref_fasta
    File ref_index

    File trgt_bed

    String docker_smrttools = "quay.io/pacbio/smrttools@sha256:01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b"
  }

  call get_reads_overlapping_repeats {
    # extract reads that overlap the trgt bed
    input:
      sample_name      = sample_name,
      mapped_bam       = mapped_bam,
      mapped_bam_bai   = mapped_bam_bai,
      trgt_bed         = trgt_bed,
      docker_smrttools = docker_smrttools
  }

  call trgt_genotype {
    # genotype the sample
    input:
      sample_name      = sample_name,
      sex              = sex,
      mapped_bam       = mapped_bam,
      mapped_bam_bai   = mapped_bam_bai,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      docker_smrttools = docker_smrttools
  }

  # call all 4 combinations of trgt plots: (repeat, methylation) x (allele, waterfall)
  call trgt_plot as trgt_plot_motifs_allele {
    input:
      sample_name      = sample_name,
      vcf              = trgt_genotype.vcf,
      spanning         = trgt_genotype.spanning_bam,
      spanning_index   = trgt_genotype.spanning_index,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      methylation      = false,
      waterfall        = false,
      docker_smrttools = docker_smrttools
  }

  call trgt_plot as trgt_plot_meth_allele {
    input:
      sample_name      = sample_name,
      vcf              = trgt_genotype.vcf,
      spanning         = trgt_genotype.spanning_bam,
      spanning_index   = trgt_genotype.spanning_index,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      methylation      = true,
      waterfall        = false,
      docker_smrttools = docker_smrttools
  }

  call trgt_plot as trgt_plot_motifs_waterfall {
    input:
      sample_name      = sample_name,
      vcf              = trgt_genotype.vcf,
      spanning         = trgt_genotype.spanning_bam,
      spanning_index   = trgt_genotype.spanning_index,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      methylation      = false,
      waterfall        = true,
      docker_smrttools = docker_smrttools
  }

  call trgt_plot as trgt_plot_meth_waterfall {
    input:
      sample_name      = sample_name,
      vcf              = trgt_genotype.vcf,
      spanning         = trgt_genotype.spanning_bam,
      spanning_index   = trgt_genotype.spanning_index,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      methylation      = true,
      waterfall        = true,
      docker_smrttools = docker_smrttools
  }

  output {
    File vcf                             = trgt_genotype.vcf
    File spanning_bam                    = trgt_genotype.spanning_bam
    File spanning_index                  = trgt_genotype.spanning_index
    File images_motifs_allele            = trgt_plot_motifs_allele.images
    File images_meth_allele              = trgt_plot_meth_allele.images
    File images_motifs_waterfall         = trgt_plot_motifs_waterfall.images
    File images_meth_waterfall           = trgt_plot_meth_waterfall.images
    File reads_overlapping_repeats       = get_reads_overlapping_repeats.repeats_reads
    File reads_overlapping_repeats_index = get_reads_overlapping_repeats.repeats_index
  }
}

task get_reads_overlapping_repeats {
  input {
    String sample_name
    File mapped_bam
    File mapped_bam_bai

    File trgt_bed

    Int expand_size = 3000

    Int threads = 1
    Int mem_gb  = 4

    String docker_smrttools
  }

  String output_name = "~{sample_name}.pbmm2.repeats.bam"
  String index_name = "~{output_name}.bai"
  String output_name_with_index = "~{output_name}##idx##~{index_name}"
  String expanded_bed_name = "~{sample_name}.expanded.bed"

  command <<<
    set -e

    awk -v sz=~{expand_size} \
      '{if ($2 < sz) { print $1"\t0\t"$3 + sz} else {print $1"\t"$2 - sz"\t"$3 + sz} }' \
      ~{trgt_bed} \
      > ~{expanded_bed_name}

    samtools view \
      --write-index \
      --use-index \
      --target-file ~{expanded_bed_name} \
      --output ~{output_name_with_index} \
      ~{mapped_bam}
  >>>

  output {
    File repeats_reads = "~{output_name}"
    File repeats_index = "~{index_name}"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}

task trgt_genotype {
  input {
    String sample_name
    String? sex
    File mapped_bam
    File mapped_bam_bai

    File ref_fasta
    File ref_index

    File trgt_bed

    Int threads = 8
    Int mem_gb  = 16

    String docker_smrttools
  }

  String prefix = "~{sample_name}.trgt"
  String karyotype = if select_first([sex, "F"]) == "M" then "XY" else "XX"

  command <<<
    set -e

    trgt --verbose genotype \
      --threads ~{threads} \
      --genome ~{ref_fasta} \
      --karyotype ~{karyotype} \
      --flank-len=300 \
      --genotyper=cluster \
      --aln-scoring=1,1,0,1,1,200 \
      --min-read-quality=-1.0 \
      --max-depth=10000 \
      --reads ~{mapped_bam} \
      --repeats ~{trgt_bed} \
      --output-prefix ~{prefix}

    gunzip -c "~{prefix}.vcf.gz" > ~{prefix}.vcf
    samtools sort -o ~{prefix}.sorted.spanning.bam ~{prefix}.spanning.bam
    samtools index ~{prefix}.sorted.spanning.bam
  >>>

  output {
    File vcf            = "~{prefix}.vcf"
    File spanning_bam   = "~{prefix}.sorted.spanning.bam"
    File spanning_index = "~{prefix}.sorted.spanning.bam.bai"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}

task trgt_plot {
  input {
    String sample_name
    File vcf
    File spanning
    File spanning_index

    File ref_fasta
    File ref_index

    File trgt_bed

    Boolean methylation = false
    Boolean waterfall = false

    Int threads = 1
    Int mem_gb  = 4

    String docker_smrttools
  }

  String show          = if (methylation) then "meth" else "motifs"
  String plot_type     = if (waterfall) then "waterfall" else "allele"
  String output_prefix = "~{sample_name}.~{show}_~{plot_type}"

  command <<<
    grep -v "^#" ~{vcf} | awk -F"[=:\t;]" '{print $9}' | while read -r i; do
      trgt plot \
        --repeat-id "$i" \
        --genome ~{ref_fasta} \
        --repeats ~{trgt_bed} \
        --spanning-reads ~{spanning} \
        --vcf ~{vcf} \
        --show ~{show} \
        --plot-type ~{plot_type} \
        --image "~{output_prefix}.$(printf '%s' "$i").trgt_plot.svg"
    done;
    python3 -m zipfile -c ~{output_prefix}.trgt_plot_alleles.zip ./*.trgt_plot.svg
  >>>

  output {
    File images = "~{output_prefix}.trgt_plot_alleles.zip"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}

task pbreports_target_enrichment_noamp {
  input {
    Array[File] demuxed_ccs

    Array[File] mapped_bams
    Array[File] mapped_bam_pbis  # !UnusedDeclaration

    Array[File] bed_files

    Int threads = 1
    Int mem_gb  = 8

    String log_level = "INFO"
    String docker_smrttools
  }

  command <<<
    set -e

    # combine all bed files
    # NOTE: make sure bed files end in a new line
    cat ~{sep=" " bed_files} > all.bed

    # FIXME hacky workaround for current report behavior
    mkdir input_reads
    for BAM in ~{sep=" " demuxed_ccs}; do
      FTMP=$(python3 -c "import uuid; print(str(uuid.uuid4()).split('-')[0])")
      ln -s "$BAM" "input_reads/${FTMP}.bam"
      if [ -f "${BAM}.pbi" ]; then
        ln -s "${BAM}.pbi" "input_reads/${FTMP}.bam.pbi"
      else
        pbindex "input_reads/${FTMP}.bam"
      fi
    done

    dataset --strict --log-level DEBUG create --type ConsensusReadSet \
      input_reads.consensusreadset.xml input_reads/*.bam
    # end of hacky workaround

    python3 -m pbreports.report.target_enrichment_noamp \
      --log-level ~{log_level} \
      --log-file pbreports.log \
      -o target_enrichment_noamp.report.json \
      --boxplots-target-order bedfile \
      --target-buffer 0 \
      --boxplots-boxes-per-panel 25 \
      --ignore-duplicates \
      --plot-read-coverage \
      input_reads.consensusreadset.xml \
      all.bed \
      ~{sep=" " mapped_bams}
  >>>

  output {
    File  report           = "target_enrichment_noamp.report.json"
    Array[File?] plot_pngs = glob("*.png")
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}

task noamp_qc {
  input {
    Array[String] sample_names
    Array[File] mapped_bams
    Array[File] vcfs
    Array[File] spanning_bams

    File ref_fasta

    File trgt_bed

    Int threads = 1
    Int mem_gb  = 8

    String docker_smrttools
  }

  String out_prefix = "puretarget_report"

  command <<<
    set -e
    mapped_bams=$(for i in ~{sep=" " mapped_bams}; do readlink -f ${i}; done | tr '\n' ',' | sed 's/,$//');
    spanning_bams=$(for i in ~{sep=" " spanning_bams}; do readlink -f ${i}; done | tr '\n' ',' | sed 's/,$//');
    vcfs=$(for i in ~{sep=" " vcfs}; do readlink -f ${i}; done | tr '\n' ',' | sed 's/,$//');
    noamp-qc --smrtlink \
      --genome ~{ref_fasta} \
      --catalog ~{trgt_bed} \
      --sample-names ~{sep="," sample_names} \
      --mapped-list "${mapped_bams}" \
      --spanning-list "${spanning_bams}" \
      --vcf-list "${vcfs}" \
      --output-prefix ~{out_prefix}
  >>>

  output {
    File report   = "~{out_prefix}.report.json"
    File genotype = "~{out_prefix}.genotype.csv"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}
