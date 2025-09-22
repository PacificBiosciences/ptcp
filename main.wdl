version 1.0

import "./common.wdl" as Common
import "./trgt.wdl" as Trgt
import "./paraphase.wdl" as Paraphase
import "./sawfish.wdl" as Sawfish

workflow ptcp {
  meta {
    title: "PureTarget tandem repeat and hard gene workflow"
    summary: "Workflow for targeted analysis of tandem repeats and hard genes"
    description: "Workflow for analyzing hard genes and repeat regions using PacBio data. Combines multiple tools including TRGT for repeat genotyping and Paraphase for haplotype phasing and (small) variant calling."
    version: "3.1.1"
  }

  parameter_meta {
    sample_sheet: {
      help: "TSV file containing sample metadata including sex information",
      label: "Sample sheet"
    }
    hifi_reads: {
      help: "Array of input BAM files",
      label: "HiFi reads"
    }
    fail_reads: {
      help: "Optional array of fail reads BAM files",
      label: "Fail reads"
    }
    ref_fasta: {
      help: "Reference genome in FASTA format",
      label: "Reference FASTA"
    }
    ref_index: {
      help: "Index file for the reference FASTA",
      label: "Reference FASTA index"
    }
    trgt_bed: {
      help: "BED file containing tandem repeat regions to analyze",
      label: "TRGT catalog BED"
    }
    paraphase_config_yaml: {
      help: "Configuration file for Paraphase analysis",
      label: "Paraphase config"
    }
    paraphase_annotation_vcf: {
      help: "Optional VCF file for havanno annotation of paraphase VCF outputs",
      label: "Paraphase annotation VCF"
    }
    genome_version: {
      help: "Reference genome version (e.g., GRCh38)",
      label: "Genome version"
    }
    ptcp_qc_bed: {
      help: "BED file defining target regions for QC",
      label: "QC Targets BED"
    }
    pt_linear_regression: {
      help: "Optional linear regression model file for SMN homology analysis.",
      label: "SMN Linear Regression Model"
    }
    log_level: {
      help: "Logging verbosity level (default: INFO)",
      label: "Log level"
    }
    # Output parameters
    sample_names: {
      help: "Array of processed sample names",
      label: "Sample names"
    }
    mapped_bam: {
      help: "Array of aligned BAM files",
      label: "Aligned BAMs"
    }
    mapped_bam_bai: {
      help: "Array of BAM index files",
      label: "BAM indices"
    }
    trgt_vcf: {
      help: "Array of VCF files containing repeat genotypes",
      label: "TRGT VCFs"
    }
    trgt_spanning_bam: {
      help: "Array of BAM files containing reads spanning repeat regions",
      label: "TRGT spanning BAMs"
    }
    images_motifs_allele: {
      help: "Array of archives containing repeat motif allele plots",
      label: "TRGT motif allele plots"
    }
    images_meth_allele: {
      help: "Array of archives containing methylation allele plots",
      label: "TRGT methylation allele plots"
    }
    images_motifs_waterfall: {
      help: "Array of archives containing repeat motif waterfall plots",
      label: "TRGT motif waterfall plots"
    }
    images_meth_waterfall: {
      help: "Array of archives containing methylation waterfall plots",
      label: "TRGT methylation waterfall plots"
    }
    paraphase_bam: {
      help: "Array of phased BAM files from Paraphase",
      label: "Paraphase phased BAMs"
    }
    paraphase_json: {
      help: "Array of JSON files from Paraphase",
      label: "Paraphase JSON"
    }
    paraphase_vcfs: {
      help: "Array of arrays containing phased VCF files",
      label: "Paraphase phased VCFs"
    }
    havanno_json: {
      help: "Array of JSON files containing havanno annotations (optional)",
      label: "Paraphase havanno annotations"
    }
    f8_vcf: {
      help: "Array of VCF files containing F8 inversion calls",
      label: "F8 VCF"
    }
    f8_json: {
      help: "Array of JSON files containing F8 inversion calls",
      label: "F8 JSON"
    }
    ptcp_qc_reports: {
      help: "Array of QC report files containing ptcp-qc statistics",
      label: "ptcp-qc stats"
    }
    sawfish_vcf: {
      help: "Array of VCF files containing structural variants",
      label: "Sawfish VCFs"
    }
    sawfish_vcf_tbi: {
      help: "Array of index files for structural variant VCFs",
      label: "Sawfish VCF indices"
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
    File sample_sheet
    Array[File] hifi_reads
    Array[File]? fail_reads

    File ref_fasta
    File ref_index

    File trgt_bed

    File paraphase_config_yaml
    File? paraphase_annotation_vcf

    String genome_version

    File ptcp_qc_bed
    File? pt_linear_regression

    String log_level        = "INFO"
    String docker_smrttools = "quay.io/pacbio/smrttools@sha256:01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b"
  }

  call Common.collect_inputs {
    input:
      sample_sheet     = sample_sheet,
      hifi_reads       = hifi_reads,
      fail_reads       = fail_reads,
      docker_smrttools = docker_smrttools
  }

  scatter (i_chunk in range(length(collect_inputs.chunks))) {
    File   chunk_bam        = collect_inputs.chunks[i_chunk]
    String sample_name      = collect_inputs.sample_names[i_chunk]
    String sample_sex       = collect_inputs.sample_sexes[i_chunk]

    File unmapped_bam = chunk_bam
    call Common.pbmm2_align_hifi {
      input:
        sample_name      = sample_name,
        unmapped_bam     = unmapped_bam,
        ref_fasta        = ref_fasta,
        out_prefix       = "~{sample_name}.mapped",
        log_level        = log_level,
        docker_smrttools = docker_smrttools
    }

    call Trgt.trgt {
      input:
        sample_name      = sample_name,
        sex              = sample_sex,
        mapped_bam       = pbmm2_align_hifi.mapped_bam,
        mapped_bam_bai   = pbmm2_align_hifi.mapped_bam_bai,
        ref_fasta        = ref_fasta,
        ref_index        = ref_index,
        trgt_bed         = trgt_bed,
        docker_smrttools = docker_smrttools
    }

    call Paraphase.paraphase {
      input:
        sample_name                = sample_name,
        sex                        = sample_sex,
        mapped_bam                 = pbmm2_align_hifi.mapped_bam,
        mapped_bam_bai             = pbmm2_align_hifi.mapped_bam_bai,
        ref_fasta                  = ref_fasta,
        ref_index                  = ref_index,
        config_file                = paraphase_config_yaml,
        genome_version             = genome_version,
        annotation_vcf             = paraphase_annotation_vcf,
        docker_smrttools           = docker_smrttools
    }

    call Sawfish.sawfish {
      input:
        sample_name      = sample_name,
        mapped_bam       = pbmm2_align_hifi.mapped_bam,
        mapped_bam_bai   = pbmm2_align_hifi.mapped_bam_bai,
        ref_fasta        = ref_fasta,
        ref_index        = ref_index,
        genome_version   = genome_version,
        docker_smrttools = docker_smrttools
    }

    call Common.tar_outputs {
      input:
        sample_name           = sample_name,
        mapped_bam            = pbmm2_align_hifi.mapped_bam,
        mapped_bam_bai        = pbmm2_align_hifi.mapped_bam_bai,
        trgt_vcf              = trgt.vcf,
        trgt_spanning_bam     = trgt.spanning_bam,
        trgt_spanning_bam_bai = trgt.spanning_bam_bai,
        paraphase_json        = paraphase.paraphase_json,
        paraphase_bam         = paraphase.paraphase_bam,
        paraphase_bam_bai     = paraphase.paraphase_bam_bai,
        f8_json               = paraphase.f8_json,
        sawfish_vcf           = sawfish.vcf,
        sawfish_vcf_tbi       = sawfish.vcf_tbi,
        docker_smrttools      = docker_smrttools
    }
  }

  call Common.ptcp_qc {
      input:
        tarballs             = tar_outputs.sample_tarball,
        qc_bed               = ptcp_qc_bed,
        pt_linear_regression = pt_linear_regression,
        docker_smrttools     = docker_smrttools
  }
  

  output {
    Array[String] sample_names                      = collect_inputs.sample_names

    Array[File] mapped_bam                          = pbmm2_align_hifi.mapped_bam
    Array[File] mapped_bam_bai                      = pbmm2_align_hifi.mapped_bam_bai

    Array[File] trgt_vcf                            = trgt.vcf
    Array[File] trgt_spanning_bam                   = trgt.spanning_bam
    Array[File] trgt_spanning_bam_bai               = trgt.spanning_bam_bai
    Array[File] images_motifs_allele                = trgt.images_motifs_allele
    Array[File] images_meth_allele                  = trgt.images_meth_allele
    Array[File] images_motifs_waterfall             = trgt.images_motifs_waterfall
    Array[File] images_meth_waterfall               = trgt.images_meth_waterfall
    Array[File] reads_overlapping_repeats           = trgt.reads_overlapping_repeats
    Array[File] reads_overlapping_repeats_bai       = trgt.reads_overlapping_repeats_bai

    Array[File] paraphase_bam                       = paraphase.paraphase_bam
    Array[File] paraphase_bam_bai                   = paraphase.paraphase_bam_bai
    Array[File] paraphase_json                      = paraphase.paraphase_json
    Array[Array[File]] paraphase_vcfs               = paraphase.paraphase_vcfs
    Array[File?] havanno_json                       = paraphase.havanno_json
    Array[File] f8_vcf                              = paraphase.f8_vcf
    Array[File] f8_json                             = paraphase.f8_json

    Array[File] ptcp_qc_reports                     = ptcp_qc.qc_reports

    Array[File] sawfish_vcf                         = sawfish.vcf
    Array[File] sawfish_vcf_tbi                     = sawfish.vcf_tbi
    Array[File] sawfish_minimap2_bam                = sawfish.minimap2_bam
    Array[File] sawfish_minimap2_bam_bai            = sawfish.minimap2_bam_bai
  }
}
