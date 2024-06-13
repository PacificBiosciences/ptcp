version 1.0

import "./common.wdl" as Common
import "./trgt.wdl" as Trgt
import "./dark_genes.wdl" as DarkGenes

workflow puretarget_darkgenes {
  input {
    File sample_sheet
    Array[File] hifi_reads
    Array[File]? fail_reads

    File ref_fasta
    File ref_index

    File pbaa_guide_bed
    File pbaa_mask_bed

    File trgt_bed

    String log_level        = "INFO"
    String docker_smrttools = "quay.io/pacbio/smrttools@sha256:01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b"
  }

  call DarkGenes.make_pbaa_guides {
    # make pbaa guide fastas for each group of guides, e.g. F8_inv1, F8_inv22, etc.
    input:
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      guide_bed        = pbaa_guide_bed,
      docker_smrttools = docker_smrttools
  }

  # mask all desired sequences in the reference
  call DarkGenes.mask_reference {
    input:
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      mask_bed         = pbaa_mask_bed,
      docker_smrttools = docker_smrttools
  }

  call Common.pbmm2_index as masked_ref_index {
    # index the masked reference
    input:
      ref_fasta        = mask_reference.masked_ref,
      log_level        = log_level,
      docker_smrttools = docker_smrttools
  }

  call Common.collect_inputs {
    # collect reads per sample
    input:
      sample_sheet     = sample_sheet,
      hifi_reads       = hifi_reads,
      fail_reads       = fail_reads,
      docker_smrttools = docker_smrttools
  }
  
  scatter (i_chunk in range(length(collect_inputs.chunks))) {
    # for each sample
    File   chunk_bam        = collect_inputs.chunks[i_chunk]
    String sample_name      = collect_inputs.sample_names[i_chunk]
    String sample_sex       = collect_inputs.sample_sexes[i_chunk]
    String requires_jasmine = collect_inputs.sample_requires_jasmine[i_chunk]

    if (requires_jasmine == "true") {
      call Common.jasmine {
        # call basemods from kinetics if necessary
        input:
          reads_bam        = chunk_bam,
          sample_name      = sample_name,
          docker_smrttools = docker_smrttools
      }
    }

    File unmapped_bam = select_first([jasmine.methyl_bam, chunk_bam])
    call Common.pbmm2_align_hifi {
      # align reads to the masked reference
      input:
        sample_name      = sample_name,
        unmapped         = unmapped_bam,
        ref_fasta        = masked_ref_index.index,
        out_prefix       = "~{sample_name}.mapped.masked",
        log_level        = log_level,
        docker_smrttools = docker_smrttools
    }

    call Common.infer_chry {
      input:
        mapped_bam       = pbmm2_align_hifi.mapped_bam,
        mapped_bam_bai   = pbmm2_align_hifi.mapped_bam_bai,
        docker_smrttools = docker_smrttools
    }

    call Trgt.trgt {
      # run tandem repeat expansion analysis subworkflow
      input:
        sample_name      = sample_name,
        sex              = sample_sex,
        mapped_bam       = pbmm2_align_hifi.mapped_bam,
        mapped_bam_bai   = pbmm2_align_hifi.mapped_bam_bai,
        ref_fasta        = mask_reference.masked_ref,
        ref_index        = mask_reference.masked_index,
        trgt_bed         = trgt_bed,
        docker_smrttools = docker_smrttools
    }

    call DarkGenes.dark_genes {
      # run dark genes analysis subworkflow
      input:
        sample_name      = sample_name,
        sex              = sample_sex,
        unmapped_bam     = unmapped_bam,
        guides           = make_pbaa_guides.guides,
        guide_indices    = make_pbaa_guides.guide_indices,
        ref_fasta        = mask_reference.masked_ref,
        ref_index        = mask_reference.masked_index,
        ref_mmi          = masked_ref_index.index,
        docker_smrttools = docker_smrttools
    }
  }

  Array[File] combined_reads = flatten(select_all([hifi_reads, fail_reads]))
  call Trgt.pbreports_target_enrichment_noamp {
    # generate target enrichment report
    input:
      demuxed_ccs      = combined_reads,
      mapped_bams      = pbmm2_align_hifi.mapped_bam,
      mapped_bam_pbis  = pbmm2_align_hifi.mapped_bam_pbi,
      bed_files        = [trgt_bed, pbaa_guide_bed],
      log_level        = log_level,
      docker_smrttools = docker_smrttools
  }

  call Trgt.noamp_qc {
    # generate noamp qc report
    input:
      sample_names     = collect_inputs.sample_names,
      mapped_bams      = pbmm2_align_hifi.mapped_bam,
      vcfs             = trgt.vcf,
      spanning_bams    = trgt.spanning_bam,
      ref_fasta        = mask_reference.masked_ref,
      trgt_bed         = trgt_bed,
      docker_smrttools = docker_smrttools
  }

  output {
    File repeats     = trgt_bed
    File pbaa_guides = pbaa_guide_bed
    File pbaa_mask   = pbaa_mask_bed

    Array[Float] chry_frequency  = infer_chry.chry_frequency

    Array[File] vcf                                 = trgt.vcf
    Array[File] spanning_bam                        = trgt.spanning_bam
    Array[File] spanning_index                      = trgt.spanning_index
    Array[File] images_motifs_allele                = trgt.images_motifs_allele
    Array[File] images_meth_allele                  = trgt.images_meth_allele
    Array[File] images_motifs_waterfall             = trgt.images_motifs_waterfall
    Array[File] images_meth_waterfall               = trgt.images_meth_waterfall
    Array[File] reads_overlapping_repeats           = trgt.reads_overlapping_repeats
    Array[File] reads_overlapping_repeats_index     = trgt.reads_overlapping_repeats_index
    Array[File?] report_target_enrichment_resources = pbreports_target_enrichment_noamp.plot_pngs
    File  genotype                                  = noamp_qc.genotype
    File  report_target_enrichment                  = pbreports_target_enrichment_noamp.report
    File  report_puretarget                         = noamp_qc.report

    Array[File] painted_bams       = dark_genes.painted_bam
    Array[File] painted_bam_bais   = dark_genes.painted_bam_bai
    Array[File] consensus_bams     = dark_genes.consensus_bam
    Array[File] consensus_bam_bais = dark_genes.consensus_bam_bai
    Array[File] minipileup_vcfs    = dark_genes.minipileup_vcf
    Array[File] f8_pbaa_vcf        = dark_genes.f8_pbaa_vcf
    Array[File] f8_pbaa_json       = dark_genes.f8_pbaa_json
  }
}
