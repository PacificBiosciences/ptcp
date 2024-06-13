version 1.0

task collect_inputs {
  meta {
    description: "Collect input reads and metadata for input reads."
  }

  input {
    File sample_sheet
    Array[File] hifi_reads
    Array[File]? fail_reads

    Int threads = 4
    Int mem_gb  = 4

    String docker_smrttools
  }

  command <<<
    set -e
    python3 <<EOF
    from pbcoretools.tasks.collect_noamp_inputs import collect_bam_files
    hifi_bams = "~{sep=';' hifi_reads}".split(";")
    fail_bams = "~{sep=';' fail_reads}".split(";") if "~{sep=';' fail_reads}" else []
    collect_bam_files(hifi_reads=hifi_bams,
                      fail_reads=fail_bams,
                      nproc=~{threads},
                      biosamples_out="biosamples.txt",
                      jasmine_tasks_out="jasmine_required.txt")
    EOF

    # if there's only one BAM for a sample, the call above symlinks rather than merging
    # TODO: temporary fix until collect_noamp_inputs copies BAMs instead of linking
    find ./ -name '*.bam' -type l \
      -exec sh -c 'for i in "$@"; do cp --preserve --remove-destination "$(readlink -f "$i")" "$i"; done' sh {} +

    # For each sample, look up the sex in the sample sheet
    while IFS= read -r i || [[ -n "$i" ]]; do
        parse_sample_sheet.py --output sex ~{sample_sheet} "$i" >> sample_sexes.txt
    done < biosamples.txt
  >>>

  output {
    Array[File] chunks                    = glob("*.unmapped.bam")
    Array[String] sample_names            = read_lines("biosamples.txt")
    Array[String] sample_sexes            = read_lines("sample_sexes.txt")
    Array[String] sample_requires_jasmine = read_lines("jasmine_required.txt")
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}

task jasmine {
  input {
    File reads_bam
    String sample_name

    Int threads = 8
    Int mem_gb  = 8

    String log_level = "DEBUG"
    String docker_smrttools
  }

  String bam_out = "~{sample_name}.unmapped.jasmine.bam"

  command <<<
    jasmine \
      --log-level ~{log_level} \
      --log-file jasmine.log \
      --alarms alarms.json \
      --num-threads ~{threads} \
      ~{reads_bam} ~{bam_out}
  >>>

  output {
    File methyl_bam = "~{bam_out}"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}

task pbmm2_align_hifi {
  input {
    String sample_name
    File unmapped

    File ref_fasta

    Int threads = 8
    Int mem_gb  = 16

    String out_prefix = "~{sample_name}.mapped"

    String log_level = "DEBUG"
    String docker_smrttools
  }

  command <<<
    set -vex

    pbmm2 \
      align \
      --log-level ~{log_level} \
      --log-file pbmm2.log \
      --alarms alarms.json \
      --num-threads ~{threads} \
      --sort \
      --preset HiFi \
      --report-json mapping_stats.report.json \
      --sample "~{sample_name}" \
      -N 1 --unmapped \
      ~{ref_fasta} \
      ~{unmapped} \
      ~{out_prefix}.bam

    pbindex --num-threads ~{threads} ~{out_prefix}.bam
  >>>

  output {
    File mapped_bam     = "~{out_prefix}.bam"
    File mapped_bam_bai = "~{out_prefix}.bam.bai"
    File mapped_bam_pbi = "~{out_prefix}.bam.pbi"
    File report         = "mapping_stats.report.json"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}

task pbmm2_index {
  input {
    File ref_fasta

    Int threads = 8
    Int mem_gb  = 16

    String log_level = "DEBUG"
    String docker_smrttools
  }

  String out_prefix = basename(ref_fasta, ".fasta")

  command <<<
    set -vex

    pbmm2 \
      index \
      --log-level ~{log_level} \
      --log-file pbmm2_index.log \
      --num-threads ~{threads} \
      --preset HiFi \
      ~{ref_fasta} \
      ~{out_prefix}.mmi
  >>>

  output {
    File index = "~{out_prefix}.mmi"
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}

task infer_chry {
  input {
    File mapped_bam
    File mapped_bam_bai

    Int threads = 1
    Int mem_gb  = 4

    String docker_smrttools
  }

  command <<<
    set -vex

    # total mapped reads
    mapped=$(samtools flagstat --output-fmt json ~{mapped_bam} | jq -r '.["QC-passed reads"]["mapped"]')
    echo "Total mapped reads: $mapped"

    chry_mapped=$(samtools idxstats ~{mapped_bam} | grep -w '^chrY' | cut -f 3)
    echo "Mapped reads on chrY: $chry_mapped"

    printf %.6f "$((10**8 * chry_mapped / mapped))e-8" > chry_frequency.txt # hack for low precision float without bc
    echo "chrY frequency: $(cat chry_frequency.txt)"
  >>>

  output {
    Float chry_frequency = read_float("chry_frequency.txt")
  }

  runtime {
    cpu: threads
    memory: "~{mem_gb} GB"
    docker: docker_smrttools
  }
}