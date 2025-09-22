version 1.0

task collect_inputs {
  meta {
    title: "Collect and process input BAMs"
    summary: "Collects (and validates) input BAM files and sample metadata"
    description: "Processes input BAM files, validates their format, merges multiple BAMs per sample if needed, and extracts sample metadata."
  }

  parameter_meta {
    sample_sheet: {
      help: "TSV file containing sample metadata",
      label: "Sample sheet"
    }
    hifi_reads: {
      help: "Array of BAM files",
      label: "Reads"
    }
    fail_reads: {
      help: "Optional array of fail reads BAM files",
      label: "Fail reads"
    }
    threads: {
      help: "Number of CPU threads to use (default: 4)",
      label: "CPU threads"
    }
    mem_gb: {
      help: "Memory allocation in gigabytes (default: 4)",
      label: "Memory (GB)"
    }
    chunks: {
      help: "Output array of processed BAM files",
      label: "Processed BAMs"
    }
    sample_names: {
      help: "Array of sample names extracted from BAMs",
      label: "Sample names"
    }
    sample_sexes: {
      help: "Array of sample sexes from sample sheet",
      label: "Sample sexes"
    }
  }

  input {
    File sample_sheet
    Array[File] hifi_reads
    Array[File]? fail_reads

    Int threads = 4
    Int mem_gb  = 4

    String docker_smrttools
  }

  Int disk_size = ceil(size(hifi_reads, 'GB') * 3 + size(select_first([fail_reads, []]), 'GB') * 3 + 20)

  # TEMP: This Python code will be moved to docker/scripts/ when stable
  command <<<
    set -e
    python3 <<EOF
    import os
    import re
    import subprocess
    import shutil
    import logging as log
    from collections import defaultdict
    from pathlib import Path
    from typing import List, Union
    from pbcore.io import BamReader

    log.basicConfig(level=log.INFO)

    def validate_bam_file(bam_file: Path):
        command = ["samtools", "quickcheck", "-u", "-vvvvv", str(bam_file)]
        try:
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"samtools quickcheck failed for {bam_file} with error:\n{result.stderr.strip()}"
                )
        except FileNotFoundError:
            raise RuntimeError("samtools is not installed or not in PATH.")

    def collect_bam_files(
        hifi_bams: List[str],
        fail_bams: List[str],
        nproc: int,
        biosamples_out: Union[Path, str],
    ) -> None:
        hifi_fail_pattern = re.compile(r"\.(hifi_reads|fail_reads)\.")
        reads_pattern = re.compile(r"\.reads$")
        bam_map = defaultdict(list)
        missing_kinetics = []

        for bam_file_str in list(hifi_bams) + list(fail_bams):
            bam_file = Path(bam_file_str)
            validate_bam_file(bam_file)

            bam_name = hifi_fail_pattern.sub(".", bam_file.name)
            bam_name = reads_pattern.sub("", bam_name)
            bam_name = bam_name.rstrip(".bam")

            sample_names = set([])
            bam = BamReader(bam_file)
            for rg in bam.readGroupTable:
                if len(rg.StringID.split("/")) < 2:
                    log.warning(f"The BAM file: {bam_file} possibly has not been through demultiplexing.")
                sample_names.add(rg.SampleName)
            if len(sample_names) > 1:
                samples_str = ", ".join(list(sample_names))
                log.warning(f"Multiple samples found in a single BAM file: {bam_file} -> {samples_str}.")

            programs = {prog["ID"] for prog in bam.peer.header["PG"]}
            if not "jasmine" in programs:
                if not bam.hasBaseFeature("Ipd"):
                    log.warning(f"Warning: Analysis output will not contain methylyation calls due to missing kinetics in {bam_name}.")
                    missing_kinetics.append(bam_name)
                else:
                    log.warning(f"Warning: Jasmine has not been run on {bam_name}.")
            bam_map[bam_name].append(bam_file)

        all_sample_names = []
        for i, (sample, bam_files) in enumerate(bam_map.items(), start=1):
            bam_name = f"sample{i:04d}_{sample}.unmapped.bam"
            if len(bam_files) == 1:
                log.info(f"Found a single bam for sample {sample}")
                shutil.copyfile(bam_files[0], bam_name)
            else:
                log.info(f"Merging {len(bam_files)} BAMs for sample {sample}")
                args = ["pbmerge", "-j", str(nproc), "-o", bam_name] + bam_files
                subprocess.check_call(args)
            all_sample_names.append(sample)

        with open(biosamples_out, "wt") as sample_out:
            sample_out.write("\n".join(all_sample_names))
        log.info(f"Wrote list of sample names to biosamples.txt")

    hifi_bams = "~{sep=';' hifi_reads}".split(";")
    fail_bams = "~{sep=';' select_first([fail_reads,[]])}".split(";") if "~{sep=';' select_first([fail_reads,[]])}" else []
    collect_bam_files(hifi_bams=hifi_bams,
                      fail_bams=fail_bams,
                      nproc=~{threads},
                      biosamples_out="biosamples.txt")
    EOF

    cat biosamples.txt | parse_sample_sheet.py --batch-input --output sex ~{sample_sheet} > sample_sexes.txt
  >>>

  output {
    Array[File] chunks                    = glob("*.unmapped.bam")
    Array[String] sample_names            = read_lines("biosamples.txt")
    Array[String] sample_sexes            = read_lines("sample_sexes.txt")
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}

task pbmm2_align_hifi {
  meta {
    title: "PBMM2"
    summary: "Align reads to reference genome using pbmm2"
    description: "Aligns reads to a reference genome using pbmm2, with specific presets for HiFi data. Generates sorted BAM output with multiple indices."
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
    }
    unmapped_bam: {
      help: "Input unmapped BAM file",
      label: "Unmapped BAM"
    }
    ref_fasta: {
      help: "Reference genome in FASTA format",
      label: "Reference FASTA"
    }
    threads: {
      help: "Number of CPU threads to use (default: 64)",
      label: "CPU threads"
    }
    mem_gb: {
      help: "Memory allocation in gigabytes (default: 64)",
      label: "Memory (GB)"
    }
    out_prefix: {
      help: "Prefix for output files",
      label: "Output prefix"
    }
    log_level: {
      help: "Logging verbosity level",
      label: "Log level"
    }
    mapped_bam: {
      help: "Output aligned BAM file",
      label: "Aligned BAM"
    }
    mapped_bam_bai: {
      help: "BAM index file",
      label: "BAM index"
    }
    mapped_bam_pbi: {
      help: "PacBio BAM index file",
      label: "PacBio index"
    }
    report: {
      help: "Alignment statistics report in JSON format",
      label: "Alignment report"
    }
  }

  input {
    String sample_name
    File unmapped_bam

    File ref_fasta

    Int threads = 64
    Int mem_gb  = 64

    String out_prefix = "~{sample_name}.mapped"

    String log_level = "DEBUG"
    String docker_smrttools
  }

  Int disk_size = ceil((size(unmapped_bam, 'GB') + size(ref_fasta, 'GB')) * 3 + 20)

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
      --unmapped \
      ~{ref_fasta} \
      ~{unmapped_bam} \
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
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}


task tar_outputs {
  meta {
    title: "Tarball sample outputs"
    summary: "Collects various output files for a sample and archives them into a single tarball."
    description: "This task takes multiple output files generated by upstream processes and combines them into a single .tar archive"
 }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
    }
    mapped_bam: {
      help: "Input aligned BAM file",
      label: "Aligned BAM"
    }
    mapped_bam_bai: {
      help: "Index for the aligned BAM file",
      label: "Aligned BAM index"
    }
    trgt_vcf: {
      help: "VCF file containing TRGT genotypes",
      label: "TRGT VCF"
    }
    trgt_spanning_bam: {
      help: "BAM file containing reads spanning TRGT regions",
      label: "TRGT spanning BAM"
    }
    trgt_spanning_bam_bai: {
      help: "Index for the TRGT spanning BAM file",
      label: "TRGT spanning BAM index"
    }
    paraphase_json: {
      help: "JSON report from Paraphase",
      label: "Paraphase report"
    }
    paraphase_bam: {
      help: "Phased BAM file from Paraphase",
      label: "Phased BAM"
    }
    paraphase_bam_bai: {
      help: "Index for the phased BAM file",
      label: "Phased BAM index"
    }
    f8_json: {
      help: "JSON report for F8 inversion calls",
      label: "F8 inversion report"
    }
    sawfish_vcf: {
      help: "VCF file containing structural variants from Sawfish",
      label: "Sawfish VCF"
    }
    sawfish_vcf_tbi: {
      help: "Index for the Sawfish VCF file",
      label: "Sawfish VCF index"
    }
    sample_tarball: {
      help: "Output tarball containing all sample results",
      label: "Sample tarball"
    }
  }

  input {
    String sample_name
    File mapped_bam
    File mapped_bam_bai
    File trgt_vcf
    File trgt_spanning_bam
    File trgt_spanning_bam_bai
    File paraphase_json
    File paraphase_bam
    File paraphase_bam_bai
    File f8_json
    File sawfish_vcf
    File sawfish_vcf_tbi
    String docker_smrttools

    Int mem_gb = 4
    Int threads = 1
  }

  Int disk_size = ceil(size(mapped_bam, "GB") + size(trgt_spanning_bam, "GB") + size(paraphase_bam, "GB") + 5)
  String output_tarball = "~{sample_name}.tar"

   command <<<
     set -e
     tar --create --verbose --file ~{output_tarball} \
      --transform 's|.*/||g' \
      ~{mapped_bam} \
      ~{mapped_bam_bai} \
      ~{trgt_vcf} \
      ~{trgt_spanning_bam} \
      ~{trgt_spanning_bam_bai} \
      ~{paraphase_json} \
      ~{paraphase_bam} \
      ~{paraphase_bam_bai} \
      ~{f8_json} \
      ~{sawfish_vcf} \
      ~{sawfish_vcf_tbi}
   >>>

   output {
     File sample_tarball = output_tarball
   }

   runtime {
     docker: docker_smrttools
     cpu: threads
     memory: "~{mem_gb} GB"
     disk: disk_size + " GB"
     disks: "local-disk " + disk_size + " SSD"
   }
 }


task ptcp_qc {
  meta {
    title: "PTCP-QC analysis"
    summary: "Runs pctcp-qc on a collection of sample outputs."
    description: "This task takes tarballs containing results for multiple samples, extracts them, generates a input JSON, and runs ptcp-qc to produce QC reports."
  }
 
  parameter_meta {
    tarballs: { 
      help: "Array of tar archives, each containing outputs for one sample.",
      label: "Input tarballs"
    }
    qc_bed: { 
      help: "BED file defining target regions for QC.",
      label: "QC Targets BED"
    }
    pt_linear_regression: {
      help: "Optional linear regression model file for SMN homology analysis.",
      label: "SMN Linear Regression Model"
    }
    qc_reports: {
      help: "Array of JSON files containing the QC reports generated by ptcp-qc.",
      label: "QC Reports" 
    }
  }
 
  input {
    Array[File] tarballs
    File qc_bed
    File? pt_linear_regression
    String docker_smrttools

    Int threads = 16
    Int mem_gb = 8
  }

  Int disk_size = length(tarballs) * 2

  String manifest_filename = "input_manifest.json"
  String output_dir = "ptcp_qc"
  String output_prefix = "qc"

  command <<<
    set -e
    mkdir -p ~{output_dir}
    mkdir -p extracted_tarballs

    for tar in ~{sep=' ' tarballs}; do
      base=$(basename ${tar} .tar)
      mkdir -p extracted_tarballs/${base}
      tar -xf ${tar} -C extracted_tarballs/${base}
    done

    python3 <<EOF
import json
import os
from collections import defaultdict

manifest = defaultdict(list)
extracted_dir = "extracted_tarballs"

for sample in os.listdir(extracted_dir):
    sample_path = os.path.join(extracted_dir, sample)
    if os.path.isdir(sample_path):
        manifest["sample_names"].append(sample)
        for root, _, files in os.walk(sample_path):
            for f in files:
                fpath = os.path.join(root, f)
                if f.endswith("mapped.bam"):
                    manifest["mapped_bams"].append(fpath)
                elif f.endswith("trgt.vcf"):
                    manifest["trgt_vcfs"].append(fpath)
                elif f.endswith("trgt.sorted.spanning.bam"):
                    manifest["trgt_spanning_bams"].append(fpath)
                elif f.endswith("paraphase.bam"):
                    manifest["paraphase_bams"].append(fpath)
                elif f.endswith("paraphase.json"):
                    manifest["paraphase_jsons"].append(fpath)
                elif f.endswith("paraphase.vcf"):
                    manifest["paraphase_vcfs"].append(fpath)
                elif f.endswith("f8inversion.json"):
                    manifest["f8_jsons"].append(fpath)
                elif f.endswith("sv.vcf.gz"):
                    manifest["sawfish_vcfs"].append(fpath)


with open("~{manifest_filename}", "w") as f:
    json.dump(manifest, f, indent=2)
print("Created manifest at ~{manifest_filename}")
EOF
    ptcp-qc -v analyze \
      --manifest ~{manifest_filename} \
      --targets-bed ~{qc_bed} \
      --output-prefix ~{output_dir}/~{output_prefix} \
      --threads ~{threads} \
      --json

    ~{if defined(pt_linear_regression) then
  "for json_file in " + output_dir + "/" + output_prefix + "*.json; do\n" +
  "  [ -f \"$json_file\" ] || continue\n" +
  "  [[ \"$json_file\" == *aggregate* ]] && continue\n" +
  "  smn_homology test --model " + pt_linear_regression +
  " --dataset \"$json_file\" > \"$json_file.tmp\"\n" +
  "  mv \"$json_file.tmp\" \"$json_file\"\n" +
  "done\n"
  else ""
}

    rm --recursive --force --verbose extracted_tarballs
  >>>

  output {
    Array[File] qc_reports = glob("~{output_dir}/~{output_prefix}*.json")
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}
