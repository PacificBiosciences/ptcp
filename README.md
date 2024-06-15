# PureTarget + Dark Genes Pipeline

- [PureTarget + Dark Genes Pipeline](#puretarget--dark-genes-pipeline)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [To run](#to-run)
  - [Makefile](#makefile)

## Inputs

[tests/inputs/templates/main.inputs.json](tests/inputs/templates/main.inputs.json)

| Type | Name | Description | Notes |
|---|---|---|---|
| `File` | `sample_sheet` | sample sheet tsv | [^1] |
| `Array[File]` | `hifi_reads` | list of `hifi_reads.bam` paths (Revio) or `reads.bam` paths (Sequel IIe) |  |
| `Array[File]?` | `fail_reads` | list of `fail_reads.bam` paths (Revio) |  |
| `File` | `ref_fasta` | reference fasta |  |
| `File` | `ref_index` | reference fasta fai |  |
| `File` | `pbaa_guide_bed` | bed file of pbaa guide regions | [^2]  |
| `File` | `pbaa_mask_bed` | bed file of reference regions to mask | [^3] |
| `File` | `trgt_bed` | bed file of repeat expansion targets in TRGT format | [^4] |
| `String?` | `log_level` | default: `'INFO'` | `['DEBUG', 'INFO', 'WARN']` |
| `String` | `docker_smrttools` | URI for Docker image or tarball | Default value will not be valid.  You must create the image yourself and override the default. |

[^1]: The sample sheet is a 4 column headered TSV with fields `movie` (e.g., `m60000e_240517_211441`), `barcode` (e.g., `bc2017--bc2017`), `sample_id` (alphanumeric characters, dashes, and underscores only), `sex` (`M`, `F`, or null string, assumed to be `F`).

[^2]: The pbaa guide bed file is a 4 column unheadered bed file with fields `chrom`, `start`, `end`, and `name` that is used to generate [pbaa guide sequences](https://github.com/PacificBiosciences/pbAA?tab=readme-ov-file#customizing-guide-sequences).  The `name` field is the guide name, following the format `sequence_name|group_name`.  Groups are similar sequences (e.g., paralogous sequences) that should be treated as a single target by pbaa.

[^3]: The mask bed file is a 3 column unheadered bed file with fields `chrom`, `start`, and `end`.  These regions will be masked in the reference fasta.

[^4]: The TRGT bed file is a 4 column unheadered bed file with fields `chrom`, `start`, `end`, and `name` following the [TRGT repeat definitions format](https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md).

## Outputs

| Type | Name | Description | Notes |
|---|---|---|---|
| `File` | `repeats` | copy of input `trgt_bed` |   |
| `File` | `pbaa_guides` | copy of input `pbaa_guide_bed` |   |
| `File` | `pbaa_mask` | copy of input `pbaa_mask_bed` |   |
| `Array[Float]` | `chry_frequency` | $\frac{reads_{chrY}}{reads_{mapped}}$ | potentially useful to infer presence of chrY |
| `Array[File]` | `trgt_vcf` | tandem repeat genotypes | [details](https://github.com/PacificBiosciences/trgt/blob/main/docs/vcf_files.md) |
| `Array[File]` | `trgt_spanning_bam` | BAM of input for trvz |   |
| `Array[File]` | `trgt_spanning_index` | index for above |   |
| `Array[File]` | `images_motifs_allele` | zip of per-region images |   |
| `Array[File]` | `images_meth_allele` | zip of per-region images |   |
| `Array[File]` | `images_motifs_waterfall` | zip of per-region images |   |
| `Array[File]` | `images_meth_waterfall` | zip of per-region images |   |
| `Array[File]` | `reads_overlapping_repeats` | BAM of reads overlapping repeat regions |   |
| `Array[File]` | `reads_overlapping_repeats_index` | index for above |   |
| `Array[File?]` | `report_target_enrichment_resources` | images of on/off-target read distribution and coverage per sample/region |   |
| `File` | `genotype` | table of repeat genotypes |   |
| `File` | `report_target_enrichment` |   |   |
| `File` | `report_puretarget` |   |   |
| `Array[File]` | `cluster_info` | list of all reads with homology to pbaa guides | [details](https://github.com/PacificBiosciences/pbaa?tab=readme-ov-file#read-information-output-file) |
| `Array[File]` | `painted_bams` | BAM of passing reads with homology to pbaa guides aligned to masked reference, with `HP` and `YC` tags |  |
| `Array[File]` | `painted_bams_bais` | index for above |   |
| `Array[File]` | `consensus_bams` | BAM of cluster consensus sequences aligned to masked reference |   |
| `Array[File]` | `consensus_bam_bais` | index for above |   |
| `Array[File]` | `minipileup_vcfs` | VCF of variants for consensus sequences against masked reference |   |
| `Array[File]` | `f8_pbaa_vcf` | VCF of *F8* inversion calls |   |
| `Array[File]` | `f8_pbaa_json` | JSON of evidence for *F8* inversion calls |   |

## To run

1) Build docker image.  Optionally tarball image or build singularity image.

    ```bash
    cd docker/smrttools/

    bash -x ./build.sh

    sha256sum=$(docker inspect --format='{{index .RepoDigests 0}}' "smrttools:${IMAGE_TAG}" 2>/dev/null | sha256sum | cut -d ' ' -f 1)

    # to save the docker tarball for DNAnexus
    # docker save "smrttools:${IMAGE_TAG}" | gzip -c > "smrttools_${IMAGE_TAG}.tar.gz"

    # to build the sif image for singularity, internal testing
    # singularity build "docker___quay.io_pacbio_smrttools@sha256_${sha256sum}.sif" docker-daemon://smrttools:${IMAGE_TAG}
    ```

2) Optional. Set up miniwdl to test pipeline.

    ```bash
    # create and source virtual environment
    python3 -m venv venv
    source venv/bin/activate

    # install miniwdl
    python3 -m pip install miniwdl
    # if submitting to slurm
    python3 -m pip install miniwdl-slurm

    # edit tests/miniwdl.cfg as desired
    ```

3) Set up inputs.json.

    ```bash
    cp tests/inputs/templates/main.inputs.json tests/inputs/main.inputs.json
    # edit with paths to test data
    ```

4) Run pipeline.

    ```bash
    miniwdl run --verbose --dir ./miniwdl_test_output --cfg ./tests/miniwdl.cfg --input ./tests/inputs/main.inputs.json main.wdl
    ```

## Makefile

The Makefile contains a few useful targets for linting the workflow, creating templates for tasks/workflows, and running tasks/workflows with miniwdl.

```bash
# lint common.wdl
make check wdl=common

# lint all wdl
make check-all

# list all tasks in common.wdl
make list-tasks wdl=common

# create template for task
make task-template wdl=common task=collect_inputs

# test a task
make task-run wdl=common task=collect_inputs

# create template for workflow
make template wdl=main

# run workflow
make run wdl=main
```
