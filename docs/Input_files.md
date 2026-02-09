# Inputs and pipeline configuration

## Table of contents

1. [Input types and descriptions](#1-input-types-and-descriptions)
    - 1.1. [PacBio read data](#11-pacbio-sequencing-data)
    - 1.2. [Sample sheet](#12-sample-sheet)
    - 1.3. [Reference genome](#13-reference-genome)
    - 1.4. [Regions and annotation files](#14-regions-and-annotation-files)
    - 1.5. [Config files](#15-configuration-files)
    - 1.6. [An image of dependencies](#16-an-image-of-dependencies)
2. [Input JSON](#2-input-json-configuration)
    - 2.1. [Example input JSON](#21-example-input-json)
    - 2.2. [Helper scripts](#22-helper-scripts)
3. [Input Definitions Table](#3-input-definitions-table)

## 1. Input types and descriptions

PTCP uses six primary input types, provided via an input JSON file (some inputs have defaults):

1. PacBio sequencing data (`.bam`)
2. Sample sheet (`.csv`)
3. Reference genome (`.fa` and `.fai`)
4. Regions and annotations (`.bed` and `.vcf`)
5. Configuration file (`.yaml`)
6. PTCP dependencies image (local `.sif` or remote/local Docker image)

In practice, only the sequencing data and sample sheet typically change between runs. The other inputs (3–6) are typically set up once when installing the pipeline and remain the same for future runs.

### 1.1. PacBio sequencing data

- HiFi reads: BAM files containing high-quality reads. Common filename patterns include `SAMPLE.hifi_reads.bam`, `SAMPLE.reads.bam`, and `SAMPLE.bam` (depending on platform).
- Fail reads: Optional BAMs of reads that failed HiFi quality filters; used by TRGT for improved tandem repeat genotyping. Common filename pattern: `SAMPLE.fail_reads.bam`. Automatically generated on Revio; available on Sequel II/IIe with SMRT Link v13.0 or later. It is highly recommended to include fail reads for PTCP.

PTCP associates BAMs to samples by *normalizing the BAM filename*:
it strips `.hifi_reads`, `.fail_reads`, `.reads`, and the `.bam` extension from the BAM basename and uses the remainder as the sample ID.
This normalized name must match `bam_id` in the sample sheet (see below).

#### Why are fail reads important for TRGT?

The PureTarget protocol produces insert sizes of about 5 kb, but large expansions of loci like *FXN*, *C9orf72*, *DMPK* and *CNBP* produce much larger molecules that may not produce reads reaching HiFi quality thresholds at typical movie times. Using all reads produced by the sequencer (including fail reads) may significantly increase the coverage of expanded alleles and prevent allelic dropouts. For detailed information on TRGT analysis of PureTarget data, see the [TRGT PureTarget documentation](https://github.com/PacificBiosciences/trgt/blob/main/docs/puretarget.md).

### 1.2. Sample sheet

A header-containing CSV file that describes each sample and its metadata.

#### Supported columns
- **`bam_id`** (required): Unique sample ID (e.g., `some_sample.bc42`). Must match the BAM filename *after PTCP normalization* (PTCP strips `.hifi_reads`, `.fail_reads`, `.reads`, and `.bam`). For example, `bam_id=SAMPLE.bc42` matches `SAMPLE.bc42.hifi_reads.bam` and `SAMPLE.bc42.fail_reads.bam`.
- **`sex`** (required): `M`, `F`, or empty string (defaults to `F`). **Case-sensitive:** only exactly `M` is treated as male in downstream karyotype logic; any other value (including `m`, `male`, etc.) will be treated as `F`. Consider validating this column (values in `{M, F, ""}`) before running to avoid silent mislabeling.
- **`bam_name`** (optional, legacy): The exact BAM filename (e.g., `some_sample_1.hifi_reads.bc42.bam`). Historically required; still accepted.

#### File naming convention

If `bam_id` is `SAMPLE`, common/recommended filenames include:
- `SAMPLE.bam`
- `SAMPLE.reads.bam`
- `SAMPLE.hifi_reads.bam`
- `SAMPLE.fail_reads.bam`

#### Example sample sheets

Minimal (recommended):

```bash
bam_id,sex
m84036_250311_030304_s2_HG00281_SMN11_SMN23.bc2081,F
m84036_250311_030304_s2_HG00324_SMN11_SMN22.bc2082,F
m84036_250311_030304_s2_HG00346_SMN11_SMN23.bc2088,F
m84036_250311_030304_s2_HG00544_CYP21A2.bc2073,F
m84036_250311_030304_s2_HG01071_CYP21A2.bc2074,F
```

Legacy three-column (still supported):

```bash
bam_name,bam_id,sex
m84036_250311_030304_s2_HG00281_SMN11_SMN23.bc2081.bam,m84036_250311_030304_s2_HG00281_SMN11_SMN23.bc2081,F
m84036_250311_030304_s2_HG00324_SMN11_SMN22.bc2082.bam,m84036_250311_030304_s2_HG00324_SMN11_SMN22.bc2082,F
m84036_250311_030304_s2_HG00346_SMN11_SMN23.bc2088.bam,m84036_250311_030304_s2_HG00346_SMN11_SMN23.bc2088,F
m84036_250311_030304_s2_HG00544_CYP21A2.bc2073.bam,m84036_250311_030304_s2_HG00544_CYP21A2.bc2073,F
m84036_250311_030304_s2_HG01071_CYP21A2.bc2074.bam,m84036_250311_030304_s2_HG01071_CYP21A2.bc2074,F
```

**Important**: Only BAM files listed in the sample sheet will be processed. Submitting any BAM not described in the sheet will cause the workflow to fail.

A [Python script](../docker/ptcp/scripts/create_sample_sheet.py) is included in the repository to help generate simple sample sheets from a given input folder containing BAM files. 
Note: this script defaults every sample's `sex` to `F`, so you must edit the output sample sheet for mixed-sex cohorts.

### 1.3. Reference genome

PTCP supports only the reference genomes listed in the [PacBio reference genomes repository](https://github.com/PacificBiosciences/reference_genomes):

- **GRCh38/hg38**: Use configuration files from the `meta` directory with the `GRCh38` substring and set the genome flag to `"38"`
- **GRCh37/hg19**: Use configuration files from the `meta` directory with the `GRCh37` substring and set the genome flag to `"37"`

#### Required files
- **FASTA file**: The reference genome sequence (`.fa`)
- **Index file**: The FASTA index (`.fai`)

### 1.4. Regions and annotation files

The latest TRGT repeat catalog, Paraphase gene configurations, `ptcp-qc` target gene configurations, and variant annotation VCF for each supported reference genome can be found in the `meta` directory. They are described in the [Input definitions table](#3-input-definitions-table) below.

⚠️ **Warning**: Mixing configuration files or using the wrong reference genome flags can lead to incorrect results or runtime errors. Always double-check that your inputs match your selected reference genome.

### 1.5. Configuration files

#### Paraphase configuration

Paraphase requires additional per-gene configurations, captured in a YAML in the `meta` folder. This is described in the [Input definitions table](#3-input-definitions-table) below.

#### Optional configuration files
- **Linear regression weights**: JSON file for SMN1/2 copy number adjustment in cases of full haplotype homology
- **Variant annotation VCF**: For annotating Paraphase small variant calls

### 1.6. An image of dependencies

PTCP is a WDL workflow that relies on many software packages. These packages are containerized in an image for PTCP to use during runtime. By default, the workflow uses `quay.io/pacbio/ptcp:X.Y` (see `docker_smrttools` in the [Input definitions table](#3-input-definitions-table)); override this if you need a different tag or a local image.

The image contains:

- [TRGT](https://github.com/PacificBiosciences/trgt)
- [Paraphase](https://github.com/PacificBiosciences/paraphase)
- [ptcp-qc](https://github.com/PacificBiosciences/ptcp-qc)
- [Sawfish](https://github.com/PacificBiosciences/sawfish)
- [SMRT Link](https://www.pacb.com/smrt-link/)
- [Rust](https://www.rust-lang.org/)

It can be invoked as a Docker container or a Singularity image file. Instructions for setting up both are provided in the [PTCP on HPC page](./PTCP_on_HPC.md#2-set-up-the-image-of-ptcp-dependencies).

## 2. Input JSON configuration

PTCP takes a single JSON file that defines all required inputs. It supports processing multiple samples in one run by linking each sample's HiFi (and optional fail) reads to shared reference and configuration files. For a template input JSON, see: [tests/inputs/templates/inputs_json_template.json](../tests/inputs/templates/inputs_json_template.json) 

### 2.1. Example input JSON

```json
{
  "ptcp.sample_sheet": "inputs/sample_sheet.csv",
  "ptcp.ref_fasta": "data/reference/GRCh38.fa",
  "ptcp.ref_index": "data/reference/GRCh38.fa.fai",
  "ptcp.trgt_bed": "meta/trgt/PureTarget_repeat_expansion_panel_2.0.repeat_definition.GRCh38.bed",
  "ptcp.paraphase_config_yaml": "meta/paraphase/paraphase_config.GRCh38.yaml",
  "ptcp.paraphase_annotation_vcf": "meta/variant_list/variant_list.GRCh38.vcf",
  "ptcp.genome_version": "38",
  "ptcp.ptcp_qc_bed": "meta/ptcp-qc/ptcp-qc.GRCh38.bed",
  "ptcp.pt_linear_regression": "regression_weights.json",
  "ptcp.docker_smrttools": "quay.io/pacbio/ptcp:3.2",
  "ptcp.hifi_reads": [
    "/data/reads/SAMPLE1.hifi_reads.bc42.bam",
    "/data/reads/SAMPLE2.hifi_reads.bc43.bam"
  ],
  "ptcp.fail_reads": [
    "/data/reads/SAMPLE1.fail_reads.bc42.bam",
    "/data/reads/SAMPLE2.fail_reads.bc43.bam"
  ]
}
```

### 2.2. Helper scripts

PTCP includes Python scripts to help generate input files:

**Sample sheet generation**: The script [`docker/ptcp/scripts/create_sample_sheet.py`](../docker/ptcp/scripts/create_sample_sheet.py) generates sample sheets from a folder containing BAM files.

**Input JSON generation**: The script [`docker/ptcp/scripts/create_input_json.py`](../docker/ptcp/scripts/create_input_json.py) generates input JSON files from BAM files, sample sheet, and template JSON. It uses a template file [`tests/inputs/templates/inputs_json_template.json`](../tests/inputs/templates/inputs_json_template.json).

## 3. Input definitions table

| Name                       | Type          | Description                                                                                               | Notes |
| -------------------------- | ------------- | --------------------------------------------------------------------------------------------------------- | ----- |
| `sample_sheet`             | `File`        | Sample sheet csv file                                                                                     | 1     |
| `ref_fasta`                | `File`        | Reference genome FASTA file                                                                               |       |
| `ref_index`                | `File`        | Reference genome index                                                                                    |       |
| `trgt_bed`                 | `File`        | BED file of tandem repeat targets in TRGT format                                                          | 2     |
| `paraphase_config_yaml`    | `File`        | YAML file defining gene families and parameters for Paraphase                                             | 3     |
| `paraphase_annotation_vcf` | `File`        | Optional VCF file with known variations for annotating small variant calls on Paraphase phased haplotypes | 4     |
| `genome_version`           | `String`      | Genome version flag for Paraphase and Sawfish                                                             | 5     |
| `ptcp_qc_bed`              | `File`        | BED file of regions used by ptcp-qc for coverage and mapping-quality metrics                              |       |
| `pt_linear_regression`     | `File`        | Optional JSON file defining linear regression weights to adjust SMN1/2 copy number                        | 6     |
| `docker_smrttools`         | `String`      | URI or tarball path for Docker image                                                                      | 7     |
| `hifi_reads`               | `Array[File]` | list of `hifi_reads.bam` paths (Revio) or `reads.bam` paths (Sequel IIe)                                  |       |
| `fail_reads`               | `Array[File]` | list of `fail_reads.bam` paths (Revio)                                                                    |       |
| `log_level`                | `String`      | default: `'INFO'`                                                                                         | 8     |

**Notes:**
1. **Sample sheet**: A header-containing CSV file with required fields `bam_id` (e.g., `some_sample.bc42`) and `sex` (`M`, `F`, or empty string, assumed to be `F`). An optional legacy `bam_name` column is still supported for backward compatibility.
2. **TRGT BED file**: A 4-column unheadered bed file with fields `chrom`, `start`, `end`, and `name` following the [TRGT repeat definitions format](https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md).
3. **Paraphase configuration**: A YAML file following the [Paraphase configuration format](https://github.com/PacificBiosciences/paraphase/blob/main/paraphase/data/38/config.yaml).
4. **Variant annotation VCF**: A VCF file that describes pathogenic or pseudogene-specific status per variant. See the [havanno README](../docker/ptcp/scripts/havanno/README.md) for more details. This is optional, and if not included, the havanno JSON will not be generated.
5. **Genome version**: Acceptable options are `'38'` and `'37'`.
6. **Linear regression weights**: Used to adjust SMN1/2 copy numbers in cases of full haplotype homology, these are written to the ptcp-qc outputs. See the [SMN homology README](../docker/ptcp/scripts/smn-homology/README.md) for details on training the regression model. This is optional, and if not included, the correction step will be skipped.
7. **Dependencies image (`docker_smrttools`)**: Defaults to `quay.io/pacbio/ptcp:X.Y`. Override to use a different tag, a local tarball, or a local image path (e.g., a `.sif` when using Apptainer/Singularity).
8. **Log level**: Acceptable options are: `'DEBUG'`, `'INFO'`, and `'WARN'`.
