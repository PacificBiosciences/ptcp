# PacBio PureTarget Carrier Pipeline (PTCP) on HPC

The PacBio PureTarget Carrier Pipeline is a WDL-based workflow for genotyping tandem repeat regions and homologous genes with segmental duplications from PacBio PureTarget HiFi data. It uses a containerized toolchain for reproducibility on workstations and HPC clusters. By cloning this repository and setting up a virtual environment you can install PTCP and run it on your HPC.

## Table of Contents

1. [Install requirements in a virtual environment](#1-install-requirements-in-a-virtual-environment)
2. [Set up the image of PTCP dependencies](#2-set-up-the-image-of-ptcp-dependencies)
    - 2.1. [Docker image of PTCP dependencies](#21-docker-image-of-ptcp-dependencies)
    - 2.2. [Singularity image file of PTCP dependencies](#22-singularity-image-file-of-ptcp-dependencies)
3. [Gather input information](#3-gather-input-information)
    - 3.1. [Create the sample-sheet.csv](#31-create-the-sample-sheetcsv)
    - 3.2. [Create the input JSON](#32-create-the-input-json)
    - 3.3. [MiniWDL configuration](#33-miniwdl-configuration)
4. [Running PTCP](#4-running-ptcp)

___

## 1. Install requirements in a virtual environment

Before you begin, make sure you have the following prerequisites available on your HPC or login node:

- Conda/Mamba (recommended on many HPCs) or `venv` for managing a Python environment
- Python 3.12 (or compatible with `miniwdl`)
- Either Docker or Singularity/Apptainer (choose one based on your environment)
- Optional: SLURM scheduler if you plan to use the `miniwdl-slurm` backend

Below are example commands used to install the requirements needed to run PTCP. The examples show `conda`/`mamba` plus either `pip` (most universal) or `uv` (faster). Pick one approach based on what your HPC supports.

> Tip: If you already have an environment named `ptcp` and want a clean rebuild, remove it first with `conda env remove -n ptcp` (or `mamba env remove -n ptcp`).

### Option A: Conda/Mamba + pip

```bash
git clone https://github.com/PacificBiosciences/ptcp.git
cd ptcp/

conda create -n ptcp python=3.12.11
conda activate ptcp
python -m pip install -r requirements.txt

which miniwdl
miniwdl --version
```

### Option B: Conda/Mamba + uv (pip-compatible, faster)

This option installs packages into the active conda environment, but uses `uv` as the installer/resolver.

```bash
git clone https://github.com/PacificBiosciences/ptcp.git
cd ptcp/

mamba create -n ptcp python=3.12.11
conda activate ptcp

python -m pip install uv
uv pip install --system -r requirements.txt

which miniwdl
miniwdl --version
```

### Option C: uv + venv (no conda)

If your HPC does not support conda/mamba, you can use a standard virtualenv instead.

```bash
git clone https://github.com/PacificBiosciences/ptcp.git
cd ptcp/

uv venv --python 3.12.11
source .venv/bin/activate
uv pip install -r requirements.txt

which miniwdl
miniwdl --version
```

## 2. Set up the image of PTCP dependencies

PTCP is a WDL workflow that relies on many software packages. These packages are containerized in an image for PTCP to use. Core tools inside the container: [TRGT](https://github.com/PacificBiosciences/trgt), [Paraphase](https://github.com/PacificBiosciences/paraphase), [Sawfish](https://github.com/PacificBiosciences/sawfish), [ptcp-qc](https://github.com/PacificBiosciences/ptcp-qc), [SMRT Link](https://www.pacb.com/smrt-link/)


It can be invoked as a Docker container or a Singularity image file. Instructions for setting up both are shown below.

### 2.1. Docker image of PTCP dependencies

A pre-built Docker image of these dependencies can be found on quay.io at: `quay.io/pacbio/ptcp:3.2`. Below is an example `docker` command used to pull the image of PTCP dependencies for version 3.2:

```bash
docker pull quay.io/pacbio/ptcp:3.2
```

### 2.2. Singularity image file of PTCP dependencies

You can build a `.sif` of PTCP dependencies for the version of PTCP you are installing by using `apptainer` (or `singularity`) to build the `.sif` from the Docker image hosted on quay.io. An example command is shown below:

```bash
apptainer pull ptcp_3.2.sif docker://quay.io/pacbio/ptcp:3.2
```

Once you have built the `.sif` file you should move it to the `ptcp` GitHub repo folder, like so:

```bash
cd /Path/to/ptcp
mkdir miniwdl_singularity_cache
cd miniwdl_singularity_cache
mv /Path/to/Downloads/ptcp_3.2.sif .
```

## 3. Gather input information

PTCP requires several inputs, all summarized in an input JSON file. For detailed information on each of the input files required, please refer to the [Input Files page](./Input_files.md). This section will briefly cover how to create the `sample-sheet.csv` and the input JSON file used to run the pipeline.

### 3.1 Create the sample-sheet.csv

The `sample-sheet.csv` contains at least two comma-separated columns: `bam_id` and `sex`. A legacy third column (`bam_name`) is also supported. For more information on how to format this file, please refer to the [Input Files page](./Input_files.md). You will need to create this file with one row per sample before you can create the input JSON in the next section.

### 3.2 Create the input JSON

The input JSON contains all of the input information required by PTCP. Because the input JSON contains sample-specific information, it needs to be generated for each run of PureTarget Carrier Panel data you want to analyze. This repository contains a script called [create_input_json.py](../docker/ptcp/scripts/create_input_json.py) to make it easier to generate the input JSON with the sample and other input information. This script requires a template JSON be passed to it with certain fields filled in.

An example template of the input JSON is provided in this repository called [inputs_json_template.json](../tests/inputs/templates/inputs_json_template.json). It is recommended that you copy this template and update all of the fields with preset example paths to point to their file locations on your server. For instance, in the example below, you should update all fields **_except for_** `ptcp.sample_sheet`, `ptcp.hifi_reads`, and `ptcp.fail_reads`:

```json
{
    "ptcp.sample_sheet": "",
    "ptcp.ref_fasta": "/path/to/reference/hg38.fa",
    "ptcp.ref_index": "/path/to/reference/hg38.fa.fai",
    "ptcp.trgt_bed": "/path/to/ptcp/meta/trgt/PureTarget_repeat_expansion_panel_2.0.repeat_definition.GRCh38.bed",
    "ptcp.paraphase_config_yaml": "/path/to/ptcp/meta/paraphase/paraphase_config.GRCh38.yaml",
    "ptcp.paraphase_annotation_vcf": "/path/to/ptcp/meta/variant_list/variant_list.GRCh38.vcf",
    "ptcp.genome_version": "38",
    "ptcp.ptcp_qc_bed": "/path/to/ptcp/meta/ptcp-qc/ptcp-qc.GRCh38.bed",
    "ptcp.hifi_reads": [],
    "ptcp.fail_reads": []
}
```

Optional inputs you may add to (or remove from) the template:

- `ptcp.paraphase_annotation_vcf` (included in the repository template; omit this key to disable `havanno_json`)
- `ptcp.pt_linear_regression` (enables SMN homology correction)
- `ptcp.docker_smrttools` (defaults to `quay.io/pacbio/ptcp:X.Y`)

Once you have your template updated, you can use it to create the input JSON with the `create_input_json.py` script like so:

```bash
python docker/ptcp/scripts/create_input_json.py \
  --data /Path/to/demux/PacBio/Run \
  --sample_sheet /Path/to/sample-sheet.csv \
  --template /Path/to/template.json \
  > ptcp_inputs.json
```

### 3.3 MiniWDL configuration

PTCP also requires a **MiniWDL configuration file** as input. This file defines important runtime settings such as image cache, task concurrency, container backends, and more. You must update the configuration file to match your local system or execution environment.

Make sure to review and customize the configuration according to your hardware, software paths, and environment preferences before running the workflow.

**Example configuration file**
You can find an example MiniWDL configuration file here at
[tests/miniwdl.cfg](../tests/miniwdl.cfg).

For more information about configuring MiniWDL, refer to the [official MiniWDL documentation](https://miniwdl.readthedocs.io/en/latest/GettingStarted/#configuration).

## 4. Running PTCP

Once you have your environment configured, the image of dependencies installed, and your inputs gathered you can run PTCP locally. Below is an example command used to run the pipeline:

```bash
conda activate ptcp
cd /Path/to/ptcp

miniwdl run \
    --verbose \
    --dir /Path/to/output_dir \
    --cfg /Path/to/miniwdl.cfg \
    --input /Path/to/ptcp_inputs.json \
    main.wdl
```

Outputs and logs will appear under /Path/to/output_dir in a timestamped run directory.
