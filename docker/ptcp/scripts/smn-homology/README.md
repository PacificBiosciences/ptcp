# SMN homology

A python script to evaluate the possibility of homology scaling in paraphase for
SMN1/2.  This happens when all copies of SMN1 are identical and all copies of
SMN2 are identical, in which case you have no discernible information to know if
the sample has twice as many copies of SMN1/2 as reported. This is important
when SMN1:SMN2 is either 2:0, 0:2 or 1:1, and in truth the number of copies is
twice as much as what was reported by paraphase.

The approach to distinguish both cases is to train a quadratic model on various
datasets with the finalized PureTarget carrier chemistry. The features are the
coverage per copy of all 11 loci (FXN, FMR1, AFF2, ARX, cyp21, gba, rpgr, hba,
hbb, f8inv1, f8inv22), plus the number of hifi reads (total) and the number of
on-target hifi reads. The output is the estimated coverage SMN1/2 is supposed to
have per copy.

The current regression looks like this.

![scatterplot of observed and predicted coverages per SMN1/2copy](puretarget-v8.json.scatter.png)


## Train a model

Included in this repo is the `ptcp-qc-outs` directory containing the JSON
files that were used to train the model. The current r-squared is `0.9112`. If
you want to train your own model, just pass a directory and an output name like
this:


```
smn_homology train \
  --train_dir ptcp-qc-outs \
  --output_file my-model.json
```


## Evaluate a file

By passing a JSON file that is an output of
[ptcp-qc](https://github.com/pacificBiosciences/ptcp-qc) you can evaluate the
probability of haplotypes having to be scaled as shown below.

```
smn_homology test \
  --model my-model.json \
  --dataset ptcp-qc-outs/sample_SMN11_SMN21.json
```

The output looks like this:
```
{
  "predicted_coverage_per_copy": 50.8511,
  "is_possible_homology": false,
  "smn1_smn1hap1": {
    "paraphase_copy_number": 1,
    "observed_coverage_per_copy": 33,
    "obs_exp_ratio": 0.649,
    "prob_copy_number_is_1": 1.0,
    "prob_copy_number_is_2": 0.0,
    "prob_copy_number_is_3": 0.0,
    "prob_copy_number_is_4": 0.0
  },
  "smn1_smn2hap1": {
    "paraphase_copy_number": 2,
    "observed_coverage_per_copy": 26,
    "obs_exp_ratio": 0.5211,
    "prob_copy_number_is_1": 0.9994,
    "prob_copy_number_is_2": 0.0006,
    "prob_copy_number_is_3": 0.0,
    "prob_copy_number_is_4": 0.0
  },
  "adjusted_smn1_cn": 1,
  "adjusted_smn2_cn": 1
}
```

In this example, paraphase compared the SMN2 coverage of 48 with the SMN1
coverage of 30 and concluded that there are two copies. However, the expected
coverage per copy is 45. The first copy is 30, which is lower than usual.
Paraphase scaling by the first copy number resultedin two copies of SMN2,
whereas the observed coverage is more in line with one copy. Note that the field
`is_possible_homology` is set to false because the genotype is not 1:1 and
contains at least one copy of each locus.

Let's look at another 0:4 example that was genotyped as 0:2
```
smn_homology \
  --model my-model.json \
  --dataset ptcp-qc-outs/sample_SMN10_SMN24.json
```

```
{
  "predicted_coverage_per_copy": 43.0051,
  "is_possible_homology": true,
  "smn1_smn2hap1": {
    "paraphase_copy_number": 2,
    "observed_coverage_per_copy": 99,
    "obs_exp_ratio": 2.3021,
    "prob_copy_number_is_1": 0.0,
    "prob_copy_number_is_2": 0.0,
    "prob_copy_number_is_3": 0.0009,
    "prob_copy_number_is_4": 0.9991
  },
  "adjusted_smn1_cn": 0,
  "adjusted_smn2_cn": 4
}
``` 
This example is a possible homology since there's only one haplotype of
SMN2. Here the number of reads per copy is 99, which is twice what the
regression estimated, so the copy number is actually twice the expected and the
adjusted CN is reported as 4.


