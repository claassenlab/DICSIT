# DICSIT

This repository contains code related to the **d**etection of **i**mportant **c**ell subsets in **si**ngle-cell
**t**ranscriptomics (DICSIT) workflow.
R markdown files for data preprocessing are contained in the [r_markdown](r_markdown) folder, together with output markdown
documents.

The original [CellCnn](https://github.com/eiriniar/CellCnn) code on which DICSIT is based is included in the [python](python) folder.
In addition, there are two python scripts for running a DICSIT/CellCnn analysis based on preprocessed data, `main.py` and `generate_param_grid.py`.
The former can be used to run a single analysis with the following command, specifying a name of the analysis.
For more detailed information about the arguments and options, call `python3 main.py -h`.

```
python3 main.py --ncell NCELL
                --nsubset NSUBSET
                --nfeatures NFEATURES
                --train_data TRAIN_DATA
                --test_data TEST_DATA
                --output_path OUTPUT_PATH
                --response_data RESPONSE_DATA
                --response RESPONSE
                --sample_col SAMPLE_COL
                --name NAME
```

The latter can be used to run analyses using different numbers of features (by default, a range from 10 to 100 in steps of 10),
for different splits of the data.
By default, it assumes that the training and test data CSV files for the different splits are named `train_data_split_{i}.csv`,
where `i` is between 1 and `n_splits` (by default, 3).
This script repeatedly runs `main.py`. By default, the parameters `ncell` and `nsubset` are set to 500 and 1000, respectively.
Here, several arguments can be provided to test different settings.
For more detailed information about the arguments and options, call `python3 generate_param_grid.py -h`.

```
python3 generate_param_grid.py --name NAME
                               --data_path DATA_PATH
                               --response_data RESPONSE_DATA
                               --output_path OUTPUT_PATH
```

Jupyter notebooks for downstream analysis of selected cell subsets are included in [notebooks](notebooks).
Additionally, there is code in [r_markdown](r_markdown) for evaluating differentially expressed genes between selected
and unselected cells.