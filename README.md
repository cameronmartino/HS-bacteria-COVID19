# Bacterial modification of the host glycosaminoglycan heparan sulfate modulates SARS-CoV-2 infectivity

This repository includes all the source code, tests and notebooks to generate the figures used in the (Martino, Kellman, Sandoval, & Clausen et al. 2020). Each analysis step is outlined here.

## Notes to the reader

We do not provide all of the data in this repository only tables. The Shen et al. data is available and Zhou et al. data is available through the National Genomics Data Center accession PRJCA002202 and NCBI’s Sequence Read Archive accession PRJNA605983 respectively. All American Gut Project sequence data and deidentified participant responses can be found in EBI under project PRJEB11419 and Qiita [https://qiita.ucsd.edu/](https://qiita.ucsd.edu/) study ID 10317.

In particular the FINRISK data that support the findings of this study are available from the THL Biobank based on a written application and following relevant Finnish legislation. Details of the application process are described in the web-site of the Biobank: [https://thl.fi/en/web/thl-biobank/for-researchers](https://thl.fi/en/web/thl-biobank/for-researchers).

## Repository Structure

* code

This directory contains the code (as notebooks) for each step in a numerically organized format. More descriptions can be found within the README file for that directory.

* data

This directory contains the input data, along with some output from computationally expensive tools.

* results

The main and extended figures and tables from the paper which can be reproduced by running the notebooks within the code directory.

# Notebooks

## Setup

All jupyter notebooks (.ipynb) were run in a [qiime2-2019.10](https://raw.githubusercontent.com/qiime2/environment-files/master/2019.10/release/qiime2-2019.10-py36-osx-conda.yml) conda environment. 
