#!/bin/bash

date

./00_make_vep_table.sh 2>&1 | tee logs/out_vep.txt

python ./01_make_score_table.py 2>&1 | tee logs/out_score.txt

python ./02_get_prioritized_varIDs.py 2>&1 | tee logs/out_prior_ids.txt

./03_filter_prioritized_vcf.sh 2>&1 | tee logs/out_prior.txt

python ./04_make_combined_table.py 2>&1 | tee logs/out_combine.txt

./05_cmd_filter_table.sh 2>&1 | tee logs/out_filter.txt

date
