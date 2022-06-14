# scripts_for_16S
Scripts and pipelines provided in this repository aid to process 16S read libraries. 

#first step. Caution: samples name can't have "_", "-" and "." .

nohup sh 16S_pipeline2.sh 2>&1 > log &

#second step

Rscript melt.r

#third step

Rscript 2021.8.30.16S.R
