# scripts_for_16S
#first step
nohup sh 16S_pipeline2.sh 2>&1 > log &

#second step
Rscript melt.r

#third step

Rscript 2021.3.16S.Alpha.r
