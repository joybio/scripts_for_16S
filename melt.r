library(reshape2)
data <- read.table("otu_table_tax.txt",header=T,row.names=1,sep="\t",quote="")
data_melt <- melt(data,id.vars="taxonomy",variable.name="sample",value.name="abundunce")
write.csv(data_melt,"data_melt.csv")




