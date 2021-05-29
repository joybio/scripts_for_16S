rm(list=ls())
setwd("D://项目/2021.5.29.hfq.2.16S/")
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)
dir.create("2.OTU table")
setwd("2.OTU table")

DATA=read.table("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
group=read.table("../2.OTU table/group.csv",sep=",",header = T)
colnames(group) <- c("Sample","Group")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#group$Group <- factor(as.character(group$Group),levels=c("A0","A5","A8","A14","B0","B5","B8","B14","C0","C5","C8","C14"))
#group$group <- substr(group$Group,1,1)
row.names(group) <- group$Sample
table(group$Sample)
#?read.table
?group_by
?dir.create
dir.create("../4.Taxonomy/")
dir.create("../4.Taxonomy/4.1 abundunce")
setwd("../4.Taxonomy/4.1 abundunce")

data <- DATA[,row.names(group)]
data_k <- DATA
data_k <- subset(data_k,data_k$k != " " & data_k$k != "")

name <- colnames(data)[1]
data_k_per <- data.frame(tapply(data_k[,c(name)], data_k$k, sum))
colnames(data_k_per) <- name

for(i in as.character(colnames(data)[2:length(colnames(data))])){
  name <- as.character(eval(parse(text="i")))
  output <- data.frame(tapply(data_k[,c(name)], data_k$k, sum))
  colnames(output) <- name
  data_k_per <- cbind(data_k_per,output)
}

# 新建data_norm数据框,赋初值
data_k_per_sum=apply(data_k_per, 2, sum)
# 计算每个样本的菌种总丰度
data_k_per_norm <- data_k_per
for(i in 1:ncol(data_k_per))
{
  for(j in 1:nrow(data_k_per)){
    data_k_per_norm[j,i]=data_k_per[j,i]/data_k_per_sum[i]
    # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
  }
}


apply(data_k_per_norm, 2, sum)
#data_ex <- subset(data_k_per_norm,select=c(A0rep3,A14rep3,A5rep3,A8rep3,B0rep3,B14rep3,B5rep3,B8rep3,C0rep4,C14rep4,C5rep4,C8rep4))

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_k_per_norm[is.na(data_k_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_k_per_norm<-data_k_per_norm[order(rowSums(data_k_per_norm),decreasing = T),]
?write.table
write.table(data_k_per_norm,"data_k_per_norm.xls",quote =F,sep="\t")
#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_ex_list<-rownames(data_k_per_norm)[1:N]
new_x<-rbind(data_k_per_norm[row.names(data_k_per_norm) %in% data_ex_list,],
             others=rowSums(data_k_per_norm[!rownames(data_k_per_norm) %in% data_ex_list,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)))#,id.vars = c('Taxonomy'))
#datm$variable <- factor(datm$variable,levels=c("A0rep1","A0rep2",  "A0rep3",  "A0rep4",  "A0rep5", "A5rep1",  "A5rep2" , "A5rep3" , "A5rep4" , "A5rep5", 
#                                               "A8rep1" , "A8rep2",  "A8rep3",  "A8rep4" , "A8rep5" ,'A14rep1', 'A14rep2', 'A14rep3', 'A14rep4' ,'A14rep5',
#                                               "B0rep1","B0rep2",  "B0rep3",  "B0rep4",  "B0rep5", "B5rep1",  "B5rep2" , "B5rep3" , "B5rep4" , "B5rep5",
#                                               "B8rep1" , "B8rep2",  "B8rep3",  "B8rep4" , "B8rep5" ,'B14rep1', 'B14rep2', 'B14rep3', 'B14rep4' ,'B14rep5',
#                                               "C0rep1","C0rep2",  "C0rep3",  "C0rep4",  "C0rep5", "C5rep1",  "C5rep2" , "C5rep3" , "C5rep4" , "C5rep5",
#                                               "C8rep1" , "C8rep2",  "C8rep3",  "C8rep4" , "C8rep5" ,'C14rep1', 'C14rep2', 'C14rep3', 'C14rep4' ,'C14rep5'))
#作图
#datm$variable <- factor(datm$variable,levels=c("A0rep3","A5rep3","A8rep3",'A14rep3',"B0rep3","B5rep3","B8rep3","B14rep3","C0rep4","C5rep4","C8rep4","C14rep4"))
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Domain.pdf",width = 10, height = 8)


data_p <- DATA
data_p <- subset(data_p,data_p$p != " " & data_p$p != "")

name <- colnames(data)[1]
data_p_per <- data.frame(tapply(data_p[,c(name)], data_p$p, sum))
colnames(data_p_per) <- name

for(i in as.character(colnames(data)[2:length(colnames(data))])){
  name <- as.character(eval(parse(text="i")))
  output <- data.frame(tapply(data_p[,c(name)], data_p$p, sum))
  colnames(output) <- name
  data_p_per <- cbind(data_p_per,output)
}

# 新建data_norm数据框,赋初值
data_p_per_sum=apply(data_p_per, 2, sum)
# 计算每个样本的菌种总丰度
data_p_per_norm <- data_p_per
for(i in 1:ncol(data_p_per))
{
  for(j in 1:nrow(data_p_per)){
    data_p_per_norm[j,i]=data_p_per[j,i]/data_p_per_sum[i]
    # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
  }
}

apply(data_p_per_norm, 2, sum)
#data_ex <- subset(data_p_per_norm,select=c(A0rep3,A14rep3,A5rep3,A8rep3,B0rep3,B14rep3,B5rep3,B8rep3,C0rep4,C14rep4,C5rep4,C8rep4))

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_p_per_norm[is.na(data_p_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_p_per_norm<-data_p_per_norm[order(rowSums(data_p_per_norm),decreasing = T),]
write.table(data_p_per_norm,"data_p_per_norm.xls",quote =F,sep="\t")

#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_ex_list<-rownames(data_p_per_norm)[1:N]
new_x<-rbind(data_p_per_norm[row.names(data_p_per_norm) %in% data_ex_list,],
             others=rowSums(data_p_per_norm[!rownames(data_p_per_norm) %in% data_ex_list,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)))#,id.vars = c('Taxonomy'))
#datm$variable <- factor(datm$variable,levels=c("A0rep1","A0rep2",  "A0rep3",  "A0rep4",  "A0rep5", "A5rep1",  "A5rep2" , "A5rep3" , "A5rep4" , "A5rep5", 
#                                               "A8rep1" , "A8rep2",  "A8rep3",  "A8rep4" , "A8rep5" ,'A14rep1', 'A14rep2', 'A14rep3', 'A14rep4' ,'A14rep5',
#                                               "B0rep1","B0rep2",  "B0rep3",  "B0rep4",  "B0rep5", "B5rep1",  "B5rep2" , "B5rep3" , "B5rep4" , "B5rep5",
#                                               "B8rep1" , "B8rep2",  "B8rep3",  "B8rep4" , "B8rep5" ,'B14rep1', 'B14rep2', 'B14rep3', 'B14rep4' ,'B14rep5',
#                                               "C0rep1","C0rep2",  "C0rep3",  "C0rep4",  "C0rep5", "C5rep1",  "C5rep2" , "C5rep3" , "C5rep4" , "C5rep5",
#                                               "C8rep1" , "C8rep2",  "C8rep3",  "C8rep4" , "C8rep5" ,'C14rep1', 'C14rep2', 'C14rep3', 'C14rep4' ,'C14rep5'))
#作图
#datm$variable <- factor(datm$variable,levels=c("A0rep3","A5rep3","A8rep3",'A14rep3',"B0rep3","B5rep3","B8rep3","B14rep3","C0rep4","C5rep4","C8rep4","C14rep4"))
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Phylum.pdf",width = 10, height = 8)

data_c <- DATA
data_c <- subset(data_c,data_c$c != " " & data_c$c != "")

name <- colnames(data)[1]
data_c_per <- data.frame(tapply(data_c[,c(name)], data_c$c, sum))
colnames(data_c_per) <- name

for(i in as.character(colnames(data)[2:length(colnames(data))])){
  name <- as.character(eval(parse(text="i")))
  output <- data.frame(tapply(data_c[,c(name)], data_c$c, sum))
  colnames(output) <- name
  data_c_per <- cbind(data_c_per,output)
}

# 新建data_norm数据框,赋初值
data_c_per_sum=apply(data_c_per, 2, sum)
# 计算每个样本的菌种总丰度
data_c_per_norm <- data_c_per
for(i in 1:ncol(data_c_per))
{
  for(j in 1:nrow(data_c_per)){
    data_c_per_norm[j,i]=data_c_per[j,i]/data_c_per_sum[i]
    # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
  }
}

apply(data_c_per_norm, 2, sum)
#data_ex <- subset(data_c_per_norm,select=c(A0rep3,A14rep3,A5rep3,A8rep3,B0rep3,B14rep3,B5rep3,B8rep3,C0rep4,C14rep4,C5rep4,C8rep4))

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_c_per_norm[is.na(data_c_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_c_per_norm<-data_c_per_norm[order(rowSums(data_c_per_norm),decreasing = T),]
write.table(data_c_per_norm,"data_c_per_norm.xls",quote =F,sep="\t")

#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_ex_list<-rownames(data_c_per_norm)[1:N]
new_x<-rbind(data_c_per_norm[row.names(data_c_per_norm) %in% data_ex_list,],
             others=rowSums(data_c_per_norm[!rownames(data_c_per_norm) %in% data_ex_list,]))


#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)))#,id.vars = c('Taxonomy'))
#datm$variable <- factor(datm$variable,levels=c("A0rep1","A0rep2",  "A0rep3",  "A0rep4",  "A0rep5", "A5rep1",  "A5rep2" , "A5rep3" , "A5rep4" , "A5rep5", 
#                                               "A8rep1" , "A8rep2",  "A8rep3",  "A8rep4" , "A8rep5" ,'A14rep1', 'A14rep2', 'A14rep3', 'A14rep4' ,'A14rep5',
#                                               "B0rep1","B0rep2",  "B0rep3",  "B0rep4",  "B0rep5", "B5rep1",  "B5rep2" , "B5rep3" , "B5rep4" , "B5rep5",
#                                               "B8rep1" , "B8rep2",  "B8rep3",  "B8rep4" , "B8rep5" ,'B14rep1', 'B14rep2', 'B14rep3', 'B14rep4' ,'B14rep5',
#                                               "C0rep1","C0rep2",  "C0rep3",  "C0rep4",  "C0rep5", "C5rep1",  "C5rep2" , "C5rep3" , "C5rep4" , "C5rep5",
#                                               "C8rep1" , "C8rep2",  "C8rep3",  "C8rep4" , "C8rep5" ,'C14rep1', 'C14rep2', 'C14rep3', 'C14rep4' ,'C14rep5'))
#作图
#datm$variable <- factor(datm$variable,levels=c("A0rep3","A5rep3","A8rep3",'A14rep3',"B0rep3","B5rep3","B8rep3","B14rep3","C0rep4","C5rep4","C8rep4","C14rep4"))
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Class.pdf",width = 10, height = 8)

data_o <- DATA
data_o <- subset(data_o,data_o$o != " " & data_o$o != "")

name <- colnames(data)[1]
data_o_per <- data.frame(tapply(data_o[,c(name)], data_o$o, sum))
colnames(data_o_per) <- name

for(i in as.character(colnames(data)[2:length(colnames(data))])){
  name <- as.character(eval(parse(text="i")))
  output <- data.frame(tapply(data_o[,c(name)], data_o$o, sum))
  colnames(output) <- name
  data_o_per <- cbind(data_o_per,output)
}

# 新建data_norm数据框,赋初值
data_o_per_sum=apply(data_o_per, 2, sum)
# 计算每个样本的菌种总丰度
data_o_per_norm <- data_o_per
for(i in 1:ncol(data_o_per))
{
  for(j in 1:nrow(data_o_per)){
    data_o_per_norm[j,i]=data_o_per[j,i]/data_o_per_sum[i]
    # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
  }
}

apply(data_o_per_norm, 2, sum)
#data_ex <- subset(data_o_per_norm,select=c(A0rep3,A14rep3,A5rep3,A8rep3,B0rep3,B14rep3,B5rep3,B8rep3,C0rep4,C14rep4,C5rep4,C8rep4))

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_o_per_norm[is.na(data_o_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_o_per_norm<-data_o_per_norm[order(rowSums(data_o_per_norm),decreasing = T),]
write.table(data_o_per_norm,"data_o_per_norm.xls",quote =F,sep="\t")

#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_ex_list<-rownames(data_o_per_norm)[1:N]
new_x<-rbind(data_o_per_norm[row.names(data_o_per_norm) %in% data_ex_list,],
             others=rowSums(data_o_per_norm[!rownames(data_o_per_norm) %in% data_ex_list,]))


#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)))#,id.vars = c('Taxonomy'))
#datm$variable <- factor(datm$variable,levels=c("A0rep1","A0rep2",  "A0rep3",  "A0rep4",  "A0rep5", "A5rep1",  "A5rep2" , "A5rep3" , "A5rep4" , "A5rep5", 
#                                               "A8rep1" , "A8rep2",  "A8rep3",  "A8rep4" , "A8rep5" ,'A14rep1', 'A14rep2', 'A14rep3', 'A14rep4' ,'A14rep5',
#                                               "B0rep1","B0rep2",  "B0rep3",  "B0rep4",  "B0rep5", "B5rep1",  "B5rep2" , "B5rep3" , "B5rep4" , "B5rep5",
#                                               "B8rep1" , "B8rep2",  "B8rep3",  "B8rep4" , "B8rep5" ,'B14rep1', 'B14rep2', 'B14rep3', 'B14rep4' ,'B14rep5',
#                                               "C0rep1","C0rep2",  "C0rep3",  "C0rep4",  "C0rep5", "C5rep1",  "C5rep2" , "C5rep3" , "C5rep4" , "C5rep5",
#                                               "C8rep1" , "C8rep2",  "C8rep3",  "C8rep4" , "C8rep5" ,'C14rep1', 'C14rep2', 'C14rep3', 'C14rep4' ,'C14rep5'))
#作图
#datm$variable <- factor(datm$variable,levels=c("A0rep3","A5rep3","A8rep3",'A14rep3',"B0rep3","B5rep3","B8rep3","B14rep3","C0rep4","C5rep4","C8rep4","C14rep4"))
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Order.pdf",width = 10, height = 8)


data_f <- DATA
data_f <- subset(data_f,data_f$f != " " & data_f$f != "")

name <- colnames(data)[1]
data_f_per <- data.frame(tapply(data_f[,c(name)], data_f$f, sum))
colnames(data_f_per) <- name

for(i in as.character(colnames(data)[2:length(colnames(data))])){
  name <- as.character(eval(parse(text="i")))
  output <- data.frame(tapply(data_f[,c(name)], data_f$f, sum))
  colnames(output) <- name
  data_f_per <- cbind(data_f_per,output)
}

# 新建data_norm数据框,赋初值
data_f_per_sum=apply(data_f_per, 2, sum)
# 计算每个样本的菌种总丰度
data_f_per_norm <- data_f_per
for(i in 1:ncol(data_f_per))
{
  for(j in 1:nrow(data_f_per)){
    data_f_per_norm[j,i]=data_f_per[j,i]/data_f_per_sum[i]
    # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
  }
}

apply(data_f_per_norm, 2, sum)
#data_ex <- subset(data_f_per_norm,select=c(A0rep3,A14rep3,A5rep3,A8rep3,B0rep3,B14rep3,B5rep3,B8rep3,C0rep4,C14rep4,C5rep4,C8rep4))

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_f_per_norm[is.na(data_f_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_f_per_norm<-data_f_per_norm[order(rowSums(data_f_per_norm),decreasing = T),]
write.table(data_f_per_norm,"data_f_per_norm.xls",quote =F,sep="\t")

#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_ex_list<-rownames(data_f_per_norm)[1:N]
new_x<-rbind(data_f_per_norm[row.names(data_f_per_norm) %in% data_ex_list,],
             others=rowSums(data_f_per_norm[!rownames(data_f_per_norm) %in% data_ex_list,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)))#,id.vars = c('Taxonomy'))
#datm$variable <- factor(datm$variable,levels=c("A0rep1","A0rep2",  "A0rep3",  "A0rep4",  "A0rep5", "A5rep1",  "A5rep2" , "A5rep3" , "A5rep4" , "A5rep5", 
#                                               "A8rep1" , "A8rep2",  "A8rep3",  "A8rep4" , "A8rep5" ,'A14rep1', 'A14rep2', 'A14rep3', 'A14rep4' ,'A14rep5',
#                                               "B0rep1","B0rep2",  "B0rep3",  "B0rep4",  "B0rep5", "B5rep1",  "B5rep2" , "B5rep3" , "B5rep4" , "B5rep5",
#                                               "B8rep1" , "B8rep2",  "B8rep3",  "B8rep4" , "B8rep5" ,'B14rep1', 'B14rep2', 'B14rep3', 'B14rep4' ,'B14rep5',
#                                               "C0rep1","C0rep2",  "C0rep3",  "C0rep4",  "C0rep5", "C5rep1",  "C5rep2" , "C5rep3" , "C5rep4" , "C5rep5",
#                                               "C8rep1" , "C8rep2",  "C8rep3",  "C8rep4" , "C8rep5" ,'C14rep1', 'C14rep2', 'C14rep3', 'C14rep4' ,'C14rep5'))
#作图
#datm$variable <- factor(datm$variable,levels=c("A0rep3","A5rep3","A8rep3",'A14rep3',"B0rep3","B5rep3","B8rep3","B14rep3","C0rep4","C5rep4","C8rep4","C14rep4"))
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Family.pdf",width = 10, height = 8)

data_g <- DATA
data_g <- subset(data_g,data_g$g != " " & data_g$g != "")

name <- colnames(data)[1]
data_g_per <- data.frame(tapply(data_g[,c(name)], data_g$g, sum))
colnames(data_g_per) <- name

for(i in as.character(colnames(data)[2:length(colnames(data))])){
  name <- as.character(eval(parse(text="i")))
  output <- data.frame(tapply(data_g[,c(name)], data_g$g, sum))
  colnames(output) <- name
  data_g_per <- cbind(data_g_per,output)
}

# 新建data_norm数据框,赋初值
data_g_per_sum=apply(data_g_per, 2, sum)
# 计算每个样本的菌种总丰度
data_g_per_norm <- data_g_per
for(i in 1:ncol(data_g_per))
{
  for(j in 1:nrow(data_g_per)){
    data_g_per_norm[j,i]=data_g_per[j,i]/data_g_per_sum[i]
    # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
  }
}

apply(data_g_per_norm, 2, sum)
#data_ex <- subset(data_g_per_norm,select=c(A0rep3,A14rep3,A5rep3,A8rep3,B0rep3,B14rep3,B5rep3,B8rep3,C0rep4,C14rep4,C5rep4,C8rep4))

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_g_per_norm[is.na(data_g_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_g_per_norm<-data_g_per_norm[order(rowSums(data_g_per_norm),decreasing = T),]
write.table(data_g_per_norm,"data_g_per_norm.xls",quote =F,sep="\t")

#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_ex_list<-rownames(data_g_per_norm)[1:N]
new_x<-rbind(data_g_per_norm[row.names(data_g_per_norm) %in% data_ex_list,],
             others=rowSums(data_g_per_norm[!rownames(data_g_per_norm) %in% data_ex_list,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)))#,id.vars = c('Taxonomy'))
#datm$variable <- factor(datm$variable,levels=c("A0rep1","A0rep2",  "A0rep3",  "A0rep4",  "A0rep5", "A5rep1",  "A5rep2" , "A5rep3" , "A5rep4" , "A5rep5", 
#                                               "A8rep1" , "A8rep2",  "A8rep3",  "A8rep4" , "A8rep5" ,'A14rep1', 'A14rep2', 'A14rep3', 'A14rep4' ,'A14rep5',
#                                               "B0rep1","B0rep2",  "B0rep3",  "B0rep4",  "B0rep5", "B5rep1",  "B5rep2" , "B5rep3" , "B5rep4" , "B5rep5",
#                                               "B8rep1" , "B8rep2",  "B8rep3",  "B8rep4" , "B8rep5" ,'B14rep1', 'B14rep2', 'B14rep3', 'B14rep4' ,'B14rep5',
#                                               "C0rep1","C0rep2",  "C0rep3",  "C0rep4",  "C0rep5", "C5rep1",  "C5rep2" , "C5rep3" , "C5rep4" , "C5rep5",
#                                               "C8rep1" , "C8rep2",  "C8rep3",  "C8rep4" , "C8rep5" ,'C14rep1', 'C14rep2', 'C14rep3', 'C14rep4' ,'C14rep5'))
#作图
#datm$variable <- factor(datm$variable,levels=c("A0rep3","A5rep3","A8rep3",'A14rep3',"B0rep3","B5rep3","B8rep3","B14rep3","C0rep4","C5rep4","C8rep4","C14rep4"))
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Genus.pdf",width = 10, height = 8)

#############################################################################
setwd("../../2.OTU table/")

DATA=read.table("otu_table_tax.txt",header = T,row.names = 1,sep = "\t",quote = "")
group=read.table("group.csv",sep=",",header = T)
row.names(group) <- group$Sample
data <- DATA[,row.names(group)]
# 给样本命名,结果如下：

colnames(group) <- c("Sample","Group")
#group$Group <- factor(as.character(group$Group),levels=c("A0","A5","A8","A14","B0","B5","B8","B14","C0","C5","C8","C14"))
?as.factor

data_norm=data
# 新建data_norm数据框,赋初值
sample_sum=apply(data, 2, sum)
# 计算每个样本的菌种总丰度

for(i in 1:ncol(data))
{
  for(j in 1:nrow(data)){
    data_norm[j,i]=data[j,i]/sample_sum[i]
    # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
  }
}

apply(data_norm, 2, sum)

#############################################################################

dir.create("../5.Alpha diversity")
dir.create("../5.Alpha diversity/5.1 Alpha指数")
dir.create("../5.Alpha diversity/5.2 Alpha多样性指数条形图")
dir.create("../5.Alpha diversity/5.3 Alpha多样性指数QQ图")
dir.create("../5.Alpha diversity/5.4 Alpha多样性指数箱型图")
dir.create("../5.Alpha diversity/5.5 Alpha多样性曲线")

setwd("../5.Alpha diversity")
library(vegan)
# 加载VEGAN包
?diversity
data_norm <- data.frame(t(data_norm))
data_norm_shannon=diversity(data_norm, "shannon")
# 使用VEGAN包中的diversity函数计算每个样品的Shannon Alpha多样性指数
data_ggplot=data.frame(data_norm_shannon)
# 多样性计算结果如下：
data_ggplot=data.frame(data_ggplot, group["Group"])
# 添加分组信息,如下：

group_EMFAC1=data_ggplot$data_norm_shannon[data_ggplot$Group=='EMFAC1']
group_SJK=data_ggplot$data_norm_shannon[data_ggplot$Group=='SJK'] 

setwd("../5.Alpha diversity/5.2 Alpha多样性指数条形图")
pdf("group_EMFAC1.hist.pdf")
hist(group_EMFAC1)
dev.off()
pdf("group_SJK.hist.pdf")
hist(group_SJK)
dev.off()
#############################################################################

# 绘制多样性指数分布直方图
# 观察形状若为倒钟形那便是接近正态分布的
setwd("../5.3 Alpha多样性指数QQ图")
pdf("group_EMFAC1.qqnorm.pdf")
qqnorm(group_A0)
dev.off()
pdf("group_SJK.qqnorm.pdf")
qqnorm(group_SJK)
dev.off()

# 绘制多样性指数QQ图
# 观察形状是一条连接主对角线的线那便是接近正态分布
shapiro.test(data_ggplot$data_norm_shannon)
# 夏皮罗-威尔克(Shapiro-Wilk)检验正态性,p>0.05,接受原假设,符合正态分布
bartlett.test(data_norm_shannon~Group,data=data_ggplot)
# 巴特利特(Bartlett)检验方差齐性,p>0.05,接受原假设,即两样本数据方差齐
with(data_ggplot,t.test(formula=data_norm_shannon~Group,conf.level=0.95))
# T检验,p>0.05,没法拒绝原假设,两组Shannon多样性指数差异不显著

setwd("../5.4 Alpha多样性指数箱型图")
library(ggplot2)
# 加载R包ggplot2
alpha_boxplot=ggplot(data_ggplot, aes(x=Group, y=data_norm_shannon, fill=Group))+
  # 添加数据、xy值、 颜色参数给画图函数ggplot
  geom_boxplot()+
  # 盒图
  labs(title="Alpha diversity", x="Group", y="Shannon index")+
  # 标题
  theme(plot.title=element_text(hjust=0.5), legend.title=element_blank())
# 标题居中

pdf('Alpha box result.pdf')
alpha_boxplot
dev.off()

#############################################################################

setwd("../5.5 Alpha多样性曲线")
data <- apply(data,2,function(x)x*50000/(sum(x)))
write.table(data,"../../2.OTU table/otu_table_norm_50000.xls",sep="\t",quote=F,row.names=T)
data<-data[order(rowSums(data),decreasing = T),]
data <- apply(data,2,round)
data <- as.data.frame(t(data))
library(vegan)    #用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样
library(picante)
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[2, ]    #Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon 指数
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
  else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}


#统计 data 丰度表中各样本的 Shannon 指数，对数底数使用 e
shannon_index <- alpha_index(data, method = 'shannon', base = exp(1))
write.csv(shannon_index,"../5.1 Alpha指数/shannon_index.csv")
richness_index <- alpha_index(data, method = 'richness', base = exp(1))
write.csv(richness_index,"../5.1 Alpha指数/richness_index.csv")
chao1_index <- alpha_index(data, method = 'chao1', base = exp(1))
write.csv(chao1_index,"../5.1 Alpha指数/chao1_index.csv")
ace_index <- alpha_index(data, method = 'ace', base = exp(1))
write.csv(ace_index,"../5.1 Alpha指数/ace_index.csv")
simpson_index <- alpha_index(data, method = 'simpson', base = exp(1))
write.csv(simpson_index,"../5.1 Alpha指数/simpson_index.csv")
pielou_index <- alpha_index(data, method = 'pielou', base = exp(1))
write.csv(pielou_index,"../5.1 Alpha指数/pielou_index.csv")
gc_index <- alpha_index(data, method = 'gc', base = exp(1))
write.csv(gc_index,"../5.1 Alpha指数/gc_index.csv")

#根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare
}

plot_richness <- data.frame()
for (n in 1:5) {
  richness_curves <- alpha_curves(data, step = 2000, method = 'richness')
  
  for (i in names(richness_curves)) {
    richness_curves_i <- (richness_curves[[i]])
    richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
    plot_richness <- rbind(plot_richness, richness_curves_i)
  }
}
library(doBy)
library(ggplot2)
#计算均值 ± 标准差（doBy 包中的 summaryBy() 函数）
plot_richness_stat <- summaryBy(alpha~sample+rare, plot_richness, FUN = c(mean, sd))
plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)
plot_richness_stat[which(plot_richness_stat$rare == 0),'alpha.sd'] <- NA
#ggplot2 作图

ggplot(plot_richness_stat, aes(rare, alpha.mean, color = sample)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = alpha.mean - alpha.sd, ymax = alpha.mean + alpha.sd), width = 500) +
  labs(x = 'number of normolised reads', y = 'Richness', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = min(rowSums(data)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 1000000, 200000), labels = as.character(seq(0, 1000000,200000)))
ggsave("richness.pdf",width = 20,height = 20)
dev.off()
#############################################################################

#Shannon指数稀释曲线等
##若简单的“geom_line()”样式波动幅度过大，不平滑等，可以尝试拟合曲线的样式
#获得作图数据。前面多生成一个点，使得 Shannon 拟合曲线更加平滑（你把 shannon_curves1 注释掉就知道我说的啥了）
shannon_curves1 <- alpha_curves(data, step = 200, rare = 200, method = 'shannon')
shannon_curves2 <- alpha_curves(data, step = 2000, method = 'shannon')
shannon_curves <- c(shannon_curves1, shannon_curves2)

plot_shannon <- data.frame()
for (i in 1:length(shannon_curves)) {
  shannon_curves_i <- shannon_curves[[i]]
  shannon_curves_i <- data.frame(rare = names(shannon_curves_i), alpha = shannon_curves_i, sample = names(shannon_curves)[i], stringsAsFactors = FALSE)
  plot_shannon <- rbind(plot_shannon, shannon_curves_i)
}

rownames(plot_shannon) <- NULL
plot_shannon$rare <- as.numeric(plot_shannon$rare)
plot_shannon$alpha <- as.numeric(plot_shannon$alpha)
plot_shannon <- plot_shannon[order(plot_shannon$sample, plot_shannon$rare), ]

#ggplot2 作图（使用到 ggalt 包的 geom_xspline() 绘制平滑拟合线）
library(ggalt)    #若未加载时先加载
pdf("Shannon Index.pdf",width = 20,height = 20)
ggplot(plot_shannon, aes(rare, alpha, color = sample)) +
  geom_xspline() +
  labs(x = 'Number of normolised reads', y = 'Shannon', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = min(rowSums(data)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 100000, 20000), labels = as.character(seq(0, 10000, 2000)))
dev.off()

#Rank-abundance曲线在R中的绘制方法
#虽然自己排序也很简单，但对于我这样的懒人来讲，还是导个包统计省事点……
library(BiodiversityR)

#统计（BiodiversityR 包 rankabundance() 实现 data 排序）
data_relative <- data/rowSums(data)
rank_dat <- data.frame()
for (i in rownames(data_relative)) {
  rank_dat_i <- data.frame(rankabundance(subset(data_relative, rownames(data_relative) == i), digits = 6))[1:2]
  rank_dat_i$sample <- i
  rank_dat <- rbind(rank_dat, rank_dat_i)
}
rank_dat <- subset(rank_dat, rank_dat$abundance != 0)

#ggplot2 作图
dev.off()
ggplot(rank_dat, aes(rank, log(abundance, 10), color = sample)) +
  geom_line() +
  labs(x = 'datas rank', y = 'Relative adundance (%)', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  scale_y_continuous(breaks = 0:-5, labels = c('100', '10', '1', '0.1', '0.01', '0.001'), limits = c(-5, 0))
ggsave("Rank-abundance.pdf",width = 20,height = 20)

#############################################################################

setwd("../../")
dir.create("6.beta diversity")
setwd("6.beta diversity")
library(ggplot2)
library(ggdendro)
library(vegan)
#读取数据
#?substring


#采用Bray Curtis方法，如需要更换其他方法，可在method参数中调整
beta_bray<-vegdist(t(data_g_per),method="bray")
#建树
hc<-hclust(beta_bray)
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type = "rectangle")
#绘图
pdf("normal.tree.pdf",width = 30, height = 16)
ggplot(dend_data$segments) + 
  theme_dendro()+
  scale_x_discrete(expand = c(0,1))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            size = 5,check_overlap = T,angle=45,vjust = 3,
            nudge_y = -0.02)
dev.off()
#PCA分析（主成分分析）
###############################################################################

library(ggbiplot)

taxa<-t(data_g_per)
taxa[is.na(taxa)]<-0
taxa <- taxa[,which(colSums(taxa) > 0)]
meta=read.table("../2.OTU table/group.csv",sep=",",header = T,row.names = 1)

#计算PCA值
pca<-prcomp(taxa,scale. = T)
#作图
pdf("PCA.pdf")
ggbiplot(pca, obs.scale = 2, var.scale = 1,var.axes = F,
         groups = meta$group, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()
###############################################################################

library(ade4)
library(ggplot2)
library(RColorBrewer)
#?substring
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table

taxa<-t(data_g_per)
taxa[is.na(taxa)]<-0
taxa <- taxa[,which(colSums(taxa) > 0)]
meta=read.csv("../2.OTU table/group.csv",sep=",",header = T)

pca<- dudi.pca(taxa, scal = T, scan = FALSE)

#坐标轴解释量（前两轴）
pca_eig <- (pca$eig)[1:2] / sum(pca$eig)

#提取样本点坐标（前两轴）
sample_site <- data.frame({pca$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCA1', 'PCA2')

#以group为分组
sample_site$level<-factor(meta$group)

library(ggplot2)
pca_plot <- ggplot(sample_site, aes(PCA1, PCA2,color=level)) +
  theme_classic()+#去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 1.5)+  #可在这里修改点的透明度、大小
  # scale_color_manual(values = brewer.pal(12,"Set2")) + #可在这里修改点的颜色
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank()
  )+
  labs(x = paste('PCA1: ', round(100 * pca_eig[1], 2), '%'), y = paste('PCA2: ', round(100 * pca_eig[2], 2), '%')) 
ggsave("PCA-2.pdf")

library(ade4)
library(ggplot2)
library(RColorBrewer)
library(vegan)
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table


taxa<-t(data_g_per)
taxa[is.na(taxa)]<-0
taxa <- taxa[,which(colSums(taxa) > 0)]
meta=read.csv("../2.OTU table/group.csv",sep=",",header = T)


taxa <- taxa[,which(colSums(taxa) > 0)]
tab.dist<-vegdist(taxa,method='euclidean')
pcoa<- dudi.pco(tab.dist, scan = T,nf=3)
####此处需要根据柱形图选择nf数目
#坐标轴解释量（前两轴）
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

#提取样本点坐标（前两轴）
sample_site <- data.frame({pcoa$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

#以group作为分组
sample_site$level<-factor(group$Group)

library(ggplot2)
pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2,color=level)) +
  theme_classic()+#去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 1.5)+  #可在这里修改点的透明度、大小
  # scale_color_manual(values = brewer.pal(8,"Set2")) + #可在这里修改点的颜色
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank()
  )+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) 
ggsave("PCoA.pdf")



