rm(list=ls())
setwd("D://16S/usearch/")

DATA=read.table("otu_table_tax.txt",header = T,row.names = 1,sep = "\t",quote = "")
DATA <- DATA[,-50]
DATA$F6_2 <- DATA$F6
DATA$I9_2 <- DATA$I9
data <- DATA
# 给样本命名,结果如下：
group=read.table("group.csv",sep=",",header = T)
colnames(group) <- c("Sample","Group")

data_norm=data
# 新建data_norm数据框,赋初值
sample_sum=apply(data, 1, sum)
# 计算每个样本的菌种总丰度

for(i in 1:nrow(data))
{
  for(j in 1:ncol(data)){
    data_norm[i,j]=data[i,j]/sample_sum[i]
    # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
  }
}

apply(data_norm, 1, sum)

library(vegan)
# 加载VEGAN包
data_norm <- data.frame(t(data_norm))
data_norm_shannon=diversity(data_norm, "shannon")
# 使用VEGAN包中的diversity函数计算每个样品的Shannon Alpha多样性指数
data_ggplot=data.frame(data_norm_shannon)
# 多样性计算结果如下：
data_ggplot=data.frame(data_ggplot, group["Group"])
# 添加分组信息,如下：

group_R=data_ggplot$data_norm_shannon[data_ggplot$Group=='R']
group_S=data_ggplot$data_norm_shannon[data_ggplot$Group=='S']
group_T=data_ggplot$data_norm_shannon[data_ggplot$Group=='T'] 
group_U=data_ggplot$data_norm_shannon[data_ggplot$Group=='U'] 
group_V=data_ggplot$data_norm_shannon[data_ggplot$Group=='V'] 
group_W=data_ggplot$data_norm_shannon[data_ggplot$Group=='W'] 
group_X=data_ggplot$data_norm_shannon[data_ggplot$Group=='X'] 
group_Y=data_ggplot$data_norm_shannon[data_ggplot$Group=='Y'] 

hist(group_R)
hist(group_S)
hist(group_T)
hist(group_U)
hist(group_V)
hist(group_W)
hist(group_X)
hist(group_Y)

# 绘制多样性指数分布直方图
# 观察形状若为倒钟形那便是接近正态分布的
qqnorm(group_R)
qqnorm(group_S)
qqnorm(group_T)
qqnorm(group_U)
qqnorm(group_V)
qqnorm(group_W)
qqnorm(group_X)
qqnorm(group_Y)
# 绘制多样性指数QQ图
# 观察形状是一条连接主对角线的线那便是接近正态分布
shapiro.test(data_ggplot$data_norm_shannon)
# 夏皮罗-威尔克(Shapiro-Wilk)检验正态性,p>0.05,接受原假设,符合正态分布
bartlett.test(data_norm_shannon~Group,data=data_ggplot)
# 巴特利特(Bartlett)检验方差齐性,p>0.05,接受原假设,即两样本数据方差齐
with(data_ggplot,t.test(formula=data_norm_shannon~Group,conf.level=0.95))
# T检验,p>0.05,没法拒绝原假设,两组Shannon多样性指数差异不显著

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

pdf('result.pdf')
alpha_boxplot
dev.off()
# 保存结果,打开result.pdf文件,结果如下：

################################################################################
library(dplyr)
library(reshape2)
library(tidyr)
#?substring
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_k <- DATA[,1:50]
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

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_k_per_norm[is.na(data_k_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_k_per_norm<-data_k_per_norm[order(rowSums(data_k_per_norm),decreasing = T),]
#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_k_per_norm_list<-rownames(data_k_per_norm)[1:N]
new_x<-rbind(data_k_per_norm[row.names(data_k_per_norm) %in% data_k_per_norm_list,],
             others=rowSums(data_k_per_norm[!rownames(data_k_per_norm) %in% data_k_per_norm,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)),id.vars = c('Taxonomy'))
datm$variable <- factor(datm$variable,levels=as.character(group$Sample))
#作图

ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Domain.pdf",width = 10, height = 8)

################################################################################

################################################################################
#?substring
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_p <- DATA[,c(1:49,51)]
data_p <- subset(data_p,data_p$p != " " & data_p$p != "")
table(data_p$p)
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

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_p_per_norm[is.na(data_p_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_p_per_norm<-data_p_per_norm[order(rowSums(data_p_per_norm),decreasing = T),]
#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_p_per_norm_list<-rownames(data_p_per_norm)[1:N]
new_x<-rbind(data_p_per_norm[row.names(data_p_per_norm) %in% data_p_per_norm_list,],
             others=rowSums(data_p_per_norm[!rownames(data_p_per_norm) %in% data_p_per_norm,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)),id.vars = c('Taxonomy'))
datm$variable <- factor(datm$variable,levels=as.character(group$Sample))
#作图

ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Phylum.pdf",width = 10, height = 8)

################################################################################
#?substring
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_c <- DATA[,c(1:49,52)]

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

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_c_per_norm[is.na(data_c_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_c_per_norm<-data_c_per_norm[order(rowSums(data_c_per_norm),decreasing = T),]
#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_c_per_norm_list<-rownames(data_c_per_norm)[1:N]
new_x<-rbind(data_c_per_norm[row.names(data_c_per_norm) %in% data_c_per_norm_list,],
             others=rowSums(data_c_per_norm[!rownames(data_c_per_norm) %in% data_c_per_norm,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)),id.vars = c('Taxonomy'))
datm$variable <- factor(datm$variable,levels=as.character(group$Sample))
#作图

ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Class.pdf",width = 10, height = 8)

################################################################################
#?substring
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_o <- DATA[,c(1:49,53)]
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

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_o_per_norm[is.na(data_o_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_o_per_norm<-data_o_per_norm[order(rowSums(data_o_per_norm),decreasing = T),]
#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_o_per_norm_list<-rownames(data_o_per_norm)[1:N]
new_x<-rbind(data_o_per_norm[row.names(data_o_per_norm) %in% data_o_per_norm_list,],
             others=rowSums(data_o_per_norm[!rownames(data_o_per_norm) %in% data_o_per_norm,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)),id.vars = c('Taxonomy'))
datm$variable <- factor(datm$variable,levels=as.character(group$Sample))
#作图
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Order.pdf",width = 10, height = 8)


################################################################################
#?substring
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_f <- DATA[,c(1:49,54)]
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

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_f_per_norm[is.na(data_f_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_f_per_norm<-data_f_per_norm[order(rowSums(data_f_per_norm),decreasing = T),]
#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_f_per_norm_list<-rownames(data_f_per_norm)[1:N]
new_x<-rbind(data_f_per_norm[row.names(data_f_per_norm) %in% data_f_per_norm_list,],
             others=rowSums(data_f_per_norm[!rownames(data_f_per_norm) %in% data_f_per_norm,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)),id.vars = c('Taxonomy'))
datm$variable <- factor(datm$variable,levels=as.character(group$Sample))
#作图
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Family.pdf",width = 10, height = 8)


################################################################################
#?substring
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_g <- DATA[,c(1:49,55)]
data_g[4,50]
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

library(reshape2)
library(ggplot2)

#将文件中的NA值改为０
data_g_per_norm[is.na(data_g_per_norm)] <- 0
#将数据中的种类根据数量的多少进行排序
data_g_per_norm<-data_g_per_norm[order(rowSums(data_g_per_norm),decreasing = T),]
#Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
N<-10
data_g_per_norm_list<-rownames(data_g_per_norm)[1:N]
new_x<-rbind(data_g_per_norm[row.names(data_g_per_norm) %in% data_g_per_norm_list,],
             others=rowSums(data_g_per_norm[!rownames(data_g_per_norm) %in% data_g_per_norm,]))

#合并数据
datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)),id.vars = c('Taxonomy'))
datm$variable <- factor(datm$variable,levels=as.character(group$Sample))
#作图
ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
  xlab("")+
  ylab("")+
  geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
ggsave("Genus.pdf",width = 10, height = 8)

###############################################################################
#聚类分析#设置工作路径
#载入工作包
library(ggplot2)
library(ggdendro)
library(vegan)
#读取数据
#?substring
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_g <- DATA[,c(1:49,55)]
data_g[4,50]
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
meta=read.table("group.csv",sep=",",header = T,row.names = 1)

#计算PCA值
pca<-prcomp(taxa,scale. = T)
#作图
pdf("PCA.pdf")
ggbiplot(pca, obs.scale = 2, var.scale = 1,var.axes = F,
         groups = meta$Group, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()
###############################################################################
setwd("D://16S/usearch/")

library(ade4)
library(ggplot2)
library(RColorBrewer)
#?substring
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_g <- DATA[,c(1:49,55)]
name <- colnames(data)[1]
data_g_per <- data.frame(tapply(data_g[,c(name)], data_g$g, sum))
colnames(data_k_per) <- name

for(i in as.character(colnames(data)[2:length(colnames(data))])){
  name <- as.character(eval(parse(text="i")))
  output <- data.frame(tapply(data_g[,c(name)], data_g$g, sum))
  colnames(output) <- name
  data_g_per <- cbind(data_g_per,output)
}
taxa<-t(data_g_per)
taxa[is.na(taxa)]<-0
taxa <- taxa[,which(colSums(taxa) > 0)]
meta=read.csv("group.csv",sep=",",header = T)
group <- meta[c(-7,-22),]
row.names(group) <- group$Sample
meta <- data.frame(group$Group)
colnames(meta) <- c("Group")
row.names(meta) <- group$Sample
pca<- dudi.pca(taxa, scal = T, scan = FALSE)

#坐标轴解释量（前两轴）
pca_eig <- (pca$eig)[1:2] / sum(pca$eig)

#提取样本点坐标（前两轴）
sample_site <- data.frame({pca$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCA1', 'PCA2')

#以group为分组
sample_site$level<-factor(meta$Group)

library(ggplot2)
pca_plot <- ggplot(sample_site, aes(PCA1, PCA2,color=level)) +
  theme_classic()+#去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 1.5)+  #可在这里修改点的透明度、大小
  scale_color_manual(values = brewer.pal(8,"Set2")) + #可在这里修改点的颜色
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
DATA=read.csv("otu_table_tax.xls",header = T,row.names = 1,sep = "\t",quote = "")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#?read.table
?group_by
data <- DATA[,1:49]
data_g <- DATA[,c(1:49,55)]
name <- colnames(data)[1]
data_g_per <- data.frame(tapply(data_g[,c(name)], data_g$g, sum))
colnames(data_k_per) <- name

for(i in as.character(colnames(data)[2:length(colnames(data))])){
  name <- as.character(eval(parse(text="i")))
  output <- data.frame(tapply(data_g[,c(name)], data_g$g, sum))
  colnames(output) <- name
  data_g_per <- cbind(data_g_per,output)
}
taxa<-t(data_g_per)
taxa[is.na(taxa)]<-0
taxa <- taxa[,which(colSums(taxa) > 0)]
meta=read.csv("group.csv",sep=",",header = T)
group <- meta[c(-7,-22),]
row.names(group) <- group$Sample
meta <- data.frame(group$Group)
colnames(meta) <- c("Group")
row.names(meta) <- group$Sample

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
  scale_color_manual(values = brewer.pal(8,"Set2")) + #可在这里修改点的颜色
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank()
  )+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) 
ggsave("PCoA.pdf")

