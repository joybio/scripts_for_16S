rm(list=ls())
setwd("~/Downloads/2021.8.30.16S/Results/")
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
library(picante)
library(doBy)
library(ggalt)    #若未加载时先加载
library(BiodiversityR)
library(ggdendro)
library(ggtree)
library(phangorn)
library(ggbiplot)
library(ade4)

dir.create("4.Taxonomy/")
dir.create("5.Alpha diversity/")
dir.create("6.Beta diversity/")


DATA=read.table("2.OTU table/otu_table_tax.txt",header = T,row.names = 1,sep = "\t",quote = "",check.names = F)
?separate
?str_replace
Taxonomy <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
DATA$taxonomy <- str_replace_all(DATA$taxonomy," .__","")
DATA$taxonomy <- str_replace_all(DATA$taxonomy,"k__","")
data <- DATA %>% separate(col = "taxonomy",into = Taxonomy,sep=";")
DATA <- subset(DATA,select = -taxonomy)
DATA_norm <- apply(DATA,2,function(x)x*50000/(sum(x)))
write.csv(DATA_norm,"2.OTU table/otu_table_norm_50000.csv")
DATA=read.table("2.OTU table/otu_table_tax.txt",header = T,row.names = 1,sep = "\t",quote = "",check.names = F)
group=read.table("2.OTU table/group.csv",sep=",",header = T)
colnames(group) <- c("Sample","Group")
#data <- separate(DATA, taxonomy, c("k","p","c","o","f","g","s"), sep = ";")
#group$Group <- factor(as.character(group$Group),levels=c("A0","A5","A8","A14","B0","B5","B8","B14","C0","C5","C8","C14"))
#group$group <- substr(group$Group,1,1)
row.names(group) <- group$Sample
table(group$Group)
#?read.table
?group_by
?dir.create
#Alpha多样性指数计算
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[2, ]    #Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
  else if (method == 'shannon') result <- vegan::diversity(x, index = 'shannon', base = base)    #Shannon 指数
  else if (method == 'simpson') result <- vegan::diversity(x, index = 'simpson')    #Gini-Simpson 指数
  else if (method == 'pielou') result <- vegan::diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
  else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

#Alpha曲线作图
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

DATA <- DATA[,row.names(group)]
index <- c("shannon","richness","chao1","ace","simpson","pielou","gc")
#otu tree
dir.create("2.OTU table/OTU_tree")
otu_DATA <- as.data.frame(t(DATA))
up=upgma(vegdist(otu_DATA, method="bray"))
pdf('2.OTU table/OTU_tree/upgma.pdf')
opar=par(no.readonly=TRUE)
# 生成图形参数列表
par(mfrow=c(3, 2), col.main="red", family="serif")
# par设置：按行填充,3行，2列，标题颜色，字体（罗马）
par(mai=c(0.2, 0.2, 0.2, 0.2))
# par设置：每个图形距边距离（英寸）
plot(up, main="by default")
plot(up, type="phylogram", main="phylogram")  # 默认
plot(up, type="cladogram", main="cladogram")
plot(up, type="fan", main="fan")
plot(up, type="unrooted", main="unrooted")
plot(up, type="radial", main="radial")
par(opar)
# 关闭par
dev.off()
# 按门水平建树并上色
## 给每个OTU按门分类分组，此处可以更改为其它分类级别，如纲、目等，即phylum替换为order或class即可
group <- read.csv("2.OTU table/group.csv",header = T,row.names = 1)
groupInfo <- split(row.names(otu_DATA), group$Group)
## 将分组信息添加到树中
tree <- groupOTU(up, groupInfo)
pdf('2.OTU table/OTU_tree/upgma.ggtree.circular.pdf')
ggtree(tree,aes(color=group),layout="circular")  + geom_tiplab2()   #+ geom_text(aes(label=node))
dev.off()
pdf('2.OTU table/OTU_tree/upgma.ggtree.radial.pdf')
ggtree(tree,aes(color=group),layout="radial")  + geom_tiplab2()   #+ geom_text(aes(label=node))
dev.off()

#Ｎ值代表选择数量排前40的物种,将剩下的物种合并成其他
N<-40
?upgma
for (i in Taxonomy){
  ID <- eval(i)
  top40_OTU <- paste0("2.OTU table/OTU_tree/",ID,"/")
  dir.create(top40_OTU)
  Taxonomy_ID <- subset(data,data[,ID] != " " & data[,ID] !=  "NA" & data[,ID] != "" & data[,ID] !=  "")
  OTU_oder_DATA <- Taxonomy_ID[order(rowSums(Taxonomy_ID[,row.names(group)]),decreasing = T),]
  new_x <- OTU_oder_DATA[1:40,row.names(group)]
  names <- row.names(new_x)
  norm = t(t(new_x)/colSums(new_x,na=T)) * 100 # normalization to total 100
  up=upgma(vegdist(new_x, method="bray"))
  data_N <- OTU_oder_DATA[1:40,]
  groupInfo <- split(row.names(new_x), data_N[,ID])
  tree <- groupOTU(up, groupInfo)
  circular_pdf <- paste0(top40_OTU,"upgma.top40.otu.ggtree.circular.pdf")
  ggtree(tree,aes(color=group),layout="circular")  + geom_tiplab2()   #+ geom_text(aes(label=node))
  ggsave(circular_pdf)
  radial_pdf <- paste0(top40_OTU,"upgma.top40.otu.ggtree.radial.pdf")
  ggtree(tree,aes(color=group),layout="radial")  + geom_tiplab2()   #+ geom_text(aes(label=node))
  ggsave(radial_pdf)
  rectangular_pdf <- paste0(top40_OTU,"upgma.top40.otu.ggtree.rectangular.pdf")
  ggtree(tree, aes(color=group))+  theme(legend.position = "right")+geom_tiplab(offset=0.1, size=1, align=TRUE, linesize=.5)
  ggsave(rectangular_pdf)
  p = ggtree(tree, aes(color=group))+  theme(legend.position = "right")+geom_tiplab(offset=0.1, size=1, align=TRUE, linesize=.5)
  tree_heatmap_pdf <- paste0(top40_OTU,"ggtree_heat_sample.pdf")
  gheatmap(p, norm, offset = .15, width=3, font.size=1, colnames_angle=-45, hjust=-.1)
  ggsave(tree_heatmap_pdf)
}
Taxonomy <- c("Phylum","Class","Order","Family","Genus","Species")
for (i in Taxonomy){
  ID <- eval(i)
  Alpha_directory_main <- paste0("5.Alpha diversity/",ID,"/")
  dir.create(Alpha_directory_main)
  Beta_directory_main <- paste0("6.beta diversity/",ID,"/")
  dir.create(Beta_directory_main)
  Alpha_directory_5.1 <- paste0(Alpha_directory_main,"5.1 Alpha_index/")
  Alpha_directory_5.2 <- paste0(Alpha_directory_main,"5.2 Alpha_hist/")
  Alpha_directory_5.3 <- paste0(Alpha_directory_main,"5.3 Alpha_qqnorm/")
  Alpha_directory_5.4 <- paste0(Alpha_directory_main,"5.4 Alpha_box/")
  Alpha_directory_5.5 <- paste0(Alpha_directory_main,"5.5 Alpha_curve/")
  dir.create(Alpha_directory_5.1)
  dir.create(Alpha_directory_5.2)
  dir.create(Alpha_directory_5.3)
  dir.create(Alpha_directory_5.4)
  dir.create(Alpha_directory_5.5)
  #物种组成丰度图
  Taxonomy_ID_per_csv <- paste0("4.Taxonomy/data_",eval(i),"_per.csv")
  Taxonomy_ID_per_pdf <- paste0("4.Taxonomy/data_",eval(i),"_per.pdf")
  #筛选有注释的Taxonomy
  Taxonomy_ID <- subset(data,data[,ID] != " " & data[,ID] !=  "NA" & data[,ID] != "")
  name <- colnames(DATA)[1]
  Taxonomy_ID_num <- data.frame(tapply(Taxonomy_ID[,c(name)], Taxonomy_ID[,ID], sum))
  colnames(Taxonomy_ID_num) <- name
  for(j in as.character(colnames(DATA)[2:length(colnames(DATA))])){
    name <- as.character(eval(parse(text="j")))
    output <- data.frame(tapply(Taxonomy_ID[,c(name)], Taxonomy_ID[,ID], sum))
    colnames(output) <- name
    Taxonomy_ID_num <- cbind(Taxonomy_ID_num,output)
  }
  Taxonomy_ID_num_sum=apply(Taxonomy_ID_num, 2, sum)
  # 计算相对丰度
  Taxonomy_ID_per_norm <- Taxonomy_ID_num
  for(m in 1:ncol(Taxonomy_ID_num))
  {
    for(n in 1:nrow(Taxonomy_ID_num)){
      Taxonomy_ID_per_norm[n,m]=Taxonomy_ID_num[n,m]/Taxonomy_ID_num_sum[m]
      # [每个样本的每个物种的丰度] / [该样本物种总丰度] = 相对丰度
    }
  }
  Taxonomy_ID_per_norm[is.na(Taxonomy_ID_per_norm)] <- 0
  #将数据中的种类根据数量的多少进行排序
  Taxonomy_ID_per_norm<-Taxonomy_ID_per_norm[order(rowSums(Taxonomy_ID_per_norm),decreasing = T),]
  ?write.table
  write.csv(Taxonomy_ID_per_norm,Taxonomy_ID_per_csv)
  #Ｎ值代表选择数量排前10的物种,将剩下的物种合并成其他
  N<-10
  data_ex_list<-rownames(Taxonomy_ID_per_norm)[1:N]
  new_x<-rbind(Taxonomy_ID_per_norm[row.names(Taxonomy_ID_per_norm) %in% data_ex_list,],
               others=rowSums(Taxonomy_ID_per_norm[!rownames(Taxonomy_ID_per_norm) %in% data_ex_list,]))
  
  #合并数据
  datm<-melt(cbind(new_x,Taxonomy=rownames(new_x)))#,id.vars = c('Taxonomy'))
  ggplot(datm,aes(x=variable,y=value,fill=Taxonomy))+
    xlab("")+
    ylab("")+
    geom_bar(position = "fill",stat = 'identity',width = 0.8)+ scale_y_continuous(expand = c(0,0)) +theme(axis.text.x=element_text(angle=45,vjust = 0.5))
  ggsave(Taxonomy_ID_per_pdf,width = 10, height = 8)
  t_Taxonomy_ID_num <- t(Taxonomy_ID_num)
  # 使用VEGAN包中的diversity函数计算每个样品的Shannon Alpha多样性指数
  Taxonomy_ID_shannon <- vegan::diversity(t_Taxonomy_ID_num, "shannon")
  # 多样性计算结果如下：
  Taxonomy_ID_shannon_ggplot=data.frame(Taxonomy_ID_shannon)
  # 添加分组信息,如下：
  Taxonomy_ID_shannon_ggplot=data.frame(Taxonomy_ID_shannon_ggplot, group["Group"])
  #绘制多样性指数hist和QQ图
  s_group <- unique(group$Group)
  for (o in s_group){
    group_hist <- paste0(Alpha_directory_5.2,"group_",i,o,".hist.pdf")
    group_qqnorm <- paste0(Alpha_directory_5.3,"group_",i,o,".qqnorm.pdf")
    group_plot = Taxonomy_ID_shannon_ggplot$Taxonomy_ID_shannon[Taxonomy_ID_shannon_ggplot$Group==o]
    pdf(group_hist)
    hist(group_plot)
    dev.off()
    pdf(group_qqnorm)
    qqnorm(group_plot)
    dev.off()
  }
  #标准化成50000条tags
  if (nrow(Taxonomy_ID_num) == 1){
    Taxonomy_ID_num_norm <- as.data.frame(t(Taxonomy_ID_num))
    }else{
      Taxonomy_ID_num_norm <- apply(Taxonomy_ID_num,2,function(x)x*50000/(sum(x)))
      Taxonomy_ID_num_norm <- Taxonomy_ID_num_norm[order(rowSums(Taxonomy_ID_num_norm),decreasing = T),]
      Taxonomy_ID_num_norm <- apply(Taxonomy_ID_num_norm,2,round)
      Taxonomy_ID_num_norm <- as.data.frame(t(Taxonomy_ID_num_norm))}
      #计算index并画boxplot
      for (p in index){
      specific_csv <- paste0(Alpha_directory_5.1,eval(p),"_index.csv")
      specific_index <- alpha_index(Taxonomy_ID_num_norm, method = eval(p), base = exp(1))
      specific_index_ggplot <- data.frame(specific_index, group["Group"])
      write.csv(specific_index_ggplot,specific_csv)
      index_pdf <- paste0(Alpha_directory_5.4,"Alpha diversity ",eval(p),".pdf")
      alpha_boxplot <- ggplot(specific_index_ggplot, aes(x=Group, y=specific_index, fill=Group))+
          # 添加数据、xy值、 颜色参数给画图函数ggplot
         geom_boxplot()+
          # 盒图
          labs(title="Alpha diversity", x="Group", y=eval(p))+
          # 标题
         theme(plot.title=element_text(hjust=0.5), legend.title=element_blank())
         # 标题居中
      ggsave(index_pdf,alpha_boxplot)
      #richness 曲线
      plot_richness <- data.frame()
      for (q in 1:5) {
        richness_curves <- alpha_curves(Taxonomy_ID_num_norm, step = 10000, method = 'richness')
        for (k in names(richness_curves)) {
          richness_curves_k <- (richness_curves[[k]])
          richness_curves_k <- data.frame(rare = names(richness_curves_k), alpha = richness_curves_k, sample = k, stringsAsFactors = FALSE)
          plot_richness <- rbind(plot_richness, richness_curves_k)
        }
      }
      plot_richness_pdf <- paste0(Alpha_directory_5.5,"Richness.pdf")
      plot_richness_stat <- summaryBy(alpha~sample+rare, plot_richness, FUN = c(mean, sd))
      plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)
      plot_richness_stat[which(plot_richness_stat$rare == 0),'alpha.sd'] <- NA
      ggplot(plot_richness_stat, aes(rare, alpha.mean, color = sample)) +
        geom_line() +
        geom_point() +
        geom_errorbar(aes(ymin = alpha.mean - alpha.sd, ymax = alpha.mean + alpha.sd), width = 500) +
        labs(x = 'number of normolised reads', y = 'Richness', color = NULL) +
        theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
        geom_vline(xintercept = min(rowSums(Taxonomy_ID_num_norm)), linetype = 2)
      ggsave(plot_richness_pdf,width = 20,height = 20)
      shannon_curves1 <- alpha_curves(Taxonomy_ID_num_norm, step = 5000, rare = 5000, method = 'shannon')
      shannon_curves2 <- alpha_curves(Taxonomy_ID_num_norm, step = 10000, method = 'shannon')
      shannon_curves <- c(shannon_curves1, shannon_curves2)
      #shannon曲线
      plot_shannon_pdf <- paste0(Alpha_directory_5.5,"Shannon.pdf")
      plot_shannon <- data.frame()
      for (k in 1:length(shannon_curves)) {
        shannon_curves_k <- shannon_curves[[k]]
        shannon_curves_k <- data.frame(rare = names(shannon_curves_k), alpha = shannon_curves_k, sample = names(shannon_curves)[k], stringsAsFactors = FALSE)
        plot_shannon <- rbind(plot_shannon, shannon_curves_k)
      }
      rownames(plot_shannon) <- NULL
      plot_shannon$rare <- as.numeric(plot_shannon$rare)
      plot_shannon$alpha <- as.numeric(plot_shannon$alpha)
      plot_shannon <- plot_shannon[order(plot_shannon$sample, plot_shannon$rare), ]
      ggplot(plot_shannon, aes(rare, alpha, color = sample)) +
        geom_xspline() +
        labs(x = 'Number of normolised reads', y = 'Shannon', color = NULL) +
        theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
        geom_vline(xintercept = min(rowSums(Taxonomy_ID_num_norm)), linetype = 2)
      ggsave(plot_shannon_pdf,width = 20,height = 20)
      #Rank-abundance
      plot_Rank_abundance.pdf <- paste0(Alpha_directory_5.5,"Rank_abundance.pdf")
      data_relative <- Taxonomy_ID_num_norm/rowSums(Taxonomy_ID_num_norm)
      rank_dat <- data.frame()
      for (l in rownames(data_relative)) {
        rank_dat_l <- data.frame(rankabundance(subset(data_relative, rownames(data_relative) == l), digits = 6))[1:2]
        rank_dat_l$sample <- l
        rank_dat <- rbind(rank_dat, rank_dat_l)
      }
      rank_dat <- subset(rank_dat, rank_dat$abundance != 0)
      ggplot(rank_dat, aes(rank, log(abundance, 10), color = sample)) +
        geom_line() +
        labs(x = 'datas rank', y = 'Relative adundance (%)', color = NULL) +
        theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent'))
      ggsave(plot_Rank_abundance.pdf,width = 20,height = 20)
      #采用Bray Curtis方法，如需要更换其他方法，可在method参数中调整
      beta_bray<-vegdist(t(Taxonomy_ID_per_norm),method="bray")
      #建树
      hc<-hclust(beta_bray)
      hcd <- as.dendrogram(hc)
      dend_data <- dendro_data(hcd, type = "rectangle")
      #tree绘图
      tree_pdf <- paste0(Beta_directory_main,"tree.pdf")
      ggplot(dend_data$segments) + 
        theme_dendro()+
        scale_x_discrete(expand = c(0,1))+
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
        geom_text(data = dend_data$labels, aes(x, y, label = label),
                size = 5,check_overlap = T,angle=45,vjust = 3)#,
      #nudge_y = -0.02)
      ggsave(tree_pdf,height = 15,width = 20)
      #PCA分析（主成分分析）
      taxa<-Taxonomy_ID_num_norm
      taxa[is.na(taxa)]<-0
      taxa <- taxa[,which(colSums(taxa) > 0)]
      meta <- read.table("2.OTU table/group.csv",sep=",",header = T,row.names = 1)
      #计算PCA值
      colnames(meta) <- c("Group")
      pca<-prcomp(taxa,scale. = T)
      #作图
      PCA_pdf <- paste0(Beta_directory_main,"PCA.pdf")
      ggbiplot(pca, obs.scale = 2, var.scale = 1,var.axes = F,
             groups = meta$Group, ellipse = TRUE, circle = TRUE) +
        scale_color_discrete(name = '') +
        theme(legend.direction = 'horizontal', legend.position = 'top')
      ggsave(PCA_pdf)
      pca<- dudi.pca(taxa, scal = T, scan = FALSE)
      #坐标轴解释量（前两轴）
      pca_eig <- (pca$eig)[1:2] / sum(pca$eig)
      #提取样本点坐标（前两轴）
      sample_site <- data.frame({pca$li})[1:2]
      sample_site$names <- rownames(sample_site)
      names(sample_site)[1:2] <- c('PCA1', 'PCA2')
      #以group为分组
      sample_site$level<-factor(meta$Group)
      PCA2_pdf <- paste0(Beta_directory_main,"PCA2.pdf")
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
      ggsave(PCA2_pdf)
      PCoA_pdf <- paste0(Beta_directory_main,"PCoA.pdf")
      tab.dist<-vegdist(taxa,method='euclidean')
      pcoa<- dudi.pco(tab.dist, scan = F,nf=3)
      pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
      sample_site <- data.frame({pcoa$li})[1:2]
      sample_site$names <- rownames(sample_site)
      names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
      sample_site$level<-factor(meta$Group)
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
      ggsave(PCoA_pdf)
    }
  }
}

