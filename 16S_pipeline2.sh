#!/bin/bash
#coding:utf-8 
###########################################################################################
#  注意命名过程不能有   -   .  等							  #
#  尤其是rep，一旦出现会覆盖掉rep							  #
###########################################################################################

### 使用 vsearch+usearch11+QIIME 处理双端测序数据
HOME=/media/ruizhi/NGS/2021.8.12.30.16S
### 操作空间
temp=${HOME}/temp2
### 临时文件存放目录
res=${HOME}/res2
### 结果文件存放目录
data=${HOME}
### 测序数据存放位置
refdb=/media/ruizhi/database/16S
### 参考测序数据库
err=${HOME}/err2
### 错误信息
log=${HOME}/log2


### 日志信息
########################################################################################### 
# 实验准备                                                                                #
# 实验数据为已合并为fasta格式的双端测序序列                                               #
# 安装usearch11、vsearch和qiime                                                           #
# 下载Greengene最新数据库，320MB                                                          #
# wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz     #
# 解压数据包后大小3.4G                                                                    #
# tar xvzf gg_13_8_otus.tar.gz                                                            #
###########################################################################################
echo "实验开始"
cd ${HOME}

mkdir -p Results/2.OTU\ table
mkdir temp2
mkdir res2
mkdir err2
mkdir log2
echo "双端测序数据合并"
#echo "Step_1.1：decompress paired reads sequences"
#gunzip ${data}/*/*.fastq.gz		### 每个样品存放于一个目录中
#prefix='ERR7765'			### 样本序号前缀
### merged.fastq 合并后的FASTQ格式reads
### merged.fasta 合并后的FASTA格式reads

echo "Step_1.2:Merge paired reads"
ls *1.fq | while read id;do(usearch -fastq_mergepairs ${data}/$id -reverse ${data}/$(basename $id '1.fq')2.fq  -relabel @ -fastqout ${temp}/$(basename $id '1.fq')merged.fastq -fastaout ${temp}/$(basename $id '1.fq')merged.fasta);done
echo "Step_2:Quality filtering"
ls ${temp}/*merged.fastq | while read id;do(usearch -fastq_filter $id -fastq_maxee 1.0 -fastaout ${temp}/$(basename $id 'merged.fastq')filtered.fasta -threads 2);done
cat ${temp}/*filtered.fasta > ${res}/merged.fasta
#cat ${temp}/*.filtered.fastq > ${res}/merged.fastq
echo ""
ls -lh ${res}/merged.fasta
echo ""

echo "Step_3:Dereplication"
## 方式1
usearch -fastx_uniques ${res}/merged.fasta -fastaout ${res}/unique.fasta --minuniquesize 2 -threads 2 -sizeout -relabel Uniq
## 方式2
# vsearch --derep_fulllength ${res}/filtered_v.fasta --output ${res}/unique_v.fa --minuniquesize 2 --threads 10 --sizeout


echo "Step_4:Discard singletons"
usearch -sortbysize ${res}/unique.fasta -fastaout ${res}/discard.fasta -minsize 2

### 聚类OTU ###
echo "Step_5:UPARSE-OTU"
## 采用步骤3中的方式1的执行结果
usearch -cluster_otus ${res}/unique.fasta -otus ${res}/otus_u51.fasta -uparseout ${res}/uparse.1.txt -relabel Otu -threads 10
## 采用步骤3中的执行结果
usearch -cluster_otus ${res}/discard.fasta -otus ${res}/otus_u52.fasta -uparseout ${res}/uparse.2.txt -relabel Otu -threads 10
## Denoise 降噪得到的OTU
usearch -unoise3 ${res}/unique.fasta -zotus ${res}/zotus_51.fasta -threads 30

## 采用QIIME的方式
echo "Picking OTUs"
#1 使用source + 绝对路径
#2 运行脚本一定使用bash，source才可使用
source ~/miniconda3/bin/activate python27
/root/miniconda3/envs/python27/bin/pick_otus.py -i ${res}/unique.fasta -o ${res} -m uclust -r ${refdb}/gg_13_8_otus/rep_set/97_otus.fasta --threads 30
echo "Picking a representative sequence set, one sequence from each OTU"
/root/miniconda3/envs/python27/bin/pick_rep_set.py -i ${res}/unique_otus.txt -f ${res}/unique.fasta -o ${res}/rep_set1.fasta
#conda activate base
### 统计fasta文件中序列的数量 ###
### 查看OTU数量  (分别为三种方式聚类得到的OTU数量以及一个降噪后得到的OTU数量) ###
echo "Number of OTU:"
grep '>' -c ${res}/otus_u51.fasta
grep '>' -c ${res}/otus_u52.fasta
grep '>' -c ${res}/zotus_51.fasta

echo "Step_6:Construct OTU table"
### Input should be reads before quality filtering and before discarding low-abundance unique sequences
### Reads must have sample indenfiers in the labels.
### The search database is the OTU (or ZOTU) sequences, which must have valid OTU identifiers in the labels.
### suggest making OTU tables for both 97% OTUs and denoised sequences (ZOTUs).
usearch -otutab ${res}/merged.fasta -otus ${res}/otus_u51.fasta -otutabout ${res}/otutab_u61.txt -mapout ${res}/map_61.txt -threads 15
usearch -otutab ${res}/merged.fasta -zotus ${res}/zotus_51.fasta -otutabout ${res}/zotutab_61.txt -mapout ${res}/zmap_61.txt -threads 15

# 基于数据库去嵌合体
echo "Step_7:Chimera detection"
otusfile='otus_u52.fasta'	### 选择上文中产生的OTU table
reffile='rdp_gold.fa'		### 选择一个参考数据库
printf "OTU表名：%s ,参考数据库名：%s \n"${otusfile} ${reffile}
### 基于所选择的参考数据库比对去除已知序列的嵌合体 ###
usearch -uchime2_ref ${res}/${otusfile} -db ${refdb}/${reffile} -chimeras ${res}/otus_chimeras.fasta -notmatched ${res}/otus_rdp.fasta -uchimeout ${res}/otus_rdp.uchime -strand plus -mode sensitive -threads 30

### 获取嵌合体的序列ID ###
echo "Get sequenceID of chimeric:"
grep '>' ${res}/otus_chimeras.fasta|sed 's/>//g' > ${res}/otus_chimeras.id

### 剔除嵌合体的序列 ###
/root/miniconda3/envs/python27/bin/filter_fasta.py -f ${res}/${otusfile} -o ${res}/otus_non_chimera.fasta -s ${res}/otus_chimeras.id -n
### 检查是否为预期的序列数量
echo "Check whether it is the expected number of sequences:"
grep '>' -c ${res}/otus_non_chimera.fasta
# 去除非细菌序列
echo "Step_8:Removal of non-bacterial sequences"
### 下载Greengene最新数据库，320MB
### wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
### 解压数据包后大小3.4G
### tar xvzf gg_13_8_otus.tar.gz
### 将OTU与97%相似聚类的代表性序列多序列比对，大约8min

time /root/miniconda3/envs/python27/bin/align_seqs.py -i ${res}/otus_non_chimera.fasta -t ${refdb}/gg_13_8_otus/rep_set_aligned/97_otus.fasta -o ${res}/aligned/
# 无法比对细菌的数量
echo "Number of incomparable bacteria:"
grep -c '>' ${res}/aligned/otus_non_chimera_failures.fasta 
# 获得不像细菌的OTU ID
echo "Obtaining OTU IDs that are not like bacteria:"
grep '>' ${res}/aligned/otus_non_chimera_failures.fasta|cut -f 1 -d ' '|sed 's/>//g' > ${res}/aligned/otus_non_chimera_failures.id
# 过滤非细菌序列

/root/miniconda3/envs/python27/bin/filter_fasta.py -f ${res}/otus_non_chimera.fasta -o ${res}/otus_rdp_align.fasta -s ${res}/aligned/otus_non_chimera_failures.id -n

# 查看我们现在还有多少OTU
echo "How much is left of OTUs?:"
grep '>' -c ${res}/otus_rdp_align.fasta
echo "Step_9:Generate representative sequences and OTU tables"
# 重命名OTU，这就是最终版的代表性序列，即Reference(可选)
echo "Step_9.1:Rename OTU"
awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' ${res}/otus_rdp_align.fasta > ${res}/rep_seqs.fasta
# 生成OTU表
echo "Step_9.2:Generate OTU tables"
usearch -usearch_global ${res}/merged.fasta -db ${res}/rep_seqs.fasta -otutabout ${res}/otu_table.txt -biomout ${res}/otu_table.biom -strand plus -id 0.97 -threads 30

#############
usearch -usearch_global res2/merged.fasta -db /media/ruizhi/database/16S/gg_13_8_otus/rep_set/97_otus.fasta -otutabout res2/gg138_otu_table.txt -biomout res2/gg138_table.biom -strand plus -id 0.97 -threads 20

mkdir Results/8.PICRUST_annotation/picrust1
source ~/miniconda3/bin/activate picrust1
#按拷备数标准化OTU表， -c指定数据库
normalize_by_copy_number.py -i ${res}/gg138_table.biom -o Results/8.PICRUST_annotation/picrust1/gg138_normalized_otus.biom
# 预测宏基因组
predict_metagenomes.py -i Results/8.PICRUST_annotation/picrust1/gg138_normalized_otus.biom -o Results/8.PICRUST_annotation/picrust1/metagenome_predictions.biom
# 转换结果为txt查看
biom convert -i Results/8.PICRUST_annotation/picrust1/metagenome_predictions.biom -o Results/8.PICRUST_annotation/picrust1/metagenome_predictions.txt --table-type="OTU table" --to-tsv
# 调整为标准表格
sed -i '/# Const/d;s/#OTU //g' Results/8.PICRUST_annotation/picrust1/metagenome_predictions.txt

# 按功能级别分类汇总
# -c指输出类型，有KEGG_Pathways, COG_Category，RFAM三种，-l是级别，分4级
categorize_by_function.py -i Results/8.PICRUST_annotation/picrust1/metagenome_predictions.biom -c KEGG_Pathways -l 2 -o Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L2.biom
# 转换结果为txt查看
biom convert -i Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L2.biom -o Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L2.txt --table-type="OTU table" --to-tsv

# -c指输出类型，有KEGG_Pathways, COG_Category，RFAM三种，-l是级别，分4级
categorize_by_function.py -i Results/8.PICRUST_annotation/picrust1/metagenome_predictions.biom -c KEGG_Pathways -l 3 -o Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L3.biom
# 转换结果为txt查看
biom convert -i Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L3.biom -o Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L3.txt --table-type="OTU table" --to-tsv

# -c指输出类型，有KEGG_Pathways, COG_Category，RFAM三种，-l是级别，分4级
categorize_by_function.py -i Results/8.PICRUST_annotation/picrust1/metagenome_predictions.biom -c KEGG_Pathways -l 4 -o Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L4.biom
# 转换结果为txt查看
biom convert -i Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L4.biom -o Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L4.txt --table-type="OTU table" --to-tsv

# 调整为标准表格
sed -i '/# Const/d;s/#OTU //g' Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L2.txt
sed -i '/# Const/d;s/#OTU //g' Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L3.txt
sed -i '/# Const/d;s/#OTU //g' Results/8.PICRUST_annotation/picrust1/metagenome_predictions.L4.txt

# 查找指定通路的丰度情况(可选，用于具体查看某个KO)
#metagenome_contributions.py -i result/normalized_otus.biom -l K00001,K00002,K00004 -o result/ko_metagenome_contributions.tab

source ~/miniconda3/bin/activate python27


#############

# OTU物种注释，表操作
### 物种注释
### 必须添加rdp_classifier.jar的路径到PATH中，否则可能报错
### 如果没有rdp_classifier,可从如下地址下载,下载后进行解压
### https://jaist.dl.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.2.zip
echo "Step_10:Species annotation"
#RDP_Classifier_Path=$HOME/miniconda3/share/rdp_classifier-2.2-1/rdp_classifier.jar
#echo "export RDP_JAR_PATH=${RDP_Classifier_Path}" >> $HOME/.bashrc
#echo "export RDP_JAR_PATH=$HOME/database/rdp_classifier_2.2/rdp_classifier-2.2.jar" >> $HOME/.bashrc
#source $HOME/.bashrc
/root/miniconda3/envs/python27/bin/assign_taxonomy.py -i ${res}/rep_seqs.fasta -r ${refdb}/gg_13_8_otus/rep_set/97_otus.fasta -t ${refdb}/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -m rdp -o ${res}
#echo "Align the representative sequence set"
#align_seqs.py  -i ${res}/rep_seqs.fasta -t ${refdb}/gg_13_8_otus/rep_set_aligned/core_set_aligned.fasta.imputed -o ${res}/pynast_aligned/
#echo "Filtering the alignment prior to tre building and removing positions which are all gaps, or not useful for phylogenetic inference"
#filter_alignment.py -i ${res}/pynast_aligned/rep_seqs_aligned.fasta -o ${res}/filtered_alignment/ --remove_outliers
#echo "Building a phylogenetic tree"
#make_phylogeny.py -i ${res}/filtered_alignment/rep_seqs_aligned_pfiltered.fasta -o ${res}/rep_phylo.tre
#echo "Building an OTU table"
#make_otu_table.py -i ${res}/unique_otus.txt -t ${refdb}/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -o ${res}/otu_table_non_chimeric.biom
### OTU表统计、格式转换、添加信息
echo "Step_11:OTU Table Statistics, Format Conversion, Adding Information"
### 文本OTU表转换为BIOM：方便操作
biom convert -i ${res}/otu_table.txt -o ${res}/otu_table.biom --table-type="OTU table" --to-json
### 添加物种信息至OTU表最后一列，命名为taxonomy
biom add-metadata -i ${res}/otu_table.biom --observation-metadata-fp ${res}/rep_seqs_tax_assignments.txt -o ${res}/otu_table_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
### 转换biom为txt格式，带有物种注释：人类可读
biom convert -i ${res}/otu_table_tax.biom -o ${res}/otu_table_tax.txt --to-tsv --header-key taxonomy
cp ${res}/otu_table_tax.txt Results/2.OTU\ table
cp ${res}/otu_table_tax.biom Results/2.OTU\ table

### 查看OTU表的基本信息：样品，OUT数量统计
biom summarize-table -i ${res}/otu_table_tax.biom -o ${res}/otu_table_tax.sum
#less ${res}/otu_table_tax.sum
# OTU表筛选
echo "Step_12:Screen of OTU table"
# 按样品数据量过滤：选择counts>3000的样品
/root/miniconda3/envs/python27/bin/filter_samples_from_otu_table.py -i ${res}/otu_table_tax.biom -o ${res}/otu_table2.biom -n 3000
# 查看过滤后结果：只有15个样品，1129个OTU
echo "The result of screen step1:"
biom summarize-table -i ${res}/otu_table2.biom
# 按OTU丰度过滤：选择相对丰度均值大于万分之一的OTU
/root/miniconda3/envs/python27/bin/filter_otus_from_otu_table.py --min_count_fraction 0.0001 -i ${res}/otu_table2.biom -o ${res}/otu_table3.biom
# 查看过滤后结果：只有15个样品，444个OTU
echo "The result of screen step2:"
biom summarize-table -i ${res}/otu_table3.biom
# 转换最终biom格式OTU表为文本OTU表格
biom convert -i ${res}/otu_table3.biom -o ${res}/otu_table3.txt --table-type="OTU table" --to-tsv
# OTU表格式调整方便R读取
sed -i '/# Const/d;s/#OTU //g;s/ID.//g' ${res}/otu_table3.txt
# 筛选最终OTU表中对应的OTU序列
/root/miniconda3/envs/python27/bin/filter_fasta.py -f ${res}/rep_seqs.fasta -b ${res}/otu_table3.biom -o ${res}/rep_seqs4.fasta
printf "最终的OTU表:%s\t最终的OTU序列:%s\n"${res}/otu_table_tax.txt ${res}/rep_seqs.fasta
cp ${res}/rep_seqs.fasta Results/2.OTU\ table
# End
##################lefse###########################
mkdir -p Results/7.Lefse\ analysis
summarize_taxa.py -i Results/2.OTU\ table/otu_table_tax.biom -o Results/7.Lefse\ analysis/taxonomy
python2 /media/ruizhi/software/lefse_custom/1-summarize_trans.py -i Results/7.Lefse\ analysis/taxonomy --prefix otu_table_tax
grep -v 'Other' Results/7.Lefse\ analysis/taxonomy/otu_table_tax_all.txt > Results/7.Lefse\ analysis/otu_table_tax_filt.txt
cp group.csv Results/7.Lefse\ analysis
mv Results/7.Lefse\ analysis/group.csv Results/7.Lefse\ analysis/group.txt
sed -i '1d' Results/7.Lefse\ analysis/group.txt
sed -i 's/,/\t/g' Results/7.Lefse\ analysis/group.txt
source ~/miniconda3/bin/activate base
python2 /media/ruizhi/software/lefse_custom/2-lefse_trans.py Results/7.Lefse\ analysis/otu_table_tax_filt.txt Results/7.Lefse\ analysis/group.txt Results/7.Lefse\ analysis/otu_table_for_lefse.txt
format_input.py Results/7.Lefse\ analysis/otu_table_for_lefse.txt Results/7.Lefse\ analysis/otu_table_tax.in -c 1 -u 2 -o 1000000
run_lefse.py Results/7.Lefse\ analysis/otu_table_tax.in Results/7.Lefse\ analysis/otu_table_tax.res -l 2
plot_res.py Results/7.Lefse\ analysis/otu_table_tax.res Results/7.Lefse\ analysis/lefse.lda.pdf --format pdf --dpi 150 --width 16
plot_cladogram.py Results/7.Lefse\ analysis/otu_table_tax.res Results/7.Lefse\ analysis/lefse.cladogram.pdf --format pdf --dpi 150
plot_features.py -f diff --archive zip Results/7.Lefse\ analysis/otu_table_tax.in Results/7.Lefse\ analysis/otu_table_tax.res  Results/7.Lefse\ analysis/lefse.biomarkers.zip

#################################################################################
source ~/miniconda3/bin/activate picrust2
picrust2_pipeline.py -s ${res}rep_seqs.fasta -i ${res}otu_table.txt -o Results/8.PICRUST_annotation/picrust2 -p 4
cp /media/ruizhi/software/16S_pipelines/picrust_readme/Readme.docx Results/8.PICRUST_annotation
#################################################################################
source ~/miniconda3/bin/activate python27
mkdir -p Results/3.Krona 
cp res2/otu_table_tax.txt Results/3.Krona
sed -i '1d' Results/3.Krona/otu_table_tax.txt
sed -i 's/#//g' Results/3.Krona/otu_table_tax.txt
cd Results/3.Krona/
Rscript /media/ruizhi/software/16S/melt.r
cd ../../
sed -i 's/,/\t/g' Results/3.Krona/data_melt.csv
sed -i 's/"//g' Results/3.Krona/data_melt.csv
awk -F '\t' '{a[$1"\t"$2]+=$3}END{for(i in a)print i"\t"a[i]}' Results/3.Krona/data_melt.csv > Results/3.Krona/data_melt.txt

awk -F '\t' '{a[$2]++}END{for(i in a)print i}' Results/3.Krona/data_melt.txt > Results/3.Krona/sample_id
sed -i '/^sample/d' Results/3.Krona/sample_id
mkdir -p Results/3.Krona/sample
while read line;do(awk -F "\t" -v var="$line" '$2==var{print $0}' Results/3.Krona/data_melt.txt > Results/3.Krona/sample/$(basename $line).txt);done < Results/3.Krona/sample_id
cd Results/3.Krona/sample
ls *.txt | while read id;do(awk -F '\t' '{print $3"\t"$1}' $id > $(basename $id)2);done
ls *.txt2 | while read id;do(sed -i 's/; /\t/g' $id);done
ls *.txt2 | while read id;do(sed -i 's/k__//g' $id);done
ls *.txt2 | while read id;do(sed -i 's/g__//g' $id);done
ls *.txt2 | while read id;do(sed -i 's/p__//g' $id);done
ls *.txt2 | while read id;do(sed -i 's/c__//g' $id);done
ls *.txt2 | while read id;do(sed -i 's/o__//g' $id);done
ls *.txt2 | while read id;do(sed -i 's/f__//g' $id);done
ls *.txt2 | while read id;do(sed -i 's/s__//g' $id);done
ls *.txt2 | while read id;do(sed -i '/^0/d' $id);done
ls *.txt2 | while read id;do(mv $id $(basename $id '2'));done

source ~/miniconda3/bin/activate base
ktImportText  `ls *.txt` -o krona.html
mv krona.html ../krona.result.html
cd ..
rm sample_id data_melt*










