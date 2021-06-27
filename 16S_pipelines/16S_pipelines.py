configfile:  "16S_pipelines.yaml"
#__date__ = "2021-6-27"
#__author__ = "Junbo Yang"
#__email__ = "yang_junbo_hi@126.com"
#__license__ = "yang_junbo_hi@126.com"

#Snakefile

#from snakemake.utils import makedirs
#from snakemake.utils import listfiles

#import numpy as np
import os
### 
################################################################################
# prepare
#
# usearch11、vsearch和qiime
# download Greengene database，320MB
# wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.t
# tar xvzf gg_13_8_otus.tar.gz
################################################################################


sample = config["samples"]
refdb = config["refdb"]
group = config["group"]

rule all:
	## LOCAL ##
#	'''
#	Defines the target rule by specifying all final output files.
#	Additionally, the cluster_logs dir is deleted if
#	snakemake was run locally.
#	'''
	input: 
		config["results_dir"] + "/rep_seqs4.fasta"
#	output:
#		expand(config["results_dir"] + "/{sample}.test.uniq.3pend.bed",sample=sample)


#-------------------------------------------------------------------------------
# Before: Dependency packages - None
#-------------------------------------------------------------------------------
#rule Prepair:
#	output:
#		config["results_dir"] + "/aligned/otus_non_chimera_failures.fasta"
#	params:
#		config["results_dir"] + "/aligned/"
#	shell:
#		"""
#		mkdir -p {params} | touch {output}
#		"""
#-------------------------------------------------------------------------------
# gunzip rule 1: Dependency packages - None
#-------------------------------------------------------------------------------
#rule gunzip:
#	input:
#		gz_raw = expand(config["results_dir"] + "/{sample}_1.fq.gz",sample=sample)	
#	output:
#		fq_raw = expand(config["results_dir"] + "/{sample}_1.fq",sample=sample)
#	message:
#		"step1：Decompress paired reads sequences"		
#	shell:
#		"""
#		gunzip {input.gz_raw} >{output.fq_raw}
#		"""
#-------------------------------------------------------------------------------
# Mergefastq rule 2: Dependency packages - usearch
#-------------------------------------------------------------------------------
rule Merge_fastq:
	input:  
		raw_R1 = config["input_dir"] + "/{sample}_1.fq",
		raw_R2 = config["input_dir"] + "/{sample}_2.fq"
	output: 
		fq = config["results_dir"] + "/{sample}.merge.fq",
		fa = config["results_dir"] + "/{sample}.merge.fa"
	message: 
		"Step2: Merge paired reads..."
	shell:
		'''
		usearch -fastq_mergepairs {input.raw_R1} -reverse {input.raw_R2} -relabel @ -fastqout {output.fq} -fastaout {output.fa} 
		'''
#-------------------------------------------------------------------------------
# Quality_filter rule 3: Dependency package - None
#-------------------------------------------------------------------------------
rule Quality_filter:
	input:
		config["results_dir"] + "/{sample}.merge.fq"
	output:
		config["results_dir"] + "/{sample}.filtered.fasta"
	message:
		"Step3: Start quality filtering..."
	shell:
		'''
		usearch -fastq_filter {input} -fastq_maxee 1.0 -fastaout {output}  -threads 2
		'''
#-------------------------------------------------------------------------------
# mergefasta rule 4: Dependency packages - None
#-------------------------------------------------------------------------------
rule Mergefasta:
	input:
		expand(config["results_dir"] + "/{sample}.filtered.fasta",sample = sample)
	output:
		config["results_dir"] + "/merge.fasta"
	message:
		"Step4: Start merging all sample..."
	shell:
		"cat {input} > {output}"

#-------------------------------------------------------------------------------
# Dereplication rule 5: Dependency packages - None
#-------------------------------------------------------------------------------
rule Dereplication:
	input:
		config["results_dir"] + "/merge.fasta"
	output:
		config["results_dir"] + "/unique.fasta"
	message:
		"Step5: Start dereplication..."
	shell:
		"usearch -fastx_uniques {input} -fastaout {output} --minuniquesize 2 -threads 2 -sizeout -relabel Uniq"

#-------------------------------------------------------------------------------
# Discard_singletons rule 6: Dependency packages - None
#-------------------------------------------------------------------------------
rule Discard_singletons:
	input:
		config["results_dir"] + "/unique.fasta"
	output:
		config["results_dir"] + "/discard.fasta"
	message:
		"Step6: Start discard singletons..."
	shell:
		"""
		usearch -sortbysize {input} -fastaout {output} -minsize 2
		"""
#-------------------------------------------------------------------------------
# OTU_cluster rule 7: Dependency packages - None
#-------------------------------------------------------------------------------
rule OTU_cluster1:
	input:
		config["results_dir"] + "/unique.fasta"
	output:
		config["results_dir"] + "/otus_unique.fasta",
		config["results_dir"] + "/otus_unique.uparse.txt"
	message:
		"Step7: Start UPARSE-OTU by unique fasta..."
	shell:
		"""
		usearch -cluster_otus {input} -otus {output[0]} -uparseout {output[1]}  -relabel Otu -threads 10
		"""
#-------------------------------------------------------------------------------
# OTU_cluster2 rule 8: Dependency packages - None
#-------------------------------------------------------------------------------
rule OTU_cluster2:
	input:
		config["results_dir"] + "/discard.fasta"
	output:
		config["results_dir"] + "/otus_discard.fasta",
		config["results_dir"] + "/otus_discard.uparse.txt"
	message:
		"Step8: Start UPARSE-OTU by non-singletons fasta..."
	shell:
		"""
		usearch -cluster_otus {input} -otus {output[0]} -uparseout {output[1]}  -relabel Otu -threads 10
		"""
#-------------------------------------------------------------------------------
# OTU_cluster3 rule 9: Dependency packages - None
#-------------------------------------------------------------------------------
rule OTU_cluster3:
	input:
		config["results_dir"] + "/unique.fasta"
	output:
		 config["results_dir"] + "/otus_denoise.fasta"
	message:
		"step9: Start denoise..."
	shell:
		"""
		usearch -unoise3 {input} -zotus {output} -threads 30
		"""
#-------------------------------------------------------------------------------
# Picking OTUs rule 10: Dependency packages python27--QIIME
#-------------------------------------------------------------------------------
rule Picking_OTUs:
	input:
		fa = config["results_dir"] + "/unique.fasta",
		refdb = config["refdb"] + '/gg_13_8_otus/rep_set/97_otus.fasta'
	output:
		config["results_dir"] + "/Picking_OTUs/unique_otus.txt"
	params:
		config["results_dir"] + "/Picking_OTUs"
	message:
		"Step10: Start picking a representative sequence set, one sequence from each OTU..."
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/pick_otus.py -i {input.fa} -o {params} -m uclust -r {input.refdb} --threads 30
		"""
#-------------------------------------------------------------------------------
# Picking OTUs rule 11: Dependency packages - None
#-------------------------------------------------------------------------------

rule Picking_OTUs2:
	input:
		config["results_dir"] + "/Picking_OTUs/unique_otus.txt",
		config["results_dir"] + "/unique.fasta"
	output:
		config["results_dir"] + "/rep_set1.fasta"
	message:
		"Step11: Continue picking...."
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/pick_rep_set.py -i {input[0]} -f {input[1]} -o {output}
		"""
#-------------------------------------------------------------------------------
# OTU number rule12: Dependency packages - None
#-------------------------------------------------------------------------------
rule OTU_number: 
	input:
		unique_fa = config["results_dir"] + "/otus_unique.fasta",
		discard_fa = config["results_dir"] + "/otus_discard.fasta",
		denoise_fa = config["results_dir"] + "/otus_denoise.fasta"
	output:
		config["results_dir"] + "/OTU_number.xls"
	message:
		"Step12: Echo the number of OTUs..."
	shell:
		"""
		echo "Number of OTU:" > {output}
		grep ">" -c {input.unique_fa} >> {output}
		grep ">" -c {input.discard_fa} >> {output}
		grep ">" -c {input.denoise_fa} >> {output}
		"""
#-------------------------------------------------------------------------------
# OTU table rule 13: Dependency packages - None
#-------------------------------------------------------------------------------

rule OTU_table1:
	input: 
		fa = config["results_dir"] + "/merged.fasta",
		otus = config["results_dir"] + "/otus_unique.fasta"
	output:
		otutable = config["results_dir"] + "/otutab_unique.txt",
		mapout = config["results_dir"] + "/map_unique.txt"
	message:
		"Step13: Making OTU tables for 97% OTUs..."
	shell:
		"""
		usearch -otutab {input.fa} -otus {input.otus} -otutabout {output.otutable} -mapout {output.mapout}
		"""
#-------------------------------------------------------------------------------
# OTU table rule 14: Dependency packages - None
#-------------------------------------------------------------------------------

rule OTU_table2:
	input:
		fa = config["results_dir"] + "/merged.fasta",
		otus = config["results_dir"] + "/otus_denoise.fasta"
	output:
		otutable = config["results_dir"] + "/otutab_denoise.txt",
		mapout = config["results_dir"] + "/map_denoise.txt"
	message:
		"Step14: Making OTU tables for 97% denoised sequences (ZOTUs)..."
	shell:
		"""
		usearch -otutab {input.fa} -otus {input.otus} -otutabout {output.otutable} -mapout {output.mapout}
		"""
#-------------------------------------------------------------------------------
# Chimera_detection rule15: Dependency packages - None
#-------------------------------------------------------------------------------

rule Chimera_detection:
	input:
		discard_fa = config["results_dir"] + "/otus_discard.fasta",
		refdb = config["refdb"] + "/rdp_gold.fa"
	output:
		Chimera = config["results_dir"] + "/otus_discard_chimeras.fasta",
		fa = config["results_dir"] + "/otus_discard_rdp.fasta",
		uchimeout = config["results_dir"] + "/otus_discard_rdp.uchime"
	message:
		"Step15: Detecting Chimera..."
	shell:
		"""
		usearch -uchime2_ref {input.discard_fa} -db {input.refdb} -chimeras {output.Chimera} -notmatched {output.fa} -uchimeout {output.uchimeout}  -strand plus -mode sensitive -threads 30
		"""
#-------------------------------------------------------------------------------
# Chimera_ID rule 16: Dependency packages - None
#-------------------------------------------------------------------------------
rule Chimera_ID:
	input:
		config["results_dir"] + "/otus_discard_chimeras.fasta"
	output:
		config["results_dir"] + "/otus_discard_chimeras.id"
	message:
		"Step16: Extract chimera_ID..."
	shell:
		"""
		grep ">" {input} | sed 's/>//g' > {output}
		"""
#-------------------------------------------------------------------------------
# Remove chimera rule 17: Dependency packages - None
#-------------------------------------------------------------------------------
rule Remove_chimera:
	input:
		config["results_dir"] + "/otus_discard.fasta",
		config["results_dir"] + "/otus_discard_chimeras.id"
	output:
		config["results_dir"] + "/otus_discard_non_chimera.fasta"
	message:
		"Step17: Removing chimera..."
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/filter_fasta.py -f {input[0]} -o {output[0]} -s {input[1]} -n
		"""
#-------------------------------------------------------------------------------
# Remove non-bacterial rule18: Dependency packages - None
#-------------------------------------------------------------------------------
rule Remove_non_bacterial:
	input:
		nc = config["results_dir"] + "/otus_discard_non_chimera.fasta",
		refdb = config["refdb"] + "/gg_13_8_otus/rep_set_aligned/97_otus.fasta"
	output:
		config["results_dir"] + "/aligned/otus_discard_non_chimera_failures.fasta"
	params:
		config["results_dir"] + "/aligned/"
	message:
		"Step18: Removing non-bacterial..."
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/align_seqs.py -i {input.nc} -t {input.refdb} -o {params}
		"""
#-------------------------------------------------------------------------------
# Obtaining_non_bac_id rule 19: Dependency packages - None
#-------------------------------------------------------------------------------
rule Obtain_non_bac_id:
	input:
		config["results_dir"] + "/aligned/otus_discard_non_chimera_failures.fasta"
	output:
		config["results_dir"] + "/aligned/otus_discard_non_chimera_failures.id"
	message:
		"Step19: Obtaining OTU IDs that are not like bacteria..."
	shell:
		"""
		python scripts/extract_non_chimera_failures_id.py -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# Remove_non_bac_id rule 20: Dependency packages - None
#-------------------------------------------------------------------------------
rule Remove_non_bac_id:
	input:
		config["results_dir"] + "/otus_discard_non_chimera.fasta",
		config["results_dir"] + "/aligned/otus_discard_non_chimera_failures.id"
	output:
		config["results_dir"] + "/otus_discard_rdp_align.fasta"
	message:
		"Step20: Removing OTU IDs that are not like bacteria..."
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/filter_fasta.py -f {input[0]} -o {output} -s {input[1]} -n
		"""
#-------------------------------------------------------------------------------
# Generate_rep_seq rule 21: Dependency packages - None
#-------------------------------------------------------------------------------
rule Generate_rep_seq:
	input:
		config["results_dir"] + "/otus_discard_rdp_align.fasta"
	output:
		config["results_dir"] + "/rep_seqs.fasta"
	message:
		"Step21: Generating representative sequences..."
	shell:
		"""
		python scripts/Generat_rep_seq.py -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# Generate_otu_tab rule 22: Dependency packages - None
#-------------------------------------------------------------------------------
rule Generate_otu_tab:
	input:
		fa = config["results_dir"] + "/merge.fasta",
		refdb = config["results_dir"] + "/rep_seqs.fasta"
	output:
		otu_tab = config["results_dir"] + "/otu_table.txt",
		otu_biom = config["results_dir"] + "/otu_table.biom"
	message:
		"Step22: Generating otu table by representative sequences..."
	shell:
		"""
		usearch -usearch_global {input.fa} -db {input.refdb} -otutabout {output.otu_tab} -biomout {output.otu_biom} -strand plus -id 0.97 -threads 30
		"""
#-------------------------------------------------------------------------------
# Generate_otu_tab_for_metagenome_annotation rule 23: Dependency packages - None
#-------------------------------------------------------------------------------
rule Generate_otu_tab2:
	input:
		fa = config["results_dir"] + "/merged.fasta",
		refdb = config["refdb"] + "/gg_13_8_otus/rep_set/97_otus.fasta"
	output:
		otu_tab = config["results_dir"] + "/gg138_otu_table.txt",
		otu_biom = config["results_dir"] + "/gg138_otu_table.biom"
	message:
		"Step23: Generating otu table by gg_13_8_otus database..."
	shell:
		"""
		usearch -usearch_global {input.fa} -db {input.refdb} -otutabout {output.otu_tab} -biomout {output.biom} -strand plus -id 0.97 -threads 20
		"""
#-------------------------------------------------------------------------------
# Normolise_otu_tab rule 24: Dependency packages - None
#-------------------------------------------------------------------------------
rule Normolise_otu_tab:
	input:
		config["results_dir"] + "/gg138_otu_table.biom"
	output:
		config["results_dir"] + "/gg138_normalized_otus.biom"
	message:
		"Step24: normalise otu table..."
	shell:
		"""
		normalize_by_copy_number.py -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# Predict_metagenomes rule 25: Dependency packages - None
#-------------------------------------------------------------------------------
rule Predict_metagenomes:
	input:
		config["results_dir"] + "/gg138_normalized_otus.biom"
	output:
		config["results_dir"] + "/metagenome_predictions.biom"
	message:
		"Step25: Predict metagenomes..."
	shell:
		"""
		predict_metagenomes.py -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# biom_convert rule 26: Dependency packages - None
#-------------------------------------------------------------------------------
rule Biom_convert:
	input:
		config["results_dir"] + "/metagenome_predictions.biom"
	output:
		config["results_dir"] + "/metagenome_predictions.txt"
	message:
		"Step26: Biom convert..."
	shell:
		"""
		biom convert -i {input} -o {output} --table-type="OTU table" --to-tsv
		"""
#-------------------------------------------------------------------------------
# format_predict_results rule 27: Dependency packages - None
#-------------------------------------------------------------------------------
rule Format_predict_results:
	input:
		config["results_dir"] + "/metagenome_predictions.txt"
	output:
		config["results_dir"] + "/metagenome_predictions.format.txt"
	message:
		"Step27: Format predict results..."
	shell:
		"""
		sed '/# Const/d;s/#OTU //g' {input} > {output}
		"""
#-------------------------------------------------------------------------------
# taxonomy annotation rule 28: Dependency packages - None
#-------------------------------------------------------------------------------
rule Taxonomy_annotation:
	input:
		config["results_dir"] + "/rep_seqs.fasta",
		config["refdb"] + "/gg_13_8_otus/rep_set/97_otus.fasta",
		config["refdb"] + "/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt"
	output:
		config["results_dir"] + "/Species_annotation"
	message:
		"Step28: Taxonomy annotation..."
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/assign_taxonomy.py -i {input[0]} -r {input[1]} -t {input[2]} -m rdp -o {output}
		"""
#-------------------------------------------------------------------------------
# otu_table_tax rule 29: Dependency packages - None
#-------------------------------------------------------------------------------
#rule Otu_table_tax:
#	input:
#		#config["results_dir"] + "/otu_table.txt"
#		config["results_dir"] + "/otu_table_pre.biom"
#	output:
#		config["results_dir"] + "/otu_table.biom"
#	message:
#		"Step29: generating OTU table and taxonomy annotation..."
#	shell:
#		""" cp {input} {output}"""
		#"""
		#biom convert -i {input} -o {output} --table-type="OTU table" --to-json
		#"""
#-------------------------------------------------------------------------------
# otu_table_tax2 rule 30: Dependency packages - None
#-------------------------------------------------------------------------------
rule Otu_table_tax2:
	input:
		config["results_dir"] + "/otu_table.biom"
	output:
		config["results_dir"] + "/rep_seqs_tax_assignments.txt",
		config["results_dir"] + "/otu_table_tax.biom"
	message:
		"Step30: generating OTU table and taxonomy annotation2..."
	shell:
		"""
		biom add-metadata -i {input} --observation-metadata-fp {output[0]} -o {output[1]} --sc-separated taxonomy --observation-header OTUID,taxonomy
		"""
#-------------------------------------------------------------------------------
# otu_table_tax3 rule 31: Dependency packages - None
#-------------------------------------------------------------------------------
rule Otu_table_tax3:
	input:
		config["results_dir"] + "/otu_table_tax.biom"
	output:
		config["results_dir"] + "/otu_table_tax.txt"
	shell:
		"""
		biom convert -i {input} -o {output} --to-tsv --header-key taxonomy
		"""
#-------------------------------------------------------------------------------
# Generate_otu_tab rule 32: Dependency packages - None
#-------------------------------------------------------------------------------
rule Generate_otu_sum:
	input:
		config["results_dir"] + "/otu_table_tax.biom"
	output:
		config["results_dir"] + "/otu_table_tax.sum"
	shell:
		"""
		biom summarize-table -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# Filter_OTU rule 33: Dependency packages - None
#-------------------------------------------------------------------------------
rule Filter_OTU1:
	input:
		config["results_dir"] + "/otu_table_tax.biom"
	output:
		config["results_dir"] + "/otu_table_filter1.biom"
	message:
		"Screen of OTU table: counts > 3000"
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/filter_samples_from_otu_table.py -i {input} -o {output} -n 3000
		"""
#-------------------------------------------------------------------------------
# Filter_OTU rule 34: Dependency packages - None
#-------------------------------------------------------------------------------
rule Filter_OTU2:
	input:
		config["results_dir"] + "/otu_table_filter1.biom"
	output:
		config["results_dir"] + "/otu_table_filter2.biom"
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/filter_otus_from_otu_table.py --min_count_fraction 0.0001 -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# Filter_OTU_sum rule 35: Dependency packages - None
#-------------------------------------------------------------------------------
rule Filter_OTU_sum:
	input:
		config["results_dir"] + "/otu_table_filter2.biom"
	output:
		config["results_dir"] + "/otu_table_filter2.sum"
	shell:
		"""
		biom summarize-table -i {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# Final_OTU_table rule 36: Dependency packages - None
#-------------------------------------------------------------------------------
rule Pre_OTU_table:
	input:
		config["results_dir"] + "/otu_table_filter2.biom"
	output:
		config["results_dir"] + "/pre_otu_table.txt"
	shell:
		"""
		biom convert -i ${res}/otu_table3.biom -o ${res}/otu_table3.txt --table-type="OTU table" --to-tsv
		"""
#-------------------------------------------------------------------------------
# Final_OTU_table rule 37: Dependency packages - None
#-------------------------------------------------------------------------------
rule Final_OTU_table:
	input:
		config["results_dir"] + "/pre_otu_table.txt"
	output:
		config["results_dir"] + "/final_otu_table.txt"
	shell:
		"""
		sed '/# Const/d;s/#OTU //g;s/ID.//g' {input} > {output}
		"""
#-------------------------------------------------------------------------------
# filter_Final_OTU_table rule 38: Dependency packages - None
#-------------------------------------------------------------------------------
rule Filter_final_OTU_table:
	input:
		config["results_dir"] + "/rep_seqs.fasta",
		config["results_dir"] + "/otu_table_filter2.biom"
	output:
		config["results_dir"] + "/rep_seqs4.fasta"
	shell:
		"""
		/root/miniconda3/envs/python27/bin/python /root/miniconda3/envs/python27/bin/filter_fasta.py -f {input[0]} -b {input[1]} -o {output}
		"""




