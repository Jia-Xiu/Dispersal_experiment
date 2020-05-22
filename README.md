# Dispersal_experiment
The dispersal microcosom experiment was performed in March-August 2018. 
The following is the scripts related to this manuscript: *Dispersal by inundation significantly alters soil bacterial communities in early but not late successional stages*, authored by Xiu Jia, Cas Cornet and Joana Falcão Salles

## Analysis 16S rRNA sequences (DNA- and RNA-based) 
Work at [QIIME2/2019.10](https://docs.qiime2.org/2019.10/) environment, using the following script
* [script_dispersal_xiu_5-11-2019.sh](https://github.com/Jia-Xiu/dispersal_experiment_2018/blob/master/script_dispersal_xiu_5-11-2019.sh)
	 
## Community and statistic analysis 

* alpha & beta-diversity analysis:
	* [alpha_diversity_dispersal_scatterplot.R](https://github.com/Jia-Xiu/dispersal_experiment_2018/blob/master/alpha_diversity_dispersal_scatterplot.R)
	* [beta_diversity_dispersal.R](https://github.com/Jia-Xiu/dispersal_experiment_2018/blob/master/beta_diversity_dispersal.R)

* Analyze the rare biopphere:
	* [define_rare_common_biospheres.R](https://github.com/Jia-Xiu/dispersal_experiment_2018/blob/master/define_rare_common_biospheres.R)
	* [rare_biosphere_analysis_dispersal.R](https://github.com/Jia-Xiu/dispersal_experiment_2018/blob/master/rare_biosphere_analysis_dispersal.R)
	
* Metagenome prediction by [PICRUSt2](https://github.com/picrust/picrust2/wiki)
	* Using [full pipeline script](https://github.com/picrust/picrust2/wiki/Full-pipeline-script)
	```
	# load module in HPC
	module load PICRUSt2/2.3.0b
	# run script with default settings (<30min)
	picrust2_pipeline.py -s rep-seqs-cleaned-0-year.fasta -i table_0_year.biom -o picrust2_out_pipeline -p 4
	picrust2_pipeline.py -s rep-seqs-cleaned-70-year.fasta -i table_70_year.biom -o picrust2_out_pipeline -p 4
	```
	* Overall pattern [PICRUST2_PCoA_dispersal.R](https://github.com/Jia-Xiu/dispersal_experiment_2018/blob/master/PICRUST2_PCoA_dispersal.R)
	* [PICRUST2_heatmap_dispersal.R](https://github.com/Jia-Xiu/dispersal_experiment_2018/blob/master/PICRUST2_heatmap_dispersal.R)


## Author
* **Jia, Xiu** 
