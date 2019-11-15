# This script is for the 16S amplicon sequence analysis (515F-806R) of dispersal experiment
# Experiment in March-August 2018
# Colleborate with Cas Cornet (Master project)
# Sequeincing on Illumina Miseq platform (250*2), Argonne
# Builted on 08-03-2019 by Jia Xiu
# Updated on 05-11-2019 by Jia Xiu

# if submit job to Peregrine HPC as .sh file, i.e. import_multiplex_seq.sh, add following header to .sh file
#!/bin/bash
#SBATCH --job-name=qiime2_dispersal
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=12GB
#SBATCH --partition=gelifes
#SBATCH -o dispersal-%j.out
#SBATCH -e dispersal-%j.error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=you-email@rug.nl

# change directory
cd $PWD


# using keemei to validate sample-metadata-dispersal-2019.tsv files


# load QIIME2
module load QIIME2/2019.10



### generate a .qza file ###
# Import paired-end sequences
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path paired-end-sequences \
  --output-path paired-end-sequences.qza



### demultiplex ###
# demultiplex the sequence reads
qiime demux emp-paired \
  --i-seqs paired-end-sequences.qza \
  --m-barcodes-file sample-metadata-dispersal-2019.tsv \
  --m-barcodes-column barcode-sequence \
  --p-no-golay-error-correction \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details demux-details.qza

# view a summary the umber of sequences were obtained per sample.
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# using command line to view .qzv on Chrome to see demutiplexed sequence counts
# qiime tools view demux.qzv



### quality control ###
# DADA2 sequence quality control
echo DADA2 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --p-chimera-method consensus \
  --o-denoising-stats stats-dada2.qza

# visualize stats
qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

# FeatureTable and FeatureData summarize
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata-dispersal-2019.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
 


### Taxonomic analysis ###

# assign taxonomy by silva database
qiime feature-classifier classify-sklearn \
  --i-classifier /data/p278113/QIIME2/silva-132-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# summerize taxonomy info
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Export taxonomy
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported



### Remove mitochondria, cloroplasti, archaea and keep sequence assigned at phyla level (D_0_ for SILVA database) ###

# filter representative sequences

# filter feature table
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include D_0__ \
  --p-exclude archaea,eukaryota,mitochondria,chloroplast \
  --o-filtered-table table-filtered.qza

# filter representative sequences
qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-include D_0__ \
  --p-exclude archaea,eukaryota,mitochondria,chloroplast \
  --o-filtered-sequences rep-seqs-filtered.qza

# summerize filtered FeatureTable and FeatureData
qiime feature-table summarize \
  --i-table table-filtered.qza \
  --o-visualization table-filtered.qzv \
  --m-sample-metadata-file sample-metadata-dispersal-2019.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-filtered.qza \
  --o-visualization rep-seqs-filtered.qzv



### generate a phylogenetic tree using the filtered represented sequences ###
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-filtered.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# or do step-by-step
# multiple sequence alignment
qiime alignment mafft \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza

# mask/filter the alignment to remove positions that are highly variable
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza

# generate a phylogenetic tree from the masked alighment by FastTree method
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza

# place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# exporting a phylogenetic tree (newick formatted file)
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path exported




### Alpha and beta diversity ###
# rarefy
# check --p-sampling-depth in the table-filtered.qzv
qiime feature-table rarefy \
    --i-table table-filtered.qza \
    --p-sampling-depth 13000 \
    --o-rarefied-table table-rarified.qza

# rarified FeatureTable and FeatureData summarize
qiime feature-table summarize \
  --i-table table-rarified.qza \
  --o-visualization table-rarified.qzv \
  --m-sample-metadata-file sample-metadata-dispersal-2019.tsv

# exporting a rarified feature table
qiime tools export \
  --input-path table-rarified.qza \
  --output-path exported/  # extracted will will do this verbosely

# convert biom to txt
biom convert -i exported/feature-table.biom -o feature-table-rarified-nontax.tsv --to-tsv


# diversity analysis:  # check --p-sampling-depth in the table-filtered.qzv
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-rarified.qza \
  --p-sampling-depth 13000 \
  --m-metadata-file sample-metadata-dispersal-2019.tsv \
  --output-dir core-metrics-results

#rarefaction curve
qiime diversity alpha-rarefaction \
  --i-table table-rarified.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 13000 \
  --m-metadata-file sample-metadata-dispersal-2019.tsv \
  --o-visualization alpha-rarefaction.qzv \
  --verbose

  
#statistics  
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata-dispersal-2019.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file sample-metadata-dispersal-2019.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata-dispersal-2019.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/unweighted-unifrac-elevation-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata-dispersal-2019.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/weighted-unifrac-elevation-significance.qzv \
  --p-pairwise

qiime taxa barplot \
  --i-table table-filtered.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata-dispersal-2019.tsv \
  --o-visualization taxa-bar-plots.qzv

qiime tools export \
  --input-path rarefied_table.qza \
  --output-path exported

biom convert -i rarefied_table_exported/feature-table.biom -o feature-table-rarefied.tsv --to-tsv

