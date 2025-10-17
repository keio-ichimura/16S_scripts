#!/bin/bash
# QIIME2 pipeline

# === 1. Import FASTQ ===
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux.qza

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# === 2. Denoise with DADA2 ===
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 20 \
  --p-trim-left-r 19 \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 210 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --p-n-threads 0

# === 3. Rename samples ===
qiime feature-table group \
  --i-table table.qza \
  --p-axis sample \
  --m-metadata-file Samplesheet/SampleSheet.txt \
  --m-metadata-column newID \
  --p-mode sum \
  --o-grouped-table table_cn.qza

# === 4. Build phylogenetic tree ===
mkdir -p phylogeny
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree phylogeny/unrooted-tree.qza \
  --o-rooted-tree phylogeny/rooted-tree.qza

# === 5. Core diversity metrics ===
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny phylogeny/rooted-tree.qza \
  --i-table table_cn.qza \
  --p-sampling-depth 18258 \
  --m-metadata-file sample-metadata_cn.txt \
  --output-dir core-metrics-results

# === 6. Export alpha/beta diversity ===
for metric in weighted_unifrac unweighted_unifrac; do
  qiime tools export \
    --input-path core-metrics-results/${metric}_pcoa_results.qza \
    --output-path core-metrics-results/
  mv core-metrics-results/ordination.txt core-metrics-results/${metric}_pcoa_results.txt

  qiime tools export \
    --input-path core-metrics-results/${metric}_distance_matrix.qza \
    --output-path core-metrics-results/
  mv core-metrics-results/distance-matrix.tsv core-metrics-results/${metric}_distance_matrix.tsv
done

for alpha in observed_features_vector chao1_index_vector shannon_vector faith_pd_vector; do
  qiime tools export \
    --input-path core-metrics-results/${alpha}.qza \
    --output-path core-metrics-results/
  mv core-metrics-results/alpha-diversity.tsv core-metrics-results/${alpha}.tsv
done

# === 7. Classify taxonomy ===
mkdir -p taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier classifier/classifier_silva99_v12.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy/rep-seqs_classified.qza

qiime metadata tabulate \
  --m-input-file taxonomy/rep-seqs_classified.qza \
  --o-visualization taxonomy/rep-seqs_classified.qzv

qiime tools export \
  --input-path taxonomy/rep-seqs_classified.qza \
  --output-path taxonomy
mv taxonomy/taxonomy.tsv taxonomy/rep-seqs_classified.tsv

# === 8. Taxa bar plot ===
qiime taxa barplot \
  --i-table table_cn.qza \
  --i-taxonomy taxonomy/rep-seqs_classified.qza \
  --m-metadata-file SampleSheet_cn.txt \
  --o-visualization taxonomy/taxa-bar-plots.qzv

# === 9. Collapse by taxonomy levels (for LEfSe etc.) ===
for level in 5 6 7; do
  qiime taxa collapse \
    --i-table table_cn.qza \
    --i-taxonomy taxonomy/rep-seqs_classified.qza \
    --p-level ${level} \
    --o-collapsed-table taxonomy/collapsed_table${level}.qza

  qiime tools export \
    --input-path taxonomy/collapsed_table${level}.qza \
    --output-path taxonomy

  mv taxonomy/feature-table.biom taxonomy/collapsed_table${level}.biom
  biom convert \
    -i taxonomy/collapsed_table${level}.biom \
    -o taxonomy/collapsed_table${level}.txt \
    --to-tsv
done