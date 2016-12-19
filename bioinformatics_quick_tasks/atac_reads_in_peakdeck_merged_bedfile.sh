# Written by Nick De Veaux on 12/19/16

# This script computes chromsome level read counts for each output ATAC bam filtered by the peak bed file
# This allows you to find the proportion of reads that fall into peaks vs. outside of peaks
# The current hardcoded data sources are pointed to ILC ATAC data and the mouse mm10 genome 

scriptDir=/mnt/xfs1/home/carriero/projects/EmilyMiraldi/gh_git/ATACSeqSF/centos7

# The following implicitly describe dependecies of this script:
export PATH=/mnt/xfs1/bioinfoCentos7/software/installs/bedtools/bedtools-2.17.0:\
/mnt/xfs1/bioinfoCentos7/software/installs/bowtie2/bowtie2-2.2.3:\
/mnt/xfs1/bioinfoCentos7/software/installs/igvtools/IGVTools_2.3.32:\
/mnt/xfs1/bioinfoCentos7/software/installs/ngsutils/git/bin:\
/mnt/xfs1/bioinfoCentos7/software/installs/PeakSplitter/cpp_v1:\
/mnt/xfs1/bioinfoCentos7/software/installs/python/2.7.9/bin:\
/mnt/xfs1/bioinfoCentos7/software/installs/samtools/0.1.19:\
/mnt/xfs1/bioinfoCentos7/software/installs/ucsc/userApps/140911:\
/mnt/xfs1/bioinfoCentos7/software/installs/gnuParallel/20150622:\
${PATH}

BedRef='/mnt/hdfs/emiraldi/GenomeInf/mm10/genebody_4col.bed'
KgRef='/mnt/hdfs/emiraldi/GenomeInf/mm10/kgXref.txt'
RSRef='/mnt/hdfs/emiraldi/GenomeInf/mm10/mm10_GRCm38_Refseq.bed'

baseDir=/mnt/ceph/users/ndeveaux/ILC_NERD_2017
peakDir=$baseDir/peakdeck/peakdeck_75bps_10kb
for i in `ls ${peakDir} | grep sorted_and_merged`; do
  noMito=${i/_sorted_and_merged.bed/_noMito.bam};
  output=${noMito/_noMito/_noMito_inside_peaks};
  statsFile=${output/.bam/_chromStats.txt}; 
  echo $baseDir/Stats/$statsFile
  echo $noMito
  echo $baseDir/$output
  bedtools intersect -abam $baseDir/$noMito -b $peakDir/$i  > $baseDir/$output
    samtools index $baseDir/$output
    python ${scriptDir}/gcdB_NJC_EM.py \
    $baseDir/$output \
    ${BedRef} \
    ${KgRef} \
    10000 \
    $baseDir/Stats/$statsFile
 done
