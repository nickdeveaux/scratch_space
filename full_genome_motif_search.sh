# This script finds all motif hits across the genome for either hg38 or mm10

fastaInput="/mnt/ceph/users/ndeveaux/reference/genome-versions/gencode24/GRCh38.primary_assembly.genome.fa"
outDir=/mnt/ceph/users/ndeveaux/reference/motif_cache_hg38
database=hg19_2016
peaks=human_bkgrnd_from_DC_ATAC_peaks
localBkgFile=/mnt/xfs1/home/ndeveaux/hiv-dc/data/atac_output/2016_10_23_qc_bfsw_out/peakdeck/peakdeck_75bps_10kb/all_peaks_merged_bkgrd_Order1.txt

# outDir=/mnt/ceph/users/ndeveaux/reference/motif_cache_mm10
# peaks=mouse_ATAC_bkgrnd_peaks
# database=mouse_2016
# fastaInput=/mnt/xfs1/bioinfo/data/Illumina/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa
# localBkgFile=/mnt/ceph/users/ndeveaux/reference/motif_cache_mm10/mouse_ATAC_bkgrnd_peaks/mouse_2016/mm10_Littman_ATAC_samples_bkgrd_Order1.txt

export PATH=/mnt/xfs1/bioinfoCentos7/software/installs/bedtools/bedtools-2.25.0/bedtools:$PATH
memehome="/mnt/xfs1/bioinfoCentos7/software/installs/meme/4.10.1/bin"
export PATH=$memehome:$PATH

function log {
    echo $(date +%F_%T) $$ ${BASHPID} $1
}

log "working with ${database} database"

# Archive the background file
bkgrdFolder=${outDir}/${peaks}
mkdir -p $bkgrdFolder
cp $localBkgFile $bkgrdFolder

motifDirectory="/mnt/hdfs/emiraldi/motif_databases"
bgOrder="1"

# Create output directory if it does not already exist.
mkdir -p $outDir || { echo "Cannot create \"$outDir\". Exiting" ; exit 1 ; }

# Where was I?
wwi=$(pwd)

# Create a temp directory dedicated to the processing of this sample.
wd=$(mktemp -d /tmp/${USER}_f4mBP_XXXXXXXXXX)

function cleanup {
    # Try to go back to where we started from.
    cd $wwi
    rm -rf $wd || echo "Check for \"$wd\" on $(hostname)."
}
# This runs the cleanup function when the script exits (normally or
# due to an error).
trap cleanup EXIT

# Switch to the temp directory.
cd $wd
log "Work in $(pwd)"

mkdir "$wd/motifs"

motifLines="motifLines.txt"
grep "MOTIF " ${motifDirectory}/${database}.meme | awk '{print $2}' > ${motifLines}

# GNU parallel's author requires attribution.
log "running FIMO on ${database}"
for motif in $(cat ${motifLines}); do         echo "fimo \
                --parse-genomic-coord \
                --thresh .0001 \
                --text \
                --bgfile ${localBkgFile} \
                --verbosity 1 \
                --motif ${motif} \
                ${motifDirectory}/${database}.meme \
                ${fastaInput} > $wd/motifs/${motif// /_}.txt  \
                2>  $wd/motifs/${motif// /_}.err"
done > $wd/job.txt

cores=$(~carriero/bin/nprocNoHT)

cat $wd/job.txt | xargs -I CMD -P $cores bash -c CMD

log "Finished running fimo search, outputting txt files"
# convert motif outputs to bedfiles
for f in $wd/motifs/*.txt; do
  motif=$(basename $f .txt)
  bed=$motif.bed
  awk -F'\t' 'NR>1 {print $2"\t"$3"\t"$4"\t"$1"\t"$7"\t"$6}' $f > $wd/motifs/$bed
done

log "Rsyncing results"      
rsync -av motifs ${outDir}/

log "Done"