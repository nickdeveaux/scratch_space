outDir=/mnt/ceph/users/ndeveaux/reference/motif_cache_mm10
peaks=mouse_ATAC_bkgrnd_peaks
database=mouse_2016
localBkgFile=${finalOutputFolder}/mm10_Littman_ATAC_samples_bkgrd_Order1.txt

fastaInput=/mnt/xfs1/bioinfo/data/Illumina/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa
export PATH=/mnt/xfs1/bioinfoCentos7/software/installs/bedtools/bedtools-2.25.0/bedtools:$PATH
memehome="/mnt/xfs1/bioinfoCentos7/software/installs/meme/4.10.1/bin"
export PATH=$memehome:$PATH

function log {
    echo $(date +%F_%T) $$ ${BASHPID} $1
}

log "working with ${database} database"

finalOutputFolder=${outDir}/${peaks}/${database}
mkdir -p ${outDir}/${peaks}
mkdir -p ${finalOutputFolder}

localBkgFile=${finalOutputFolder}/mm10_Littman_ATAC_samples_bkgrd_Order1.txt

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


mv "$wd/motifs" $outDir