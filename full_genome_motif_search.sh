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
fimo \
                --parse-genomic-coord \
                --thresh .0001 \
                —-text \
                --bgfile ${localBkgFile} \
                --verbosity 1 \
                ${motifDirectory}/${database}.meme \
                ${fastaInput} > $wd/output.txt  \
                2>  $wd/log.err


