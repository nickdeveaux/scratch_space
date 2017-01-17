baseChipDir=/mnt/ceph/users/ndeveaux/atac_chipseq_concordance_test/chipseq_output
baseGenomeInfDir=/mnt/hdfs/emiraldi/GenomeInf/mm10
outDir=/mnt/ceph/users/ndeveaux/chipseq_mm10_th17_prior_gen_2017_01
declare -a macsDirs=("2016_11_28_MACS_out" "2016_11_15_MACS_out") 
declare -a tssDistances=("5000" "10000")
export PATH=/mnt/xfs1/bioinfoCentos7/software/installs/bedtools/bedtools-2.17.0:${PATH}

function log {
    echo $(date +%F_%T) $$ ${BASHPID} $1
}

# Where was I?
wwi=$(pwd)

# Create a temp directory dedicated to the processing of this sample.
wd=$(mktemp -d /tmp/${USER}_chipseq_prior_XXXXXXXXXX)

function cleanup {
    # Try to go back to where we started from.
    cd ${wwi}
    rm -rf ${wd} || echo "Check for \"${wd}\" on $(hostname)."
}
# This runs the cleanup function when the script exits (normally or
# due to an error).
trap cleanup EXIT

# Switch to the temp directory.
cd ${wd}
log "Working in $(pwd) on $(hostname):"

mkdir -p $outDir
mkdir outBound
mkdir references

for tssDistance in ${tssDistances[@]}; do
    bedtools slop -b $tssDistance -i $baseGenomeInfDir/tss_mm10.bed -g $baseGenomeInfDir/mm10.chrom.sizes > references/${tssDistance}_tss_mm10.bed
done
cp $baseGenomeInfDir/genebodies_mm10_bp10000.bed references

for ref in `ls references`; do
    ref_dir=${ref/.bed/}
    mkdir -p outBound/$ref_dir
    for macDir in ${macsDirs[@]}; do
        mkdir -p outBound/$ref_dir/$macDir
        for sample in `ls $baseChipDir/$macDir`; do
            log "references/$ref_dir $sample $macDir"
            bedtools intersect -a references/$ref -b $baseChipDir/$macDir/$sample/${sample}_peaks.bed -u > outBound/$ref_dir/$macDir/${sample}_intersection.bed
        done
    done
done

rsync -av outBound $outDir


