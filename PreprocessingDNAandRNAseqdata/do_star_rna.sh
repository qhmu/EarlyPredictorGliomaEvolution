#! /bin/sh

echo [Directory] `pwd`
echo [Machine] `uname -n`
echo [Start] `date`
echo [args] $*
time1=$( date "+%s" )

export LD_LIBRARY_PATH=/home/qmu/tools/lib:/opt/intel/advisor_xe_2016.1.40.463413/lib64:/home/qmu/tools/xz-5.2.3/tmp/lib

lfq=$1
rfq=$2
pid=$3

echo Processing $pid
mkdir -p ${pid}.star
ulimit -n 4096
STAR  	--runThreadN 24   --genomeDir /home/share/jgwang/references/hg19_STAR \
	--outFileNamePrefix ${pid}.star/${pid} \
	--readFilesIn ${lfq} ${rfq} \
	--readFilesCommand zcat      --limitBAMsortRAM 0   \
	--outSAMtype BAM   SortedByCoordinate      --outSAMstrandField intronMotif   \
	--outSAMattributes NH   HI   NM   MD   AS   XS      --outSAMunmapped Within  \
	--outSAMheaderHD @HD   VN:1.4      \
	--outFilterMultimapNmax 20   --outFilterMultimapScoreRange 1   \
	--outFilterScoreMinOverLread 0.33   --outFilterMatchNminOverLread 0.33   \
	--outFilterMismatchNmax 10   --alignIntronMax 500000   \
	--alignMatesGapMax 1000000   --alignSJDBoverhangMin 1   --sjdbOverhang 100   --sjdbScore 2

samtools index ${pid}.star/${pid}Aligned.sortedByCoord.out.bam

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))
echo [End]
