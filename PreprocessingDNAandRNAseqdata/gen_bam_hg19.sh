#! /bin/bash
set -euxo pipefail
export LD_LIBRARY_PATH=/home/qmu/tools/lib:/opt/intel/advisor_xe_2016.1.40.463413/lib64:/home/qmu/tools/xz-5.2.3/tmp/lib 
picard=/home/share/jgwang/softwares/picard/build/libs/picard.jar

fq1=$1
fq2=$2
P=$3

REF=/home/qmu/reference/hg19bwa/hg19.fa

R1=${fq1}
R2=${fq2}

echo Mapping ${P} to human genome hg19...

bwa mem -t 24 -T 0 -M -R "@RG\tID:$P" ${REF} ${R1} ${R2} |samtools view -Shub - > ${P}.bam

samtools sort -@ 24 -m 1G ${P}.bam ${P}.sort

java -Djava.io.tmpdir=/home/qmu/tmp \
	-Xmx8g -jar ${picard} MarkDuplicates \
	INPUT=${P}.sort.bam \
	OUTPUT=${P}.sort.MD.bam \
	METRICS_FILE=${P}.sort.MD.bam.txt \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=false \
	VALIDATION_STRINGENCY=LENIENT

samtools index ${P}.sort.MD.bam
rm ${P}.bam ${P}.sort.bam #${P}.sort.MD.bam.txt
