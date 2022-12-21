#! /bin/bash
REFGTF=/home/share/jgwang/references/Homo_sapiens.GRCh37.75.gtf
while read bam
do
        sid=${bam%%.*}
	echo $sid
        featureCounts -T 24 -p -s 0 -B -t exon -g gene_name -a $REFGTF -o ${sid}.count.fc.txt ${bam}
done <$1
