#!/bin/bash
set -x
while read -r fq1 fq2 prefix
do
  fastp -w 16 -i $fq1 -o ${prefix}_1.fastq.gz \
 	-I $fq2 -O ${prefix}_2.fastq.gz \
 	-q 15 -u 40 --length_required 45 \
 	--detect_adapter_for_pe \
 	-h fastp_$prefix\.html -j fastp_$prefix\.json -R report_$prefix\.txt
done <$1


