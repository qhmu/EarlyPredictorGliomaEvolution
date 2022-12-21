#!/bin/bash

source activate episome 

normal=$1 #the normal sample, bam file
tumor=$2 #the tumor, bam file
pid=$3 #the patient id, used as prefix 
outdir=$4 #output directory
mkdir -p ${outdir}

/home/mhluiaa/Genomics/cnvkit/kit/cnvkit.py batch ${tumor} \
        --normal ${normal} \
        -p 24 \
        --targets /home/mhluiaa/Genomics/cnvkit/ref/S04380110_Regions.target.bed \
        --antitargets /home/mhluiaa/Genomics/cnvkit/ref/S04380110_Regions.antitarget.bed \
        --annotate /home/mhluiaa/Genomics/cnvkit/ref/refFlat.txt \
        --fasta /home/qmu/reference/hg19.fa \
        --output-dir ${outdir}/${pid}

        #--access /home/mhluiaa/Genomics/cnvkit/data/access-5kb-mappable.hg19.bed \
        # --output-reference ${outdir}/${pid}_reference.cnn \

/home/mhluiaa/Genomics/cnvkit/kit/cnvkit.py call ${outdir}/${pid}/*.cns -o ${outdir}/${pid}/${pid}.call.cns -y
/home/mhluiaa/Genomics/cnvkit/kit/cnvkit.py export seg ${outdir}/${pid}/*.call.cns -o ${outdir}/${pid}/${pid}.call.seg


echo "Run complete. Files can be found at ${outdir}/${pid}"
source deactivate
