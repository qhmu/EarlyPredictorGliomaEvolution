#!/bin/bash

set -e
source activate py27_savi
#cd $HOME/TCC
N=${1} #the 1st input is the normal controli, bam file
P=${2} #the 2nd input is the plasma, bam file
L=${3} #the 2nd input is the lung tumor, bam file
nN=${4} #the name of the normal sample to appear in th result
nP=${5} #the name of the tumor to be appear in the result
nL=${6} #the name of the tumor to be appear in the result
OUTPUT=${7} #the folder to store the results
export LD_LIBRARY_PATH=/home/qmu/tools/lib:/opt/intel/advisor_xe_2016.1.40.463413/lib64:/home/qmu/tools/xz-5.2.3/tmp/lib
#INPUT=/home/qmu/projects/pangliomaevolution/data/SMC/
#OUTPUT=/home/qmu/projects/pangliomaevolution/results/savi/${P%.*}.triplet
V=/home/share/jgwang/softwares/vcfs
S=/home/share/jgwang/softwares/SAVI
#R=/home/share/jgwang/references/WGS_Tao_hg19_broad/Homo_sapiens_assembly19.fasta
R=/home/qmu/reference/hg19.fa

for i in {1..12}; do
	mkdir -p ${OUTPUT}/${i}
	${S}/savi.py --bams ${N},${P},${L} \
		--names ${nN},${nP},${nL} \
		--ref ${R} \
		--region chr${i} -v \
		--outputdir ${OUTPUT}/${i} \
		--annvcf ${V}/cbio.fix.sort.vcf,${V}/meganormal186TCGA.fix.sort.vcf,${V}/CosmicCodingMuts.v72.May52015.jw.vcf,${V}/CosmicNonCodingVariants.v72.May52015vcf,${V}/All_20150605.vcf,${V}/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf 1> ${OUTPUT}/${i}/err.out 2>${OUTPUT}/${i}/err.log & 
done
wait
for i in {13..22} {X,Y}; do
	mkdir -p ${OUTPUT}/${i}
	${S}/savi.py --bams ${N},${P},${L} \
		--names ${nN},${nP},${nL} \
		--ref ${R} \
		--region chr${i} -v \
		--outputdir ${OUTPUT}/${i} \
		--annvcf ${V}/cbio.fix.sort.vcf,${V}/meganormal186TCGA.fix.sort.vcf,${V}/CosmicCodingMuts.v72.May52015.jw.vcf,${V}/CosmicNonCodingVariants.v72.May52015vcf,${V}/All_20150605.vcf,${V}/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf 1> ${OUTPUT}/${i}/err.out 2>${OUTPUT}/${i}/err.log & 
done
wait

echo SAVI done!

echo Start merging...

O=${OUTPUT}/${nN}.report.coding.PDfilter.txt

head -n 1 ${OUTPUT}/1/report/report.coding.PDfilter.txt > ${O}

for i in {1..22} {X,Y}; do
        F=${OUTPUT}/${i}/report/report.coding.PDfilter.txt
        if [ -e ${F} ]; then
                tail -n +2 ${F} >> ${O}
        else
                echo ${F} does not exist
        fi
done

source deactivate
