#! /bin/bash
lfq=$1
rtq=$2
otd=${3}.starfusion

STARFusion=/home/qmu/tools/STAR-Fusion-v1.5.0/STAR-Fusion
$STARFusion --examine_coding_effect --left_fq ${lfq}  --right_fq ${rtq} --CPU 24 --genome_lib_dir /home/qmu/reference/GRCh37_gencode_v19_CTAT_lib_Nov012017/ctat_genome_lib_build_dir -O ${otd}
mv ${otd}/star-fusion.fusion_predictions.abridged.coding_effect.tsv $PWD/${3}.star-fusion.fusion_predictions.abridged.coding_effect.tsv
#mv ${otd}/star-fusion.fusion_predictions.tsv $PWD/${3}.star-fusion.fusion_predictions.tsv
#rm -rf ${otd}/*
mv ${3}.star-fusion.fusion_predictions.abridged.coding_effect.tsv ${otd} 
