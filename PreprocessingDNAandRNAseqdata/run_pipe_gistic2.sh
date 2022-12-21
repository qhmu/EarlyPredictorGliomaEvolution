#!/bin/bash
mkdir -p I_NonHM_0.1_0.1
mkdir -p R_NonHM_0.1_0.1
mkdir -p IR_NonHM_0.1_0.1

mkdir -p I_NonHM_0.321_0.415
mkdir -p R_NonHM_0.321_0.415
mkdir -p IR_NonHM_0.321_0.415

mkdir -p log

nohup sh do_gistic2.sh ./I_seg_NonHM_20220908.seg 0.1 0.1 ./I_NonHM_0.1_0.1 > ./log/I_NonHM_0.1.out &
nohup sh do_gistic2.sh ./R_seg_NonHM_20220908.seg 0.1 0.1 ./R_NonHM_0.1_0.1 > ./log/R_NonHM_0.1.out &
nohup sh do_gistic2.sh ./IR_seg_NonHM_20220908.seg 0.1 0.1 ./IR_NonHM_0.1_0.1 > ./log/IR_NonHM_0.1.out &

nohup sh do_gistic2.sh ./I_seg_NonHM_20220908.seg 0.321 0.415 ./I_NonHM_0.321_0.415 > ./log/I_NonHM_0.321.out &
nohup sh do_gistic2.sh ./R_seg_NonHM_20220908.seg 0.321 0.415 ./R_NonHM_0.321_0.415 > ./log/R_NonHM_0.321.out &
nohup sh do_gistic2.sh ./IR_seg_NonHM_20220908.seg 0.321 0.415 ./IR_NonHM_0.321_0.415 > ./log/IR_NonHM_0.321.out &