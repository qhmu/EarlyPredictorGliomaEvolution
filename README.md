# EarlyPredictorGliomaEvolution

## Overview

This is the custom code for the analysis presented in the manuscript "**Tracing early predictors of glioma evolution under therapy**", which is currently under review. In this study we aimed to identify early predictors of glioma evolution by learning from longitudinally sampled glioma pairs.


## Analysis contributors
[Wang Lab at HKUST](http://wang-lab.ust.hk/)


## Structure of the scripts and folders

Each folder correspond to one part of the analysis. For example, *PreprocessingDNAandRNAseqdata* folder contains the scripts to preprocess the DNA and RNA sequencing data, which include trimming the raw reads, read mapping, somatic mutation detection, gene fusion detection, copy number alteration detection, etc.

## Software and package versions

FastQC v0.11.5, fastp v0.20.1, BWA v0.7.15-r1140, samtools 1.2, picard MarkDuplicates 2.9.2, CNVkit 0.9.546,
SAVI 2.0, CELLO 1.0, STAR 2.6.1d, featureCounts 1.5.1, STAR-Fusion 1.5.0, R 3.6.3, ComplexHeatmap 1.2.0, 
survival 2.44-1.1, ggplot2 3.2.1, survminer 0.4.6, ape 5.4-1, xgboost1.2.0.1, SHAPforxgboost 0.0.4,  DEseq 1.26.0, fgsea 1.12.0


## Additional resources

The processed DNA and RNA sequencing data of the longitudinal glioma pairs have been uploaded to the CELLO2 website (www.wang-lab-hkust.com:3838/cello2 ) for visualization and exploration.

## Contact
For technical issues please send an email to qmu@connect.ust.hk or jgwang@ust.hk.
