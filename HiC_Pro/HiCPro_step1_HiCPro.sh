#!/bin/bash
#BSUB -J HiCpro_s1_HiCPro[1-10]
#BSUB -e HiCpro_s1_HiCPro.%J.e
#BSUB -o HiCpro_s1_HiCPro.%J.o
#BSUB -q ser 
#BSUB -n 12
#BSUB -J HiCpro_s1_HiCPro[1-10]
#BSUB -R "span[ptile=12]"

FASTQFILE=inputfiles_HiCPro.txt; export FASTQFILE
make --file /data/med-wulf/software/HiC-Pro/build/HiC-Pro-master/scripts/Makefile CONFIG_FILE=/scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/config-hicpro.txt CONFIG_SYS=/data/med-wulf/software/HiC-Pro/build/HiC-Pro-master//config-system.txt all_sub 2>&1
