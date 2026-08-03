#!/bin/bash
#BSUB -J HiCpro_s2_
#BSUB -e HiCpro_s2_.%J.e
#BSUB -o HiCpro_s2_.%J.o
#BSUB -q ser 
#BSUB -n 12
#BSUB -R "span[ptile=12]"
make --file /data/med-wulf/software/HiC-Pro/build/HiC-Pro-master/scripts/Makefile CONFIG_FILE=/scratch/2024-05-06/med-wulf/Hi-C-Result/HiC/config-hicpro.txt CONFIG_SYS=/data/med-wulf/software/HiC-Pro/build/HiC-Pro-master//config-system.txt all_persample 2>&1
