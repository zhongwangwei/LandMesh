#!/bin/bash
#SBATCH -J make
#SBATCH -p first
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
# #SBATCH --mem=100G
#SBATCH -o make.o%j
#SBATCH -e make.e%j
# #SBATCH --exclusive
#SBATCH -t 0-1:00
ulimit -s unlimited 

source /home/zhangr23/.bashrc_CoLM202X_intel
make

