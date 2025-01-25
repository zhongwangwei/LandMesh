#!/bin/bash
#SBATCH -J landmeshv7
#SBATCH -p first
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=80
##SBATCH --nodelist=tms10
##SBATCH --mem=500G
#SBATCH --exclusive
##SBATCH -t 0-8:00

# 生成时间戳
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# 使用时间戳生成输出文件名
#SBATCH -o mkgrd_${TIMESTAMP}.o
#SBATCH -e mkgrd_${TIMESTAMP}.e

ulimit -s unlimited
source ~/.bashrc_CoLM202X_intel
echo "Script started at $(date)" > log.landmeshv7-2880x1440-20250122-usgs-mode4
./mkgrd.x landmeshv7.nml &>> log.landmeshv7-2880x1440-20250122-usgs-mode4
echo "Script finished at $(date)" >> log.landmeshv7-2880x1440-20250122-usgs-mode4

#echo "Script started at $(date)" > log.landmeshv7-case1-20250112
#./mkgrd.x landmeshv7-case1.nml &>> log.landmeshv7-case1-20250112
#echo "Script finished at $(date)" >> log.landmeshv7-case1-20250112
