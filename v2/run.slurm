#!/bin/bash
#SBATCH -J mesh-teach 
#SBATCH -p second
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
source /home/zhangr23/.bashrc_CoLM202X_intel
#echo "Script started at $(date)" > log.landmeshv4-nxp36
#./mkgrd.x nxp36.nml &>> log.landmeshv4-nxp36
#echo "Script finished at $(date)" >> log.landmeshv4-nxp36
#wait

#echo "Script started at $(date)" > log.landmeshv4-nxp72
#./mkgrd.x nxp72.nml &>> log.landmeshv4-nxp72
#echo "Script finished at $(date)" >> log.landmeshv4-nxp72
#wait
 
#echo "Script started at $(date)" > log.landmeshv4-nxp144
#./mkgrd.x nxp144.nml &>> log.landmeshv4-nxp144
#echo "Script finished at $(date)" >> log.landmeshv4-nxp144
#wait

