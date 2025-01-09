#!/bin/bash

ulimit -s unlimited
source /home/zhangr23/.bashrc_CoLM202X_intel
make 
wait
echo "make finish"
