#!/bin/bash

ulimit -s unlimited
source ~/.bashrc_CoLM202X_intel
make 
wait
echo "make finish"
