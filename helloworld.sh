#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/%j.out
echo "Hello, World! From $HOSTNAME"
sleep 15
echo "Goodbye, World! From $HOSTNAME"
