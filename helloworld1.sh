#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/my_directory/%j.out
#SBATCH --job-name=hello
echo "Hello, World! From $HOSTNAME"
sleep 15
echo "Goodbye, World! From $HOSTNAME"
