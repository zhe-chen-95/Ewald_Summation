#!/bin/bash
#
##SBATCH --nodes=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --mem=2GB
#SBATCH --gres=gpu:1
#SBATCH --job-name=myHPC_final
#SBATCH --mail-type=END
##SBATCH --mail-user=gl1705@nyu.edu
#SBATCH --output=slurm_%j.out

module load cuda/10.1.105 gcc/6.3.0 fftw/intel/3.3.6-pl2

make

./main 8 64 2048 24 2

