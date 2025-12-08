#!/bin/bash
#SBATCH --job-name=ucug30
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=23:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aw1907@princeton.edu
#SBATCH --output=ucug30_%j.out

module purge
module load gcc/11
module load openmpi/gcc/4.1.6

cd $HOME/software/RNA_lammps/mainSimulations

srun $HOME/.local/bin/lmp_rna -in lmp_continue1_ucug30.in