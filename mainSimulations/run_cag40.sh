#!/bin/bash
#SBATCH --job-name=cag40
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aw1907@princeton.edu
#SBATCH --output=cag40_%j.out

module purge
module load gcc/11
module load openmpi/gcc/4.1.6

cd $HOME/software/RNA_lammps/mainSimulations

srun $HOME/.local/bin/lmp_rna -in lmp_continue3_cag40.in