#!/bin/bash
#SBATCH --job-name=vrp_mpi
#SBATCH --output=vrp_mpi_output.txt
#SBATCH --ntasks=4             # Número total de processos MPI
#SBATCH --cpus-per-task=1      # Número de threads por processo MPI
#SBATCH --time=01:00:00        # Tempo máximo de execução (HH:MM:SS)
#SBATCH --partition=compute    # Nome da partição

mpirun -np 4 ./q3MPI