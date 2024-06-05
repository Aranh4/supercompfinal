#!/bin/bash
#SBATCH --job-name=q3mpi
#SBATCH --output=q3mpi_output.txt
#SBATCH --ntasks=4             # Número total de processos MPI
#SBATCH --cpus-per-task=1      # Número de threads por processo MPI
#SBATCH --time=01:00:00        # Tempo máximo de execução (HH:MM:SS)
#SBATCH --partition=normal    # Nome da partição
#SBATCH --mem=500M            # Memória por processo MPI

mpirun -np 4 ./q3MPI