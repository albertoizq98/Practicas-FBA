#!/bin/bash
#
#SBATCH -p hpc-bio-ochoa
#SBATCH --chdir=/home/alumno07/Practicas_FBA/
#SBATCH -J TEST
#SBATCH --mail-type=END
#SBATCH --mail-user=alberto.izquierdom@um.es


python3.6 EntregableFBA_Ejemplo3.py
