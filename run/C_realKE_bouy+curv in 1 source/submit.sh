#!/bin/bash
#SBATCH --account=def-cdegroot
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=1024M
#SBATCH --time=4-10:00

module load openfoam/11

decomposePar



srun foamRun -parallel 
