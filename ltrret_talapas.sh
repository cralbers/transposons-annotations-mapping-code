#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=kimltrret
#SBATCH --output=kimltrret.out
#SBATCH --error=kimltrret.err
#SBATCH --time=0-23:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mail-user=calbers@uoregon.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=libudalab

conda activate LTR_retriever

cd /home/calbers/libudalab/kim_ltrretriever

/home/calbers/libudalab/kim_ltrretriever/LTR_retriever/LTR_retriever -genome sequence.fasta -inharvest genome.fa.rawLTR.scn -verbose -u .0000002808 > kim_LTR_retriever.out