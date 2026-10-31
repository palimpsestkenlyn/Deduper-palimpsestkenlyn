#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=deduper
#SBATCH --output=deduper_%j.out

/usr/bin/time -v python hodapp_deduper.py \
    -f /scratch/bgmp/hodapp/deduper/C1_SE_uniqAlign.sorted.sam \
    -o /scratch/bgmp/hodapp/deduper/C1_SE_uniqAlign.deduped.sam \
    -u /projects/bgmp/hodapp/bioinfo/Bi624/Deduper-palimpsestkenlyn/STL96.txt