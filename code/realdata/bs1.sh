{\rtf1\ansi\ansicpg1252\cocoartf2639
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/env bash\
#SBATCH --job-name=covidMRP\
#SBATCH --time=05:00:00\
#SBATCH --ntasks=1\
#SBATCH --nodes=1\
#SBATCH --cpus-per-task=1\
#SBATCH --output=/dev/null\
#SBATCH --mem-per-cpu=8000\
#SBATCH --array=1:5\
\
module load R/4.4.1-gfbf-2023a\
\
R CMD BATCH --no-save --no-restore SBC_posterior.R script_$SLURM_ARRAY_TASK_ID }