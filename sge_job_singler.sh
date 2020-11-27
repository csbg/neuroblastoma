#!/usr/bin/qsub

### max runtime: HH:MM:SS
#$ -l h_rt=5:00:00

### estimated memory usage
#$ -l mt=5000m

### request number of CPUs:
###   -pe smp <x>  for a fixed number of cpus that you need
###   -pe smp <a>-<b> any number of CPUs between <a> and <b>
### #$ -pe smp 2-

### priority; 0 is highest priority
#$ -p -500

### send email for (b)eginning, (e)nd, (a)bort and (s)uspension ...
#$ -m beas

### ... to the following address
#$ -M wolfgang.esser-skala@sbg.ac.at

### run in current working directory
#$ -cwd

### clone your environment variables
#$ -V

### jobname
#$ -N neuroblastoma_singler

Rscript cell_types_singler.R
