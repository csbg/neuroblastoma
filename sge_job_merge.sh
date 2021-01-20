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

### send email for (b)eginning, (e)nd, (a)bort and (s)uspension
#$ -m beas

### run in current working directory
#$ -cwd

### clone your environment variables
#$ -V

### merge stdout and stderr
#$ -j yes

### set stdout filename
#$ -o $JOB_ID_$JOB_NAME_out.txt

### jobname
#$ -N neuroblastoma_merge_all_samples

Rscript merge_samples.R
