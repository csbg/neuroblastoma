#!/usr/bin/qsub

### max runtime: HH:MM:SS
#$ -l h_rt=24:00:00

### estimated memory usage
#$ -l mt=32g

### request number of CPUs:
###   -pe smp <x>  for a fixed number of cpus that you need
###   -pe smp <a>-<b> any number of CPUs between <a> and <b>
### #$ -pe smp 1-

### priority; 0 is highest priority
#$ -p -500

### send email for (b)eginning, (e)nd, (a)bort and (s)uspension ...
#$ -m beas

### run in current working directory
#$ -cwd

### clone your environment variables
#$ -V

### merge stdout and stderr
#$ -j yes

### set stdout filename
#$ -o $JOB_ID_$SGE_TASK_ID_$JOB_NAME_out.txt

### create an array job - read sample names from these lines in samples.txt
#$ -t 1-6

### jobname
#$ -N velocyto

./run_velocyto_sample.sh `sed "${SGE_TASK_ID}q;d" samples.txt`
