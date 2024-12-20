#!/usr/bin/bash

##############################
#       Job blueprint        #
##############################

# echo "First arg: $1"

# Give your job a name, so you can recognize it in the queue overview
#SBATCH --job-name=llFBA

# Define, how many nodes you need. Here, we ask for 1 node.
# Each node has 16 or 20 CPU cores.
#SBATCH --nodes=1
# You can further define the number of tasks with --ntasks-per-*
# See "man sbatch" for details. e.g. --ntasks=4 will ask for 4 cpus.

# Define, how long the job will run in real time. This is a hard cap meaning
# that if the job runs longer than what is written here, it will be
# force-stopped by the server. If you make the expected time too long, it will
# take longer for the job to start. Here, we say the job will take 5 minutes.
#              d-hh:mm:ss
#SBATCH --time=0-23:40:00

# Define the partition on which the job shall run. May be omitted.
#SBATCH --partition small

# How much memory you need.
# --mem will define memory per node and
# --mem-per-cpu will define memory per CPU/core. Choose one of those.
##SBATCH --mem-per-cpu=1500MB
#SBATCH --mem=30GB    # this one is not in effect, due to the double hash
##SBATCH --exclusive

# Turn on mail notification. There are many possible self-explaining values:
# NONE, BEGIN, END, FAIL, ALL (including all aforementioned)
# For more values, check "man sbatch"
#SBATCH --mail-type=END

# You may not place any commands before the last SBATCH directive
julia --project loopless_fba.jl $1 $2 $3 $4 $5 &> ll_fba_$1_$SLURM_JOB_ID.txt

# Finish the script
exit 0
