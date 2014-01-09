#!/usr/local/bin/bash

#This script allow to cycle through you running jobs and move them from one queue to another...
#ARGS:
#$1 = n cycles
#$2 = sleep time
#$3 = starting queue
#$4 = destination queue

if [ $# -lt 4 ]
then
        echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
        echo -e "jobs_switcher.sh <n_cycles> <sleep time> <starting queue> <destination queue> \n"

        exit 1
fi

cycles=$1
s_time=$2
s_queue=$3
e_queue=$4

let i=1
while [ $i -le $cycles ]
do
	sleep ${s_time}
	jobs=$(bjobs -r -q ${s_queue} | tail -n+2| cut -f 1 -d " ")
	n_jobs=`bjobs -r -q ${s_queue} | fgrep -v JOBID | wc -l`
	
	if [ ${n_jobs} -gt 0 ]
	then
		for job in $jobs
		do
			bswitch ${e_queue} $job
		done
		echo "${i}st batch moved"
	else
		echo "NO running jobs to switch"
	fi
	
	
	let i=$i+1
done
