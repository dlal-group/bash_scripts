set -e

in_file=$1
out_file=$2
in_jobid=$3

bsub -o log/call.%J.o -e log/call.%J.e -w "done(${in_jobid})" -q normal -M29000 -n 32 -R'span[hosts=1] select[mem>29000] rusage[mem=29000]' -- ./call_inner.sh $in_file $out_file
