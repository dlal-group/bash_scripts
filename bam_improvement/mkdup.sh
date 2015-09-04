set -e

in_file=$1
out_file=$2

module add hgi/picard-tools/1.127

bsub -o log/mkdup.%J.o -e log/mkdup.%J.e -q normal -M2000 -R'select[mem>2000] rusage[mem=2000]' -- picard-tools MarkDuplicates INPUT=$in_file OUTPUT=$out_file METRICS=${out_file}.metrics | sed -n 's/Job <\([0-9]*\)>.*/\1/p'
