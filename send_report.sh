#!/usr/local/bin/bash

#This bit should go in the script that will call the sender
# PID=$!
# wait $!
# status=$?
# wdir=`pwd -p`
# cmd=`history | tail -n2| head -1| cut -f 2- -d " "`
# email=mc14@sanger.ac.uk
# send_report.sh ${status} ${email} ${wdir} ${cmd}


#small script to generate a mail report
status=$1
email=$2
wdir=$3
cmd=$4
info=$5

echo -e "The runner has finished, all done!\nWorking directory: ${wdir}\nCommand line: ${cmd}\nAdditional info:\n ${info}" | mutt -s "Script report: ${status}" ${email}