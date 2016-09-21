#!/usr/local/bin/bash

# Download file to ftp site.
#set base remote dir
#REMOTE_DIR="esgi-vbseq/REL-2012-06-07/v1"
#set connection parameters
#SERVER="sftpsrv.sanger.ac.uk"
file=$1
SERVER="storage.hsr.it"
Password="ingivb2012"   # Change above to suit.

#now for each remote/local directory do the same thing...
ftp -i -n $SERVER <<End-Of-Session
        user ingivb "$Password"
        interactive mode off
 #       cd $REMOTE_DIR
        pwd
        get "$file"
        bye
End-Of-Session
#done
exit 0
