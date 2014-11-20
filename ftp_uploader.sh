#!/usr/local/bin/bash

# Upload file to ftp site.
#set connection parameters
if [ $# -lt 4 ]
then
	echo "ERROR!!Missing arguments!!"
	echo "USAGE:"
	echo "ftp_uploader.sh <ftpserver> <username> <password> <filename/dirname> [dir_mode]"
exit 1
fi

SERVER=$1
USER=$2
Password=$3   # Change above to suit.
FILENAME=$4

if [ $# -eq 5 ] 
then
DIRNAME=$4
#upload all gz content of a folder

cd $DIRNAME
ftp -v -i -n $SERVER <<End-Of-Session
        user "$USER" "$Password"
        pwd
        mput "*.*"
        bye
End-Of-Session
exit 0
else

#now for each remote/local directory do the same thing...
#REMOTE_DIR=$4
ftp -v -i -n $SERVER <<End-Of-Session
        user "$USER" "$Password"
        pwd
        put "${FILENAME}"
        bye
End-Of-Session
exit 0

fi
