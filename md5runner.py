#Script to create a table to perform annotation with ANCESTRAL alleles on VCF data
import sys 
import os 
import subprocess
import threading
from threading import Thread

try: from queue import Queue
except ImportError:
	from Queue import Queue # Python 2.x

def worker(queue):
	for cmd in iter(queue.get, None):
		# print cmd
		subprocess.check_call(cmd, stdout=None, stderr=subprocess.STDOUT)

file_list=sys.argv[1]
# file_list="/lacie/01/ftp/files/WES_FVG/presbyacusis_fastqfiles/file_list.txt"

#read the file list
files_to_md=open(file_list,'r')

fastqs=files_to_md.read().split()

q = Queue()
threads = [Thread(target=worker, args=(q,)) for _ in range(40)]
for t in threads: # start workers
	t.daemon = True
	t.start()

for i in range(0,len(fastqs)):
	out_file = fastqs[i]+".md5"
	# command=["md5sum",fastqs[i],">",out_file]
	command="md5sum " + fastqs[i] + " > "+ out_file
	q.put_nowait(command)

for _ in threads: q.put(None) # signal no more commands
for t in threads: t.join()    # wait for completion


# def worker():
# 	"""thread worker function"""
# 	print 'Worker'
# 	return

# threads = []
# for i in range(5):
# 	t = threading.Thread(target=worker)
# 	threads.append(t)
# 	t.start()