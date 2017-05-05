import os
import glob

#did this directory have JEC wiggled runs in it or not?
includeJEC = len(glob.glob('*JES_up*'))>0 and len(glob.glob('*JES_down*'))>0 and len(glob.glob('*JER_up*'))>0 and len(glob.glob('*JER_down*'))>0  

#first get the list of all the output files and root files, and the number of original jobs
outputfilelist = glob.glob('output_JID_*.log')
rootfilelist = glob.glob('*_tree.root')
nJobs = int(os.popen('cat ana.listOfJobs_all | wc -l').read()) if len(glob.glob('ana.listOfJobs_all'))!=0 else int(os.popen('cat ana.listOfJobs | wc -l').read())
#make a list of the failed job numbers
failedjobnumbers = []
#first look for jobs that didn't return anything
print 'len(rootfilelist)=%d, nJobs=%d'%(len(rootfilelist),nJobs)
if len(rootfilelist) < nJobs :
	#if there were JEC files run 
	if includeJEC :
		#for each job number
		for i in range(nJobs/5) :
			theseRootFiles = glob.glob('*_'+str(i)+'_tree.root')
			#there should be five files per job
			if len(theseRootFiles) < 5 :
				print 'Missing some output from job number '+str(i)+', checking which of the JEC wiggles it is'
				newfailedjobnumbers = [5*i,5*i+1,5*i+2,5*i+3,5*i+4]
				for rfilename in theseRootFiles :
					if rfilename.find('JES')==-1 and rfilename.find('JER')==-1: newfailedjobnumbers.pop(newfailedjobnumbers.index(5*i))
					elif rfilename.find('JES_up')!=-1 : newfailedjobnumbers.pop(newfailedjobnumbers.index(5*i+1))
					elif rfilename.find('JES_down')!=-1 : newfailedjobnumbers.pop(newfailedjobnumbers.index(5*i+2))
					elif rfilename.find('JER_up')!=-1 : newfailedjobnumbers.pop(newfailedjobnumbers.index(5*i+3))
					elif rfilename.find('JER_down')!=-1 : newfailedjobnumbers.pop(newfailedjobnumbers.index(5*i+4))
				failedjobnumbers+=newfailedjobnumbers
	#otherwise it's simpler
	else :
		#for each job number
		for i in range(nJobs) :
			#check if the there's an outputted file
			if len(glob.glob('*_'+str(i)+'_tree.root'))==0 :
				print 'Job number '+str(i)+' had no output!'
				failedjobnumbers.append(i)
#now check the file sizes to find any that are abnormally small
totalsize = 0.
for rootfile in rootfilelist :
	totalsize+=os.path.getsize(rootfile)
expected_contribution = totalsize/len(rootfilelist)
for rootfile in rootfilelist :
	filesize = os.path.getsize(rootfile)
	if filesize/expected_contribution<0.80 :
		print 'File '+rootfile+' is too small, its size is '+str(filesize)+' bytes, contributing '+str(filesize/expected_contribution)+' of its expectation'
		jobnumber = int(rootfile.rstrip('_tree.root').split('_')[len(rootfile.rstrip('_tree.root').split('_'))-1])
		if includeJEC and (not 'singleel' in rootfile.lower() and not 'singlemu' in rootfile.lower()) :
			jobnumber *= 5
			if rootfile.find('JES_up')!=-1 : jobnumber+=1
			elif rootfile.find('JES_down')!=-1 : jobnumber+=2
			elif rootfile.find('JER_up')!=-1 : jobnumber+=3
			elif rootfile.find('JER_down')!=-1 : jobnumber+=4
		if not jobnumber in failedjobnumbers : failedjobnumbers.append(jobnumber)
#sort the list of failed job numbers
failedjobnumbers.sort()
#open the list of all the jobs and add the failed ones to the new file
linecount = 0
if not os.path.isfile('ana.listOfJobs_all') :
	print 'TOTAL LIST OF JOBS DOES NOT EXIST YET, COPYING CURRENT LIST OF JOBS!!'
	os.system('mv ana.listOfJobs ana.listOfJobs_all')
os.system('rm -rf ana.listOfJobs')
joblist = open('ana.listOfJobs_all','r')
for job in joblist.readlines() :
	if linecount in failedjobnumbers : os.system('echo "'+job+'" >> ana.listOfJobs')
	linecount+=1
print 'Total new list of jobs: '
os.system('cat ana.listOfJobs')
os.system('bash cleanup.bash')
