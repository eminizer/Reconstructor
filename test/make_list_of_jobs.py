from optparse import OptionParser
import os

parser = OptionParser()
#Run options to change
parser.add_option('--n_jobs',   type='int',    action='store', default=1,  dest='n_jobs',  
		  help='Number of total grid jobs')
#Sample options to change
parser.add_option('--name',  type='string', action='store',   dest='name',     
		  help='Name of sample or process (used to name output files, etc.)')
parser.add_option('--xSec',  type='string', action='store', default='1.', dest='xSec', 
		  help='Cross section for this process')
#Probably won't change each time
parser.add_option('--input',   type='string', action='store', default='input', dest='input',
                  help='Path to input file holding list of files to run on')
parser.add_option('--on_grid',   type='string', action='store', default='no',  dest='on_grid',
                  help='Changes everything to relative paths if running on the grid, default is "no"')
parser.add_option('--max_events', type='int',    action='store', default=-1,  dest='max_events',
                  help='Maximum number of events to process (default is -1 for "all")')
parser.add_option('--print_every',type='int',    action='store', default=1000,  dest='print_every',
                  help='Print progress after how many events?')

(options, args) = parser.parse_args()

for i in range(options.n_jobs) :
	cmd = 'echo "python ./tardir/run_reconstructor.py --name '+options.name+' --xSec '+options.xSec
	if options.input!='input' :
		cmd+= '--input '+options.input
	if options.on_grid!='no' :
		cmd+= ' --on_grid '+options.on_grid
	if options.n_jobs!=1 :
		cmd = cmd+' --n_jobs '+str(options.n_jobs)+' --i_job '+str(i)
	if options.max_events!=-1 :
		cmd = cmd+' --max_events '+str(options.max_events)
	if options.print_every!=1000 :
		cmd+= ' --print_every '+str(options.print_every)
	print cmd
	os.system(cmd+'" >> ana.listOfJobs')
	if options.name.lower().find('singleel')==-1 and options.name.lower().find('singlemu')==-1 :
		#and also the JES and JER up/down jobs
		#os.system(cmd+' --JES up" >> ana.listOfJobs')
		#os.system(cmd+' --JES down" >> ana.listOfJobs')
		#os.system(cmd+' --JER up" >> ana.listOfJobs')
		#os.system(cmd+' --JER down" >> ana.listOfJobs')
		continue

print 'done. Completed file: '
os.system('cat ana.listOfJobs')
