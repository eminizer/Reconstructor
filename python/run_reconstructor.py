#Imports
import os
import sys
from DataFormats.FWLite import Events, Handle
from ROOT import TChain, TFile
from optparse import OptionParser
from reconstructor import Reconstructor

##########								Parser Options								##########

parser = OptionParser()
#Run options
parser.add_option('--input', 	  type='string', action='store', default='input', dest='input',	   	  
	help='Path to input file holding list of files to run on (default: "input")')
parser.add_option('--ttree_dir_name', type='string', action='store', default='B2GTTreeMaker', dest='ttree_dir_name',	   	  
	help='Name of directory that holds trees in B2GTTree file (default: "B2GTTreeMaker")')
parser.add_option('--ttree_name', type='string', action='store', default='B2GTree', dest='ttree_name',	   	  
	help='Name of tree in B2GTTree file (default: "B2GTree")')
parser.add_option('--on_grid', 	  type='string', action='store', default='no',	  dest='on_grid',	  
	help='Changes everything to relative paths if running on the grid, default is "no"')
parser.add_option('--max_events', type='int',    action='store', default=-1,	  dest='max_events',  
	help='Maximum number of events to process (default is -1 for "all")')
parser.add_option('--print_every',type='int',    action='store', default=1000,	  dest='print_every', 
	help='Print progress after how many events?')
parser.add_option('--n_jobs', 	  type='int',    action='store', default=1,		  dest='n_jobs',	  
	help='Number of total grid jobs')
parser.add_option('--i_job', 	  type='int',    action='store', default=0,		  dest='i_job',	   	  
	help='Which job is this in the sequence?')
#Sample options
parser.add_option('--name', 		 type='string', action='store', 			  	dest='name', 		    
	help='Name of sample or process (used to name output files, etc.)')
parser.add_option('--generator', 	 type='string', action='store', default='none', dest='generator', 		
	help='Monte Carlo generator for this file (powheg, madgraph, pythia8, mcatnlo); default is "none"')
parser.add_option('--JEC', type='string', action='store', default='nominal',  dest='JEC',  
	help='JEC systematics: shift JEC up or down (default is nominal)')
(options, args) = parser.parse_args()

##########							Set Up Event Loop								##########

print 'Opening files for sample '+options.name+' . . .'  
#Build path to input file
input_files_list = ''
if options.on_grid == 'yes' :
	input_files_list += 'tardir/'
else :
	input_files_list += './'
input_files_list += options.input
if not options.input.endswith('.txt') :
	input_files_list += '.txt'
print 'Using input file '+input_files_list+''
#open input file in read only mode 
input_files_list = open(input_files_list,'r')
chain = TChain(options.ttree_name)
print 'Getting these files: '
#Read files in line by line and get the tree and total weight value from each
totweight = 0.
for input_file in input_files_list :
	print '	'+input_file.rstrip()+''
	chain.AddFile(input_file.rstrip()+'/'+options.ttree_dir_name+'/'+options.ttree_name)
	f = TFile.Open(input_file.rstrip()) 
	histo=f.Get('EventCounter/totweight')
	newweight=histo.GetBinContent(1)
	print '		Added %.2f to total weight'%(newweight)
	totweight+=newweight
print 'TOTAL SUM OF EVENT WEIGHTS = '+str(totweight)
ntotalevents = chain.GetEntries()
#Set filename for analyzer from sample name
filename = options.name
if options.JEC.lower() != 'nominal' :
	filename+='_JEC_'+options.JEC.lower()
if options.n_jobs>1 :
	filename+='_'+str(options.i_job)
filename+='_tree.root'
#Initialize analyzer
data=False
if options.name.lower().find('singlemu')!=-1 or options.name.lower().find('singleel')!=-1 :
	data=True
analyzer = Reconstructor(filename, chain, data, options.generator, options.JEC.lower(), options.on_grid, totweight)

#Counters
real_count = 0
count = 0

##########								Main Event Loop								##########

print 'Files opened, starting event loop'
for event in range(ntotalevents) :
	#increment the "real" counter
	real_count+=1
	#check the grid split
	if ((real_count-1)-options.i_job) % options.n_jobs != 0 :
		continue
	count+=1
	#check the max events 
	if count == options.max_events+1 :
		print 'Processed event number '+str(count-1)+', exiting'
		break
	#print progress
	if count % options.print_every == 0 or count == 1:
			print ( 'Count at '+str(count)+' out of '+str(ntotalevents/options.n_jobs)+', (%.4f%% complete)'
				%(float(count) / float(ntotalevents/options.n_jobs) * 100.0) )
	#analyze event and add to TTree
	analyzer.analyze(event)
	#reset analyzer
	analyzer.reset()
#clean up after yourself
del analyzer