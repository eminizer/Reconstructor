#Imports
import os
import sys
from ROOT import TChain, TFile, TH2D, TTree
from array import array
from math import *
from time import *
from datetime import timedelta
from optparse import OptionParser
from branch import Branch

################################   addBranch function  #################################
def AddBranch(readname=None,writename=None,ttreetype='F',inival=-900.,size='1',dictlist=None) :
	newBranch = Branch(readname=readname,writename=writename,ttreetype=ttreetype,inival=inival,size=size)
	for d in dictlist :
		if writename!=None :
			d[writename] = newBranch
		else :
			d[readname] = newBranch
	return newBranch

##########								Parser Options								##########

startTime = time()

parser = OptionParser()
#Run options
parser.add_option('--input', 	  type='string', action='store', default='input', dest='input',	   	  
	help='Path to input file holding list of files to run on (default: "input")')
parser.add_option('--ttree_dir_name', type='string', action='store', default='B2GTTreeMaker', dest='ttree_dir_name',	   	  
	help='Name of directory that holds trees in B2GTTree file (default: "B2GTTreeMaker")')
parser.add_option('--ttree_name', type='string', action='store', default='B2GTree', dest='ttree_name',	   	  
	help='Name of tree in B2GTTree file (default: "B2GTree")')
parser.add_option('--max_events', type='int',    action='store', default=-1,	  dest='max_events',  
	help='Maximum number of events to process (default is -1 for "all")')
parser.add_option('--print_every',type='int',    action='store', default=100,	  dest='print_every', 
	help='Print progress after how many events?')
parser.add_option('--n_jobs', 	  type='int',    action='store', default=1,		  dest='n_jobs',	  
	help='Number of total grid jobs')
parser.add_option('--i_job', 	  type='int',    action='store', default=0,		  dest='i_job',	   	  
	help='Which job is this in the sequence?')
(options, args) = parser.parse_args()

##########							Set Up Event Loop								##########

#Build path to input file
input_files_list = './'+options.input
if not options.input.endswith('.txt') :
	input_files_list += '.txt'
print 'Using input file '+input_files_list+''
#open input file in read only mode 
input_files_list = open(input_files_list,'r')
#Set up the garbage output file and the chain
garbageFileName = 'garbage'
if options.n_jobs!=1 :
	garbageFileName+='_'+str(options.i_job)
garbageFileName+='.root' 
garbageFile = TFile(garbageFileName,'recreate') 
chain = TChain(options.ttree_name)
print 'Getting these files: '
for input_file in input_files_list :
	print '	'+input_file.rstrip()+''
	chain.Add(input_file.rstrip()+'/'+options.ttree_dir_name+'/'+options.ttree_name)
garbageTree = TTree('garbageTree','recreate')
garbageTree.SetDirectory(garbageFile)
#Get the total number of events
ntotalevents = chain.GetEntries()
print 'number of total events = %d'%(ntotalevents) 
#how many events should be in this job?
nanalysisevents = abs(int(ntotalevents/options.n_jobs))
if options.max_events!=-1 :
	nanalysisevents = min(options.max_events, nanalysisevents)
if options.n_jobs>1 :
	#adjust so they don't all pile up in the last job
	while ntotalevents-options.n_jobs*nanalysisevents > nanalysisevents :
		nanalysisevents+=1
analysisTree = chain
if options.n_jobs>1 :
	analysisTree = chain.CopyTree('','',nanalysisevents,options.i_job*nanalysisevents); analysisTree.SetDirectory(garbageFile)
#define the output file and the histograms in it
outputfilename = 'MC_btagging_efficiency'
if options.n_jobs!=1 :
	outputfilename+='_'+str(options.i_job)
outputfilename+='.root'
outputfile = TFile(outputfilename,'recreate')
aetabins = array('d',[0.,0.2,0.3,0.9,1.2,1.6,2.1,2.4,3.0])
ptbins  = array('d',[0.,15.,30.,60.,100.,140.,180.,240.,300.,500.,1000.])
naetabins = len(aetabins)-1
nptbins = len(ptbins)-1
bjet_pass = TH2D('bjet_pass','b-jet passing CSVv2L pt vs |eta|; |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
bjet_fail = TH2D('bjet_fail','b-jet failing CSVv2L pt vs |eta|; |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
bjet_ratio = TH2D('bjet_ratio','b-jet CSVv2L efficiency (vs. pt vs. |eta|); |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
cjet_pass = TH2D('cjet_pass','c-jet passing CSVv2L pt vs |eta|; |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
cjet_fail = TH2D('cjet_fail','c-jet failing CSVv2L pt vs |eta|; |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
cjet_ratio = TH2D('cjet_ratio','c-jet CSVv2L efficiency (vs. pt vs. |eta|); |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
udsg_pass = TH2D('udsg_pass','(u,d,s,g)-jet passing CSVv2L pt vs |eta|; |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
udsg_fail = TH2D('udsg_fail','(u,d,s,g)-jet failing CSVv2L pt vs |eta|; |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
udsg_ratio = TH2D('udsg_ratio','(u,d,s,g)-jet CSVv2L efficiency (vs. pt vs. |eta|); |#eta|; p_{T}',naetabins,aetabins,nptbins,ptbins)
#add the branches we'll need
allbranches = {}
size 					= AddBranch(readname='jetAK4CHS_size',ttreetype='i',dictlist=[allbranches])
pts 					= AddBranch(readname='jetAK4CHS_Pt',size='jetAK4CHS_size',dictlist=[allbranches])
etas 					= AddBranch(readname='jetAK4CHS_Eta',size='jetAK4CHS_size',dictlist=[allbranches])
csvv2s 					= AddBranch(readname='jetAK4CHS_CSVv2',size='jetAK4CHS_size',dictlist=[allbranches])
neutralMultiplicity 	= AddBranch(readname='jetAK4CHS_neutralMultiplicity',size='jetAK4CHS_size',dictlist=[allbranches])
neutralHadronEnergyFrac = AddBranch(readname='jetAK4CHS_neutralHadronEnergyFrac',size='jetAK4CHS_size',dictlist=[allbranches])
neutralEmEnergyFrac 	= AddBranch(readname='jetAK4CHS_neutralEmEnergyFrac',size='jetAK4CHS_size',dictlist=[allbranches])
chargedHadronEnergyFrac = AddBranch(readname='jetAK4CHS_chargedHadronEnergyFrac',size='jetAK4CHS_size',dictlist=[allbranches])
chargedEmEnergyFrac 	= AddBranch(readname='jetAK4CHS_chargedEmEnergyFrac',size='jetAK4CHS_size',dictlist=[allbranches])
chargedMultiplicity 	= AddBranch(readname='jetAK4CHS_chargedMultiplicity',size='jetAK4CHS_size',dictlist=[allbranches])
flavors 				= AddBranch(readname='jetAK4CHS_HadronFlavour',size='jetAK4CHS_size',dictlist=[allbranches])
for branch in allbranches.values() :
	branch.initialize(analysisTree,garbageTree)

#count
count=0
#If there is no max_events argument, set it to the total number of analyzed events
maxEvents = options.max_events if options.max_events!=-1 else nanalysisevents

setupDoneTime = time(); setupTime = timedelta(seconds=setupDoneTime-startTime)
print 'Done Setting up; time taken: %02d:%02d:%02d'%(setupTime.seconds/3600,(setupTime.seconds%3600)/60,(setupTime.seconds%60))

##########								Main Event Loop								##########
print 'Files opened, starting event loop'
for event in range(ntotalevents) : 
	count+=1
	#check the max events 
	if count == maxEvents+1 :
		print 'Processed event number '+str(count-1)+', exiting'
		break
	#print progress
	if count % options.print_every == 0 and count!=0 :
			timeSinceSetup = timedelta(seconds=time()-setupDoneTime)
			percentDone = float(count) / float(min(nanalysisevents,maxEvents)) * 100.0
			timeLeft = timedelta(seconds=(100.-percentDone)*timeSinceSetup.total_seconds()/percentDone)
			print ( 'Count at %d out of %d, (%.4f%% complete; time elapsed = %02d:%02d:%02d, approx. time left = %02d:%02d:%02d)'
				   %(count,min(nanalysisevents,maxEvents),percentDone,timeSinceSetup.seconds/3600,(timeSinceSetup.seconds%3600)/60,(timeSinceSetup.seconds%60),timeLeft.seconds/3600,(timeLeft.seconds%3600)/60,(timeLeft.seconds%60))  )
	
	#get the event from the chain
	check = analysisTree.GetEntry(event)
	#print '%d jets in this event'%(size.getReadValue()) #DEBUG
	#for every jet in the event
	for i in range(size.getReadValue()) :
		#make cuts
		jetisValid = False
		aeta = abs(etas.getReadValue(i))
		pt = pts.getReadValue(i)
	#	print 'pt = %.1f, |eta| = %.2f'%(pt,aeta) #DEBUG
		if pt>30. and aeta<2.4 :
			if aeta<2.7 :
				if neutralHadronEnergyFrac.getReadValue(i)<0.99 :
					if neutralEmEnergyFrac.getReadValue(i)<0.99 :
						if neutralMultiplicity.getReadValue(i)+chargedMultiplicity.getReadValue(i)>1 :
							if aeta<2.4 :
								if chargedHadronEnergyFrac.getReadValue(i)>0. and chargedMultiplicity.getReadValue(i)>0 and chargedEmEnergyFrac.getReadValue(i)<0.99 :
									jetisValid=True
							else :
								jetisValid=True
			elif aeta<3.0 :
				if neutralEmEnergyFrac.getReadValue(i)>0.1 and neutralHadronEnergyFrac.getReadValue(i)<0.98 and neutralMultiplicity.getReadValue(i)>2 :
					jetisValid=True
			else :
				if neutralEmEnergyFrac.getReadValue(i)<0.9 and neutralMultiplicity.getReadValue(i)>10 :
					jetisValid=True
		if not jetisValid :
			continue
	#	print 'valid jet!' #DEBUG
		#add to histograms
		flav = flavors.getReadValue(i)
		passtag = csvv2s.getReadValue(i)>0.5426
		if flav==5 :
			if passtag :
				bjet_pass.Fill(aeta,pt)
			else :
				bjet_fail.Fill(aeta,pt)
		elif flav==4 :
			if passtag :
				cjet_pass.Fill(aeta,pt)
			else :
				cjet_fail.Fill(aeta,pt)
		elif flav==0 :
			if passtag :
				udsg_pass.Fill(aeta,pt)
			else :
				udsg_fail.Fill(aeta,pt)
		else :
			print 'WARNING: Flavor %s not recognized!!!'%(flav)
	#reset branches
	for branch in allbranches.values() :
		branch.reset()
#set the values in the ratio histograms
for i in range(bjet_pass.GetNbinsX()*bjet_pass.GetNbinsY()) :
	if not (bjet_pass.IsBinOverflow(i) or bjet_pass.IsBinUnderflow(i)) :
		bp = bjet_pass.GetBinContent(i)
		ball = bp+bjet_fail.GetBinContent(i)
		if ball!=0. :
			bjet_ratio.SetBinContent(i,bp/ball)
		cp = cjet_pass.GetBinContent(i)
		call = cp+cjet_fail.GetBinContent(i)
		if call!=0. :
			cjet_ratio.SetBinContent(i,cp/call)
		lp = udsg_pass.GetBinContent(i)
		lall = lp+udsg_fail.GetBinContent(i)
		if lall!=0. :
			udsg_ratio.SetBinContent(i,lp/lall)

#closeout
garbageFile.Close()
os.system('rm -rf '+garbageFileName)
outputfile.cd()
outputfile.Write()
outputfile.Close()
