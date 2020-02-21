#Imports
import os
import sys
from ROOT import TMinuit, Long, Double, TChain, TFile, TH1D, TH2D, TF2
from array import array
from math import *
from time import *
from datetime import timedelta
from optparse import OptionParser
from reconstructor import Reconstructor

#global 2D histograms needed for alpha and epsilon fits
csvb_qq_global_hist = TH2D(); csvb_gg_global_hist = TH2D()
qq_global_hist_projection = TH1D(); gg_global_hist_projection = TH1D()
BETA = [-1.]

#function to symmetrize the distributions
def symmetrize(hist) :
	newhist = hist.Clone(); newhist.Reset()
	#loop over the bins and add half the content in each to the bins with the same and opposite c* values
	for i in range(1,hist.GetXaxis().GetNbins()+1) :
		bincenter = hist.GetXaxis().GetBinCenter(i)
		content = hist.GetBinContent(i)
		add_err = hist.GetBinError(i)
		old_err = newhist.GetBinError(i)
		newhist.Fill(bincenter,0.5*content)
		newhist.Fill(-1.*bincenter,0.5*content)
		newhist.SetBinError(newhist.FindFixBin(bincenter),sqrt(old_err**2+add_err**2))
		newhist.SetBinError(newhist.FindFixBin(-1.*bincenter),sqrt(old_err**2+add_err**2))
	return newhist

#minuit fitting function for alpha
def alpha_fcn(npar, deriv, f, par, flag) :
	alpha = par[0]
	lnL=0.
	#loop over the costheta bins
	for i in range(1,qq_global_hist_projection.GetNbinsX()+1) :
		cstar = qq_global_hist_projection.GetXaxis().GetBinCenter(i)
		#increment the values
		cont=qq_global_hist_projection.GetBinContent(i)
		f_L = (2.+(BETA[0]*BETA[0])*(cstar*cstar)-(BETA[0]*BETA[0])+alpha*(1.-(BETA[0]*BETA[0])*(cstar*cstar)))/(2.*(2-2.*(BETA[0]*BETA[0])/3.+alpha*(1.-(BETA[0]*BETA[0])/3.)))
		if f_L>0 :
			lnL+=-2.0*cont*log(f_L)
		else :
			lnL+=10000000000.
	#print '	lnL = %.8f, alpha = %.8f'%(lnL,alpha) #DEBUG
	f[0] = lnL

#minuit fitting function for epsilon
def epsilon_fcn(npar, deriv, f, par, flag) :
	epsilon = par[0]
	lnL = 0.0
	#loop over the costheta bins
	for i in range(1,gg_global_hist_projection.GetXaxis().GetNbins()+1) :
		cstar = gg_global_hist_projection.GetXaxis().GetBinCenter(i)
		#increment the values
		cont=gg_global_hist_projection.GetBinContent(i)
		beta2 = (BETA[0])*(BETA[0])
		cstar2 = cstar*cstar
		bc2 = beta2*cstar2
		omb2 = 1. - beta2
		beta4 = beta2*beta2
		rnorm = ( (epsilon*(34.*beta4-100.*beta2+98.)+2.*beta4-36.*beta2+66.)*atanh(BETA[0])/(BETA[0])
					-9./5.*epsilon*beta4+6.*beta2*(epsilon*(beta2-2.38889)-0.5)-16.*(1.+epsilon)*(beta4-2.*beta2+1.)/(1.-beta2)
					-epsilon*(18.*beta4-68.*beta2+82.)+18.*beta2-43. )
		num = (7.+9.*bc2)/(1.-bc2)*((1.+bc2)/2.+omb2*beta2*(1.-cstar2)/(1.-bc2))*(1.+epsilon*bc2)
		f_L = num/rnorm
		lnL+=-2.0*cont*log(f_L)
	#print '	lnL = %.8f, epsilon = %.8f'%(lnL,epsilon) #DEBUG
	f[0] = lnL

#Helper function for making the renormalization dictionary
def make_renormalization_dict(name,alpha,epsilon,muRup,muRdown,muFup,muFdown,scup,scdown,pdfas,pdfasup,pdfasdown,topptrwhist) :
	global qq_global_hist_projection
	global gg_global_hist_projection
	returndict = {}
	#start with the simple values
	returndict['muRup']=muRup.GetMean() 
	returndict['muRdown']=muRdown.GetMean() 
	returndict['muFup']=muFup.GetMean() 
	returndict['muFdown']=muFdown.GetMean() 
	returndict['scup']=scup.GetMean() 
	returndict['scdown']=scdown.GetMean() 
	returndict['pdfas']=pdfas.GetMean() 
	returndict['pdfasup']=pdfasup.GetMean() 
	returndict['pdfasdown']=pdfasdown.GetMean() 
	returndict['topptrwmean']=topptrwhist.GetMean()
	printlines = []
	if alpha==1.0 and epsilon==1.0 and (name.find('_TT')!=-1 or name.find('_TTJets')!=-1) :
		#now we've got to do the fits for the "deweighting" alpha and epsilon values
		print 'fitting for alpha/epsilon values... '
		#Define the Minuit instance we'll use
		alpha_minuit = TMinuit(1); alpha_minuit.SetFCN(alpha_fcn)
		epsilon_minuit = TMinuit(1); epsilon_minuit.SetFCN(epsilon_fcn)
		#miscellaneous minuit stuff
		ierflag = Long(1); arglist = array('d',[-1])
		alpha_minuit.mnexcm('SET PRINT', arglist, 1,ierflag); alpha_minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
		epsilon_minuit.mnexcm('SET PRINT', arglist, 1,ierflag); epsilon_minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
		arglist[0]=100000.
		#for each bin in beta
		for i in range(1,csvb_qq_global_hist.GetXaxis().GetNbins()+1) :
			BETA[0] = csvb_qq_global_hist.GetXaxis().GetBinCenter(i)
			#qq_global_hist_projection = symmetrize(csvb_qq_global_hist.ProjectionY('qq_proj_'+str(i),i,i))
			qq_global_hist_projection = csvb_qq_global_hist.ProjectionY('qq_proj_'+str(i),i,i)
			#add parameter
			alpha_minuit.mnparm(0,'alpha',0.1,0.5,0.,0.,ierflag)
			#minimize
			alpha_minuit.mnexcm('MIGRAD',arglist,1,ierflag)
			#get back the fitted alpha value
			fitted_alpha=Double(0.0); fitted_alpha_err=Double(0.0)
			alpha_minuit.GetParameter(0,fitted_alpha,fitted_alpha_err)
			#print 'NEW NAME = alpha_'+str(csvb_qq_global_hist.GetXaxis().GetBinLowEdge(i))
			returndict['alpha_'+str(csvb_qq_global_hist.GetXaxis().GetBinLowEdge(i))]=fitted_alpha
			printlines.append('fitted alpha value in bin %d = %.4f +/- %.4f'%(i,fitted_alpha,fitted_alpha_err))
		#epsilon stuff
		for i in range(1,csvb_gg_global_hist.GetXaxis().GetNbins()+1) :
			BETA[0] = csvb_gg_global_hist.GetXaxis().GetBinCenter(i)
			#gg_global_hist_projection = symmetrize(csvb_gg_global_hist.ProjectionY('gg_proj_'+str(i),i,i))
			gg_global_hist_projection = csvb_gg_global_hist.ProjectionY('gg_proj_'+str(i),i,i)
			epsilon_minuit.mnparm(0,'epsilon',0.1,0.5,0.,0.,ierflag)
			epsilon_minuit.mnexcm('MIGRAD',arglist,1,ierflag)
			fitted_epsilon=Double(0.0); fitted_epsilon_err=Double(0.0)
			epsilon_minuit.GetParameter(0,fitted_epsilon,fitted_epsilon_err)
			returndict['epsilon_'+str(csvb_gg_global_hist.GetXaxis().GetBinLowEdge(i))]=fitted_epsilon
			printlines.append('fitted epsilon value in bin %d = %.4f +/- %.4f'%(i,fitted_epsilon,fitted_epsilon_err))
	else :
		returndict['alpha']=alpha; returndict['epsilon']=epsilon
		fitted_alpha_err = 1.0; fitted_epsilon_err = 1.0
	#print the values
	print 'muRup value = %.4f'%(returndict['muRup'])
	print 'muRdown value = %.4f'%(returndict['muRdown'])
	print 'muFup value = %.4f'%(returndict['muFup'])
	print 'muFdown value = %.4f'%(returndict['muFdown'])
	print 'scup value = %.4f'%(returndict['scup'])
	print 'scdown value = %.4f'%(returndict['scdown'])
	print 'pdfas value = %.4f'%(returndict['pdfas'])
	print 'pdfasup value = %.4f'%(returndict['pdfasup'])
	print 'pdfasdown value = %.4f'%(returndict['pdfasdown'])
	print 'top pT reweight (v1) mean value = %.4f'%(returndict['topptrwmean'])
	for line in printlines :
		print line
	#return the dictionary
	return returndict

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
parser.add_option('--on_grid', 	  type='string', action='store', default='no',	  dest='on_grid',	  
	help='Changes everything to relative paths if running on the grid, default is "no"')
parser.add_option('--max_events', type='int',    action='store', default=-1,	  dest='max_events',  
	help='Maximum number of events to process (default is -1 for "all")')
parser.add_option('--print_every',type='int',    action='store', default=100,	  dest='print_every', 
	help='Print progress after how many events?')
parser.add_option('--n_jobs', 	  type='int',    action='store', default=1,		  dest='n_jobs',	  
	help='Number of total grid jobs')
parser.add_option('--i_job', 	  type='int',    action='store', default=0,		  dest='i_job',	   	  
	help='Which job is this in the sequence?')
#Sample options
parser.add_option('--name', 		 type='string', action='store', 			  	dest='name', 		    
	help='Name of sample or process (used to name output files, etc.)')
parser.add_option('--xSec', 		 type='float', action='store', default=1.0,	dest='xSec', 		    
	help="Cross section of sample's process")
parser.add_option('--kFac', 		 type='float', action='store', default=1.0,	dest='kFac', 		    
	help="NNLO/NLO xSec correction factor")
parser.add_option('--alpha', 		 type='float', action='store', default=1.0,	dest='alphaOverride', 		    
	help="Cross section of sample's process")
parser.add_option('--epsilon', 		 type='float', action='store', default=1.0,	dest='epsilonOverride', 		    
	help="Cross section of sample's process")
parser.add_option('--JEC', type='choice', action='store', default='nominal',  dest='jec', 
	choices=['nominal',
			 'JES_up','JES_dn','JER_up','JER_dn',
			 #'AK4JESPU_up','AK4JESEta_up','AK4JESPt_up','AK4JESScale_up','AK4JESTime_up','AK4JESFlav_up','AK4JERStat_up','AK4JERSys_up',
			 #'AK4JESPU_dn','AK4JESEta_dn','AK4JESPt_dn','AK4JESScale_dn','AK4JESTime_dn','AK4JESFlav_dn','AK4JERStat_dn','AK4JERSys_dn',
			 #'AK8JESPU_up','AK8JESEta_up','AK8JESPt_up','AK8JESScale_up','AK8JESTime_up','AK8JESFlav_up','AK8JERStat_up','AK8JERSys_up',
			 #'AK8JESPU_dn','AK8JESEta_dn','AK8JESPt_dn','AK8JESScale_dn','AK8JESTime_dn','AK8JESFlav_dn','AK8JERStat_dn','AK8JERSys_dn',
			 ],  
	help='JEC systematics variations')
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
#Set up the garbage output file and the chain
garbageFileName = 'garbage_'+options.name 
if options.n_jobs!=1 : 
	garbageFileName+='_'+str(options.i_job) 
garbageFileName+='.root' 
garbageFile = TFile(garbageFileName,'recreate') 
chain = TChain(options.ttree_name)
#is it real data?
data=False
if options.name.lower().find('singlemu')!=-1 or options.name.lower().find('singleel')!=-1 :
	data=True
#Read files in line by line and get the tree, pileup,cstar vs. beta,pdf_alphas/F_up/down,scale_comb_sf_up/down,pdf_alphas_sf_up/down histograms, and total weight value from each
print 'Getting these files: '
total_pileup_histo               = TH1D('total_pileup_histo','total MC pileup; pileup; events',100,0.,100.); total_pileup_histo.SetDirectory(0)
total_cstar_vs_beta_qqbar_histo  = TH2D('total_cstar_vs_beta_qqbar_histo','MC truth c* vs. #beta (q#bar{q} events); #beta; c*; events',10,0.,1.,20,-1.,1.); total_cstar_vs_beta_qqbar_histo.SetDirectory(0)
total_cstar_vs_beta_gg_histo     = TH2D('total_cstar_vs_beta_gg_histo','MC truth c* vs. #beta (qg/gg events); #beta; c*; events',10,0.,1.,20,-1.,1.); total_cstar_vs_beta_gg_histo.SetDirectory(0)
total_mu_R_sf_up_histo           = TH1D('total_mu_R_sf_up_histo','total #mu_{R} sf (up); #mu_{R} sf up; events',100,-1.,3.); total_mu_R_sf_up_histo.SetDirectory(0)
total_mu_R_sf_down_histo         = TH1D('total_mu_R_sf_down_histo','total #mu_{R} sf (down); #mu_{R} sf down; events',100,-1.,3.); total_mu_R_sf_down_histo.SetDirectory(0)
total_mu_F_sf_up_histo           = TH1D('total_mu_F_sf_up_histo','total #mu_{F} sf (up); #mu_{F} sf up; events',100,-1.,3.); total_mu_F_sf_up_histo.SetDirectory(0)
total_mu_F_sf_down_histo         = TH1D('total_mu_F_sf_down_histo','total #mu_{F} sf (down); #mu_{F} sf down; events',100,-1.,3.); total_mu_F_sf_down_histo.SetDirectory(0)
total_scale_comb_sf_up_histo     = TH1D('total_scale_comb_sf_up_histo','total comb. scale sf (up); comb. scale sf up; events',100,-1.,3.); total_scale_comb_sf_up_histo.SetDirectory(0)
total_scale_comb_sf_down_histo   = TH1D('total_scale_comb_sf_down_histo','total comb. scale sf (down); comb. scale sf down; events',100,-1.,3.); total_scale_comb_sf_down_histo.SetDirectory(0)
total_pdf_alphas_sf_histo        = TH1D('total_pdf_alphas_sf_histo','total pdf/alpha_{s} sf; pdf/alpha_{s} sf; events',100,-1.,3.); total_pdf_alphas_sf_histo.SetDirectory(0)
total_pdf_alphas_sf_up_histo     = TH1D('total_pdf_alphas_sf_up_histo','total pdf/alpha_{s} sf (up); pdf/alpha_{s} sf up; events',100,-1.,3.); total_pdf_alphas_sf_up_histo.SetDirectory(0)
total_pdf_alphas_sf_down_histo   = TH1D('total_pdf_alphas_sf_down_histo','total pdf/alpha_{s} sf (down); pdf/alpha_{s} sf down; events',100,-1.,3.); total_pdf_alphas_sf_down_histo.SetDirectory(0)
total_top_pt_reweight_v1_histo 	 = TH1D('total_top_pt_reweight_v1_histo','total top p_{T} reweighting factors; top p_{T} reweight; events',150,-3.,3.); total_top_pt_reweight_v1_histo.SetDirectory(0)
totweight = 0.
for input_file in input_files_list :
	print '	'+input_file.rstrip()+''
	chain.Add(input_file.rstrip()+'/'+options.ttree_dir_name+'/'+options.ttree_name)
	f = TFile.Open(input_file.rstrip()) 
	pu_histo = f.Get('EventCounter/pileup')
	total_pileup_histo.Add(pu_histo)
	csvbqq_histo = f.Get('EventCounter/cstar_vs_beta_qqbar')
	total_cstar_vs_beta_qqbar_histo.Add(csvbqq_histo)
	csvbgg_histo = f.Get('EventCounter/cstar_vs_beta_gg')
	total_cstar_vs_beta_gg_histo.Add(csvbgg_histo)
	muRup_histo = f.Get('EventCounter/mu_R_sf_up')
	total_mu_R_sf_up_histo.Add(muRup_histo)
	muRdown_histo = f.Get('EventCounter/mu_R_sf_down')
	total_mu_R_sf_down_histo.Add(muRdown_histo)
	muFup_histo = f.Get('EventCounter/mu_F_sf_up')
	total_mu_F_sf_up_histo.Add(muFup_histo)
	muFdown_histo = f.Get('EventCounter/mu_F_sf_down')
	total_mu_F_sf_down_histo.Add(muFdown_histo)
	scup_histo = f.Get('EventCounter/scale_comb_sf_up')
	total_scale_comb_sf_up_histo.Add(scup_histo)
	scdown_histo = f.Get('EventCounter/scale_comb_sf_down')
	total_scale_comb_sf_down_histo.Add(scdown_histo)
	pdfas_histo = f.Get('EventCounter/pdf_alphas_sf')
	total_pdf_alphas_sf_histo.Add(pdfas_histo)
	pdfasup_histo = f.Get('EventCounter/pdf_alphas_sf_up')
	total_pdf_alphas_sf_up_histo.Add(pdfasup_histo)
	pdfasdown_histo = f.Get('EventCounter/pdf_alphas_sf_down')
	total_pdf_alphas_sf_down_histo.Add(pdfasdown_histo)
	if not data :
		top_pt_reweight_v1_histo = f.Get('EventCounter/top_pt_rw_v1')
		total_top_pt_reweight_v1_histo.Add(top_pt_reweight_v1_histo)
	histo=f.Get('EventCounter/totweight')
	newweight=histo.GetBinContent(1)
	print '		Added %.2f to total weight'%(newweight)
	totweight+=newweight
	f.Close()
print 'TOTAL SUM OF EVENT WEIGHTS = '+str(totweight)
#set the global histograms (normalized)
if options.name.lower().find('mcatnlo')!=-1 :
	total_cstar_vs_beta_qqbar_histo.Scale(1./totweight) 
	total_cstar_vs_beta_gg_histo.Scale(1./totweight)
csvb_qq_global_hist=total_cstar_vs_beta_qqbar_histo
csvb_gg_global_hist=total_cstar_vs_beta_gg_histo
renormalization_dict = make_renormalization_dict(options.name,options.alphaOverride,options.epsilonOverride,total_mu_R_sf_up_histo,total_mu_R_sf_down_histo,total_mu_F_sf_up_histo,total_mu_F_sf_down_histo,
												 total_scale_comb_sf_up_histo,total_scale_comb_sf_down_histo,total_pdf_alphas_sf_histo,total_pdf_alphas_sf_up_histo,total_pdf_alphas_sf_down_histo,total_top_pt_reweight_v1_histo)
#Get the total number of events
ntotalevents = chain.GetEntries()
print 'number of total events = %d'%(ntotalevents) 
#find how many events should be in this tree 
nanalysisevents = abs(int(ntotalevents/options.n_jobs))
if options.max_events!=-1 :
	nanalysisevents = min(options.max_events, nanalysisevents)
if options.n_jobs>1 :
	#adjust so they don't all pile up too badly in the last job but every job runs on some events
	while ntotalevents-options.n_jobs*nanalysisevents > nanalysisevents and (options.n_jobs-1)*(nanalysisevents+1)<ntotalevents :
		nanalysisevents+=1
#copy the subset tree to be analyzed 
garbageFile.cd()
analysisTree = chain.CopyTree('','',nanalysisevents,options.i_job*nanalysisevents); analysisTree.SetDirectory(garbageFile) 
nanalysisevents = analysisTree.GetEntries() 
print 'number of analyzed events for this job = %d'%(nanalysisevents) 
#Set filename for analyzer from sample name
filename = options.name
if options.jec != 'nominal' :
	filename+='_'+options.jec
if options.n_jobs>1 :
	filename+='_'+str(options.i_job)
filename+='_tree.root'
#Initialize analyzer
analyzer = Reconstructor(filename, analysisTree, data, options.xSec, options.kFac, options.jec, options.on_grid, total_pileup_histo, totweight, renormalization_dict) 

#Counter 
count = 0

#If there is no max_events argument, set it to the total number of analyzed events
maxEvents = options.max_events if options.max_events!=-1 else nanalysisevents
setupDoneTime = time(); setupTime = timedelta(seconds=setupDoneTime-startTime)
print 'Done Setting up; time taken: %02d:%02d:%02d'%(setupTime.seconds/3600,(setupTime.seconds%3600)/60,(setupTime.seconds%60))

##########								Main Event Loop								##########
print 'Files opened, starting event loop'
for event in range(nanalysisevents) : 
	count+=1
	#check the max events 
	if count == maxEvents+1 :
		print 'Processed event number '+str(count-1)+', exiting'
		break
	#print progress
	if count % options.print_every == 0 or count == 1:
			timeSinceSetup = timedelta(seconds=time()-setupDoneTime)
			percentDone = float(count) / float(min(nanalysisevents,maxEvents)) * 100.0
			timeLeft = timedelta(seconds=(100.-percentDone)*timeSinceSetup.total_seconds()/percentDone)
			print ( 'Count at %d of %d, (%.4f%% complete; elapsed = %02d:%02d:%02d, remaining = %02d:%02d:%02d)'
				   %(count,min(nanalysisevents,maxEvents),percentDone,timeSinceSetup.seconds/3600,(timeSinceSetup.seconds%3600)/60,(timeSinceSetup.seconds%60),timeLeft.seconds/3600,(timeLeft.seconds%3600)/60,(timeLeft.seconds%60))  )
	#a couple vetoes for bugged events
	if options.name=='DYJets_M-50_HT-2500toInf' :
		if options.n_jobs==40 :
			if (options.i_job==11 and count==4839) or (options.i_job==32 and count==8219) or (options.i_job==34 and count==4157) :
				continue
		elif options.n_jobs==400 :
			if (options.i_job==115  and count==293) or (options.i_job==329 and count==146) or (options.i_job==344 and count==747) :
				continue
	elif options.name=='WJets_HT-2500toInf' :
		if options.n_jobs==10 :
			if (options.i_job==4 and count==17253) or (options.i_job==5 and count==16924) :
				continue
		elif options.n_jobs==100 :
			if (options.i_job==47 and count==1829) or (options.i_job==57 and count==1508) :
				continue
	elif options.name=='QCD_HT-2000toInf' :
		if options.n_jobs==50 :
			if (options.i_job==20 and count==4468) or (options.i_job==30 and count==19666) or (options.i_job==47 and count==21334) :
				continue
	#analyze event and add to TTree
	analyzer.analyze(event)
	#reset analyzer
	analyzer.reset()
#clean up after yourself
del analyzer
garbageFile.Close()
os.system('rm -rf '+garbageFileName)
