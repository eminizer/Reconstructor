#Imports
import os
import sys
from DataFormats.FWLite import Events, Handle
from ROOT import TMinuit, Long, Double, TChain, TFile, TH1D, TH2D, TF2
from array import array
from math import *
from optparse import OptionParser
from reconstructor import Reconstructor

#global 2D histograms needed for alpha and epsilon fits
csvb_qq_global_hist = TH2D(); csvb_gg_global_hist = TH2D()

#minuit fitting function for alpha
def alpha_fcn(npar, deriv, f, par, flag) :
	alpha = par[0]
	lnL=0.
	#loop over the costheta bins
	for i in range(1,csvb_qq_global_hist.GetYaxis().GetNbins()+1) :
		cstar = csvb_qq_global_hist.GetYaxis().GetBinCenter(i)
		#loop over the beta bins
		nbetabins = csvb_qq_global_hist.GetXaxis().GetNbins()
		for j in range(1,nbetabins+1) :
			beta = csvb_qq_global_hist.GetXaxis().GetBinCenter(j)
			#increment the values
			globalBin = csvb_qq_global_hist.FindFixBin(beta,cstar)
			cont=csvb_qq_global_hist.GetBinContent(globalBin)
			f_L = (2.+(beta*beta)*(cstar*cstar)-(beta*beta)+alpha*(1.-(beta*beta)*(cstar*cstar)))/(2.*(2-2.*(beta*beta)/3.+alpha*(1.-(beta*beta)/3.)))
			lnL+=-2.0*cont*log(f_L)
	print '	lnL = %.8f, alpha = %.8f'%(lnL,alpha) #DEBUG
	f[0] = lnL

#minuit fitting function for epsilon
def epsilon_fcn(npar, deriv, f, par, flag) :
	epsilon = par[0]
	lnL = 0.0
	#loop over the costheta bins
	for i in range(1,csvb_gg_global_hist.GetYaxis().GetNbins()+1) :
		cstar = csvb_gg_global_hist.GetYaxis().GetBinCenter(i)
		#loop over the beta bins
		nbetabins = csvb_gg_global_hist.GetXaxis().GetNbins()
		for j in range(1,nbetabins+1) :
			beta = csvb_gg_global_hist.GetXaxis().GetBinCenter(j)
			#increment the values
			globalBin = csvb_gg_global_hist.FindFixBin(beta,cstar)
			cont=csvb_gg_global_hist.GetBinContent(globalBin)
			beta2 = beta*beta
			cstar2 = cstar*cstar
			bc2 = beta2*cstar2
			omb2 = 1. - beta2
			beta4 = beta2*beta2
			rnorm = ( (epsilon*(34.*beta4-100.*beta2+98.)+2.*beta4-36.*beta2+66.)*atanh(beta)/beta
						-9./5.*epsilon*beta4+6.*beta2*(epsilon*(beta2-2.38889)-0.5)-16.*(1.+epsilon)*(beta4-2.*beta2+1.)/(1.-beta2)
						-epsilon*(18.*beta4-68.*beta2+82.)+18.*beta2-43. )
			f_L = (7.+9.*bc2)/(1.-bc2)*((1.+bc2)/2.+omb2*beta2*(1.-cstar2)/(1.-bc2))*(1.+epsilon*bc2)/rnorm
			lnL+=-2.0*cont*log(f_L)
	print '	lnL = %.8f, epsilon = %.8f'%(lnL,epsilon) #DEBUG
	f[0] = lnL

#Helper function for making the renormalization dictionary
def make_renormalization_dict(muRup,muRdown,muFup,muFdown,scup,scdown,pdfas,pdfasup,pdfasdown) :
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
	#now we've got to do the fits for the "deweighting" alpha and epsilon values
	print 'fitting for alpha... '
	#Define the Minuit instance we'll use
	alpha_minuit = TMinuit(1); alpha_minuit.SetFCN(alpha_fcn)
	#miscellaneous minuit stuff
	ierflag = Long(1); arglist = array('d',[100000.])
	#add parameter
	alpha_minuit.mnparm(0,'alpha',0.1,0.5,0.,0.,ierflag)
	#minimize
	alpha_minuit.mnexcm('MIGRAD',arglist,1,ierflag)
	#get back the fitted alpha value
	fitted_alpha=Double(0.0); fitted_alpha_err=Double(0.0)
	alpha_minuit.GetParameter(0,fitted_alpha,fitted_alpha_err)	
	returndict['alpha']=fitted_alpha
	#Do the same thing except for epsilon
	print 'fitting for epsilon...'
	epsilon_minuit = TMinuit(1); epsilon_minuit.SetFCN(epsilon_fcn)
	ierflag = Long(1); arglist = array('d',[100000.])
	epsilon_minuit.mnparm(0,'epsilon',0.1,0.5,0.,0.,ierflag)
	epsilon_minuit.mnexcm('MIGRAD',arglist,1,ierflag)
	fitted_epsilon=Double(0.0); fitted_epsilon_err=Double(0.0)
	epsilon_minuit.GetParameter(0,fitted_epsilon,fitted_epsilon_err)
	returndict['epsilon']=fitted_epsilon
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
	print 'Fitted alpha value = %.4f +/- %.4f'%(fitted_alpha, fitted_alpha_err)
	print 'Fitted epsilon value = %.4f +/- %.4f'%(fitted_epsilon, fitted_epsilon_err)
	#return the dictionary
	return returndict

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
parser.add_option('--xSec', 		 type='float', action='store', 			  	dest='xSec', 		    
	help="Cross section of sample's process")
parser.add_option('--JES', type='string', action='store', default='nominal',  dest='JES',  
	help='JES systematics: shift JES up or down (default is nominal)')
parser.add_option('--JER', type='string', action='store', default='nominal',  dest='JER',  
	help='JER systematics: shift JER up or down (default is nominal)')
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
#Read files in line by line and get the tree, pileup,cstar vs. beta,pdf_alphas/F_up/down,scale_comb_sf_up/down,pdf_alphas_sf_up/down histograms, and total weight value from each
total_pileup_histo               = TH1D('total_pileup_histo','total MC pileup; pileup; events',100,0.,100.); total_pileup_histo.SetDirectory(0)
total_cstar_vs_beta_qqbar_histo  = TH2D('total_cstar_vs_beta_qqbar_histo','MC truth c* vs. #beta (q#bar{q} events); #beta; c*; events',10,0.,1.,20,0.,1.); total_cstar_vs_beta_qqbar_histo.SetDirectory(0)
total_cstar_vs_beta_gg_histo     = TH2D('total_cstar_vs_beta_gg_histo','MC truth c* vs. #beta (qg/gg events); #beta; c*; events',10,0.,1.,20,0.,1.); total_cstar_vs_beta_gg_histo.SetDirectory(0)
total_mu_R_sf_up_histo           = TH1D('total_mu_R_sf_up_histo','total #mu_{R} sf (up); #mu_{R} sf up; events',4,-1.,3.); total_mu_R_sf_up_histo.SetDirectory(0);
total_mu_R_sf_down_histo         = TH1D('total_mu_R_sf_down_histo','total #mu_{R} sf (down); #mu_{R} sf down; events',4,-1.,3.); total_mu_R_sf_down_histo.SetDirectory(0);
total_mu_F_sf_up_histo           = TH1D('total_mu_F_sf_up_histo','total #mu_{F} sf (up); #mu_{F} sf up; events',4,-1.,3.); total_mu_F_sf_up_histo.SetDirectory(0);
total_mu_F_sf_down_histo         = TH1D('total_mu_F_sf_down_histo','total #mu_{F} sf (down); #mu_{F} sf down; events',4,-1.,3.); total_mu_F_sf_down_histo.SetDirectory(0);
total_scale_comb_sf_up_histo     = TH1D('total_scale_comb_sf_up_histo','total comb. scale sf (up); comb. scale sf up; events',4,-1.,3.); total_scale_comb_sf_up_histo.SetDirectory(0);
total_scale_comb_sf_down_histo   = TH1D('total_scale_comb_sf_down_histo','total comb. scale sf (down); comb. scale sf down; events',4,-1.,3.); total_scale_comb_sf_down_histo.SetDirectory(0);
total_pdf_alphas_sf_histo        = TH1D('total_pdf_alphas_sf_histo','total pdf/alpha_{s} sf; pdf/alpha_{s} sf; events',4,-1.,3.); total_pdf_alphas_sf_histo.SetDirectory(0);
total_pdf_alphas_sf_up_histo     = TH1D('total_pdf_alphas_sf_up_histo','total pdf/alpha_{s} sf (up); pdf/alpha_{s} sf up; events',4,-1.,3.); total_pdf_alphas_sf_up_histo.SetDirectory(0);
total_pdf_alphas_sf_down_histo   = TH1D('total_pdf_alphas_sf_down_histo','total pdf/alpha_{s} sf (down); pdf/alpha_{s} sf down; events',4,-1.,3.); total_pdf_alphas_sf_down_histo.SetDirectory(0);
totweight = 0.
for input_file in input_files_list :
	print '	'+input_file.rstrip()+''
	chain.AddFile(input_file.rstrip()+'/'+options.ttree_dir_name+'/'+options.ttree_name)
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
	histo=f.Get('EventCounter/totweight')
	newweight=histo.GetBinContent(1)
	print '		Added %.2f to total weight'%(newweight)
	totweight+=newweight
print 'TOTAL SUM OF EVENT WEIGHTS = '+str(totweight)
#set the global histograms (normalized)
total_cstar_vs_beta_qqbar_histo.Sumw2(); total_cstar_vs_beta_qqbar_histo.Scale(1./total_cstar_vs_beta_qqbar_histo.Integral()) 
total_cstar_vs_beta_gg_histo.Sumw2(); total_cstar_vs_beta_gg_histo.Scale(1./total_cstar_vs_beta_gg_histo.Integral())
csvb_qq_global_hist=total_cstar_vs_beta_qqbar_histo
csvb_gg_global_hist=total_cstar_vs_beta_gg_histo
renormalization_dict = make_renormalization_dict(total_mu_R_sf_up_histo,total_mu_R_sf_down_histo,total_mu_F_sf_up_histo,total_mu_F_sf_down_histo,
												 total_scale_comb_sf_up_histo,total_scale_comb_sf_down_histo,total_pdf_alphas_sf_histo,total_pdf_alphas_sf_up_histo,total_pdf_alphas_sf_down_histo)
ntotalevents = chain.GetEntries()
#Set filename for analyzer from sample name
filename = options.name
if options.JES.lower() != 'nominal' :
	filename+='_JES_'+options.JES.lower()
if options.JER.lower() != 'nominal' :
	filename+='_JER_'+options.JER.lower()
if options.n_jobs>1 :
	filename+='_'+str(options.i_job)
filename+='_tree.root'
#Initialize analyzer
data=False
if options.name.lower().find('singlemu')!=-1 or options.name.lower().find('singleel')!=-1 :
	data=True
analyzer = Reconstructor(filename, chain, data, options.xSec, options.JES.lower(), options.JER.lower(), options.on_grid, total_pileup_histo, totweight, renormalization_dict)

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
