#Imports
import os
import sys
from ROOT import *
from array import array
from math import *
from optparse import OptionParser

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
def get_alpha_epsilon() :
	#do the fits for the "deweighting" alpha and epsilon values
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
	#Do the same thing except for epsilon
	print 'fitting for epsilon...'
	epsilon_minuit = TMinuit(1); epsilon_minuit.SetFCN(epsilon_fcn)
	ierflag = Long(1); arglist = array('d',[100000.])
	epsilon_minuit.mnparm(0,'epsilon',0.1,0.5,0.,0.,ierflag)
	epsilon_minuit.mnexcm('MIGRAD',arglist,1,ierflag)
	fitted_epsilon=Double(0.0); fitted_epsilon_err=Double(0.0)
	epsilon_minuit.GetParameter(0,fitted_epsilon,fitted_epsilon_err)
	#print the values
	print 'Fitted alpha value = %.4f +/- %.4f'%(fitted_alpha, fitted_alpha_err)
	print 'Fitted epsilon value = %.4f +/- %.4f'%(fitted_epsilon, fitted_epsilon_err)
	#return the dictionary
	return fitted_alpha, fitted_alpha_err, fitted_epsilon, fitted_epsilon_err

##########								Parser Options								##########

parser = OptionParser()
#Run options
parser.add_option('--input', 	  type='string', action='store', default='input', dest='input',	   	  
	help='Path to input file holding list of files to run on (default: "input")')
(options, args) = parser.parse_args()

##########							Set Up Event Loop								##########

print 'Opening files . . .'  
#Build path to input file
input_files_list = options.input
if not options.input.endswith('.txt') :
	input_files_list += '.txt'
print 'Using input file '+input_files_list+''
#open input file in read only mode 
input_files_list = open(input_files_list,'r')
#Set up the garbage output file and the chain
outFileName = 'alpha_epsilon_plots.root' 
outFile = TFile(outFileName,'recreate') 
chain = TChain('B2GTree')
print 'Getting these files: '
#Read files in line by line and get the tree, pileup,cstar vs. beta,pdf_alphas/F_up/down,scale_comb_sf_up/down,pdf_alphas_sf_up/down histograms, and total weight value from each
total_cstar_vs_beta_qqbar_histo  = TH2D('total_cstar_vs_beta_qqbar_histo','MC truth c* vs. #beta (q#bar{q} events); #beta; c*; events',10,0.,1.,20,-1.,1.); total_cstar_vs_beta_qqbar_histo.SetDirectory(0)
total_cstar_vs_beta_gg_histo     = TH2D('total_cstar_vs_beta_gg_histo','MC truth c* vs. #beta (qg/gg events); #beta; c*; events',10,0.,1.,20,-1.,1.); total_cstar_vs_beta_gg_histo.SetDirectory(0)
totweight = 0.
for input_file in input_files_list :
	print '	'+input_file.rstrip()+''
	chain.AddFile(input_file.rstrip()+'/B2GTTreeMaker/B2GTree')
	f = TFile.Open(input_file.rstrip()) 
	csvbqq_histo = f.Get('EventCounter/cstar_vs_beta_qqbar')
	total_cstar_vs_beta_qqbar_histo.Add(csvbqq_histo)
	csvbgg_histo = f.Get('EventCounter/cstar_vs_beta_gg')
	total_cstar_vs_beta_gg_histo.Add(csvbgg_histo)
	histo=f.Get('EventCounter/totweight')
	newweight=histo.GetBinContent(1)
	print '		Added %.2f to total weight'%(newweight)
	totweight+=newweight
	f.Close()
print 'TOTAL SUM OF EVENT WEIGHTS = '+str(totweight)
#set the global histograms (normalized)
total_cstar_vs_beta_qqbar_histo.Sumw2(); #total_cstar_vs_beta_qqbar_histo.Scale(1./total_cstar_vs_beta_qqbar_histo.Integral()) 
total_cstar_vs_beta_gg_histo.Sumw2(); #total_cstar_vs_beta_gg_histo.Scale(1./total_cstar_vs_beta_gg_histo.Integral())
csvb_qq_global_hist=total_cstar_vs_beta_qqbar_histo
csvb_gg_global_hist=total_cstar_vs_beta_gg_histo
fit_alpha, fit_alpha_err, fit_epsilon, fit_epsilon_err = get_alpha_epsilon()
#alpha fitting function
f_alpha_num = '(2.+(x*x)*(y*y)-(x*x)+%f*(1.-(x*x)*(y*y)))'%(fit_alpha)
f_alpha_den = '(2.*(2-2.*(x*x)/3.+(%f)*(1.-(x*x)/3.)))'%(fit_alpha)
f_alpha = TF2('f_alpha',f_alpha_num+'/'+f_alpha_den,0.,1.,0.,1.)
#epsilon fitting function
rnorm = '(%f*(34.*x*x*x*x-100.*x*x+98.)+2.*x*x*x*x-36.*x*x+66.)'%(fit_epsilon)
rnorm+= '*atanh(x)/x-9./5.*%f*x*x*x*x+6.*x*x'%(fit_epsilon)
rnorm+= '*(%f*(x*x-2.38889)-0.5)-16.*(1.+%f)'%(fit_epsilon,fit_epsilon)
rnorm+= '*(x*x*x*x-2.*x*x+1.)/(1.-x*x)-%f'%(fit_epsilon)
rnorm+= '*(18.*x*x*x*x-68.*x*x+82.)+18.*x*x-43.'
f_epsilon = TF2('f_epsilon','(7.+9.*x*x*y*y)/(1.-x*x*y*y)*((1.+x*x*y*y)/2.+(1.-x*x)*x*x*(1.-y*y)/(1.-x*x*y*y))*(1.+%f*x*x*y*y)/%s'%(fit_epsilon,rnorm),0.,1.,0.,1.)
#plot alpha 2D plot and total function
canv_alpha_hist = TCanvas('canv_alpha_hist','canv_alpha_hist',1100,900)
csvb_qq_global_hist.Draw('COLZ')
canv_alpha_func = TCanvas('canv_alpha_func','canv_alpha_func',1100,900)
f_alpha.Draw('COLZ')
#plot epsilon 2D plot and total function
canv_epsilon_hist = TCanvas('canv_epsilon_hist','canv_epsilon_hist',1100,900)
csvb_gg_global_hist.Draw('COLZ')
canv_epsilon_func = TCanvas('canv_epsilon_func','canv_epsilon_func',1100,900)
f_epsilon.Draw('COLZ')
#plot alpha projected plots
canv_alpha_proj = TCanvas('canv_alpha_proj','canv_alpha_proj',1100,900)
alpha_proj_histo = csvb_qq_global_hist.ProjectionY()
alpha_proj_histo.Scale(1./alpha_proj_histo.Integral())
alpha_proj_histo.GetYaxis().SetRangeUser(0.,0.07)
npoints_alpha = 10*csvb_qq_global_hist.GetYaxis().GetNbins()
xpoints_alpha = array('d',npoints_alpha*[0.])
ypoints_alpha = array('d',npoints_alpha*[0.])
for i in range(npoints_alpha) :
	pointspacing = 1./(npoints_alpha)
	cstar_low = i*pointspacing
	cstar = (i+0.5)*pointspacing
	cstar_hi  = (i+1)*pointspacing
	xpoints_alpha[i]=cstar
	ypoints_alpha[i]=f_alpha.Integral(0.,1.,cstar_low,cstar_hi)
sumypoints_alpha = 0.
for i in range(npoints_alpha) :
	sumypoints_alpha+=ypoints_alpha[i]
for i in range(npoints_alpha) :
	fac = npoints_alpha/csvb_qq_global_hist.GetYaxis().GetNbins()
	ypoints_alpha[i]  = fac*ypoints_alpha[i]/sumypoints_alpha
alpha_graph = TGraph(npoints_alpha,xpoints_alpha,ypoints_alpha)
alpha_graph.SetMarkerStyle(20)
alpha_graph.SetLineWidth(3)
alpha_graph.SetLineColor(kRed)
alpha_proj_histo.SetMarkerStyle(20)
alpha_proj_histo.Draw('E0')
alpha_graph.Draw('L SAME')
leg_alpha = TLegend(0.1,0.7,0.48,0.9)
leg_alpha.AddEntry(alpha_proj_histo,"MC Events","PE0")
leg_alpha.AddEntry(alpha_graph,"Fit (#alpha="+str(fitted_alpha)+")","L")
leg_alpha.Draw('SAME')
#plot epsilon projected plots
canv_epsilon_proj = TCanvas('canv_epsilon_proj','canv_epsilon_proj',1100,900)
epsilon_proj_histo = csvb_gg_global_hist.ProjectionY()
epsilon_proj_histo.Scale(1./epsilon_proj_histo.Integral())
epsilon_proj_histo.GetYaxis().SetRangeUser(0.,0.15)
npoints_epsilon = 10*csvb_gg_global_hist.GetYaxis().GetNbins()
xpoints_epsilon = array('d',npoints_epsilon*[0.])
ypoints_epsilon = array('d',npoints_epsilon*[0.])
for i in range(npoints_epsilon) :
	pointspacing = 1./(npoints_epsilon)
	cstar_low = i*pointspacing
	cstar = (i+0.5)*pointspacing
	cstar_hi  = (i+1)*pointspacing
	xpoints_epsilon[i]=cstar
	ypoints_epsilon[i]=f_epsilon.Integral(0.,1.,cstar_low,cstar_hi)
sumypoints_epsilon = 0.
for i in range(npoints_epsilon) :
	sumypoints_epsilon+=ypoints_epsilon[i]
for i in range(npoints_epsilon) :
	fac = npoints_epsilon/csvb_gg_global_hist.GetYaxis().GetNbins()
	ypoints_epsilon[i]  = fac*ypoints_epsilon[i]/sumypoints_epsilon
epsilon_graph = TGraph(npoints_epsilon,xpoints_epsilon,ypoints_epsilon)
epsilon_graph.SetMarkerStyle(20)
epsilon_graph.SetLineWidth(3)
epsilon_graph.SetLineColor(kRed)
epsilon_proj_histo.SetMarkerStyle(20)
epsilon_proj_histo.Draw('E0')
epsilon_graph.Draw('L SAME')
leg_epsilon = TLegend(0.1,0.7,0.48,0.9)
leg_epsilon.AddEntry(epsilon_proj_histo,"MC Events","PE0")
leg_epsilon.AddEntry(epsilon_graph,"Fit (#epsilon="+str(fitted_epsilon)+")","L")
leg_epsilon.Draw('SAME')

#write plots to file
outFile.cd()
csvb_qq_global_hist.Write()
csvb_gg_global_hist.Write()
f_alpha.Write()
f_epsilon.Write()
canv_alpha_hist.Write()
canv_alpha_func.Write()
canv_epsilon_hist.Write()
canv_epsilon_func.Write()
alpha_proj_histo.Write()
alpha_graph.Write()
canv_alpha_proj.Write()
epsilon_proj_histo.Write()
epsilon_graph.Write()
canv_epsilon_proj.Write()
outFile.Write()
outFile.Close()


