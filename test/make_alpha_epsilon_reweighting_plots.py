#Imports
import os
import sys
from ROOT import *
from array import array
from math import *
from optparse import OptionParser

#global 2D histograms needed for alpha and epsilon fits
csvb_qq_global_hist = TH2D(); csvb_gg_global_hist = TH2D()
qq_global_hist_projection = TH1D(); gg_global_hist_projection = TH1D()
BETA = [-1.]
#hope you put enough colors here lol
colors = [kBlue,kGreen+2,kRed+2,kMagenta,kCyan,kYellow+2,kOrange+7,kViolet-6,kGreen+4,kRed-7]

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
	print '	lnL = %.8f, alpha = %.8f'%(lnL,alpha) #DEBUG
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
	print '	lnL = %.8f, epsilon = %.8f'%(lnL,epsilon) #DEBUG
	f[0] = lnL

#Helper function for making the renormalization dictionary
def get_alpha_epsilon_list() :
	global qq_global_hist_projection
	global gg_global_hist_projection
	alphalist = []; epsilonlist = []
	#do the fits for the "deweighting" alpha and epsilon values
	print 'fitting for alpha and epsilon values... '
	#Define the Minuit instances we'll use
	alpha_minuit = TMinuit(1); alpha_minuit.SetFCN(alpha_fcn)
	epsilon_minuit = TMinuit(1); epsilon_minuit.SetFCN(epsilon_fcn)
	#miscellaneous minuit stuff
	ierflag = Long(1); arglist = array('d',[100000.])
	for i in range(1,csvb_qq_global_hist.GetXaxis().GetNbins()+1) :
		BETA[0] = csvb_qq_global_hist.GetXaxis().GetBinCenter(i)
		beta_low = csvb_qq_global_hist.GetXaxis().GetBinLowEdge(i)
		beta_hi = beta_low+csvb_qq_global_hist.GetXaxis().GetBinWidth(i)
		qq_global_hist_projection = symmetrize(csvb_qq_global_hist.ProjectionY('qq_proj_'+str(i),i,i))
		#add parameter
		alpha_minuit.mnparm(0,'alpha',0.1,0.5,0.,0.,ierflag)
		#minimize
		alpha_minuit.mnexcm('MIGRAD',arglist,1,ierflag)
		#get back the fitted alpha value
		fitted_alpha=Double(0.0); fitted_alpha_err=Double(0.0)
		alpha_minuit.GetParameter(0,fitted_alpha,fitted_alpha_err)	
		alphalist.append((float(BETA[0]),float(fitted_alpha),float(fitted_alpha_err),float(beta_low),float(beta_hi)))
	#Do the same thing except for epsilon
	for i in range(1,csvb_gg_global_hist.GetXaxis().GetNbins()+1) :
		BETA[0] = csvb_gg_global_hist.GetXaxis().GetBinCenter(i)
		beta_low = csvb_gg_global_hist.GetXaxis().GetBinLowEdge(i)
		beta_hi = beta_low+csvb_gg_global_hist.GetXaxis().GetBinWidth(i)
		gg_global_hist_projection = symmetrize(csvb_gg_global_hist.ProjectionY('gg_proj_'+str(i),i,i))	
		epsilon_minuit.mnparm(0,'epsilon',0.1,0.5,0.,0.,ierflag)
		epsilon_minuit.mnexcm('MIGRAD',arglist,1,ierflag)
		fitted_epsilon=Double(0.0); fitted_epsilon_err=Double(0.0)
		epsilon_minuit.GetParameter(0,fitted_epsilon,fitted_epsilon_err)
		epsilonlist.append((float(BETA[0]),float(fitted_epsilon),float(fitted_epsilon_err),float(beta_low),float(beta_hi)))
	#print the values
	for item in alphalist :
		print 'alpha=%.4f +/- %.4f in bin with center at beta=%.3f'%(item[1],item[2],item[0])
	for item in epsilonlist :
		print 'epsilon=%.4f +/- %.4f in bin with center at beta=%.3f'%(item[1],item[2],item[0])
	#return the lists
	return alphalist, epsilonlist

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
#Set up the output file
outFileName = 'alpha_epsilon_plots.root' 
outFile = TFile(outFileName,'recreate') 
print 'Getting these files: '
#Read files in line by line and get the tree, pileup,cstar vs. beta,pdf_alphas/F_up/down,scale_comb_sf_up/down,pdf_alphas_sf_up/down histograms, and total weight value from each
total_cstar_vs_beta_qqbar_histo  = TH2D('total_cstar_vs_beta_qqbar_histo','MC truth c* vs. #beta (q#bar{q} events); #beta; c*; events',10,0.,1.,20,-1.,1.); total_cstar_vs_beta_qqbar_histo.SetDirectory(0)
total_cstar_vs_beta_gg_histo     = TH2D('total_cstar_vs_beta_gg_histo','MC truth c* vs. #beta (qg/gg events); #beta; c*; events',10,0.,1.,20,-1.,1.); total_cstar_vs_beta_gg_histo.SetDirectory(0)
totweight = 0.
for input_file in input_files_list :
	print '	'+input_file.rstrip()+''
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
#set the global histograms
total_cstar_vs_beta_qqbar_histo.Scale(1./totweight); total_cstar_vs_beta_gg_histo.Scale(1./totweight)
csvb_qq_global_hist=total_cstar_vs_beta_qqbar_histo
csvb_gg_global_hist=total_cstar_vs_beta_gg_histo
alphalist, epsilonlist = get_alpha_epsilon_list()

#setup for plots
total_alpha_canv = TCanvas('total_alpha_canv','total_alpha_canv',1100,900)
total_epsilon_canv = TCanvas('total_epsilon_canv','total_epsilon_canv',1100,900)
ind_alpha_canvs = []
ind_epsilon_canvs = []
ind_alpha_histos = []
ind_alpha_funcs = []
ind_epsilon_histos = []
ind_epsilon_funcs = []
leg_alpha = TLegend(0.1,0.7,0.48,0.9)
leg_epsilon = TLegend(0.1,0.7,0.48,0.9)

#loop over the different bins for alpha and add functions and histograms
for i in range(len(alphalist)) :
	#get the beta value
	beta = alphalist[i][0]
	#add the histogram
	ind_alpha_histos.append(symmetrize(csvb_qq_global_hist.ProjectionY('qq_proj_'+str(i),i,i+1)))
	ind_alpha_histos[i].SetTitle('c* for qqbar events (beta bin %d); c*'%(i))
	#add the function
	thisalpha = alphalist[i][1]
	f_alpha_num = '(2.+(%f*%f)*(x*x)-(%f*%f)+%f*(1.-(%f*%f)*(x*x)))'%(beta,beta,beta,beta,thisalpha,beta,beta)
	f_alpha_den = '(2.*(2-2.*(%f*%f)/3.+(%f)*(1.-(%f*%f)/3.)))'%(beta,beta,thisalpha,beta,beta)
	fake_f = TF1('f_alpha_fake_'+str(i),'(%s)/(%s)'%(f_alpha_num,f_alpha_den),-1.,1.)
	ind_alpha_funcs.append(TF1('f_alpha_'+str(i),'(%s)/(%f*%f*(%s))'%(f_alpha_num,len(alphalist),fake_f.Integral(-1.,1.),f_alpha_den),-1.,1.))
	#append a new canvas
	ind_alpha_canvs.append(TCanvas('alpha_canv_'+str(i),'alpha_canv_'+str(i),1100,900))
#same, but for epsilon
for i in range(len(epsilonlist)) :
	#get the beta value
	beta = epsilonlist[i][0]
	#add the histogram
	ind_epsilon_histos.append(symmetrize(csvb_gg_global_hist.ProjectionY('gg_proj_'+str(i),i,i+1)))
	ind_epsilon_histos[i].SetTitle('c* for gg/qg events (beta bin %d); c*'%(i))
	#add the function
	thisepsilon = epsilonlist[i][1]
	rnorm = '((%f*(34.*%f*%f*%f*%f-100.*%f*%f+98.)+2.*%f*%f*%f*%f-36.*%f*%f+66.)'%(thisepsilon,beta,beta,beta,beta,beta,beta,beta,beta,beta,beta,beta,beta)
	rnorm+= '*atanh(%f)/%f-9./5.*%f*%f*%f*%f*%f+6.*%f*%f'%(beta,beta,thisepsilon,beta,beta,beta,beta,beta,beta)
	rnorm+= '*(%f*(%f*%f-2.38889)-0.5)-16.*(1.+%f)'%(thisepsilon,beta,beta,thisepsilon)
	rnorm+= '*(%f*%f*%f*%f-2.*%f*%f+1.)/(1.-%f*%f)-%f'%(beta,beta,beta,beta,beta,beta,beta,beta,thisepsilon)
	rnorm+= '*(18.*%f*%f*%f*%f-68.*%f*%f+82.)+18.*%f*%f-43.)'%(beta,beta,beta,beta,beta,beta,beta,beta)
	numerator = '(7.+9.*%f*%f*x*x)/(1.-%f*%f*x*x)*((1.+%f*%f*x*x)/2.+(1.-%f*%f)*%f*%f*(1.-x*x)/(1.-%f*%f*x*x))*(1.+%f*%f*%f*x*x)'%(beta,beta,beta,beta,beta,beta,beta,beta,beta,beta,beta,beta,thisepsilon,beta,beta)
	fake_f = TF1('f_epsilon_'+str(i),'(%s)/(%s)'%(numerator,rnorm),-1.,1.)
	ind_epsilon_funcs.append(TF1('f_epsilon_'+str(i),'(%s)/(%f*%f*(%s))'%(numerator,len(epsilonlist),fake_f.Integral(-1.,1.),rnorm),-1.,1.))
	#append a new canvas
	ind_epsilon_canvs.append(TCanvas('epsilon_canv_'+str(i),'epsilon_canv_'+str(i),1100,900))

#write out raw histograms and functions, setting properties after
outFile.cd()
for i in range(len(ind_alpha_histos)) :
	ind_alpha_histos[i].Write()
	ind_alpha_funcs[i].Write()
	ind_alpha_funcs[i].SetLineWidth(4);
	ind_alpha_funcs[i].SetLineColor(colors[i])
	ind_alpha_histos[i].SetMarkerStyle(20)
	ind_alpha_histos[i].SetMarkerColor(colors[i])
	ind_alpha_histos[i].SetLineColor(colors[i])
	ind_alpha_histos[i].SetLineWidth(2)
	ind_alpha_histos[i].GetYaxis().SetRangeUser(0.,0.1)
for i in range(len(ind_epsilon_histos)) :
	ind_epsilon_histos[i].Write()
	ind_epsilon_funcs[i].Write()
	ind_epsilon_funcs[i].SetLineWidth(4);
	ind_epsilon_funcs[i].SetLineColor(colors[i])
	ind_epsilon_histos[i].SetMarkerStyle(20)
	ind_epsilon_histos[i].SetMarkerColor(colors[i])
	ind_epsilon_histos[i].SetLineColor(colors[i])
	ind_epsilon_histos[i].SetLineWidth(2)
	ind_epsilon_histos[i].GetYaxis().SetRangeUser(0.,0.1)

#renormalize histograms and plot on canvases with their functions, also add to the total canvas
for i in range(len(ind_alpha_histos)) :
	ind_alpha_histos[i].Scale(1./ind_alpha_histos[i].Integral())
	ind_alpha_canvs[i].cd()
	ind_alpha_histos[i].Draw()
	ind_alpha_funcs[i].Draw('L SAME')
	outFile.cd()
	ind_alpha_canvs[i].Write()
	total_alpha_canv.cd()
	ind_alpha_histos[i].Draw() if i==0 else ind_alpha_histos[i].Draw('SAME')
	ind_alpha_funcs[i].Draw('L SAME')
	leg_alpha.AddEntry(ind_alpha_histos[i],'%.2f < #beta < %.2f, #alpha=%.2f'%(alphalist[i][3],alphalist[i][4],alphalist[i][1]),'L')
for i in range(len(ind_epsilon_histos)) :
	ind_epsilon_histos[i].Scale(1./ind_epsilon_histos[i].Integral())
	ind_epsilon_canvs[i].cd()
	ind_epsilon_histos[i].Draw()
	ind_epsilon_funcs[i].Draw('L SAME')
	outFile.cd()
	ind_epsilon_canvs[i].Write()
	total_epsilon_canv.cd()
	ind_epsilon_histos[i].Draw() if i==0 else ind_epsilon_histos[i].Draw('SAME')
	ind_epsilon_funcs[i].Draw('L SAME')
	leg_epsilon.AddEntry(ind_epsilon_histos[i],'%.2f < #beta < %.2f, #epsilon=%.2f'%(epsilonlist[i][3],epsilonlist[i][4],epsilonlist[i][1]),'L')

#plot the legends on the total canvases
total_alpha_canv.cd()
leg_alpha.Draw('SAME')
total_epsilon_canv.cd()
leg_epsilon.Draw('SAME')

#write out the total canvases
outFile.cd()
total_alpha_canv.Write()
total_epsilon_canv.Write()

##plot alpha 2D plot and total function
#canv_alpha_hist = TCanvas('canv_alpha_hist','canv_alpha_hist',1100,900)
#csvb_qq_global_hist.Draw('COLZ')
#canv_alpha_func = TCanvas('canv_alpha_func','canv_alpha_func',1100,900)
#f_alpha.Draw('COLZ')
##plot epsilon 2D plot and total function
#canv_epsilon_hist = TCanvas('canv_epsilon_hist','canv_epsilon_hist',1100,900)
#csvb_gg_global_hist.Draw('COLZ')
#canv_epsilon_func = TCanvas('canv_epsilon_func','canv_epsilon_func',1100,900)
#f_epsilon.Draw('COLZ')
##plot alpha projected plots
#canv_alpha_proj = TCanvas('canv_alpha_proj','canv_alpha_proj',1100,900)
#alpha_proj_histo = symmetrize(csvb_qq_global_hist.ProjectionY())
#alpha_proj_histo.Scale(1./alpha_proj_histo.Integral()) 
#alpha_proj_histo.GetYaxis().SetRangeUser(0.,0.07) 
#npoints_alpha = 10*csvb_qq_global_hist.GetYaxis().GetNbins()
#xpoints_alpha = array('d',npoints_alpha*[0.])
#ypoints_alpha = array('d',npoints_alpha*[0.])
#for i in range(npoints_alpha) :
#	pointspacing = 2./(npoints_alpha)
#	cstar_low = -1.0+i*pointspacing
#	cstar = -1.0+(i+0.5)*pointspacing
#	cstar_hi  = -1.0+(i+1)*pointspacing
#	xpoints_alpha[i]=cstar
#	ypoints_alpha[i]=f_alpha.Integral(0.,1.,cstar_low,cstar_hi)
#sumypoints_alpha = 0. 
#for i in range(npoints_alpha) : 
#	sumypoints_alpha+=ypoints_alpha[i] 
#for i in range(npoints_alpha) : 
#	fac = npoints_alpha/csvb_qq_global_hist.GetYaxis().GetNbins() 
#	ypoints_alpha[i]  = fac*ypoints_alpha[i]/sumypoints_alpha 
#alpha_graph = TGraph(npoints_alpha,xpoints_alpha,ypoints_alpha)
#alpha_graph.SetMarkerStyle(20)
#alpha_graph.SetLineWidth(3)
#alpha_graph.SetLineColor(kRed)
#alpha_proj_histo.SetMarkerStyle(20)
#alpha_proj_histo.Draw('E0')
#alpha_graph.Draw('L SAME')
#leg_alpha = TLegend(0.1,0.7,0.48,0.9)
#leg_alpha.AddEntry(alpha_proj_histo,"MC Events","PE0")
#leg_alpha.AddEntry(alpha_graph,"Fit (#alpha="+str(fit_alpha)+")","L")
#leg_alpha.Draw('SAME')
##plot epsilon projected plots
#canv_epsilon_proj = TCanvas('canv_epsilon_proj','canv_epsilon_proj',1100,900)
#epsilon_proj_histo = symmetrize(csvb_gg_global_hist.ProjectionY())
#epsilon_proj_histo.Scale(1./epsilon_proj_histo.Integral())
#epsilon_proj_histo.GetYaxis().SetRangeUser(0.,0.15)
#npoints_epsilon = 10*csvb_gg_global_hist.GetYaxis().GetNbins()
#xpoints_epsilon = array('d',npoints_epsilon*[0.])
#ypoints_epsilon = array('d',npoints_epsilon*[0.])
#for i in range(npoints_epsilon) :
#	pointspacing = 2./(npoints_epsilon)
#	cstar_low = -1.0+i*pointspacing
#	cstar = -1.0+(i+0.5)*pointspacing
#	cstar_hi  = -1.0+(i+1)*pointspacing
#	xpoints_epsilon[i]=cstar
#	ypoints_epsilon[i]=f_epsilon.Integral(0.,1.,cstar_low,cstar_hi)
#sumypoints_epsilon = 0.
#for i in range(npoints_epsilon) :
#	sumypoints_epsilon+=ypoints_epsilon[i]
#for i in range(npoints_epsilon) :
#	fac = npoints_epsilon/csvb_gg_global_hist.GetYaxis().GetNbins()
#	ypoints_epsilon[i]  = fac*ypoints_epsilon[i]/sumypoints_epsilon
#epsilon_graph = TGraph(npoints_epsilon,xpoints_epsilon,ypoints_epsilon)
#epsilon_graph.SetMarkerStyle(20)
#epsilon_graph.SetLineWidth(3)
#epsilon_graph.SetLineColor(kRed)
#epsilon_proj_histo.SetMarkerStyle(20)
#epsilon_proj_histo.Draw('E0')
#epsilon_graph.Draw('L SAME')
#leg_epsilon = TLegend(0.1,0.7,0.48,0.9)
#leg_epsilon.AddEntry(epsilon_proj_histo,"MC Events","PE0")
#leg_epsilon.AddEntry(epsilon_graph,"Fit (#epsilon="+str(fit_epsilon)+")","L")
#leg_epsilon.Draw('SAME')

#write plots to file
outFile.cd()
outFile.Write()
outFile.Close()


