#This script plots the differential top and antitop pT cross sections in powheg MC truth and in NNLO QCD + EW theory (from arXiv:1705.04105)
#It takes their ratio, and fits the result with an exponential to give top-by-top SFs that can be combined into event-by-event top pT reweights
#It calculates the average resulting top pT reweight over the whole generated sample to find the correction factor for the inidividual reweights

#imports
from ROOT import gROOT, TH2D, TFile, TH1D, kRed, kBlue, kBlack, TCanvas, TLegend, TLatex, TF2
import CMS_lumi, tdrstyle
from array import array
from math import exp,sqrt

#batch mode/TDR style
gROOT.SetBatch()
tdrstyle.setTDRStyle()
iPeriod = 4 #13TeV iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"

#pT histogram bins
nbins_MC = 300
pT_low_MC = 0.; pT_high_MC = 1500.
pT_bins_theory = array('d',[0., 65., 125., 200., 290., 400., 550.])
nbins_theory = len(pT_bins_theory)-1

#start with the normalized cross section histograms from powheg MC
print 'getting top/antitop pT spectra from MC'
ntuple_filenames = ['root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_0.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_1.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_10.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_11.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_12.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_13.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_14.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_15.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_16.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_17.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_18.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_19.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_2.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_20.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_21.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_22.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_3.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_4.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_5.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_6.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_7.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_8.root',
					'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_9.root',
					]
total_MC_atop_vs_top_pT_histo = TH2D('MC_atop_vs_top_pT','antitop vs. top pT spectrum from powheg+pythia MC; top p_{T} [GeV]; antitop p_{T} [GeV]',nbins_MC,pT_low_MC,pT_high_MC,nbins_MC,pT_low_MC,pT_high_MC)
#sum up the histograms from each file
for ntfn in ntuple_filenames :
	print '  adding histogram from file '+ntfn
	ntf = TFile.Open(ntfn)
	total_MC_atop_vs_top_pT_histo.Add(ntf.Get('EventCounter/gen_atop_vs_top_pt'))
	ntf.Close()
#project onto each axis
top_pT_from_MC_finer_bins  = total_MC_atop_vs_top_pT_histo.ProjectionX()	
atop_pT_from_MC_finer_bins = total_MC_atop_vs_top_pT_histo.ProjectionY()
#open the output file to save the new histograms
outfile = TFile('theory_driven_top_pT_reweighting.root','recreate')
top_pT_from_MC = TH1D('top_pT_from_MC','top p_{T} spectrum from powheg+pythia MC; top p_{T} [GeV]; event fraction',nbins_theory,pT_bins_theory)
atop_pT_from_MC = TH1D('atop_pT_from_MC','antitop p_{T} spectrum from powheg+pythia MC; antitop p_{T} [GeV]; event fraction',nbins_theory,pT_bins_theory)
#populate histograms with theory binning
for i in range(1,nbins_MC+1) :
	thispt = top_pT_from_MC_finer_bins.GetBinCenter(i)
	top_pT_from_MC.Fill(thispt,top_pT_from_MC_finer_bins.GetBinContent(i))
	atop_pT_from_MC.Fill(thispt,atop_pT_from_MC_finer_bins.GetBinContent(i))
#reset errors to sqrt(N)
for i in range(1,nbins_theory) :
	top_pT_from_MC.SetBinError(i,sqrt(top_pT_from_MC.GetBinContent(i)))
	atop_pT_from_MC.SetBinError(i,sqrt(atop_pT_from_MC.GetBinContent(i)))
#renormalize
top_pT_from_MC.Scale(1./top_pT_from_MC.Integral())
atop_pT_from_MC.Scale(1./atop_pT_from_MC.Integral())
print 'Done.'

#normalized cross section histograms from theory
print 'getting top/antitop pT spectra from theory'
top_cross_sections_from_theory  = [2.1120110700000001e+02,3.0176205599999997e+02,2.0840831900000001e+02,7.8436658899999998e+01,2.2455439800000001e+01,5.6133840499999996e+00]
top_cross_section_errors_from_theory  = [(2.1592937692000137e+02-2.0122322881064878e+02)/2.,(3.0677389105744965e+02-2.8934788160854748e+02)/2.,(2.1203843533810590e+02-1.9940882396706320e+02)/2.,
										 (7.9868717936622332e+01-7.4972824354074049e+01)/2.,(2.2842218642561143e+01-2.1387898129119471e+01)/2.,(5.7111894562548029e+00-5.3592205745254828e+00)/2.]
atop_cross_sections_from_theory = [2.0986704700000001e+02,3.0143660799999998e+02,2.0896114900000001e+02,7.9168844500000006e+01,2.2596528200000002e+01,5.7615639600000002e+00]
atop_cross_section_errors_from_theory = [(2.1462265652322915e+02-2.0005969468086477e+02)/2.,(3.0654099929526978e+02-2.8887777690470381e+02)/2.,(2.1260090839067880e+02-1.9998117899166979e+02)/2.,
										 (8.0636868200388975e+01-7.5572772635720739e+01)/2.,(2.3004256478968600e+01-2.1547841550995489e+01)/2.,(5.8554252898981005e+00-5.4984672066943716e+00)/2.]
top_pT_from_theory  = TH1D('top_pT_from_theory','top p_{T} spectrum from NNLO QCD + EW theory; top p_{T} [GeV]; event fraction',nbins_theory,pT_bins_theory)
atop_pT_from_theory = TH1D('atop_pT_from_theory','antitop p_{T} spectrum from NNLO QCD + EW theory; antitop p_{T} [GeV]; event fraction',nbins_theory,pT_bins_theory)
for i in range(len(top_cross_sections_from_theory)) :
	top_pT_from_theory.SetBinContent(i+1,top_cross_sections_from_theory[i])
	top_pT_from_theory.SetBinError(i+1,top_cross_section_errors_from_theory[i])
	atop_pT_from_theory.SetBinContent(i+1,atop_cross_sections_from_theory[i])
	atop_pT_from_theory.SetBinError(i+1,atop_cross_section_errors_from_theory[i])
top_pT_from_theory.Scale(1./top_pT_from_theory.Integral())
atop_pT_from_theory.Scale(1./atop_pT_from_theory.Integral())
print 'Done.'

#ratio plots
print 'Plotting theory/MC ratios'
top_pT_ratios  = TH1D('top_pT_ratios','; top p_{T} [GeV]; NNLO QCD + EW theory/MC',nbins_theory,pT_bins_theory)
atop_pT_ratios = TH1D('atop_pT_ratios','; antitop p_{T} [GeV]; NNLO QCD + EW theory/MC',nbins_theory,pT_bins_theory)
for i in range(1,len(top_cross_sections_from_theory)+1) :
	top_cont_theory = top_pT_from_theory.GetBinContent(i)
	top_cont_MC = top_pT_from_MC.GetBinContent(i)
	top_err_theory = top_pT_from_theory.GetBinError(i)
	top_err_MC = top_pT_from_MC.GetBinError(i)
	top_ratio = top_cont_theory/top_cont_MC
	top_ratio_err = top_ratio*sqrt((top_err_theory/top_cont_theory)**2+(top_err_MC/top_cont_MC)**2)
	top_pT_ratios.SetBinContent(i,top_ratio); top_pT_ratios.SetBinError(i,top_ratio_err)
	atop_cont_theory = atop_pT_from_theory.GetBinContent(i)
	atop_cont_MC = atop_pT_from_MC.GetBinContent(i)
	atop_err_theory = atop_pT_from_theory.GetBinError(i)
	atop_err_MC = atop_pT_from_MC.GetBinError(i)
	atop_ratio = atop_cont_theory/atop_cont_MC
	atop_ratio_err = atop_ratio*sqrt((atop_err_theory/atop_cont_theory)**2+(atop_err_MC/atop_cont_MC)**2)
	atop_pT_ratios.SetBinContent(i,atop_ratio); atop_pT_ratios.SetBinError(i,atop_ratio_err)
print 'Done.'

#fit with exponentials
print 'Fitting ratios with exponential functions'
top_pT_ratios.Fit('expo')
top_fit = top_pT_ratios.GetFunction('expo')
top_p0 = top_fit.GetParameter(0)
top_p0_err = top_fit.GetParError(0)
top_p1 = top_fit.GetParameter(1)
top_p1_err = top_fit.GetParError(1)
print '  top spectrum fit: ratio = exp(p0+p1*pT) where p0=%.6f+/-%.6f and p1=%.6f+/-%.6f'%(top_p0,top_p0_err,top_p1,top_p1_err)
atop_pT_ratios.Fit('expo')
atop_fit = atop_pT_ratios.GetFunction('expo')
atop_p0 = atop_fit.GetParameter(0)
atop_p0_err = atop_fit.GetParError(0)
atop_p1 = atop_fit.GetParameter(1)
atop_p1_err = atop_fit.GetParError(1)
print '  antitop spectrum fit: ratio = exp(p0+p1*pT) where p0=%.6f+/-%.6f and p1=%.6f+/-%.6f'%(atop_p0,atop_p0_err,atop_p1,atop_p1_err)
print 'Done.'

#calculate weighted average of factors over the whole sample and print it out
avgeventweight = 0.; weightsum = 0.
#loop over all the bins in the finely-binned 2D histogram
nglobalbins = total_MC_atop_vs_top_pT_histo.GetSize()
for i in range(nglobalbins) :
	if not total_MC_atop_vs_top_pT_histo.IsBinOverflow(i) and not total_MC_atop_vs_top_pT_histo.IsBinUnderflow(i) :
		#calculate the "event weight" based on the central top and antitop pTs
		binx = array('i',[0]); biny = array('i',[0]); binz = array('i',[0])
		total_MC_atop_vs_top_pT_histo.GetBinXYZ(i,binx,biny,binz)
		thisbincontent = total_MC_atop_vs_top_pT_histo.GetBinContent(i)
		thisbintoppT = total_MC_atop_vs_top_pT_histo.GetXaxis().GetBinCenter(binx[0])
		thisbinatoppT = total_MC_atop_vs_top_pT_histo.GetYaxis().GetBinCenter(biny[0])
		#bring the pTs back in range of the calculations
		if thisbintoppT>pT_bins_theory[-1] : thisbintoppT=pT_bins_theory[-1]
		if thisbinatoppT>pT_bins_theory[-1] : thisbinatoppT=pT_bins_theory[-1]
		topfac = exp(top_p0+top_p1*thisbintoppT)
		atopfac = exp(atop_p0+atop_p1*thisbinatoppT)
		avgeventweight+=(thisbincontent*sqrt(topfac*atopfac)); weightsum+=thisbincontent
avgeventweight/=weightsum
print 'average event weight over entire MC sample = %.8f'%(avgeventweight)

#make some nice plots of the top/antitop cross section ratios with fits
print 'Plotting top and antitop cross section ratios'
nfitbins=200
#for top fit
top_fit_hist = TH1D('top_fit_hist',';top p_{T} [GeV];theory/MC ratio',nfitbins,pT_bins_theory[0],pT_bins_theory[-1])
top_fit_line = TH1D('top_fit_line','',nfitbins,pT_bins_theory[0],pT_bins_theory[-1])
for i in range(1,nfitbins+1) :
	thispt = top_fit_hist.GetBinCenter(i)
	thiscont = exp(top_p0+top_p1*thispt)
	thiserr = thiscont*sqrt(top_p0_err**2+(thispt*top_p1_err)**2)
	top_fit_hist.SetBinContent(i,thiscont)
	top_fit_hist.SetBinError(i,thiserr)
	top_fit_line.SetBinContent(i,thiscont)
top_ratio_canv = TCanvas('top_ratio_canv','top_ratio_canv',1100,900)
top_ratio_canv.SetTopMargin(0.07)
top_fit_hist.Draw('E3')
top_fit_hist.GetYaxis().SetRangeUser(0.55,1.22)
top_fit_line.Draw('SAME L HIST')
top_pT_ratios.Draw('HIST SAME PE1')
top_pT_ratios.SetLineColor(kRed+2)
top_pT_ratios.SetLineWidth(3)
top_fit_hist.SetFillColor(kBlack)
top_fit_hist.SetFillStyle(3013)
top_fit_line.SetLineWidth(2)
top_fit_hist.SetMarkerStyle(1)
top_fit_hist.SetStats(0)
top_leg = TLegend(0.5,0.7,0.78,0.88)
top_leg.SetBorderSize(0)
top_leg.AddEntry(top_pT_ratios,'measured ratios','PE')
top_leg.AddEntry(top_fit_hist,'fit','LF')
top_leg.Draw()
CMS_lumi.CMS_lumi(top_ratio_canv, iPeriod, 11)
text_theory = 'theory: NNLO QCD + EW (arXiv:1705.04105)'
textheory = TLatex()
textheory.SetNDC()
textheory.SetTextAngle(0)
textheory.SetTextColor(kBlack)
textheory.SetTextFont(42)
textheory.SetTextAlign(11) 
textheory.SetTextSize(0.04)
textheory.DrawLatex(0.20,0.24,text_theory)
text_MC = 'MC: powheg+pythia8'
texMC = TLatex()
texMC.SetNDC()
texMC.SetTextAngle(0)
texMC.SetTextColor(kBlack)
texMC.SetTextFont(42)
texMC.SetTextAlign(11) 
texMC.SetTextSize(0.04)
texMC.DrawLatex(0.2,0.19,text_MC)
text_func = 'ratio = e^{%.4f - %.4f x p_{T}}'%(top_p0,abs(top_p1))
texfunc = TLatex()
texfunc.SetNDC()
texfunc.SetTextAngle(0)
texfunc.SetTextColor(kBlack)
texfunc.SetTextFont(42)
texfunc.SetTextAlign(11) 
texfunc.SetTextSize(0.04)
texfunc.DrawLatex(0.2,0.29,text_func)
top_ratio_canv.Write()
#for antitop fit
atop_fit_hist = TH1D('atop_fit_hist',';antitop p_{T} [GeV];theory/MC ratio',nfitbins,pT_bins_theory[0],pT_bins_theory[-1])
atop_fit_line = TH1D('atop_fit_line','',nfitbins,pT_bins_theory[0],pT_bins_theory[-1])
for i in range(1,nfitbins+1) :
	thispt = atop_fit_hist.GetBinCenter(i)
	thiscont = exp(atop_p0+atop_p1*thispt)
	thiserr = thiscont*sqrt(atop_p0_err**2+(thispt*atop_p1_err)**2)
	atop_fit_hist.SetBinContent(i,thiscont)
	atop_fit_hist.SetBinError(i,thiserr)
	atop_fit_line.SetBinContent(i,thiscont)
atop_ratio_canv = TCanvas('atop_ratio_canv','atop_ratio_canv',1100,900)
atop_ratio_canv.SetTopMargin(0.07)
atop_fit_hist.Draw('E3')
atop_fit_hist.GetYaxis().SetRangeUser(0.55,1.22)
atop_fit_line.Draw('SAME L HIST')
atop_pT_ratios.Draw('HIST SAME PE1')
atop_pT_ratios.SetLineColor(kBlue+2)
atop_pT_ratios.SetLineWidth(3)
atop_fit_hist.SetFillColor(kBlack)
atop_fit_hist.SetFillStyle(3013)
atop_fit_line.SetLineWidth(2)
atop_fit_hist.SetMarkerStyle(1)
atop_fit_hist.SetStats(0)
atop_leg = TLegend(0.5,0.7,0.78,0.88)
atop_leg.SetBorderSize(0)
atop_leg.AddEntry(atop_pT_ratios,'measured ratios','PE')
atop_leg.AddEntry(atop_fit_hist,'fit','LF')
atop_leg.Draw()
CMS_lumi.CMS_lumi(atop_ratio_canv, iPeriod, 11)
text_theory = 'theory: NNLO QCD + EW (arXiv:1705.04105)'
textheory = TLatex()
textheory.SetNDC()
textheory.SetTextAngle(0)
textheory.SetTextColor(kBlack)
textheory.SetTextFont(42)
textheory.SetTextAlign(11) 
textheory.SetTextSize(0.04)
textheory.DrawLatex(0.20,0.24,text_theory)
text_MC = 'MC: powheg+pythia8'
texMC = TLatex()
texMC.SetNDC()
texMC.SetTextAngle(0)
texMC.SetTextColor(kBlack)
texMC.SetTextFont(42)
texMC.SetTextAlign(11) 
texMC.SetTextSize(0.04)
texMC.DrawLatex(0.2,0.19,text_MC)
text_func = 'ratio = e^{%.4f - %.4f x p_{T}}'%(atop_p0,abs(atop_p1))
texfunc = TLatex()
texfunc.SetNDC()
texfunc.SetTextAngle(0)
texfunc.SetTextColor(kBlack)
texfunc.SetTextFont(42)
texfunc.SetTextAlign(11) 
texfunc.SetTextSize(0.04)
texfunc.DrawLatex(0.2,0.29,text_func)
atop_ratio_canv.Write()
print 'Done.'

#make a nice plot of the reweighting function
print 'Plotting event weights'
event_weight_canv = TCanvas('event_weight_canv','event_weight_canv',1300,1000)
event_weight_canv.SetRightMargin(0.10)
weightfunc = TF2('weightfunc','sqrt(exp(%.8f+%.8f*x)*exp(%.8f+%.8f*y))'%(top_p0,top_p1,atop_p0,atop_p1),0.,550.,0.,550.)
weightfunc.Draw('COLZ')
weightfunc.SetTitle('top pT event weights')
weightfunc.GetYaxis().SetTitle('antitop p_{T} [GeV]')
weightfunc.GetXaxis().SetTitle('top p_{T} [GeV]')
event_weight_canv.Update()
event_weight_canv.Write()
print 'Done.'

#write and close the output file
outfile.Write()
outfile.Close()
