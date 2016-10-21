#imports
from ROOT import *

#output file 
outfile = TFile('chi2_comp_plot.root','recreate')

#prepend for all file names
pp = '../total_ttree_files/'

#list of semileptonic ttbar filenames
semilep_TT_filenames = []
semilep_TT_filenames.append(pp+'mcatnlo_semilep_TT_all.root')

#list of background filenames
background_filenames = []
background_filenames.append(pp+'mcatnlo_dilep_TT_all.root')
background_filenames.append(pp+'mcatnlo_had_TT_all.root')

#chain up the files
sig_chain = TChain('tree')
bg_chain  = TChain('tree')
for fn in semilep_TT_filenames :
	sig_chain.Add(fn)
for fn in background_filenames :
	bg_chain.Add(fn)

#declare the histograms
all_histos = []
chi2_low = -50.; chi2_hi = 75.; nbins = 25
sig_fs_histo = TH1F('sig_fs_histo','Kinematic Fit #chi^{2} Values; #chi^{2}; fraction',nbins,chi2_low,chi2_hi); sig_fs_histo.SetLineColor(kRed);  all_histos.append(sig_fs_histo)
sig_sr_histo = TH1F('sig_sr_histo','',nbins,chi2_low,chi2_hi); 													sig_sr_histo.SetLineColor(kBlue); all_histos.append(sig_sr_histo)
bg_fs_histo = TH1F('bg_fs_histo','',nbins,chi2_low,chi2_hi); 													bg_fs_histo.SetLineColor(kBlack); all_histos.append(bg_fs_histo)
bg_sr_histo = TH1F('bg_sr_histo','',nbins,chi2_low,chi2_hi); 			bg_sr_histo.SetLineStyle(2); 			bg_sr_histo.SetLineColor(kBlack); all_histos.append(bg_sr_histo)
for histo in all_histos :
	histo.SetDirectory(0)

#declare the canvas
canv = TCanvas('canv','canv',1100,900)

#draw chains into histos
weights = '12917.*weight'
com_cuts = 'fullselection==1 && validminimization==1'
sig_chain.Draw('chi2>>'+sig_fs_histo.GetName()+'('+str(nbins)+','+str(chi2_low)+','+str(chi2_hi)+')','('+weights+')*('+com_cuts+')')
sig_chain.Draw('chi2>>'+sig_sr_histo.GetName()+'('+str(nbins)+','+str(chi2_low)+','+str(chi2_hi)+')','('+weights+')*('+com_cuts+' && hadt_isttagged)')
bg_chain.Draw('chi2>>'+bg_fs_histo.GetName()+'('+str(nbins)+','+str(chi2_low)+','+str(chi2_hi)+')','('+weights+')*('+com_cuts+')')
bg_chain.Draw('chi2>>'+bg_sr_histo.GetName()+'('+str(nbins)+','+str(chi2_low)+','+str(chi2_hi)+')','('+weights+')*('+com_cuts+' && hadt_isttagged)')
for histo in all_histos :
	histo.Add(gROOT.FindObject(histo.GetName()))
	histo.Scale(1./histo.Integral())
	histo.GetYaxis().SetTitleOffset(1.2)
	histo.SetLineWidth(4)

#set up the legend
leg = TLegend(0.1,0.7,0.48,0.9)
leg.AddEntry(sig_fs_histo,'fully selected semileptonic ttbar','L')
leg.AddEntry(sig_sr_histo,'signal region semileptonic ttbar','L')
leg.AddEntry(bg_fs_histo,'fully selected dileptonic/hadronic ttbar','L')
leg.AddEntry(bg_sr_histo,'signal region dileptonic/hadronic ttbar','L')

#draw histos onto canvas
canv.cd()
all_histos[0].Draw('HIST')
for i in range(1,len(all_histos)) :
	all_histos[i].Draw('HIST SAME')
#draw the legend
leg.Draw()

#save histo and Canvas to file
outfile.cd()
for histo in all_histos :
	histo.Write()
canv.Write()

#close outputfile
outfile.Close()