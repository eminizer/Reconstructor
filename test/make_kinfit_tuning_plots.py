#imports
from ROOT import *

#name of file to run on
inputfilename = '../total_ttree_files/mcatnlo_semilep_TT_all.root'

#open the input file
infile = TFile(inputfilename)

#star up the output file
outfile = TFile('kinfit_tuning_plots.root','recreate')

#get the tree from the file
tree = infile.Get('tree')

#declare the plots
tmass_low=100.; tmass_high=350.; wmass_low=50.; wmass_high=110.
allhistos = []
t1leptM = TH1F('t1leptM','Leptonic top mass, correct hypothesis, fully merged; M_{top}^{lep} (GeV); fraction',25,tmass_low,tmass_high); allhistos.append(t1leptM)
t1hadtM = TH1F('t1hadtM','Hadronic top mass, correct hypothesis, fully merged; M_{top}^{had} (GeV); fraction',25,tmass_low,tmass_high); allhistos.append(t1hadtM)
t2leptM = TH1F('t2leptM','Leptonic top mass, correct hypothesis, partially merged; M_{top}^{lep} (GeV); fraction',25,tmass_low,tmass_high); allhistos.append(t2leptM)
t2hadtM = TH1F('t2hadtM','Hadronic top mass, correct hypothesis, partially merged; M_{top}^{had} (GeV); fraction',25,tmass_low,tmass_high); allhistos.append(t2hadtM)
t2hadWM = TH1F('t2hadWM','Hadronic W mass, correct hypothesis, partially merged; M_{W}^{had} (GeV); fraction',30,wmass_low,wmass_high); allhistos.append(t2hadWM)
for histo in allhistos :
	histo.SetDirectory(0)

#plot plots
weights = '12900.*weight'
com_cuts = 'fullselection==1 && ismatchable==1'
tree.Draw('leptcorprefitM>>t1leptM(25,'+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*(ttype==1 && ('+com_cuts+'))')
tree.Draw('hadtcorprefitM>>t1hadtM(25,'+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*(ttype==1 && ('+com_cuts+'))')
tree.Draw('leptcorprefitM>>t2leptM(25,'+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*(ttype==2 && ('+com_cuts+'))')
tree.Draw('hadtcorprefitM>>t2hadtM(25,'+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*(ttype==2 && ('+com_cuts+'))')
tree.Draw('hadWcorprefitM>>t2hadWM(30,'+str(wmass_low)+','+str(wmass_high)+')','('+weights+')*(ttype==2 && ('+com_cuts+'))')
for histo in allhistos :
	histo.Add(gROOT.FindObject(histo.GetName()))

#declare canvases
allcanvs = []
t1leptM_canv = TCanvas('t1leptM_canv','t1leptM_canv',1100,900); allcanvs.append(t1leptM_canv)
t1hadtM_canv = TCanvas('t1hadtM_canv','t1hadtM_canv',1100,900); allcanvs.append(t1hadtM_canv)
t2leptM_canv = TCanvas('t2leptM_canv','t2leptM_canv',1100,900); allcanvs.append(t2leptM_canv)
t2hadtM_canv = TCanvas('t2hadtM_canv','t2hadtM_canv',1100,900); allcanvs.append(t2hadtM_canv)
t2hadWM_canv = TCanvas('t2hadWM_canv','t2hadWM_canv',1100,900); allcanvs.append(t2hadWM_canv)

#draw plots on canvases and fit with gaussians
for i in range(len(allhistos)) :
	print 'Doing '+str(allhistos[i].GetName())
	allcanvs[i].cd()
	allhistos[i].Draw()
	allhistos[i].Fit('gaus')

#save plots and canvases
outfile.cd()
for histo in allhistos :
	histo.Write()
for canv in allcanvs :
	canv.Write()