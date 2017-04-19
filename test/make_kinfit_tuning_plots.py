#imports
from ROOT import *

#name of file to run on
inputfilename = '../total_ttree_files/powheg_TT_skim_all_4_18_2017.root'

#open the input file
infile = TFile(inputfilename)

#get the tree from the file
fullTree = infile.Get('tree')
#fullTree.SetDirectory(0)

#declare the plots
#tmass_low=100.; tmass_high=350.; ntbins=75
tmass_low=100.; tmass_high=300.; ntbins=100
allhistos = []
t1leptM = TH1F('t1leptM','type-1 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t1leptM)
t1hadtM = TH1F('t1hadtM','type-1 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t1hadtM)
t2leptM = TH1F('t2leptM','type-2 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t2leptM)
t2hadtM = TH1F('t2hadtM','type-2 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t2hadtM)
t3leptM = TH1F('t3leptM','type-3 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t3leptM)
t3hadtM = TH1F('t3hadtM','type-3 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t3hadtM)
for histo in allhistos :
	histo.SetDirectory(0)

#start up the output file
outfile = TFile('kinfit_tuning_plots_4_18_2017.root','recreate')

#set up tree
com_cuts = 'fullselection==1 && ismatchable==1'
tree=fullTree.CopyTree(com_cuts)
#tree.SetDirectory(0)

#plot plots
weights = '35867.*weight*sf_pileup*sf_lep_ID*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'
print 'Drawing type-1 leptonic top mass...'
tree.Draw('leptcorprefitM>>t1leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==1)')
print 'Drawing type-1 hadronic top mass...'
tree.Draw('hadtcorprefitM>>t1hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==1)')
print 'Drawing type-2 leptonic top mass...'
tree.Draw('leptcorprefitM>>t2leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==2)')
print 'Drawing type-2 hadronic top mass...'
tree.Draw('hadtcorprefitM>>t2hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==2)')
print 'Drawing type-3 leptonic top mass...'
tree.Draw('leptcorprefitM>>t3leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3)')
print 'Drawing type-3 hadronic top mass...'
tree.Draw('hadtcorprefitM>>t3hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3)')
print 'Done.'
for histo in allhistos :
	thisName = histo.GetName()
	#print 'histo name = %s'%(thisName) #DEBUG
	toAdd = gROOT.FindObject(thisName)
	#print '	object = %s'%(toAdd) #DEBUG
	histo.Add(toAdd)
	if histo.Integral()==0 :
		print thisName+' has an integral of zero'
	histo.Scale(1./histo.Integral())
	histo.SetMarkerStyle(20)
	histo.GetYaxis().SetTitleOffset(1.2)

#declare canvases
allcanvs = []
t1leptM_canv = TCanvas('t1leptM_canv','t1leptM_canv',1100,900); allcanvs.append(t1leptM_canv)
t1hadtM_canv = TCanvas('t1hadtM_canv','t1hadtM_canv',1100,900); allcanvs.append(t1hadtM_canv)
t2leptM_canv = TCanvas('t2leptM_canv','t2leptM_canv',1100,900); allcanvs.append(t2leptM_canv)
t2hadtM_canv = TCanvas('t2hadtM_canv','t2hadtM_canv',1100,900); allcanvs.append(t2hadtM_canv)
t3leptM_canv = TCanvas('t3leptM_canv','t3leptM_canv',1100,900); allcanvs.append(t3leptM_canv)
t3hadtM_canv = TCanvas('t3hadtM_canv','t3hadtM_canv',1100,900); allcanvs.append(t3hadtM_canv)
total_canv = TCanvas('total_canv','total_canv',1100,900)

colors = [kRed,kBlue,kGreen,kOrange,kMagenta,kCyan+2]
legEntries = ['type-1, leptonic','type-1, hadronic','type-2, leptonic','type-2, hadronic','type-3, leptonic','type-3, hadronic']

#draw plots on canvases and fit with gaussians, print results
lines = []
total_leg = TLegend(0.1,0.7,0.48,0.9)
for i in range(len(allhistos)) :
	histo = allhistos[i]
	print 'Doing '+str(allhistos[i].GetName())
	allcanvs[i].cd()
	histo.Draw()
	gausfunc1 = TF1('gausfunc1','gaus',150.,195.)
	histo.Fit(gausfunc1,"","",150.,195.)
	gmean = gausfunc1.GetParameter(1)
	gausfunc = TF1('gausfunc','gaus',gmean-20.,gmean+20.)
	#gausfunc = TF1('gausfunc','[0]*ROOT::Math::lognormal_pdf(x,[1],[2])',tmass_low,tmass_high)
	#gausfunc.SetParameters(10000000000.,172.5,0.5)
	gausfunc.SetLineColor(kRed)
	gausfunc.SetLineWidth(4)
	histo.Fit(gausfunc,"","",gmean-20.,gmean+20.)
	gausfunc.Draw('SAME')
	gmean = gausfunc.GetParameter(1); gwidth = gausfunc.GetParameter(2)
	leg = TLegend(0.1,0.7,0.48,0.9)
	leg.AddEntry(gausfunc,"Gaussian Fit","L")
	leg.Draw('SAME')
	lines.append('%s, Gaussian: %.1f +/- %.1f'%(allhistos[i].GetName(),gmean,gwidth))
	total_canv.cd()
	histo.SetMarkerColor(colors[i])
	histo.SetFillColor(colors[i])
	gausfunc.SetLineColor(colors[i])
	histo.Draw() if i==0 else histo.Draw('SAME')
	gausfunc.Draw('SAME')
	total_leg.AddEntry(histo,'#splitline{%s}{M_{t}=%.1f #pm %.1f}'%(legEntries[i],gmean,gwidth),'F')
total_canv.cd()
total_leg.Draw('SAME')

print '\n\n\nFit Results:'	
for line in lines :
	print line

#save plots and canvases
outfile.cd()
for histo in allhistos :
	histo.Write()
for canv in allcanvs :
	canv.Write()
total_canv.Write()