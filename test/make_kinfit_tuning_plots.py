#imports
from ROOT import *

#name of file to run on
inputfilename = '../total_ttree_files/powheg_TT_skim_all_4_13_2017.root'

#open the input file
infile = TFile(inputfilename)

#star up the output file
outfile = TFile('kinfit_tuning_plots_4_13_2017.root','recreate')

#get the tree from the file
fullTree = infile.Get('tree')

#declare the plots
tmass_low=100.; tmass_high=350.; ntbins=25
allhistos = []
t1leptM = TH1F('t1leptM','type-1 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t1leptM)
t1hadtM = TH1F('t1hadtM','type-1 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t1hadtM)
t2leptM = TH1F('t2leptM','type-2 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t2leptM)
t2hadtM = TH1F('t2hadtM','type-2 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t2hadtM)
t31leptM = TH1F('t31leptM','type-3 (1 hadronic-side jet) leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t31leptM)
t31hadtM = TH1F('t31hadtM','type-3 (1 hadronic-side jet) hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t31hadtM)
t32leptM = TH1F('t32leptM','type-3 (2 hadronic-side jets) leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t32leptM)
t32hadtM = TH1F('t32hadtM','type-3 (2 hadronic-side jets) hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t32hadtM)
t33leptM = TH1F('t33leptM','type-3 (3 hadronic-side jets) leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t33leptM)
t33hadtM = TH1F('t33hadtM','type-3 (3 hadronic-side jets) hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t33hadtM)
t4leptM = TH1F('t4leptM','type-4 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t4leptM)
t4hadtM = TH1F('t4hadtM','type-4 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t4hadtM)
for histo in allhistos :
	histo.SetDirectory(0)

#declare the Lorentzian functions
allfuncs = []
t1leptfunc = TF1('t1leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[3]*[2])',tmass_low,tmass_high); allfuncs.append(t1leptfunc)
t1hadtfunc = TF1('t1hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t1hadtfunc)
t2leptfunc = TF1('t2leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t2leptfunc)
t2hadtfunc = TF1('t2hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t2hadtfunc)
t31leptfunc = TF1('t31leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t31leptfunc)
t31hadtfunc = TF1('t31hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t31hadtfunc)
t32leptfunc = TF1('t32leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t32leptfunc)
t32hadtfunc = TF1('t32hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t32hadtfunc)
t33leptfunc = TF1('t33leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t33leptfunc)
t33hadtfunc = TF1('t33hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t33hadtfunc)
t4leptfunc = TF1('t4leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t4leptfunc)
t4hadtfunc = TF1('t4hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t4hadtfunc)
for func in allfuncs :
	func.SetParName(0,'const'); func.SetParName(1,'mean'); func.SetParName(2,'width')
	func.SetParameter(0,3500000.); func.SetParameter(1,173.); func.SetParameter(2,20.)

#plot plots
weights = '35867.*weight*sf_pileup*sf_lep_ID*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'
com_cuts = 'eventType<2 && fullselection==1 && ismatchable==1'
tree=fullTree.CopyTree(com_cuts)
print 'Drawing type-1 leptonic top mass...'
tree.Draw('leptcorprefitM>>t1leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==1)')
print 'Drawing type-1 hadronic top mass...'
tree.Draw('hadtcorprefitM>>t1hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==1)')
print 'Drawing type-2 leptonic top mass...'
tree.Draw('leptcorprefitM>>t2leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==2)')
print 'Drawing type-2 hadronic top mass...'
tree.Draw('hadtcorprefitM>>t2hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==2)')
print 'Drawing type-3 (1 had-side jet) leptonic top mass...'
tree.Draw('leptcorprefitM>>t31leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3 && nHadSideJetsCorrect==1)')
print 'Drawing type-3 (1 had-side jet) hadronic top mass...'
tree.Draw('hadtcorprefitM>>t31hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3 && nHadSideJetsCorrect==1)')
print 'Drawing type-3 (2 had-side jets) leptonic top mass...'
tree.Draw('leptcorprefitM>>t32leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3 && nHadSideJetsCorrect==2)')
print 'Drawing type-3 (2 had-side jets) hadronic top mass...'
tree.Draw('hadtcorprefitM>>t32hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3 && nHadSideJetsCorrect==2)')
print 'Drawing type-3 (3 had-side jets) leptonic top mass...'
tree.Draw('leptcorprefitM>>t33leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3 && nHadSideJetsCorrect==3)')
print 'Drawing type-3 (3 had-side jets) hadronic top mass...'
tree.Draw('hadtcorprefitM>>t33hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3 && nHadSideJetsCorrect==3)')
print 'Drawing type-4 leptonic top mass...'
tree.Draw('leptcorprefitM>>t4leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==4)')
print 'Drawing type-4 hadronic top mass...'
tree.Draw('hadtcorprefitM>>t4hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==4)')
print 'Done.'
for histo in allhistos :
	thisName = histo.GetName()
	#print 'histo name = %s'%(thisName) #DEBUG
	toAdd = gROOT.FindObject(thisName)
	#print '	object = %s'%(toAdd) #DEBUG
	histo.Add(toAdd)
	histo.Scale(1./histo.Integral())
	histo.SetMarkerStyle(20)
	histo.GetYaxis().SetTitleOffset(1.2)

#declare canvases
allcanvs = []
t1leptM_canv = TCanvas('t1leptM_canv','t1leptM_canv',1100,900); allcanvs.append(t1leptM_canv)
t1hadtM_canv = TCanvas('t1hadtM_canv','t1hadtM_canv',1100,900); allcanvs.append(t1hadtM_canv)
t2leptM_canv = TCanvas('t2leptM_canv','t2leptM_canv',1100,900); allcanvs.append(t2leptM_canv)
t2hadtM_canv = TCanvas('t2hadtM_canv','t2hadtM_canv',1100,900); allcanvs.append(t2hadtM_canv)
t31leptM_canv = TCanvas('t31leptM_canv','t31leptM_canv',1100,900); allcanvs.append(t31leptM_canv)
t31hadtM_canv = TCanvas('t31hadtM_canv','t31hadtM_canv',1100,900); allcanvs.append(t31hadtM_canv)
t32leptM_canv = TCanvas('t32leptM_canv','t32leptM_canv',1100,900); allcanvs.append(t32leptM_canv)
t32hadtM_canv = TCanvas('t32hadtM_canv','t32hadtM_canv',1100,900); allcanvs.append(t32hadtM_canv)
t33leptM_canv = TCanvas('t33leptM_canv','t33leptM_canv',1100,900); allcanvs.append(t33leptM_canv)
t33hadtM_canv = TCanvas('t33hadtM_canv','t33hadtM_canv',1100,900); allcanvs.append(t33hadtM_canv)
t4leptM_canv = TCanvas('t4leptM_canv','t4leptM_canv',1100,900); allcanvs.append(t4leptM_canv)
t4hadtM_canv = TCanvas('t4hadtM_canv','t4hadtM_canv',1100,900); allcanvs.append(t4hadtM_canv)

#draw plots on canvases and fit with gaussians, print results
lines = []
for i in range(len(allhistos)) :
	histo = allhistos[i]
	func = allfuncs[i]
	print 'Doing '+str(allhistos[i].GetName())
	allcanvs[i].cd()
	histo.Draw()
	histo.Fit(func)
	func.SetLineColor(kRed)
	func.SetLineWidth(4)
	func.Draw('SAME')
	lmean = func.GetParameter(1); lwidth = func.GetParameter(2)
	gausfunc = TF1('gausfunc','gaus',tmass_low,tmass_high)
	histo.Fit(gausfunc)
	gausfunc.SetLineColor(kBlue)
	gausfunc.SetLineWidth(4)
	gausfunc.Draw('SAME')
	gmean = gausfunc.GetParameter(1); gwidth = gausfunc.GetParameter(2)
	leg = TLegend(0.1,0.7,0.48,0.9)
	leg.AddEntry(func,"Lorentzian","L")
	leg.AddEntry(gausfunc,"Gaussian","L")
	leg.Draw('SAME')
	lines.append('%s, Lorentzian: %.1f +/- %.1f, Gaussian: %.1f +/- %.1f'%(allhistos[i].GetName(),lmean,lwidth,gmean,gwidth))

print '\n\n\nFit Results:'	
for line in lines :
	print line

#save plots and canvases
outfile.cd()
for histo in allhistos :
	histo.Write()
for canv in allcanvs :
	canv.Write()