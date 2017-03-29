#imports
from ROOT import *

#name of file to run on
inputfilename = '../total_ttree_files/powheg_TT_skim_all.root'

#open the input file
infile = TFile(inputfilename)

#star up the output file
outfile = TFile('kinfit_tuning_plots.root','recreate')

#get the tree from the file
fullTree = infile.Get('tree')

#declare the plots
tmass_low=0.; tmass_high=350.; ntbins=35
Wmass_low=0.; Wmass_high=350.; nWbins=35
allhistos = []
t1leptM = TH1F('t1leptM','type-1 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t1leptM)
t1hadtM = TH1F('t1hadtM','type-1 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t1hadtM)
t2leptM = TH1F('t2leptM','type-2 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t2leptM)
t2hadtM = TH1F('t2hadtM','type-2 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t2hadtM)
t3leptM = TH1F('t3leptM','type-3 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t3leptM)
t3hadtM = TH1F('t3hadtM','type-3 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t3hadtM)
t3hadWM = TH1F('t3hadWM','type-3 hadronic W mass, correct hypothesis; M_{W}^{had} (GeV); fraction',nWbins,Wmass_low,Wmass_high); allhistos.append(t3hadWM)
for histo in allhistos :
	histo.SetDirectory(0)

#declare the Lorentzian functions
allfuncs = []
t1leptfunc = TF1('t1leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t1leptfunc)
t1hadtfunc = TF1('t1hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t1hadtfunc)
t2leptfunc = TF1('t2leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t2leptfunc)
t2hadtfunc = TF1('t2hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t2hadtfunc)
t3leptfunc = TF1('t3leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t3leptfunc)
t3hadtfunc = TF1('t3hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(t3hadtfunc)
t3hadWfunc = TF1('t3hadWfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',Wmass_low,Wmass_high); allfuncs.append(t3hadWfunc)
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
print 'Drawing type-3 leptonic top mass...'
tree.Draw('leptcorprefitM>>t3leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3)')
print 'Drawing type-3 hadronic top mass...'
tree.Draw('hadtcorprefitM>>t3hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3)')
print 'Drawing type-3 hadronic W mass...'
tree.Draw('hadWcorprefitM>>t3hadWM('+str(nWbins)+','+str(Wmass_low)+','+str(Wmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3)')
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
t3leptM_canv = TCanvas('t3leptM_canv','t3leptM_canv',1100,900); allcanvs.append(t3leptM_canv)
t3hadtM_canv = TCanvas('t3hadtM_canv','t3hadtM_canv',1100,900); allcanvs.append(t3hadtM_canv)
t3hadWM_canv = TCanvas('t3hadWM_canv','t3hadWM_canv',1100,900); allcanvs.append(t3hadWM_canv)

#draw plots on canvases and fit with gaussians
for i in range(len(allhistos)) :
	leg = TLegend(0.1,0.7,0.48,0.9)
	histo = allhistos[i]
	func = allfuncs[i]
	print 'Doing '+str(allhistos[i].GetName())
	allcanvs[i].cd()
	histo.Draw()
	histo.Fit(func)
	func.SetLineColor(kRed)
	func.SetLineWidth(4)
	func.Draw('SAME')
	gausfunc = TF1('gausfunc','gaus',tmass_low,tmass_high)
	histo.Fit(gausfunc)
	gausfunc.SetLineColor(kBlue)
	gausfunc.SetLineWidth(4)
	gausfunc.Draw('SAME')
	leg.AddEntry(func,"Lorentzian","L")
	leg.AddEntry(gausfunc,"Gaussian","L")
	leg.Draw('SAME')


#save plots and canvases
outfile.cd()
for histo in allhistos :
	histo.Write()
for canv in allcanvs :
	canv.Write()