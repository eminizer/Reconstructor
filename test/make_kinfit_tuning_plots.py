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
tmass_low=100.; tmass_high=350.; ntbins=25
allhistos = []
leptM = TH1F('leptM','Leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(leptM)
hadtM = TH1F('hadtM','Hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(hadtM)
for histo in allhistos :
	histo.SetDirectory(0)

#declare the Lorentzian functions
allfuncs = []
leptfunc = TF1('leptfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(leptfunc)
hadtfunc = TF1('hadtfunc','[0]/((x*x-[1]*[1])*(x*x-[1]*[1])+[1]*[1]*[2]*[2])',tmass_low,tmass_high); allfuncs.append(hadtfunc)
for func in allfuncs :
	func.SetParName(0,'const'); func.SetParName(1,'mean'); func.SetParName(2,'width')
	func.SetParameter(0,3500000.); func.SetParameter(1,173.); func.SetParameter(2,20.)

#plot plots
weights = '12917.*weight'
com_cuts = 'fullselection==1 && ismatchable==1 && hadt_isttagged==1'
tree.Draw('leptcorprefitM>>leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+')')
tree.Draw('hadtcorprefitM>>hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+')')
for histo in allhistos :
	histo.Add(gROOT.FindObject(histo.GetName()))
	histo.Scale(1./histo.Integral())
	histo.SetMarkerStyle(20)
	histo.GetYaxis().SetTitleOffset(1.2)

#declare canvases
allcanvs = []
leptM_canv = TCanvas('leptM_canv','leptM_canv',1100,900); allcanvs.append(leptM_canv)
hadtM_canv = TCanvas('hadtM_canv','hadtM_canv',1100,900); allcanvs.append(hadtM_canv)

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