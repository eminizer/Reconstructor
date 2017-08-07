from ROOT import *
import CMS_lumi, tdrstyle
from glob import glob
from datetime import date
import os
from optparse import OptionParser

tdrstyle.setTDRStyle()

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--sample',  metavar='F', type='string', action='store', 
                              dest='sample',  help='')
parser.add_option('--outtag',  metavar='F', type='string', action='store', 
                              default='', 
                              dest='outtag') ## name for output file
(options, args) = parser.parse_args()
sample = options.sample
outtag = options.outtag
if sample==None :
	optstring = raw_input('samplename (outfilename_append) >  ')
	optlist = optstring.split()
	if len(optlist)<1 :
		print 'must input sample name'
		exit()
	sample = optlist[0]
	if len(optlist)>1 :
		outtag = optlist[1]

#open the input files in a chain
if not os.path.isdir('../'+sample) :
	print 'sample name %s is not valid!'%(sample)
	exit()
infileslist = glob('../'+sample+'/aggregated_*.root')
fullchain = TChain('tree')
for f in infileslist :
	fullchain.Add(f)

#set up output file
outname = 'reconstruction_comparison_plots_'+sample+'_'+str(date.today())
if outtag!='' :
	outname+='_'+outtag
outname+='.root'
outfile = TFile(outname,'recreate')

#skim the chain
com_cuts = '(fullselection==1 && eventType<2 && ((eventTopology<3 && ((lepflavor==1 && (lep_relPt>30. || lep_dR>0.4)) || (lepflavor==2 && lep_relPt>30. && lep_dR>0.4))) || (eventTopology==3 && ((lepflavor==1 && (lep_relPt>30. || lep_dR>0.4)) || (lepflavor==2 && lep_relPt>20. && lep_dR>0.4)))))'
chain = fullchain.CopyTree(com_cuts)

#Define cuts and draw into the histograms
metfilters = 'metfilters==1'
onelepton = 'onelepton==1'
isolepton = 'isolepton==1'
jetcuts = 'jetcuts==1'
fullselection = 'fullselection==1'

weights = '(((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_lep_trk*sf_btag_eff*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'

chain.Draw("cstar:cstar_MC>>cstar_comp_int_1(40,-1.,1.,40,-1.,1.)",weights+"*("+com_cuts+' && eventType==0 && eventTopology==1)',"COLZ")
chain.Draw("(cstar-cstar_MC)/cstar_MC>>cstar_res_int_1(100,-2.,2.)",weights+"*("+com_cuts+' && eventType==0 && eventTopology==1)',"COLZ")
chain.Draw("x_F:x_F_MC>>x_F_comp_int_1(30,0.,0.6,30,0.,0.6)",weights+"*("+com_cuts+' && eventTopology==1)',"COLZ")
chain.Draw("(abs(x_F)-abs(x_F_MC))/abs(x_F_MC)>>x_F_res_int_1(100,-2.,2.)",weights+"*("+com_cuts+' && eventTopology==1)',"COLZ")
chain.Draw("M:M_MC>>M_comp_int_1(40,750.,2750.,40,750.,2750.)",weights+"*("+com_cuts+' && eventTopology==1)',"COLZ")
chain.Draw("(M-M_MC)/M_MC>>M_res_int_1(100,-2.,2.)",weights+"*("+com_cuts+' && eventTopology==1)')

chain.Draw("cstar:cstar_MC>>cstar_comp_int_2(40,-1.,1.,40,-1.,1.)",weights+"*("+com_cuts+' && eventType==0 && eventTopology==2)',"COLZ")
chain.Draw("(cstar-cstar_MC)/cstar_MC>>cstar_res_int_2(100,-2.,2.)",weights+"*("+com_cuts+' && eventType==0 && eventTopology==2)',"COLZ")
chain.Draw("x_F:x_F_MC>>x_F_comp_int_2(30,0.,0.6,30,0.,0.6)",weights+"*("+com_cuts+' && eventTopology==2)',"COLZ")
chain.Draw("(abs(x_F)-abs(x_F_MC))/abs(x_F_MC)>>x_F_res_int_2(100,-2.,2.)",weights+"*("+com_cuts+' && eventTopology==2)',"COLZ")
chain.Draw("M:M_MC>>M_comp_int_2(40,750.,2750.,40,750.,2750.)",weights+"*("+com_cuts+' && eventTopology==2)',"COLZ")
chain.Draw("(M-M_MC)/M_MC>>M_res_int_2(100,-2.,2.)",weights+"*("+com_cuts+' && eventTopology==2)')

chain.Draw("cstar:cstar_MC>>cstar_comp_int_3(40,-1.,1.,40,-1.,1.)",weights+"*("+com_cuts+' && eventType==0 && eventTopology==3)',"COLZ")
chain.Draw("(cstar-cstar_MC)/cstar_MC>>cstar_res_int_3(100,-2.,2.)",weights+"*("+com_cuts+' && eventType==0 && eventTopology==3)',"COLZ")
chain.Draw("x_F:x_F_MC>>x_F_comp_int_3(30,0.,0.6,30,0.,0.6)",weights+"*("+com_cuts+' && eventTopology==3)',"COLZ")
chain.Draw("(abs(x_F)-abs(x_F_MC))/abs(x_F_MC)>>x_F_res_int_3(100,-2.,2.)",weights+"*("+com_cuts+' && eventTopology==3)',"COLZ")
chain.Draw("M:M_MC>>M_comp_int_3(40,400.,1200.,40,400.,1200.)",weights+"*("+com_cuts+' && eventTopology==3)',"COLZ")
chain.Draw("(M-M_MC)/M_MC>>M_res_int_3(100,-2.,2.)",weights+"*("+com_cuts+' && eventTopology==3)')


#get the histograms and set titles
all_histos = []
cstar_comp_1 = TH2D('cstar_comp_1','Reconstructed vs. Generated c* (type-1 events); c* (generated); c_{r} (reconstructed)',40,-1.,1.,40,-1.,1.) ; all_histos.append(cstar_comp_1)
cstar_res_1 = TH1D('cstar_res_1','observable resolution (type-1 events); (reconstructed-generated)/generated; fraction',100,-2.,2.); all_histos.append(cstar_res_1)
x_F_comp_1 = TH2D('x_F_comp_1','Reconstructed vs. Generated |x_{F}| (type-1 events); |x_{F}| (generated); |x_{r}| (reconstructed)',30,0.,0.6,30,0.,0.6) ; all_histos.append(x_F_comp_1)
x_F_res_1 = TH1D('x_F_res_1','|x_{F}| resolution (type-1 events); (|x_{r}|-|x_{F}|)/|x_{F}|',100,-2.,2.); all_histos.append(x_F_res_1)
M_comp_1 = TH2D('M_comp_1','Reconstructed vs. Generated M (type-1 events); M (generated) (GeV); M_{r} (reconstructed) (GeV)',40,750.,2750.,40,750.,2750.); all_histos.append(M_comp_1)
M_res_1 = TH1D('M_res_1','M resolution (type-1 events); (M_{r} - M)/M',100,-2.,2.); all_histos.append(M_res_1)
cstar_comp_2 = TH2D('cstar_comp_2','Reconstructed vs. Generated c* (type-2 events); c* (generated); c_{r} (reconstructed)',40,-1.,1.,40,-1.,1.) ; all_histos.append(cstar_comp_2)
cstar_res_2 = TH1D('cstar_res_2','observable resolution (type-2 events); (reconstructed-generated)/generated; fraction',100,-2.,2.); all_histos.append(cstar_res_2)
x_F_comp_2 = TH2D('x_F_comp_2','Reconstructed vs. Generated |x_{F}| (type-2 events); |x_{F}| (generated); |x_{r}| (reconstructed)',30,0.,0.6,30,0.,0.6) ; all_histos.append(x_F_comp_2)
x_F_res_2 = TH1D('x_F_res_2','|x_{F}| resolution (type-2 events); (|x_{r}|-|x_{F}|)/|x_{F}|',100,-2.,2.); all_histos.append(x_F_res_2)
M_comp_2 = TH2D('M_comp_2','Reconstructed vs. Generated M (type-2 events); M (generated) (GeV); M_{r} (reconstructed) (GeV)',40,750.,2750.,40,750.,2750.); all_histos.append(M_comp_2)
M_res_2 = TH1D('M_res_2','M resolution (type-2 events); (M_{r} - M)/M',100,-2.,2.); all_histos.append(M_res_2)
cstar_comp_3 = TH2D('cstar_comp_3','Reconstructed vs. Generated c* (type-3 events); c* (generated); c_{r} (reconstructed)',40,-1.,1.,40,-1.,1.) ; all_histos.append(cstar_comp_3)
cstar_res_3 = TH1D('cstar_res_3','observable resolution (type-3 events); (reconstructed-generated)/generated; fraction',100,-2.,2.); all_histos.append(cstar_res_3)
x_F_comp_3 = TH2D('x_F_comp_3','Reconstructed vs. Generated |x_{F}| (type-3 events); |x_{F}| (generated); |x_{r}| (reconstructed)',30,0.,0.6,30,0.,0.6) ; all_histos.append(x_F_comp_3)
x_F_res_3 = TH1D('x_F_res_3','|x_{F}| resolution (type-3 events); (|x_{r}|-|x_{F}|)/|x_{F}|',100,-2.,2.); all_histos.append(x_F_res_3)
M_comp_3 = TH2D('M_comp_3','Reconstructed vs. Generated M (type-3 events); M (generated) (GeV); M_{r} (reconstructed) (GeV)',40,400.,1200.,40,400.,1200.); all_histos.append(M_comp_3)
M_res_3 = TH1D('M_res_3','M resolution (type-3 events); (M_{r} - M)/M',100,-2.,2.); all_histos.append(M_res_3)

for histo in all_histos :
	histo.SetDirectory(0)

cstar_res_1.SetLineWidth(3); cstar_res_1.SetMarkerStyle(20); cstar_res_1.SetLineColor(kRed); cstar_res_1.SetMarkerColor(kRed)
cstar_res_2.SetLineWidth(3); cstar_res_2.SetMarkerStyle(20); cstar_res_2.SetLineColor(kRed); cstar_res_2.SetMarkerColor(kRed)
cstar_res_3.SetLineWidth(3); cstar_res_3.SetMarkerStyle(20); cstar_res_3.SetLineColor(kRed); cstar_res_3.SetMarkerColor(kRed)
x_F_res_1.SetLineWidth(3); x_F_res_1.SetMarkerStyle(20); x_F_res_1.SetLineColor(kBlue); x_F_res_1.SetMarkerColor(kBlue)
x_F_res_2.SetLineWidth(3); x_F_res_2.SetMarkerStyle(20); x_F_res_2.SetLineColor(kBlue); x_F_res_2.SetMarkerColor(kBlue)
x_F_res_3.SetLineWidth(3); x_F_res_3.SetMarkerStyle(20); x_F_res_3.SetLineColor(kBlue); x_F_res_3.SetMarkerColor(kBlue)
M_res_1.SetLineWidth(3); M_res_1.SetMarkerStyle(20); M_res_1.SetLineColor(kGreen+2); M_res_1.SetMarkerColor(kGreen+2)
M_res_2.SetLineWidth(3); M_res_2.SetMarkerStyle(20); M_res_2.SetLineColor(kGreen+2); M_res_2.SetMarkerColor(kGreen+2)
M_res_3.SetLineWidth(3); M_res_3.SetMarkerStyle(20); M_res_3.SetLineColor(kGreen+2); M_res_3.SetMarkerColor(kGreen+2)
cstar_comp_1.Add(gROOT.FindObject('cstar_comp_int_1'))
cstar_res_1.Add(gROOT.FindObject('cstar_res_int_1'))
x_F_comp_1.Add(gROOT.FindObject('x_F_comp_int_1'))
x_F_res_1.Add(gROOT.FindObject('x_F_res_int_1'))
M_comp_1.Add(gROOT.FindObject('M_comp_int_1'))
M_res_1.Add(gROOT.FindObject('M_res_int_1'))
cstar_comp_2.Add(gROOT.FindObject('cstar_comp_int_2'))
cstar_res_2.Add(gROOT.FindObject('cstar_res_int_2'))
x_F_comp_2.Add(gROOT.FindObject('x_F_comp_int_2'))
x_F_res_2.Add(gROOT.FindObject('x_F_res_int_2'))
M_comp_2.Add(gROOT.FindObject('M_comp_int_2'))
M_res_2.Add(gROOT.FindObject('M_res_int_2'))
cstar_comp_3.Add(gROOT.FindObject('cstar_comp_int_3'))
cstar_res_3.Add(gROOT.FindObject('cstar_res_int_3'))
x_F_comp_3.Add(gROOT.FindObject('x_F_comp_int_3'))
x_F_res_3.Add(gROOT.FindObject('x_F_res_int_3'))
M_comp_3.Add(gROOT.FindObject('M_comp_int_3'))
M_res_3.Add(gROOT.FindObject('M_res_int_3'))

cstar_res_1.Scale(1./cstar_res_1.Integral())
x_F_res_1.Scale(1./x_F_res_1.Integral())
M_res_1.Scale(1./M_res_1.Integral())
cstar_res_2.Scale(1./cstar_res_2.Integral())
x_F_res_2.Scale(1./x_F_res_2.Integral())
M_res_2.Scale(1./M_res_2.Integral())
cstar_res_3.Scale(1./cstar_res_3.Integral())
x_F_res_3.Scale(1./x_F_res_3.Integral())
M_res_3.Scale(1./M_res_3.Integral())

#make the TCanvases
all_canvs = []
cstar_comp_canv_1 = TCanvas('cstar_comp_canv_1','cstar_comp_canv_1',1100,900); all_canvs.append(cstar_comp_canv_1)
x_F_comp_canv_1 = TCanvas('x_F_comp_canv_1','x_F_comp_canv_1',1100,900); all_canvs.append(x_F_comp_canv_1)
M_comp_canv_1 = TCanvas('M_comp_canv_1','M_comp_canv_1',1100,900); all_canvs.append(M_comp_canv_1)
res_canv_1 = TCanvas('res_canv_1','res_canv_1',1100,900); all_canvs.append(res_canv_1)
cstar_comp_canv_2 = TCanvas('cstar_comp_canv_2','cstar_comp_canv_2',1100,900); all_canvs.append(cstar_comp_canv_2)
x_F_comp_canv_2 = TCanvas('x_F_comp_canv_2','x_F_comp_canv_2',1100,900); all_canvs.append(x_F_comp_canv_2)
M_comp_canv_2 = TCanvas('M_comp_canv_2','M_comp_canv_2',1100,900); all_canvs.append(M_comp_canv_2)
res_canv_2 = TCanvas('res_canv_2','res_canv_2',1100,900); all_canvs.append(res_canv_2)
cstar_comp_canv_3 = TCanvas('cstar_comp_canv_3','cstar_comp_canv_3',1100,900); all_canvs.append(cstar_comp_canv_3)
x_F_comp_canv_3 = TCanvas('x_F_comp_canv_3','x_F_comp_canv_3',1100,900); all_canvs.append(x_F_comp_canv_3)
M_comp_canv_3 = TCanvas('M_comp_canv_3','M_comp_canv_3',1100,900); all_canvs.append(M_comp_canv_3)
res_canv_3 = TCanvas('res_canv_3','res_canv_3',1100,900); all_canvs.append(res_canv_3)
for canv in all_canvs :
	canv.SetLeftMargin(0.14); canv.SetRightMargin(0.05) 
	canv.SetTopMargin(0.10); canv.SetBottomMargin(0.12)

#make the legend
leg = TLegend(0.1,0.7,0.48,0.9)
leg.AddEntry(cstar_res_1,'c*','PE')
leg.AddEntry(x_F_res_1,'|x_F|','PE')
leg.AddEntry(M_res_1,'M','PE')

#set up the CMS Lumi objects
iPeriod = 0 #free form since it's only simulation
iPos = 0 #out of frame by default (iPos = 10*(alignment) + position (1/2/3 = left/center/right))
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

#plot the plots
cstar_comp_canv_1.cd()
cstar_comp_1.Draw("COLZ")
CMS_lumi.CMS_lumi(cstar_comp_canv_1, iPeriod, iPos)
x_F_comp_canv_1.cd()
x_F_comp_1.Draw("COLZ")
CMS_lumi.CMS_lumi(x_F_comp_canv_1, iPeriod, iPos)
M_comp_canv_1.cd()
M_comp_1.Draw("COLZ")
CMS_lumi.CMS_lumi(M_comp_canv_1, iPeriod, iPos)
res_canv_1.cd()
cstar_res_1.Draw("PE")
x_F_res_1.Draw("PE SAME")
M_res_1.Draw("PE SAME")
leg.Draw("SAME")
CMS_lumi.CMS_lumi(res_canv_1, iPeriod, 11)

cstar_comp_canv_2.cd()
cstar_comp_2.Draw("COLZ")
CMS_lumi.CMS_lumi(cstar_comp_canv_2, iPeriod, iPos)
x_F_comp_canv_2.cd()
x_F_comp_2.Draw("COLZ")
CMS_lumi.CMS_lumi(x_F_comp_canv_2, iPeriod, iPos)
M_comp_canv_2.cd()
M_comp_2.Draw("COLZ")
CMS_lumi.CMS_lumi(M_comp_canv_2, iPeriod, iPos)
res_canv_2.cd()
cstar_res_2.Draw("PE")
x_F_res_2.Draw("PE SAME")
M_res_2.Draw("PE SAME")
leg.Draw("SAME")
CMS_lumi.CMS_lumi(res_canv_2, iPeriod, 11)

cstar_comp_canv_3.cd()
cstar_comp_3.Draw("COLZ")
CMS_lumi.CMS_lumi(cstar_comp_canv_3, iPeriod, iPos)
x_F_comp_canv_3.cd()
x_F_comp_3.Draw("COLZ")
CMS_lumi.CMS_lumi(x_F_comp_canv_3, iPeriod, iPos)
M_comp_canv_3.cd()
M_comp_3.Draw("COLZ")
CMS_lumi.CMS_lumi(M_comp_canv_3, iPeriod, iPos)
res_canv_3.cd()
cstar_res_3.Draw("PE")
x_F_res_3.Draw("PE SAME")
M_res_3.Draw("PE SAME")
leg.Draw("SAME")
CMS_lumi.CMS_lumi(res_canv_3, iPeriod, 11)

#write the canvas
outfile.cd()
for canv in all_canvs :
	canv.Write()
outfile.Write()
outfile.Close()
