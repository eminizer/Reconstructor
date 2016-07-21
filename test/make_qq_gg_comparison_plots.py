from ROOT import *
from optparse import OptionParser
import glob

#output file
output_filename = 'qq_gg_comparison_plots.root'
outputfile = TFile(output_filename,'recreate')

#directory names
prepend = '/uscms_data/d3/eminizer/ttbar_run2/CMSSW_7_2_0/src/Analysis/TTree_Maker/test/'
qq_dirs = ['Powheg_qq_semilep_TT_Mtt_700_to_1000','Powheg_qq_semilep_TT_Mtt_1000_to_Inf']
gg_dirs = ['Powheg_gg_semilep_TT_Mtt_700_to_1000','Powheg_gg_semilep_TT_Mtt_1000_to_Inf']

#set up chains
qq_chain = TChain('tree') 
gg_chain = TChain('tree')
qq_filenamelist = []
for qqdir in qq_dirs :
	qq_filenamelist += glob.glob(prepend+qqdir+'/*_tree.root')
for qqfilename in qq_filenamelist :
	if qqfilename.find('skim')==-1 and qqfilename.find('JES')==-1 and qqfilename.find('JER')==-1 :
		qq_chain.Add(qqfilename)
gg_filenamelist = []
for ggdir in gg_dirs :
	gg_filenamelist += glob.glob(prepend+ggdir+'/*_tree.root')
for ggfilename in gg_filenamelist :
	if ggfilename.find('skim')==-1 and ggfilename.find('JES')==-1 and ggfilename.find('JER')==-1 :
		gg_chain.Add(ggfilename)

#set up plots
qq_plots_ct10 = []
qg_plots_ct10 = []
gg_plots_ct10 = []
qq_c_ct10 = TH1D('qq_c_ct10','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); qq_plots_ct10.append(qq_c_ct10)
qg_c_ct10 = TH1D('qg_c_ct10','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); qg_plots_ct10.append(qg_c_ct10)
gg_c_ct10 = TH1D('gg_c_ct10','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); gg_plots_ct10.append(gg_c_ct10)
qq_x_ct10 = TH1D('qq_x_ct10','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); qq_plots_ct10.append(qq_x_ct10)
qg_x_ct10 = TH1D('qg_x_ct10','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); qg_plots_ct10.append(qg_x_ct10)
gg_x_ct10 = TH1D('gg_x_ct10','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); gg_plots_ct10.append(gg_x_ct10)
qq_M_ct10 = TH1D('qq_M_ct10','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); qq_plots_ct10.append(qq_M_ct10)
qg_M_ct10 = TH1D('qg_M_ct10','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); qg_plots_ct10.append(qg_M_ct10)
gg_M_ct10 = TH1D('gg_M_ct10','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); gg_plots_ct10.append(gg_M_ct10)
qq_plots_cteq = []
qg_plots_cteq = []
gg_plots_cteq = []
qq_c_cteq = TH1D('qq_c_cteq','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); qq_plots_cteq.append(qq_c_cteq)
qg_c_cteq = TH1D('qg_c_cteq','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); qg_plots_cteq.append(qg_c_cteq)
gg_c_cteq = TH1D('gg_c_cteq','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); gg_plots_cteq.append(gg_c_cteq)
qq_x_cteq = TH1D('qq_x_cteq','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); qq_plots_cteq.append(qq_x_cteq)
qg_x_cteq = TH1D('qg_x_cteq','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); qg_plots_cteq.append(qg_x_cteq)
gg_x_cteq = TH1D('gg_x_cteq','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); gg_plots_cteq.append(gg_x_cteq)
qq_M_cteq = TH1D('qq_M_cteq','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); qq_plots_cteq.append(qq_M_cteq)
qg_M_cteq = TH1D('qg_M_cteq','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); qg_plots_cteq.append(qg_M_cteq)
gg_M_cteq = TH1D('gg_M_cteq','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); gg_plots_cteq.append(gg_M_cteq)
qq_plots_gjr = []
qg_plots_gjr = []
gg_plots_gjr = []
qq_c_gjr = TH1D('qq_c_gjr','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); qq_plots_gjr.append(qq_c_gjr)
qg_c_gjr = TH1D('qg_c_gjr','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); qg_plots_gjr.append(qg_c_gjr)
gg_c_gjr = TH1D('gg_c_gjr','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); gg_plots_gjr.append(gg_c_gjr)
qq_x_gjr = TH1D('qq_x_gjr','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); qq_plots_gjr.append(qq_x_gjr)
qg_x_gjr = TH1D('qg_x_gjr','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); qg_plots_gjr.append(qg_x_gjr)
gg_x_gjr = TH1D('gg_x_gjr','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); gg_plots_gjr.append(gg_x_gjr)
qq_M_gjr = TH1D('qq_M_gjr','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); qq_plots_gjr.append(qq_M_gjr)
qg_M_gjr = TH1D('qg_M_gjr','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); qg_plots_gjr.append(qg_M_gjr)
gg_M_gjr = TH1D('gg_M_gjr','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); gg_plots_gjr.append(gg_M_gjr)
all_plots = [qq_plots_ct10,qg_plots_ct10,gg_plots_ct10,qq_plots_cteq,qg_plots_cteq,gg_plots_cteq,qq_plots_gjr,qg_plots_gjr,gg_plots_gjr]
#plot attributes
for plot in qq_plots_ct10 :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlue)
for plot in qg_plots_ct10 :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kRed)
for plot in gg_plots_ct10 :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlack)
for plot in qq_plots_cteq :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlue); plot.SetLineStyle(2)
for plot in qg_plots_cteq :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kRed); plot.SetLineStyle(2)
for plot in gg_plots_cteq :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlack); plot.SetLineStyle(2)
for plot in qq_plots_gjr :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlue); plot.SetLineStyle(7)
for plot in qg_plots_gjr :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kRed); plot.SetLineStyle(7)
for plot in gg_plots_gjr :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlack); plot.SetLineStyle(7)

#Make trees and draw into plots
qq_tree = qq_chain.CloneTree()
gg_tree = gg_chain.CloneTree()
qq_tree.Draw('cstar_MC>>qq_c_p_ct10(40,-1.0,1.0)','19748*weight*hadt_pt>0.'); 					qq_c_ct10.Add(gROOT.FindObject('qq_c_p_ct10'))
qq_tree.Draw('abs(x_F_MC)>>qq_x_p_ct10(40,0.0,1.0)','19748*weight*hadt_pt>0.'); 					qq_x_ct10.Add(gROOT.FindObject('qq_x_p_ct10'))
qq_tree.Draw('M_MC>>qq_M_p_ct10(80,750.,2750.)','19748*weight*hadt_pt>0.'); 						qq_M_ct10.Add(gROOT.FindObject('qq_M_p_ct10'))
gg_tree.Draw('cstar_MC>>qg_c_p_ct10(40,-1.0,1.0)','19748*weight*(hadt_pt>0. && addTwice==0)'); 	qg_c_ct10.Add(gROOT.FindObject('qg_c_p_ct10'))
gg_tree.Draw('abs(x_F_MC)>>qg_x_p_ct10(40,0.0,1.0)','19748*weight*(hadt_pt>0. && addTwice==0)'); qg_x_ct10.Add(gROOT.FindObject('qg_x_p_ct10'))
gg_tree.Draw('M_MC>>qg_M_p_ct10(80,750.,2750.)','19748*weight*(hadt_pt>0. && addTwice==0)'); 	qg_M_ct10.Add(gROOT.FindObject('qg_M_p_ct10'))
gg_tree.Draw('cstar_MC>>gg_c_p_ct10(40,-1.0,1.0)','19748*weight*(hadt_pt>0. && addTwice==1)'); 	gg_c_ct10.Add(gROOT.FindObject('gg_c_p_ct10'))
gg_tree.Draw('abs(x_F_MC)>>gg_x_p_ct10(40,0.0,1.0)','19748*weight*(hadt_pt>0. && addTwice==1)'); gg_x_ct10.Add(gROOT.FindObject('gg_x_p_ct10'))
gg_tree.Draw('M_MC>>gg_M_p_ct10(80,750.,2750.)','19748*weight*(hadt_pt>0. && addTwice==1)'); 	gg_M_ct10.Add(gROOT.FindObject('gg_M_p_ct10'))
qq_tree.Draw('cstar_MC>>qq_c_p_cteq(40,-1.0,1.0)','19748*weight*(cteq_weights[0]/CT10_weights[0])*hadt_pt>0.'); 					qq_c_cteq.Add(gROOT.FindObject('qq_c_p_cteq'))
qq_tree.Draw('abs(x_F_MC)>>qq_x_p_cteq(40,0.0,1.0)','19748*weight*(cteq_weights[0]/CT10_weights[0])*hadt_pt>0.'); 					qq_x_cteq.Add(gROOT.FindObject('qq_x_p_cteq'))
qq_tree.Draw('M_MC>>qq_M_p_cteq(80,750.,2750.)','19748*weight*(cteq_weights[0]/CT10_weights[0])*hadt_pt>0.'); 						qq_M_cteq.Add(gROOT.FindObject('qq_M_p_cteq'))
gg_tree.Draw('cstar_MC>>qg_c_p_cteq(40,-1.0,1.0)','19748*weight*(cteq_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==0)'); 	qg_c_cteq.Add(gROOT.FindObject('qg_c_p_cteq'))
gg_tree.Draw('abs(x_F_MC)>>qg_x_p_cteq(40,0.0,1.0)','19748*weight*(cteq_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==0)'); qg_x_cteq.Add(gROOT.FindObject('qg_x_p_cteq'))
gg_tree.Draw('M_MC>>qg_M_p_cteq(80,750.,2750.)','19748*weight*(cteq_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==0)'); 	qg_M_cteq.Add(gROOT.FindObject('qg_M_p_cteq'))
gg_tree.Draw('cstar_MC>>gg_c_p_cteq(40,-1.0,1.0)','19748*weight*(cteq_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==1)'); 	gg_c_cteq.Add(gROOT.FindObject('gg_c_p_cteq'))
gg_tree.Draw('abs(x_F_MC)>>gg_x_p_cteq(40,0.0,1.0)','19748*weight*(cteq_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==1)'); gg_x_cteq.Add(gROOT.FindObject('gg_x_p_cteq'))
gg_tree.Draw('M_MC>>gg_M_p_cteq(80,750.,2750.)','19748*weight*(cteq_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==1)'); 	gg_M_cteq.Add(gROOT.FindObject('gg_M_p_cteq'))
qq_tree.Draw('cstar_MC>>qq_c_p_gjr(40,-1.0,1.0)','19748*weight*(GJR_weights[0]/CT10_weights[0])*hadt_pt>0.'); 					qq_c_gjr.Add(gROOT.FindObject('qq_c_p_gjr'))
qq_tree.Draw('abs(x_F_MC)>>qq_x_p_gjr(40,0.0,1.0)','19748*weight*(GJR_weights[0]/CT10_weights[0])*hadt_pt>0.'); 					qq_x_gjr.Add(gROOT.FindObject('qq_x_p_gjr'))
qq_tree.Draw('M_MC>>qq_M_p_gjr(80,750.,2750.)','19748*weight*(GJR_weights[0]/CT10_weights[0])*hadt_pt>0.'); 						qq_M_gjr.Add(gROOT.FindObject('qq_M_p_gjr'))
gg_tree.Draw('cstar_MC>>qg_c_p_gjr(40,-1.0,1.0)','19748*weight*(GJR_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==0)'); 	qg_c_gjr.Add(gROOT.FindObject('qg_c_p_gjr'))
gg_tree.Draw('abs(x_F_MC)>>qg_x_p_gjr(40,0.0,1.0)','19748*weight*(GJR_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==0)'); qg_x_gjr.Add(gROOT.FindObject('qg_x_p_gjr'))
gg_tree.Draw('M_MC>>qg_M_p_gjr(80,750.,2750.)','19748*weight*(GJR_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==0)'); 	qg_M_gjr.Add(gROOT.FindObject('qg_M_p_gjr'))
gg_tree.Draw('cstar_MC>>gg_c_p_gjr(40,-1.0,1.0)','19748*weight*(GJR_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==1)'); 	gg_c_gjr.Add(gROOT.FindObject('gg_c_p_gjr'))
gg_tree.Draw('abs(x_F_MC)>>gg_x_p_gjr(40,0.0,1.0)','19748*weight*(GJR_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==1)'); gg_x_gjr.Add(gROOT.FindObject('gg_x_p_gjr'))
gg_tree.Draw('M_MC>>gg_M_p_gjr(80,750.,2750.)','19748*weight*(GJR_weights[0]/CT10_weights[0])*(hadt_pt>0. && addTwice==1)'); 	gg_M_gjr.Add(gROOT.FindObject('gg_M_p_gjr'))

#save plots
outputfile.cd()
for plotgroup in all_plots :
	for plot in plotgroup :
		plot.Write()

#normalize plots
for plotgroup in all_plots :
	for plot in plotgroup :
		plot.Scale(1./plot.Integral())

#Make Canvases
c_canv = TCanvas('c_canv','c_canv',1100,900)
x_canv = TCanvas('x_canv','x_canv',1100,900)
M_canv = TCanvas('M_canv','M_canv',1100,900)

#Make legend
leg = TLegend(0.62,0.67,0.9,0.9)
leg.AddEntry(qq_plots_ct10[0],'q#bar{q}, CT10 PDFs','L')
leg.AddEntry(qq_plots_cteq[0],'q#bar{q}, cteq66 PDFs','L')
leg.AddEntry(qq_plots_gjr[0],'q#bar{q}, GJR08VFnloE PDFs','L')
leg.AddEntry(qg_plots_ct10[0],'qg, CT10 PDFs','L')
leg.AddEntry(qg_plots_cteq[0],'qg, cteq66 PDFs','L')
leg.AddEntry(qg_plots_gjr[0],'qg, GJR08VFnloE PDFs','L')
leg.AddEntry(gg_plots_ct10[0],'gg, CT10 PDFs','L')
leg.AddEntry(gg_plots_cteq[0],'gg, cteq66 PDFs','L')
leg.AddEntry(gg_plots_gjr[0],'gg, GJR08VFnloE PDFs','L')

#Plot plots
c_canv.cd()
qq_c_ct10.Draw("HIST")
qg_c_ct10.Draw("SAME HIST")
gg_c_ct10.Draw("SAME HIST")
qq_c_cteq.Draw("SAME HIST")
qg_c_cteq.Draw("SAME HIST")
gg_c_cteq.Draw("SAME HIST")
qq_c_gjr.Draw("SAME HIST")
qg_c_gjr.Draw("SAME HIST")
gg_c_gjr.Draw("SAME HIST")
leg.Draw()
x_canv.cd()
qq_x_ct10.Draw("HIST")
qg_x_ct10.Draw("SAME HIST")
gg_x_ct10.Draw("SAME HIST")
qq_x_cteq.Draw("SAME HIST")
qg_x_cteq.Draw("SAME HIST")
gg_x_cteq.Draw("SAME HIST")
qq_x_gjr.Draw("SAME HIST")
qg_x_gjr.Draw("SAME HIST")
gg_x_gjr.Draw("SAME HIST")
leg.Draw()
M_canv.cd()
qq_M_ct10.Draw("HIST")
qg_M_ct10.Draw("SAME HIST")
gg_M_ct10.Draw("SAME HIST")
qq_M_cteq.Draw("SAME HIST")
qg_M_cteq.Draw("SAME HIST")
gg_M_cteq.Draw("SAME HIST")
qq_M_gjr.Draw("SAME HIST")
qg_M_gjr.Draw("SAME HIST")
gg_M_gjr.Draw("SAME HIST")
leg.Draw()

#save canvases
outputfile.cd()
c_canv.Write()
x_canv.Write()
M_canv.Write()