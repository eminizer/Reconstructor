from ROOT import *
from optparse import OptionParser
import glob

#output file
output_filename = 'qq_gg_comparison_plots.root'
outputfile = TFile(output_filename,'recreate')

#directory names
prepend = '/uscms_data/d3/eminizer/ttbar_13TeV/CMSSW_8_0_20/src/Analysis/Reconstructor/test/'
qq_dirs = ['mcatnlo_qq_semilep_TT_just_MC']
gg_dirs = ['mcatnlo_gg_semilep_TT_just_MC']

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
qq_plots = []
qg_plots = []
gg_plots = []
qq_c = TH1D('qq_c','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); qq_plots.append(qq_c)
qg_c = TH1D('qg_c','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); qg_plots.append(qg_c)
gg_c = TH1D('gg_c','cos(#theta *) (c*) distributions; c*; fraction of events/0.05',40,-1.0,1.0); gg_plots.append(gg_c)
qq_x = TH1D('qq_x','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); qq_plots.append(qq_x)
qg_x = TH1D('qg_x','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); qg_plots.append(qg_x)
gg_x = TH1D('gg_x','Feynman x (|x_{F}|) distributions; |x_{F}|; fraction of events/0.025',40,0.0,1.0); gg_plots.append(gg_x)
qq_M = TH1D('qq_M','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); qq_plots.append(qq_M)
qg_M = TH1D('qg_M','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); qg_plots.append(qg_M)
gg_M = TH1D('gg_M','Combined t#bar{t} mass (M) distributions; M (GeV); fraction of events/25 GeV',80,750.,2750.); gg_plots.append(gg_M)
all_plots = [qq_plots,qg_plots,gg_plots]
#plot attributes
for plot in qq_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlue)
for plot in qg_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kRed)
for plot in gg_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlack)

#Make trees and draw into plots
qq_tree = qq_chain.CloneTree()
gg_tree = gg_chain.CloneTree()
weightstring = '12917.*weight'; cutstring = 'weight!=0'
qq_tree.Draw('cstar_MC>>qq_c_p(40,-1.0,1.0)','('+weightstring+')*('+cutstring+')'); 					qq_c.Add(gROOT.FindObject('qq_c_p'))
qq_tree.Draw('abs(x_F_MC)>>qq_x_p(40,0.0,1.0)','('+weightstring+')*('+cutstring+')'); 					qq_x.Add(gROOT.FindObject('qq_x_p'))
qq_tree.Draw('M_MC>>qq_M_p(80,750.,2750.)','('+weightstring+')*('+cutstring+')'); 						qq_M.Add(gROOT.FindObject('qq_M_p'))
gg_tree.Draw('cstar_MC>>qg_c_p(40,-1.0,1.0)','('+weightstring+')*('+cutstring+' && addTwice==0)'); 	qg_c.Add(gROOT.FindObject('qg_c_p'))
gg_tree.Draw('abs(x_F_MC)>>qg_x_p(40,0.0,1.0)','('+weightstring+')*('+cutstring+' && addTwice==0)'); qg_x.Add(gROOT.FindObject('qg_x_p'))
gg_tree.Draw('M_MC>>qg_M_p(80,750.,2750.)','('+weightstring+')*('+cutstring+' && addTwice==0)'); 	qg_M.Add(gROOT.FindObject('qg_M_p'))
gg_tree.Draw('cstar_MC>>gg_c_p(40,-1.0,1.0)','('+weightstring+')*('+cutstring+' && addTwice==1)'); 	gg_c.Add(gROOT.FindObject('gg_c_p'))
gg_tree.Draw('abs(x_F_MC)>>gg_x_p(40,0.0,1.0)','('+weightstring+')*('+cutstring+' && addTwice==1)'); gg_x.Add(gROOT.FindObject('gg_x_p'))
gg_tree.Draw('M_MC>>gg_M_p(80,750.,2750.)','('+weightstring+')*('+cutstring+' && addTwice==1)'); 	gg_M.Add(gROOT.FindObject('gg_M_p'))

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
leg.AddEntry(qq_plots[0],'q#bar{q}','L')
leg.AddEntry(qg_plots[0],'qg','L')
leg.AddEntry(gg_plots[0],'gg','L')

#Plot plots
c_canv.cd()
qq_c.Draw("HIST")
qg_c.Draw("SAME HIST")
gg_c.Draw("SAME HIST")
leg.Draw()
x_canv.cd()
qq_x.Draw("HIST")
qg_x.Draw("SAME HIST")
gg_x.Draw("SAME HIST")
leg.Draw()
M_canv.cd()
qq_M.Draw("HIST")
qg_M.Draw("SAME HIST")
gg_M.Draw("SAME HIST")
leg.Draw()

#save canvases
outputfile.cd()
c_canv.Write()
x_canv.Write()
M_canv.Write()