from ROOT import *
import CMS_lumi, tdrstyle
import glob
from optparse import OptionParser

tdrstyle.setTDRStyle()

parser = OptionParser()
#Run options
parser.add_option('--input', 	  type='string', action='store', default='input', dest='input',	   	  
	help='Path to input file holding list of files to run on (default: "input")')
(options, args) = parser.parse_args()

#output file
output_filename = 'qq_gg_comparison_plots.root'
outputfile = TFile(output_filename,'recreate')

#list of input files
inputfileslist = options.input
if not options.input.endswith('.txt') :
	inputfileslist += '.txt'
print 'Using input file '+inputfileslist+''
inputfileslist = open(inputfileslist,'r')

#set up plots
qq_plots = []
qg_plots = []
gg_plots = []
qq_c = TH1D('qq_c','; c*; Fraction/0.02',100,-1.0,1.0); qq_plots.append(qq_c)
qg_c = TH1D('qg_c','; c*; Fraction/0.02',100,-1.0,1.0); qg_plots.append(qg_c)
gg_c = TH1D('gg_c','; c*; Fraction/0.02',100,-1.0,1.0); gg_plots.append(gg_c)
qq_x = TH1D('qq_x','; x_{F}; Fraction/0.02',100,-1.0,1.0); qq_plots.append(qq_x)
qg_x = TH1D('qg_x','; x_{F}; Fraction/0.02',100,-1.0,1.0); qg_plots.append(qg_x)
gg_x = TH1D('gg_x','; x_{F}; Fraction/0.02',100,-1.0,1.0); gg_plots.append(gg_x)
qq_x_us = TH1D('qq_x_us','; |x_{F}|; Fraction/0.02',50,0.0,1.0); qq_plots.append(qq_x_us)
qg_x_us = TH1D('qg_x_us','; |x_{F}|; Fraction/0.02',50,0.0,1.0); qg_plots.append(qg_x_us)
gg_x_us = TH1D('gg_x_us','; |x_{F}|; Fraction/0.02',50,0.0,1.0); gg_plots.append(gg_x_us)
qq_M = TH1D('qq_M','; m_{t#bar{t}} [GeV]; Fraction/10 GeV',500,0.,5000.); qq_plots.append(qq_M)
qg_M = TH1D('qg_M','; m_{t#bar{t}} [GeV]; Fraction/10 GeV',500,0.,5000.); qg_plots.append(qg_M)
gg_M = TH1D('gg_M','; m_{t#bar{t}} [GeV]; Fraction/10 GeV',500,0.,5000.); gg_plots.append(gg_M)
all_plots = [qq_plots,qg_plots,gg_plots]
#plot attributes
for plot in qq_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlue); plot.SetLineStyle(1)
for plot in qg_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kRed); plot.SetLineStyle(7)
for plot in gg_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlack); plot.SetLineStyle(10)
for plotgroup in all_plots :
	for plot in plotgroup :
		plot.SetStats(0)

#open files and add to plots
for inputfile in inputfileslist :
	thisfile = TFile.Open(inputfile.rstrip())
	#print inputfile #DEBUG
	#print thisfile #DEBUG
	#print thisfile.Get('EventCounter/cstar') #DEBUG
	new_qq_c = (thisfile.Get('EventCounter/cstar')).ProjectionY('new_qq_c',1,1); qq_c.Add(new_qq_c)
	new_qg_c = (thisfile.Get('EventCounter/cstar')).ProjectionY('new_qg_c',2,2); qg_c.Add(new_qg_c)
	new_gg_c = (thisfile.Get('EventCounter/cstar')).ProjectionY('new_gg_c',3,3); gg_c.Add(new_gg_c)
	new_qq_x = (thisfile.Get('EventCounter/x_F')).ProjectionY('new_qq_x',1,1); qq_x.Add(new_qq_x)
	new_qg_x = (thisfile.Get('EventCounter/x_F')).ProjectionY('new_qg_x',2,2); qg_x.Add(new_qg_x)
	new_gg_x = (thisfile.Get('EventCounter/x_F')).ProjectionY('new_gg_x',3,3); gg_x.Add(new_gg_x)
	new_qq_M = (thisfile.Get('EventCounter/M')).ProjectionY('new_qq_M',1,1); qq_M.Add(new_qq_M)
	new_qg_M = (thisfile.Get('EventCounter/M')).ProjectionY('new_qg_M',2,2); qg_M.Add(new_qg_M)
	new_gg_M = (thisfile.Get('EventCounter/M')).ProjectionY('new_gg_M',3,3); gg_M.Add(new_gg_M)
	thisfile.Close()

#add to the unsigned Feyman x plots by copyig bin contents
for bin in range(1,qq_x.GetNbinsX()+1) :
	signedbincenter = qq_x.GetXaxis().GetBinCenter(bin)
	unsignedbin = qq_x_us.GetXaxis().FindFixBin(abs(signedbincenter))
	unsignedbincontent_qq = qq_x_us.GetBinContent(unsignedbin)
	unsignedbincontent_qg = qg_x_us.GetBinContent(unsignedbin)
	unsignedbincontent_gg = gg_x_us.GetBinContent(unsignedbin)
	to_add_qq = qq_x.GetBinContent(bin)
	to_add_qg = qg_x.GetBinContent(bin)
	to_add_gg = gg_x.GetBinContent(bin)
	qq_x_us.SetBinContent(unsignedbin,unsignedbincontent_qq+to_add_qq)
	qg_x_us.SetBinContent(unsignedbin,unsignedbincontent_qg+to_add_qg)
	gg_x_us.SetBinContent(unsignedbin,unsignedbincontent_gg+to_add_gg)

#save plots
outputfile.cd()
for plotgroup in all_plots :
	for plot in plotgroup :
		plot.Write()

#curtail the range of the M plots
qq_M.GetXaxis().SetRangeUser(300.,1300.)
qg_M.GetXaxis().SetRangeUser(300.,1300.)
gg_M.GetXaxis().SetRangeUser(300.,1300.)

#normalize plots
for plotgroup in all_plots :
	for plot in plotgroup :
		plot.Scale(1./plot.Integral())

#Make Canvases
c_canv = TCanvas('c_canv','c_canv',1100,900)
x_canv = TCanvas('x_canv','x_canv',1100,900)
x_canv_unsigned = TCanvas('x_canv_unsigned','x_canv_unsigned',1100,900)
M_canv = TCanvas('M_canv','M_canv',1100,900)

#Make legend
leg = TLegend(0.65,0.53,0.90,0.78)
leg.AddEntry(qq_plots[0],'q#bar{q}','L')
leg.AddEntry(qg_plots[0],'qg','L')
leg.AddEntry(gg_plots[0],'gg','L')

#set axis limits
qq_c.SetMinimum(0.); qq_c.SetMaximum(0.032)
qq_x.SetMinimum(0.); qq_x.SetMaximum(0.100)
qq_x_us.SetMinimum(0.); qq_x_us.SetMaximum(0.155)
qq_M.SetMinimum(0.); qq_M.SetMaximum(0.058)

#Plot plots
c_canv.cd()
qq_c.Draw("HIST")
qg_c.Draw("SAME HIST")
gg_c.Draw("SAME HIST")
leg.Draw("SAME")
x_canv.cd()
qq_x.Draw("HIST")
qg_x.Draw("SAME HIST")
gg_x.Draw("SAME HIST")
leg.Draw("SAME")
x_canv_unsigned.cd()
qq_x_us.Draw("HIST")
qg_x_us.Draw("SAME HIST")
gg_x_us.Draw("SAME HIST")
leg.Draw("SAME")
M_canv.cd()
qq_M.Draw("HIST")
qg_M.Draw("SAME HIST")
gg_M.Draw("SAME HIST")
leg.Draw("SAME")

#plot the CMS_Lumi lines on the canvases
iPeriod = 0 #free form since it's only simulation
iPos = 11 #iPos = 10*(alignment) + position (1/2/3 = left/center/right)
if iPos==0 : CMS_lumi.relPosX = 0.12
#CMS_lumi.cmsText = "POWHEG"
CMS_lumi.cmsText = "CMS"
#CMS_lumi.cmsTextSize = 0.85
CMS_lumi.cmsTextSize = 1.1
CMS_lumi.writeExtraText = 1
#CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.extraText = "Simulation"
#CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.CMS_lumi(c_canv, iPeriod, iPos)
CMS_lumi.CMS_lumi(x_canv, iPeriod, iPos)
CMS_lumi.CMS_lumi(x_canv_unsigned, iPeriod, iPos)
CMS_lumi.CMS_lumi(M_canv, iPeriod, 33)

#update canvases
c_canv.Update()
x_canv.Update()
x_canv_unsigned.Update()
M_canv.Update()

#save canvases
outputfile.cd()
c_canv.Write()
x_canv.Write()
x_canv_unsigned.Write()
M_canv.Write()

#save and close output file
outputfile.Write()
outputfile.Close()