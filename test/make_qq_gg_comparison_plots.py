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

#directory name
inputdir = '/eos/uscms/store/user/eminizer/B2GTTress_Mar16/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/B2GAnaFW_80X_V2p4_PR66_Mar16_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170316_154342/0000'

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
qq_c = TH1D('qq_c','; c*; fraction/0.02',100,-1.0,1.0); qq_plots.append(qq_c)
qg_c = TH1D('qg_c','; c*; fraction/0.02',100,-1.0,1.0); qg_plots.append(qg_c)
gg_c = TH1D('gg_c','; c*; fraction/0.02',100,-1.0,1.0); gg_plots.append(gg_c)
qq_x = TH1D('qq_x','; x_{F}; fraction/0.02',100,-1.0,1.0); qq_plots.append(qq_x)
qg_x = TH1D('qg_x','; x_{F}; fraction/0.02',100,-1.0,1.0); qg_plots.append(qg_x)
gg_x = TH1D('gg_x','; x_{F}; fraction/0.02',100,-1.0,1.0); gg_plots.append(gg_x)
qq_M = TH1D('qq_M','; M (GeV); fraction/25 GeV',500,0.,5000.); qq_plots.append(qq_M)
qg_M = TH1D('qg_M','; M (GeV); fraction/25 GeV',500,0.,5000.); qg_plots.append(qg_M)
gg_M = TH1D('gg_M','; M (GeV); fraction/25 GeV',500,0.,5000.); gg_plots.append(gg_M)
all_plots = [qq_plots,qg_plots,gg_plots]
#plot attributes
for plot in qq_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlue)
for plot in qg_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kRed)
for plot in gg_plots :
	plot.SetLineWidth(3); plot.SetMarkerStyle(25); plot.SetLineColor(kBlack)
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
M_canv = TCanvas('M_canv','M_canv',1100,900)

#Make legend
leg = TLegend(0.62,0.67,0.75,0.8)
leg.AddEntry(qq_plots[0],'q#bar{q}','L')
leg.AddEntry(qg_plots[0],'qg','L')
leg.AddEntry(gg_plots[0],'gg','L')

#set axis limits
qq_c.SetMinimum(0.); qq_c.SetMaximum(0.032)
qq_x.SetMinimum(0.); qq_x.SetMaximum(0.100)
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
M_canv.cd()
qq_M.Draw("HIST")
qg_M.Draw("SAME HIST")
gg_M.Draw("SAME HIST")
leg.Draw("SAME")

#plot the CMS_Lumi lines on the canvases
iPeriod = 0 #free form since it's only simulation
iPos = 11 #iPos = 10*(alignment) + position (1/2/3 = left/center/right)
if iPos==0 : CMS_lumi.relPosX = 0.12
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"
#CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.CMS_lumi(c_canv, iPeriod, iPos)
CMS_lumi.CMS_lumi(x_canv, iPeriod, iPos)
CMS_lumi.CMS_lumi(M_canv, iPeriod, 33)

#update canvases
c_canv.Update()
x_canv.Update()
M_canv.Update()

#save canvases
outputfile.cd()
c_canv.Write()
x_canv.Write()
M_canv.Write()

#save and close output file
outputfile.Write()
outputfile.Close()