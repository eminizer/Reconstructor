from ROOT import *
import CMS_lumi, tdrstyle
from glob import glob
from datetime import date
import os
from optparse import OptionParser

tdrstyle.setTDRStyle()

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--ttbar_file',  metavar='F', type='string', action='store', 
                              dest='ttbar_file',  help='')
parser.add_option('--wjets_file',  metavar='F', type='string', action='store', 
                              dest='wjets_file',  help='')
parser.add_option('--outtag',  metavar='F', type='string', action='store', 
                              default='', 
                              dest='outtag') ## name for output file
(options, args) = parser.parse_args()
ttbar_file = options.ttbar_file
wjets_file = options.wjets_file
outtag = options.outtag

#open the input files in chains
if not (os.path.isfile(ttbar_file) and os.path.isfile(wjets_file)) :
	print 'sample name %s is not valid!'%(sample)
	exit()
fullttbarchain = TChain('tree')
fullwjetschain = TChain('tree')
fullttbarchain.Add(ttbar_file)
fullwjetschain.Add(wjets_file)

#set up output file
outname = 'wjets_2d_region_plots_'+str(date.today())
if outtag!='' :
	outname+='_'+outtag
outname+='.root'
outfile = TFile(outname,'recreate')

#skim the chains
com_cuts = '(fullselection==1 || wjets_cr_selection==1) && eventTopology<3'
print 'skimming ttbar chain...'
ttbarchain = fullttbarchain.CopyTree(com_cuts)
print 'Done.\nskimming W+Jets chain...'
wjetschain = fullwjetschain.CopyTree(com_cuts)
print 'Done.'

#Draw into the histograms
weights = '(((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_btag_eff*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'

ttbarchain.Draw("scaled_lept_M:chi2>>t1_sig_int(50,-50.,200.,40,50.,450.)",weights+"*("+com_cuts+' && eventTopology==1)',"COLZ")
wjetschain.Draw("scaled_lept_M:chi2>>t1_bg_int(50,-50.,200.,40,50.,450.)",weights+"*("+com_cuts+' && eventTopology==1)',"COLZ")
ttbarchain.Draw("scaled_lept_M:chi2>>t2_sig_int(50,-50.,200.,40,50.,450.)",weights+"*("+com_cuts+' && eventTopology==2)',"COLZ")
wjetschain.Draw("scaled_lept_M:chi2>>t2_bg_int(50,-50.,200.,40,50.,450.)",weights+"*("+com_cuts+' && eventTopology==2)',"COLZ")

#get the histograms and set titles
all_histos = []
t1_sig = TH2D('t1_sig','; #chi^{2}; reconstructed t_{lep} mass [GeV]',50,-50.,200.,40,50.,450.) ; all_histos.append(t1_sig)
t1_bg = TH2D('t1_bg','; #chi^{2}; reconstructed t_{lep} mass [GeV]',50,-50.,200.,40,50.,450.) ; all_histos.append(t1_bg)
t2_sig = TH2D('t2_sig','; #chi^{2}; reconstructed t_{lep} mass [GeV]',50,-50.,200.,40,50.,450.) ; all_histos.append(t2_sig)
t2_bg = TH2D('t2_bg','; #chi^{2}; reconstructed t_{lep} mass [GeV]',50,-50.,200.,40,50.,450.) ; all_histos.append(t2_bg)

for histo in all_histos :
	histo.SetDirectory(0)
	histo.SetStats(0)
	histo.GetYaxis().SetTitleOffset(1.1)

t1_sig.Add(gROOT.FindObject('t1_sig_int'))
t1_bg.Add(gROOT.FindObject('t1_bg_int'))
t2_sig.Add(gROOT.FindObject('t2_sig_int'))
t2_bg.Add(gROOT.FindObject('t2_bg_int'))

#make the TCanvases
all_canvs = []
t1_sig_2D = TCanvas('t1_sig_2D','t1_sig_2D',1100,900); all_canvs.append(t1_sig_2D)
t1_bg_2D = TCanvas('t1_bg_2D','t1_bg_2D',1100,900); all_canvs.append(t1_bg_2D)
t2_sig_2D = TCanvas('t2_sig_2D','t2_sig_2D',1100,900); all_canvs.append(t2_sig_2D)
t2_bg_2D = TCanvas('t2_bg_2D','t2_bg_2D',1100,900); all_canvs.append(t2_bg_2D)

for canv in all_canvs :
	canv.SetLeftMargin(0.14); canv.SetRightMargin(0.05) 
	canv.SetTopMargin(0.10); canv.SetBottomMargin(0.12)

#set up the CMS Lumi objects
iPeriod = 0 #free form since it's only simulation
iPos = 11 #on the left by default (iPos = 10*(alignment) + position (1/2/3 = left/center/right))
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

#plot the plots
t1_sig_2D.cd()
t1_sig.Draw("COLZ")
t1_sig_chantxt = TLatex()
t1_sig_chantxt.SetNDC()
t1_sig_chantxt.SetTextAngle(0)
t1_sig_chantxt.SetTextColor(kBlack)
t1_sig_chantxt.SetTextFont(42)
t1_sig_chantxt.SetTextAlign(11) 
t1_sig_chantxt.SetTextSize(0.6*0.11)
t1_sig_chantxt.DrawLatex(0.16,1-0.11+0.2*0.11,'type-1 t#bar{t} events')
CMS_lumi.CMS_lumi(t1_sig_2D, iPeriod, iPos)
t1_sig_2D.Update()

t1_bg_2D.cd()
t1_bg.Draw("COLZ")
t1_bg_chantxt = TLatex()
t1_bg_chantxt.SetNDC()
t1_bg_chantxt.SetTextAngle(0)
t1_bg_chantxt.SetTextColor(kBlack)
t1_bg_chantxt.SetTextFont(42)
t1_bg_chantxt.SetTextAlign(11) 
t1_bg_chantxt.SetTextSize(0.6*0.11)
t1_bg_chantxt.DrawLatex(0.16,1-0.11+0.2*0.11,'type-1 W+Jets events')
CMS_lumi.CMS_lumi(t1_bg_2D, iPeriod, iPos)
t1_bg_2D.Update()

t2_sig_2D.cd()
t2_sig.Draw("COLZ")
t2_sig_chantxt = TLatex()
t2_sig_chantxt.SetNDC()
t2_sig_chantxt.SetTextAngle(0)
t2_sig_chantxt.SetTextColor(kBlack)
t2_sig_chantxt.SetTextFont(42)
t2_sig_chantxt.SetTextAlign(11) 
t2_sig_chantxt.SetTextSize(0.6*0.11)
t2_sig_chantxt.DrawLatex(0.16,1-0.11+0.2*0.11,'type-2 t#bar{t} events')
CMS_lumi.CMS_lumi(t2_sig_2D, iPeriod, iPos)
t2_sig_2D.Update()

t2_bg_2D.cd()
t2_bg.Draw("COLZ")
t2_bg_chantxt = TLatex()
t2_bg_chantxt.SetNDC()
t2_bg_chantxt.SetTextAngle(0)
t2_bg_chantxt.SetTextColor(kBlack)
t2_bg_chantxt.SetTextFont(42)
t2_bg_chantxt.SetTextAlign(11) 
t2_bg_chantxt.SetTextSize(0.6*0.11)
t2_bg_chantxt.DrawLatex(0.16,1-0.11+0.2*0.11,'type-2 W+Jets events')
CMS_lumi.CMS_lumi(t2_bg_2D, iPeriod, iPos)
t2_bg_2D.Update()




#write the canvas
outfile.cd()
for canv in all_canvs :
	canv.Write()
outfile.Write()
outfile.Close()
