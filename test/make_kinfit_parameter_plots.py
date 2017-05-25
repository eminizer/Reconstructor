from ROOT import *
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--file',  metavar='F', type='string', action='store', 
                              default='../total_ttree_files/mcatnlo_TT_skim_all.root', 
                              dest='file',  help='') ## TTree file path
parser.add_option('--outname',  metavar='F', type='string', action='store', 
                              default='kinfit_parameter_plots', 
                              dest='outname',  help='') ## name for output file
(options, args) = parser.parse_args()

#set up output file
outfile = TFile(options.outname+'.root','recreate')

#open the input file
infile = TFile(options.file)

#get the tree from the file
tree = infile.Get('tree')

#Define cuts and draw into the histograms
metfilters = 'metfilters==1'
onelepton = 'onelepton==1'
isolepton = 'isolepton==1'
jetcuts = 'jetcuts==1'
fullselection = 'fullselection==1'

all_cuts = '('+fullselection+' && eventType<2)'

print 'all_cuts = '+all_cuts

weights = '(35867.*weight*sf_pileup*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas)'

tree.Draw("par_0>>vpz_1_int(100,-500.,500.)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("par_1>>lsf_1_int(100,0.95,1.05)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("par_2>>lbsf_1_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("par_3>>hsf1_1_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("par_4>>hsf2_1_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("chi2>>chi2_1_int(100,-75.,125.)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")

tree.Draw("par_0>>vpz_2_int(100,-500.,500.)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("par_1>>lsf_2_int(100,0.95,1.05)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("par_2>>lbsf_2_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("par_3>>hsf1_2_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("par_4>>hsf2_2_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("par_5>>hsf3_2_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("chi2>>chi2_2_int(100,-75.,125.)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")

tree.Draw("par_0>>vpz_3_int(100,-500.,500.)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("par_1>>lsf_3_int(100,0.95,1.05)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("par_2>>lbsf_3_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("par_3>>hsf1_3_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("par_4>>hsf2_3_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("par_5>>hsf3_3_int(60,0.6,1.2)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("chi2>>chi2_3_int(100,-75.,125.)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")


#get the histograms and set titles
vpz_1 = TH1D('vpz_1','neutrino p_{Z} (type-1 events); p_{Z} (GeV); events',100,-500.,500.)
lsf_1 = TH1D('lsf_1','lepton momentum scalefactor (type-1 events); #lambda; events',100,0.95,1.05)
lbsf_1 = TH1D('lbsf_1','Jet momentum scalefactors (type-1 events); #lambda; events',60,0.6,1.2)
hsf1_1 = TH1D('hsf1_1','Jet momentum scalefactors (type-1 events); #lambda; events',60,0.6,1.2)
hsf2_1 = TH1D('hsf2_1','Jet momentum scalefactors (type-1 events); #lambda; events',60,0.6,1.2)
chi2_1 = TH1D('chi2_1','Kinematic fit #chi^{2} (type-1 events); #chi^{2}; events',100,-75.,125.)

vpz_2 = TH1D('vpz_2','neutrino p_{Z} (type-2 events); p_{Z} (GeV); events',100,-500.,500.)
lsf_2 = TH1D('lsf_2','lepton momentum scalefactor (type-2 events); #lambda; events',100,0.95,1.05)
lbsf_2 = TH1D('lbsf_2','Jet momentum scalefactors (type-2 events); #lambda; events',60,0.6,1.2)
hsf1_2 = TH1D('hsf1_2','Jet momentum scalefactors (type-2 events); #lambda; events',60,0.6,1.2)
hsf2_2 = TH1D('hsf2_2','Jet momentum scalefactors (type-2 events); #lambda; events',60,0.6,1.2)
hsf3_2 = TH1D('hsf3_2','Jet momentum scalefactors (type-2 events); #lambda; events',60,0.6,1.2)
chi2_2 = TH1D('chi2_2','Kinematic fit #chi^{2} (type-2 events); #chi^{2}; events',100,-75.,125.)

vpz_3 = TH1D('vpz_3','neutrino p_{Z} (type-3 events); p_{Z} (GeV); events',100,-500.,500.)
lsf_3 = TH1D('lsf_3','lepton momentum scalefactor (type-3 events); #lambda; events',100,0.95,1.05)
lbsf_3 = TH1D('lbsf_3','Jet momentum scalefactors (type-3 events); #lambda; events',60,0.6,1.2)
hsf1_3 = TH1D('hsf1_3','Jet momentum scalefactors (type-3 events); #lambda; events',60,0.6,1.2)
hsf2_3 = TH1D('hsf2_3','Jet momentum scalefactors (type-3 events); #lambda; events',60,0.6,1.2)
hsf3_3 = TH1D('hsf3_3','Jet momentum scalefactors (type-3 events); #lambda; events',60,0.6,1.2)
chi2_3 = TH1D('chi2_3','Kinematic fit #chi^{2} (type-3 events); #chi^{2}; events',100,-75.,125.)

vpz_1.SetLineWidth(3); vpz_1.SetMarkerStyle(20);
lsf_1.SetLineWidth(3); lsf_1.SetMarkerStyle(20);
lbsf_1.SetLineWidth(3); lbsf_1.SetMarkerStyle(20); lbsf_1.SetLineColor(kRed); lbsf_1.SetMarkerColor(kRed)
hsf1_1.SetLineWidth(3); hsf1_1.SetMarkerStyle(20); hsf1_1.SetLineColor(kBlue); hsf1_1.SetMarkerColor(kBlue)
hsf2_1.SetLineWidth(3); hsf2_1.SetMarkerStyle(20); hsf2_1.SetLineColor(kGreen+2); hsf2_1.SetMarkerColor(kGreen+2)
chi2_1.SetLineWidth(3); chi2_1.SetMarkerStyle(20);

vpz_2.SetLineWidth(3); vpz_2.SetMarkerStyle(20);
lsf_2.SetLineWidth(3); lsf_2.SetMarkerStyle(20);
lbsf_2.SetLineWidth(3); lbsf_2.SetMarkerStyle(20); lbsf_2.SetLineColor(kRed); lbsf_2.SetMarkerColor(kRed)
hsf1_2.SetLineWidth(3); hsf1_2.SetMarkerStyle(20); hsf1_2.SetLineColor(kBlue); hsf1_2.SetMarkerColor(kBlue)
hsf2_2.SetLineWidth(3); hsf2_2.SetMarkerStyle(20); hsf2_2.SetLineColor(kGreen+2); hsf2_2.SetMarkerColor(kGreen+2)
hsf3_2.SetLineWidth(3); hsf3_2.SetMarkerStyle(20); hsf3_2.SetLineColor(kOrange); hsf3_2.SetMarkerColor(kOrange)
chi2_2.SetLineWidth(3); chi2_2.SetMarkerStyle(20);

vpz_3.SetLineWidth(3); vpz_3.SetMarkerStyle(20);
lsf_3.SetLineWidth(3); lsf_3.SetMarkerStyle(20);
lbsf_3.SetLineWidth(3); lbsf_3.SetMarkerStyle(20); lbsf_3.SetLineColor(kRed); lbsf_3.SetMarkerColor(kRed)
hsf1_3.SetLineWidth(3); hsf1_3.SetMarkerStyle(20); hsf1_3.SetLineColor(kBlue); hsf1_3.SetMarkerColor(kBlue)
hsf2_3.SetLineWidth(3); hsf2_3.SetMarkerStyle(20); hsf2_3.SetLineColor(kGreen+2); hsf2_3.SetMarkerColor(kGreen+2)
hsf3_3.SetLineWidth(3); hsf3_3.SetMarkerStyle(20); hsf3_3.SetLineColor(kOrange); hsf3_3.SetMarkerColor(kOrange)
chi2_3.SetLineWidth(3); chi2_3.SetMarkerStyle(20);

vpz_1.Add(gROOT.FindObject('vpz_1_int'))
lsf_1.Add(gROOT.FindObject('lsf_1_int'))
lbsf_1.Add(gROOT.FindObject('lbsf_1_int'))
hsf1_1.Add(gROOT.FindObject('hsf1_1_int'))
hsf2_1.Add(gROOT.FindObject('hsf2_1_int'))
chi2_1.Add(gROOT.FindObject('chi2_1_int'))

vpz_2.Add(gROOT.FindObject('vpz_2_int'))
lsf_2.Add(gROOT.FindObject('lsf_2_int'))
lbsf_2.Add(gROOT.FindObject('lbsf_2_int'))
hsf1_2.Add(gROOT.FindObject('hsf1_2_int'))
hsf2_2.Add(gROOT.FindObject('hsf2_2_int'))
hsf3_2.Add(gROOT.FindObject('hsf3_2_int'))
chi2_2.Add(gROOT.FindObject('chi2_2_int'))

vpz_3.Add(gROOT.FindObject('vpz_3_int'))
lsf_3.Add(gROOT.FindObject('lsf_3_int'))
lbsf_3.Add(gROOT.FindObject('lbsf_3_int'))
hsf1_3.Add(gROOT.FindObject('hsf1_3_int'))
hsf2_3.Add(gROOT.FindObject('hsf2_3_int'))
hsf3_3.Add(gROOT.FindObject('hsf3_3_int'))
chi2_3.Add(gROOT.FindObject('chi2_3_int'))

#make the TCanvases
vpz_canv_1 = TCanvas('vpz_canv_1','vpz_canv_1',1100,900)
lsf_canv_1 = TCanvas('lsf_canv_1','lsf_canv_1',1100,900)
hsf_canv_1 = TCanvas('hsf_canv_1','hsf_canv_1',1100,900)
chi2_canv_1 = TCanvas('chi2_canv_1','chi2_canv_1',1100,900)

vpz_canv_2 = TCanvas('vpz_canv_2','vpz_canv_2',1100,900)
lsf_canv_2 = TCanvas('lsf_canv_2','lsf_canv_2',1100,900)
hsf_canv_2 = TCanvas('hsf_canv_2','hsf_canv_2',1100,900)
chi2_canv_2 = TCanvas('chi2_canv_2','chi2_canv_2',1100,900)

vpz_canv_3 = TCanvas('vpz_canv_3','vpz_canv_3',1100,900)
lsf_canv_3 = TCanvas('lsf_canv_3','lsf_canv_3',1100,900)
hsf_canv_3 = TCanvas('hsf_canv_3','hsf_canv_3',1100,900)
chi2_canv_3 = TCanvas('chi2_canv_3','chi2_canv_3',1100,900)

#make the legends
leg_1 = TLegend(0.1,0.7,0.48,0.9)
leg_1.AddEntry(lbsf_1,'#lambda_{1} (leptonic-side b-jet)','PE')
leg_1.AddEntry(hsf1_1,'#lambda_{2} (hadronic subjet 1)','PE')
leg_1.AddEntry(hsf2_1,'#lambda_{3} (hadronic subjet 2)','PE')

leg_2 = TLegend(0.1,0.7,0.48,0.9)
leg_2.AddEntry(lbsf_2,'#lambda_{1} (leptonic-side b-jet)','PE')
leg_2.AddEntry(hsf1_2,'#lambda_{2} (hadronic-side jet 1)','PE')
leg_2.AddEntry(hsf2_2,'#lambda_{3} (hadronic-side jet 2)','PE')
leg_2.AddEntry(hsf3_2,'#lambda_{4} (hadronic-side jet 3)','PE')

#plot the plots
vpz_canv_1.cd()
vpz_1.Draw("PE")
lsf_canv_1.cd()
lsf_1.Draw("PE")
hsf_canv_1.cd()
lbsf_1.Draw("PE")
hsf1_1.Draw("PE SAME")
hsf2_1.Draw("PE SAME")
leg_1.Draw("SAME")
chi2_canv_1.cd()
chi2_1.Draw("PE")

vpz_canv_2.cd()
vpz_2.Draw("PE")
lsf_canv_2.cd()
lsf_2.Draw("PE")
hsf_canv_2.cd()
lbsf_2.Draw("PE")
hsf1_2.Draw("PE SAME")
hsf2_2.Draw("PE SAME")
hsf3_2.Draw("PE SAME")
leg_2.Draw("SAME")
chi2_canv_2.cd()
chi2_2.Draw("PE")

vpz_canv_3.cd()
vpz_3.Draw("PE")
lsf_canv_3.cd()
lsf_3.Draw("PE")
hsf_canv_3.cd()
lbsf_3.Draw("PE")
hsf1_3.Draw("PE SAME")
hsf2_3.Draw("PE SAME")
hsf3_3.Draw("PE SAME")
leg_2.Draw("SAME")
chi2_canv_3.cd()
chi2_3.Draw("PE")

#write the canvas
outfile.cd()
vpz_canv_1.Write()
lsf_canv_1.Write()
hsf_canv_1.Write()
chi2_canv_1.Write()
vpz_canv_2.Write()
lsf_canv_2.Write()
hsf_canv_2.Write()
chi2_canv_2.Write()
vpz_canv_3.Write()
lsf_canv_3.Write()
hsf_canv_3.Write()
chi2_canv_3.Write()
outfile.Write()
outfile.Close()
