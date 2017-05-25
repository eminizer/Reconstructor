from ROOT import *
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--file',  metavar='F', type='string', action='store', 
                              default='../total_ttree_files/mcatnlo_TT_skim_all.root', 
                              dest='file',  help='') ## TTree file path
parser.add_option('--outname',  metavar='F', type='string', action='store', 
                              default='reconstruction_comparison_plots', 
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

all_cuts = '('+fullselection+' && eventType==0)'

print 'all_cuts = '+all_cuts

weights = '(35867.*weight*sf_pileup*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas)'

tree.Draw("cstar:cstar_MC>>cstar_comp_int_1(40,-1.,1.,40,-1.,1.)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("(cstar-cstar_MC)/cstar_MC>>cstar_res_int_1(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("x_F:x_F_MC>>x_F_comp_int_1(30,0.,0.6,30,0.,0.6)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("(abs(x_F)-abs(x_F_MC))/abs(x_F_MC)>>x_F_res_int_1(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("M:M_MC>>M_comp_int_1(40,750.,2750.,40,750.,2750.)",weights+"*("+all_cuts+' && eventTopology==1)',"COLZ")
tree.Draw("(M-M_MC)/M_MC>>M_res_int_1(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==1)')

tree.Draw("cstar:cstar_MC>>cstar_comp_int_2(40,-1.,1.,40,-1.,1.)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("(cstar-cstar_MC)/cstar_MC>>cstar_res_int_2(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("x_F:x_F_MC>>x_F_comp_int_2(30,0.,0.6,30,0.,0.6)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("(abs(x_F)-abs(x_F_MC))/abs(x_F_MC)>>x_F_res_int_2(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("M:M_MC>>M_comp_int_2(40,750.,2750.,40,750.,2750.)",weights+"*("+all_cuts+' && eventTopology==2)',"COLZ")
tree.Draw("(M-M_MC)/M_MC>>M_res_int_2(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==2)')

tree.Draw("cstar:cstar_MC>>cstar_comp_int_3(40,-1.,1.,40,-1.,1.)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("(cstar-cstar_MC)/cstar_MC>>cstar_res_int_3(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("x_F:x_F_MC>>x_F_comp_int_3(30,0.,0.6,30,0.,0.6)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("(abs(x_F)-abs(x_F_MC))/abs(x_F_MC)>>x_F_res_int_3(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("M:M_MC>>M_comp_int_3(40,750.,2750.,40,750.,2750.)",weights+"*("+all_cuts+' && eventTopology==3)',"COLZ")
tree.Draw("(M-M_MC)/M_MC>>M_res_int_3(100,-2.,2.)",weights+"*("+all_cuts+' && eventTopology==3)')


#get the histograms and set titles
cstar_comp_1 = TH2D('cstar_comp_1','Reconstructed vs. Generated c* (type-1 events); c* (generated); c_{r} (reconstructed)',40,-1.,1.,40,-1.,1.) 
cstar_res_1 = TH1D('cstar_res_1','observable resolution (type-1 events); (reconstructed-generated)/generated; fraction',100,-2.,2.)
x_F_comp_1 = TH2D('x_F_comp_1','Reconstructed vs. Generated |x_{F}| (type-1 events); |x_{F}| (generated); |x_{r}| (reconstructed)',30,0.,0.6,30,0.,0.6) 
x_F_res_1 = TH1D('x_F_res_1','|x_{F}| resolution (type-1 events); (|x_{r}|-|x_{F}|)/|x_{F}|',100,-2.,2.)
M_comp_1 = TH2D('M_comp_1','Reconstructed vs. Generated M (type-1 events); M (generated) (GeV); M_{r} (reconstructed) (GeV)',40,750.,2750.,40,750.,2750.)
M_res_1 = TH1D('M_res_1','M resolution (type-1 events); (M_{r} - M)/M',100,-2.,2.)
cstar_comp_2 = TH2D('cstar_comp_2','Reconstructed vs. Generated c* (type-2 events); c* (generated); c_{r} (reconstructed)',40,-1.,1.,40,-1.,1.) 
cstar_res_2 = TH1D('cstar_res_2','observable resolution (type-2 events); (reconstructed-generated)/generated; fraction',100,-2.,2.)
x_F_comp_2 = TH2D('x_F_comp_2','Reconstructed vs. Generated |x_{F}| (type-2 events); |x_{F}| (generated); |x_{r}| (reconstructed)',30,0.,0.6,30,0.,0.6) 
x_F_res_2 = TH1D('x_F_res_2','|x_{F}| resolution (type-2 events); (|x_{r}|-|x_{F}|)/|x_{F}|',100,-2.,2.)
M_comp_2 = TH2D('M_comp_2','Reconstructed vs. Generated M (type-2 events); M (generated) (GeV); M_{r} (reconstructed) (GeV)',40,750.,2750.,40,750.,2750.)
M_res_2 = TH1D('M_res_2','M resolution (type-2 events); (M_{r} - M)/M',100,-2.,2.)
cstar_comp_3 = TH2D('cstar_comp_3','Reconstructed vs. Generated c* (type-3 events); c* (generated); c_{r} (reconstructed)',40,-1.,1.,40,-1.,1.) 
cstar_res_3 = TH1D('cstar_res_3','observable resolution (type-3 events); (reconstructed-generated)/generated; fraction',100,-2.,2.)
x_F_comp_3 = TH2D('x_F_comp_3','Reconstructed vs. Generated |x_{F}| (type-3 events); |x_{F}| (generated); |x_{r}| (reconstructed)',30,0.,0.6,30,0.,0.6) 
x_F_res_3 = TH1D('x_F_res_3','|x_{F}| resolution (type-3 events); (|x_{r}|-|x_{F}|)/|x_{F}|',100,-2.,2.)
M_comp_3 = TH2D('M_comp_3','Reconstructed vs. Generated M (type-3 events); M (generated) (GeV); M_{r} (reconstructed) (GeV)',40,750.,2750.,40,750.,2750.)
M_res_3 = TH1D('M_res_3','M resolution (type-3 events); (M_{r} - M)/M',100,-2.,2.)

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

#Save plots
outfile.cd()
cstar_comp_1.Write()
cstar_res_1.Write()
x_F_comp_1.Write()
x_F_res_1.Write()
M_comp_1.Write()
M_res_1.Write()
cstar_comp_2.Write()
cstar_res_2.Write()
x_F_comp_2.Write()
x_F_res_2.Write()
M_comp_2.Write()
M_res_2.Write()
cstar_comp_3.Write()
cstar_res_3.Write()
x_F_comp_3.Write()
x_F_res_3.Write()
M_comp_3.Write()
M_res_3.Write()

#make the TCanvases
cstar_comp_canv_1 = TCanvas('cstar_comp_canv_1','cstar_comp_canv_1',1100,900)
x_F_comp_canv_1 = TCanvas('x_F_comp_canv_1','x_F_comp_canv_1',1100,900)
M_comp_canv_1 = TCanvas('M_comp_canv_1','M_comp_canv_1',1100,900)
res_canv_1 = TCanvas('res_canv_1','res_canv_1',1100,900)
cstar_comp_canv_2 = TCanvas('cstar_comp_canv_2','cstar_comp_canv_2',1100,900)
x_F_comp_canv_2 = TCanvas('x_F_comp_canv_2','x_F_comp_canv_2',1100,900)
M_comp_canv_2 = TCanvas('M_comp_canv_2','M_comp_canv_2',1100,900)
res_canv_2 = TCanvas('res_canv_2','res_canv_2',1100,900)
cstar_comp_canv_3 = TCanvas('cstar_comp_canv_3','cstar_comp_canv_3',1100,900)
x_F_comp_canv_3 = TCanvas('x_F_comp_canv_3','x_F_comp_canv_3',1100,900)
M_comp_canv_3 = TCanvas('M_comp_canv_3','M_comp_canv_3',1100,900)
res_canv_3 = TCanvas('res_canv_3','res_canv_3',1100,900)

#make the legend
leg = TLegend(0.1,0.7,0.48,0.9)
leg.AddEntry(cstar_res_1,'c*','PE')
leg.AddEntry(x_F_res_1,'|x_F|','PE')
leg.AddEntry(M_res_1,'M','PE')

#plot the plots
cstar_comp_canv_1.cd()
cstar_comp_1.Draw("COLZ")
x_F_comp_canv_1.cd()
x_F_comp_1.Draw("COLZ")
M_comp_canv_1.cd()
M_comp_1.Draw("COLZ")
res_canv_1.cd()
cstar_res_1.Draw("PE")
x_F_res_1.Draw("PE SAME")
M_res_1.Draw("PE SAME")
leg.Draw("SAME")
cstar_comp_canv_2.cd()
cstar_comp_2.Draw("COLZ")
x_F_comp_canv_2.cd()
x_F_comp_2.Draw("COLZ")
M_comp_canv_2.cd()
M_comp_2.Draw("COLZ")
res_canv_2.cd()
cstar_res_2.Draw("PE")
x_F_res_2.Draw("PE SAME")
M_res_2.Draw("PE SAME")
leg.Draw("SAME")
cstar_comp_canv_3.cd()
cstar_comp_3.Draw("COLZ")
x_F_comp_canv_3.cd()
x_F_comp_3.Draw("COLZ")
M_comp_canv_3.cd()
M_comp_3.Draw("COLZ")
res_canv_3.cd()
cstar_res_3.Draw("PE")
x_F_res_3.Draw("PE SAME")
M_res_3.Draw("PE SAME")
leg.Draw("SAME")

#write the canvas
outfile.cd()
cstar_comp_canv_1.Write()
x_F_comp_canv_1.Write()
M_comp_canv_1.Write()
res_canv_1.Write()
cstar_comp_canv_2.Write()
x_F_comp_canv_2.Write()
M_comp_canv_2.Write()
res_canv_2.Write()
cstar_comp_canv_3.Write()
x_F_comp_canv_3.Write()
M_comp_canv_3.Write()
res_canv_3.Write()
outfile.Write()
outfile.Close()
