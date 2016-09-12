from ROOT import *
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--file',  metavar='F', type='string', action='store', 
                              default='/uscms_data/d3/eminizer/ttbar_13TeV/CMSSW_8_0_16/src/Analysis/Reconstructor/test/total_ttree_files/mcatnlo_qq_semilep_TT_skim_all.root', 
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

#get the total weight from the file
weighthisto = infile.Get('totweight_histo')
totweight = weighthisto.GetBinContent(weighthisto.FindFixBin(0.5))

#Define cuts and draw into the histograms
preselection = '(hadt_pt>300. && scaled_hadt_M>50.)'

#mu_selection = '(muTrig==1 && muon1_ID==1 && muon1_pt>40. && muon1_pt>ele1_pt && abs(muon1_eta)<2.4)'
mu_selection = '(muon1_ID==1 && muon1_pt>40. && muon1_pt>ele1_pt && abs(muon1_eta)<2.4)'

#el_selection = '(elTrig==1 && ele1_ID==1 && ele1_pt>40. && ele1_pt>muon1_pt && abs(ele1_eta)<2.4)'
el_selection = '(ele1_ID==1 && ele1_pt>40. && ele1_pt>muon1_pt && abs(ele1_eta)<2.4)'

addl_lep_side_cuts = '(scaled_lept_M>140. && scaled_lept_M<250.)'

addl_had_cuts = '(scaled_hadt_M>50. && M_scaled>750. && hadt_tau21>0.1)'

all_cuts = '('+preselection+' && ('+mu_selection+' || '+el_selection+') && '+addl_lep_side_cuts+' && '+addl_had_cuts+')'

print 'all_cuts = '+all_cuts

weights = '((12900.*weight)/'+str(totweight)+')'

tree.Draw("cstar_scaled:cstar_MC>>cstar_comp_int(40,-1.,1.,40,-1.,1.)",weights+"*"+all_cuts,"COLZ")
tree.Draw("x_F_scaled:x_F_MC>>x_F_comp_int(30,0.,0.3,30,0.,0.3)",weights+"*"+all_cuts,"COLZ")
tree.Draw("M_scaled:M_MC>>M_comp_int(30,750.,1350.,30,750.,1350.)",weights+"*"+all_cuts,"COLZ")
tree.Draw("M_scaled-M_MC>>M_res_int(50,-500.,500.)",weights+"*"+all_cuts)

#get the histograms and set titles
cstar_comp = TH2D('cstar_comp','Reconstructed vs. Generated c*; c* (generated); c_{r} (reconstructed)',40,-1.,1.,40,-1.,1.) 
x_F_comp = TH2D('x_F_comp','Reconstructed vs. Generated |x_{F}|; |x_{F}| (generated); |x_{r}| (reconstructed)',30,0.,0.3,30,0.,0.3) 
M_comp = TH2D('M_comp','Reconstructed vs. Generated M; M (generated) (GeV); M_{r} (reconstructed) (GeV)',30,750.,1350.,30,750.,1350.)
M_res = TH1D('M_res','M resolution; M_{r} - M (Reconstructed - Generated) (GeV)',50,-500.,500.)
M_res.SetLineWidth(4); M_res.SetMarkerStyle(25)
cstar_comp.Add(gROOT.FindObject('cstar_comp_int'))
x_F_comp.Add(gROOT.FindObject('x_F_comp_int'))
M_comp.Add(gROOT.FindObject('M_comp_int'))
M_res.Add(gROOT.FindObject('M_res_int'))

#Save plots
outfile.cd()
cstar_comp.Write()
x_F_comp.Write()
M_comp.Write()
M_res.Write()

#make the TCanvases
cstar_canv = TCanvas('cstar_canv','cstar_canv',1100,900)
x_F_canv = TCanvas('x_F_canv','x_F_canv',1100,900)
M_canv = TCanvas('M_canv','M_canv',1100,900)
M_res_canv = TCanvas('M_res_canv','M_res_canv',1100,900)

#plot the plots
cstar_canv.cd()
cstar_comp.Draw("COLZ")
x_F_canv.cd()
x_F_comp.Draw("COLZ")
M_canv.cd()
M_comp.Draw("COLZ")
M_res_canv.cd()
M_res.Draw("HIST")

#write the canvas
outfile.cd()
cstar_canv.Write()
x_F_canv.Write()
M_canv.Write()
M_res_canv.Write()