from ROOT import *
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--file',  metavar='F', type='string', action='store', 
                              default='../total_ttree_files/mcatnlo_semilep_TT_all.root', 
                              dest='file',  help='') ## TTree file path
parser.add_option('--outname',  metavar='F', type='string', action='store', 
                              default='reconstruction_comparison_plots_semilep_sr', 
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
signalregion = fullselection+' && hadt_isttagged==1 '

#all_cuts = '('+fullselection+')'
all_cuts = '('+signalregion+')'

print 'all_cuts = '+all_cuts

weights = '(12900.*weight)'

tree.Draw("cstar:cstar_MC>>cstar_comp_int(40,-1.,1.,40,-1.,1.)",weights+"*"+all_cuts,"COLZ")
tree.Draw("x_F:x_F_MC>>x_F_comp_int(30,0.,0.6,30,0.,0.6)",weights+"*"+all_cuts,"COLZ")
tree.Draw("M:M_MC>>M_comp_int(40,750.,2750.,40,750.,2750.)",weights+"*"+all_cuts,"COLZ")
tree.Draw("(M-M_MC)/M_MC>>M_res_int(30,-1.5,1.5)",weights+"*"+all_cuts)

#get the histograms and set titles
cstar_comp = TH2D('cstar_comp','Reconstructed vs. Generated c*; c* (generated); c_{r} (reconstructed)',40,-1.,1.,40,-1.,1.) 
x_F_comp = TH2D('x_F_comp','Reconstructed vs. Generated |x_{F}|; |x_{F}| (generated); |x_{r}| (reconstructed)',30,0.,0.6,30,0.,0.6) 
M_comp = TH2D('M_comp','Reconstructed vs. Generated M; M (generated) (GeV); M_{r} (reconstructed) (GeV)',40,750.,2750.,40,750.,2750.)
M_res = TH1D('M_res','M resolution; M_{r} - M (Reconstructed - Generated) (GeV)',30,-1.5,1.5)
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