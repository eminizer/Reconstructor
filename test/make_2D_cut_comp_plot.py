from ROOT import *
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--file', metavar='F', type='string', action='store', 
                            default='/uscms_data/d3/eminizer/ttbar_13TeV/CMSSW_8_0_16/src/Analysis/Reconstructor/test/total_ttree_files/mcatnlo_semilep_TT_skim_all.root', 
                            dest='file',  help='') ## TTree file path
parser.add_option('--outname',  metavar='F', type='string', action='store', 
                              default='2D_cut_comp_plot', 
                              dest='outname',  help='') ## name for output file
(options, args) = parser.parse_args()

#set up output file
outfile = TFile(options.outname+'.root','recreate')

#set up the plot
histo = TH2D('histo','lepton 2D isolation cut; #Delta R(lepton,nearest jet); p_{T,rel}(lepton,nearest jet)',30,0.05,1.55,30,5.,80.)

#open the input file
infile = TFile(options.file)

#get the tree from the file
tree = infile.Get('tree')

#get the total weight from the file
weighthisto = infile.Get('totweight_histo')
totweight = weighthisto.GetBinContent(weighthisto.FindFixBin(0.5))

#Define cuts and draw into the histograms
preselection = '(hadt_pt>300. && hadt_M>50.)'

#singlemu_selection = '((mu_trigger==1 && muon1_ID==1 && muon1_pt>40. && muon1_pt>ele1_pt && abs(muon1_eta)<2.4) && (muon2_ID!=1 || muon2_pt<40. || abs(muon2_eta)>2.4))'
singlemu_selection = '((muon1_ID==1 && muon1_pt>40. && muon1_pt>ele1_pt && abs(muon1_eta)<2.4) && (muon2_ID!=1 || muon2_pt<40. || abs(muon2_eta)>2.4))'
el_rejection = '((ele1_ID!=1 || ele1_pt<40. || abs(ele1_eta)>2.4) && (ele2_ID!=1 || ele2_pt<40. || abs(ele2_eta)>2.4))'

#singleel_selection = '((el_trigger==1 && ele1_ID==1 && ele1_pt>40. && ele1_pt>muon1_pt && abs(ele1_eta)<2.4) && (ele2_ID!=1 || ele2_pt<40. || abs(ele2_eta)>2.4))'
singleel_selection = '((ele1_ID==1 && ele1_pt>40. && ele1_pt>muon1_pt && abs(ele1_eta)<2.4) && (ele2_ID!=1 || ele2_pt<40. || abs(ele2_eta)>2.4))'
mu_rejection = '((muon1_ID!=1 || muon1_pt<40. || abs(muon1_eta)>2.4) && (muon2_ID!=1 || muon2_pt<40. || abs(muon2_eta)>2.4))'

addl_had_cuts = '(hadt_M>140. && M_scaled>750.)'

all_mu_cuts = '('+preselection+' && '+singlemu_selection+' && '+el_rejection+' && '+addl_had_cuts+')'
all_el_cuts = '('+preselection+' && '+singleel_selection+' && '+mu_rejection+' && '+addl_had_cuts+')'

print 'all_mu_cuts = '+all_mu_cuts
print 'all_el_cuts = '+all_el_cuts

weights = '((12900.*weight)/'+str(totweight)+')'

tree.Draw("muon1_relPt:muon1_dR>>mu_histo(30,0.05,1.55,30,5.,80.)",weights+"*"+all_mu_cuts,"COLZ")
tree.Draw("ele1_relPt:ele1_dR>>el_histo(30,0.05,1.55,30,5.,80.)",weights+"*"+all_el_cuts,"COLZ")

#sum the histograms and save
mu_histo = gROOT.FindObject('mu_histo')
el_histo = gROOT.FindObject('el_histo')
histo.Add(mu_histo); histo.Add(el_histo)
outfile.cd()
histo.Write()
mu_histo.Write()
el_histo.Write()

#make the TCanvas
canv = TCanvas('canv','canv',1100,900)

#make the red box to go on the plot
box = TBox(0.05,5.,0.5,25.)
box.SetLineColor(kRed)
box.SetFillColor(kRed)

#plot the plot
canv.cd()
histo.Draw("COLZ")
box.Draw("SAME")

#write the canvas
outfile.cd()
canv.Write()
