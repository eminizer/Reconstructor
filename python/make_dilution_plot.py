#this script makes the dilution plot, showing for all skimmed qqbar -> ttbar events the 
#ratio of the difference to the sum of correct/incorrect assignments of the initial parton directions

#imports
from ROOT import gROOT, TFile, TChain, TH1F, TLorentzVector, kBlue, TCanvas, TLine, kBlack
import CMS_lumi, tdrstyle
from math import sqrt
import os
from array import array
from time import *
from datetime import timedelta
from branch import Branch
from angleReconstructor import getObservables

#constants
SQRT_S=13000.0

#some setup stuff
gROOT.SetBatch()
tdrstyle.setTDRStyle()

#list of input files (the nominal ttbar ntuple files)
ifns = [
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_0.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_1.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_10.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_11.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_12.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_13.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_14.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_15.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_16.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_17.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_18.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_19.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_2.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_20.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_21.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_22.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_3.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_4.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_5.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_6.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_7.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_8.root',
'root://131.225.204.161:1094//store/user/eminizer/B2GTTrees_Dec17_2018/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/aggregated_powheg_TT_9.root',
]

#open a garbage file
gfp = TFile('dilution_plot_garbage.root','recreate')

#chain up the files
chain = TChain('B2GTree')
for ifn in ifns :
	chain.Add(os.path.join(ifn,'B2GTTreeMaker','B2GTree'))

#make the dict of branches that we'll need
branches = {}
branches['part_1_fac'] = Branch('MC_part1_factor')
branches['part_1_ID'] = Branch('MC_part1_ID')
branches['part_2_fac'] = Branch('MC_part2_factor')
branches['part_2_ID'] = Branch('MC_part2_ID')
branches['t_pt'] = Branch('MC_t_pt')
branches['t_eta'] = Branch('MC_t_eta')
branches['t_phi'] = Branch('MC_t_phi')
branches['t_E'] = Branch('MC_t_E')
branches['tbar_pt'] = Branch('MC_tbar_pt')
branches['tbar_eta'] = Branch('MC_tbar_eta')
branches['tbar_phi'] = Branch('MC_tbar_phi')
branches['tbar_E'] = Branch('MC_tbar_E')
for b in branches.values() :
	b.initialize(chain,None)

#set up the plot file
pfp = TFile('dilution_plots.root','recreate')

#set up the histograms to fill
#bins = array('d',[0.000,0.025,0.050,0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,
#				  0.525,0.550,0.575,0.600,1.00])
bins = array('d',[0.000,0.025,0.050,0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,
				  0.525,0.550,0.575,0.600,0.650,0.700,0.800,1.00])
correct = TH1F('correct','correct assignments as |x_{F}|; |x_{F}|; # correct assignments',28,bins)
incorrect = TH1F('incorrect','incorrect assignments as |x_{F}|; |x_{F}|; # incorrect assignments',28,bins)
dilution = TH1F('dilution',';|x_{F}|;Dilution factor',28,bins)

setupDoneTime = time()

#loop over all the events and fill the (in)correct histograms
nTotalEvents = chain.GetEntries()
count=0
for ievent in range(nTotalEvents)  :
	count+=1
	#if count>100: #DEBUG
	#	break #DEBUG
	if count % 100000 == 0 or count == 1:
		timeSinceSetup = timedelta(seconds=time()-setupDoneTime)
		percentDone = float(count) / float(nTotalEvents) * 100.0
		timeLeft = timedelta(seconds=(100.-percentDone)*timeSinceSetup.total_seconds()/percentDone)
		print ( 'Count at %d of %d, (%.4f%% complete; elapsed = %02d:%02d:%02d, remaining = %02d:%02d:%02d)'
			   %(count,nTotalEvents,percentDone,timeSinceSetup.seconds/3600,(timeSinceSetup.seconds%3600)/60,(timeSinceSetup.seconds%60),timeLeft.seconds/3600,(timeLeft.seconds%3600)/60,(timeLeft.seconds%60))  )
	chain.GetEntry(ievent)
	#get the parton IDs and check for the correct event type
	p1id = branches['part_1_ID'].getReadValue()
	p2id = branches['part_2_ID'].getReadValue()
	#print('p1ID = {}, p2ID = {}'.format(p1id,p2id)) #DEBUG
	if not p1id+p2id==0 :
		continue
	#print('-------------------------------------------\np1ID = {}, p2ID = {}'.format(p1id,p2id)) #DEBUG
	#make the fourvector of the top/antitop system
	tvec = TLorentzVector()
	tvec.SetPtEtaPhiE(branches['t_pt'].getReadValue(),branches['t_eta'].getReadValue(),
						branches['t_phi'].getReadValue(),branches['t_E'].getReadValue())
	tbarvec = TLorentzVector()
	tbarvec.SetPtEtaPhiE(branches['tbar_pt'].getReadValue(),branches['tbar_eta'].getReadValue(),
							branches['tbar_phi'].getReadValue(),branches['tbar_E'].getReadValue())
	Q = tvec+tbarvec
	#get the feynman x value
	xF = 2*Q.Pz()/SQRT_S
	#print('x_F = {}'.format(xF)) #DEBUG
	#check if the assignment would be correct, and then fill the corresponding histogram
	qforward = branches['part_1_fac'].getReadValue()>0 if p1id>0 else branches['part_2_fac'].getReadValue()>0
	Qforward = Q.Eta()>0.
	#print('p1 fac = {}, p2 fac = {}, Q eta = {}, qforward = {}, Qforward = {}'.format(branches['part_1_fac'].getReadValue(), #DEBUG
	#																				  branches['part_2_fac'].getReadValue(),Q.Eta(),qforward,Qforward)) #DEBUG
	if qforward==Qforward :
		correct.Fill(abs(xF))
	else :
		incorrect.Fill(abs(xF))

print 'Done.'

#build the dilution plot histogram
for i in range(1,dilution.GetNbinsX()+1) :
	corcont = correct.GetBinContent(i)
	corerr = correct.GetBinError(i)
	icorcont = incorrect.GetBinContent(i)
	icorerr = incorrect.GetBinError(i)
	print('bin {}, corcont={}, corerr={}, icorcont={}, icorerr={}'.format(i,corcont,corerr,icorcont,icorerr)) #DEBUG
	dcont = (corcont-icorcont)/(corcont+icorcont)
	derr = (2./((corcont+icorcont)**2))*sqrt((icorcont*corerr)**2+(corcont*icorerr)**2)
	dilution.SetBinContent(i,dcont)
	dilution.SetBinError(i,derr)
dilution.SetStats(0)
dilution.SetLineColor(kBlue)
dilution.SetLineWidth(2)
dilution.SetMarkerColor(kBlue)

#plot on canvas
pfp.cd()
dilution_plot = TCanvas('dilution_plot','dilution_plot',1100,900)
#dilution_plot.SetTopMargin(0.11)
dilution.Draw('PE')
dilution.GetYaxis().SetRangeUser(0.,1.20)
#plot the line at one on the canvas
oneline=TLine(0.0,1.0,1.0,1.0)
oneline.SetLineWidth(2)
oneline.SetLineColor(kBlack)
oneline.SetLineStyle(7)
oneline.Draw('SAME')
#plot the CMS_Lumi lines on the canvas
iPeriod = 0 #free form since it's only simulation
iPos = 33#1 #iPos = 10*(alignment) + position (1/2/3 = left/center/right)
#CMS_lumi.cmsText = "POWHEG"
CMS_lumi.cmsText = "CMS"
#CMS_lumi.cmsTextSize = 0.85
CMS_lumi.cmsTextSize = 1.1
CMS_lumi.writeExtraText = 1
#CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.extraText = "Simulation"
#CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.CMS_lumi(dilution_plot, iPeriod, iPos)
dilution_plot.Update()
dilution_plot.Write()

#handling remaining file I/O
pfp.Write()
pfp.Close()
gfp.Close()
print 'Deleting garbage file'
os.system('rm -f dilution_plot_garbage.root')
print 'Done!'


