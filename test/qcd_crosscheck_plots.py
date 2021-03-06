import os
from ROOT import *
import CMS_lumi, tdrstyle
from glob import glob
from math import *
from datetime import date
from optparse import OptionParser
import multiprocessing
import copy

def_weights  = '(((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_ttag_eff_merged*sf_ttag_eff_semimerged*sf_ttag_eff_notmerged*sf_btag_eff_heavy*sf_btag_eff_light*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'
qcda_weights = '(((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH))*weight*sf_pileup*sf_ttag_eff_merged*sf_ttag_eff_semimerged*sf_ttag_eff_notmerged*sf_btag_eff_heavy*sf_btag_eff_light*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'
qcdb_weights = def_weights
qcdc_weights = '(((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH))*weight*sf_pileup*sf_ttag_eff_merged*sf_ttag_eff_semimerged*sf_ttag_eff_notmerged*sf_btag_eff_heavy*sf_btag_eff_light*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'

gROOT.SetBatch()

#treeskimmer helper function
def treeskimmer(thischain,thisname,thiscuts,lock,qcdsb=False) :
	print 'Skimming chains for '+thisname+' samples'
	if qcdsb :
		skimmedtreefile = TFile(outname+'_'+thisname.replace(' ','_').replace('#','')+'_qcd_sb_skimmed_tree.root','recreate')
	else :	
		skimmedtreefile = TFile(outname+'_'+thisname.replace(' ','_').replace('#','')+'_sr_skimmed_tree.root','recreate')
	#print 'thischain = %s, thiscuts = %s'%(thischain,thiscuts) #DEBUG
	newtree = thischain.CopyTree(thiscuts)
	lock.acquire()
	skimmedtreefile.cd()
	newtree.Write()
	skimmedtreefile.Write()
	skimmedtreefile.Close()
	lock.release()
	print 'Done skimming chains for '+thisname+' samples'

tdrstyle.setTDRStyle()
iPeriod = 4 #13TeV iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"

# COMMAND LINE OPTIONS
parser = OptionParser()
#Run options
parser.add_option('--leptype', 	  type='string', action='store', default='all_leptons', dest='leptype',	   	  
	help='Use SingleMu or SingleEl data? ("muons" of "electrons")')
parser.add_option('--ttbar_generator', 	  type='string', action='store', default='powheg', dest='ttbar_generator',	   	  
	help='Use "powheg" or "mcatnlo" ttbar MC?')
parser.add_option('--mode', 	  type='string', action='store', default='', dest='mode')
parser.add_option('--n_threads', 	  type='int', action='store', default=7, dest='n_threads')
parser.add_option('--outtag', 	  type='string', action='store', default='', dest='outtag')
parser.add_option('--skimcut', 	  type='string', action='store', 
						default='def', 
						dest='skimcut')
parser.add_option('--use_preskimmed_files', type='string', action='store', default='no', dest='usepreskim')
parser.add_option('--use_top_pt_reweighting', type='string', action='store', default='yes', dest='usetopptrw')
parser.add_option('--MC_weights', 	  type='string', action='store', 
						default='def', 
						dest='MC_weights')
(options, args) = parser.parse_args()

leptype=options.leptype.lower()
generator = options.ttbar_generator.lower()
mode = options.mode.lower()
if mode in ['wjets_cr','wjcr','wjetscr'] :
	mode = 'wjetscr'
	if options.outtag=='' :
		options.outtag='wjets_cr'
	else :
		options.outtag=options.outtag+'_wjets_cr'
	options.skimcut='wjets_cr_selection==1'
elif mode in ['qcd_a_sr','qcdasr','qcd_a_sr'] :
	mode = 'qcd_a_sr'
	if options.outtag=='' :
		options.outtag='qcd_a_sr_sb'
	else :
		options.outtag=options.outtag+'_qcd_a_sr_sb'
	options.skimcut='qcd_A_SR_selection==1'
	options.MC_weights = qcda_weights
elif mode in ['qcd_b_sr','qcdbsr','qcd_b_sr'] :
	mode = 'qcd_b_sr'
	if options.outtag=='' :
		options.outtag='qcd_b_sr_sb'
	else :
		options.outtag=options.outtag+'_qcd_b_sr_sb'
	options.skimcut='qcd_B_SR_selection==1'
elif mode in ['qcd_c_sr','qcdcsr','qcd_c_sr'] :
	mode = 'qcd_c_sr'
	if options.outtag=='' :
		options.outtag='qcd_c_sr_sb'
	else :
		options.outtag=options.outtag+'_qcd_c_sr_sb'
	options.skimcut='qcd_C_SR_selection==1'
	options.MC_weights = qcdc_weights
elif mode in ['qcd_a_cr','qcdacr','qcd_a_cr'] :
	mode = 'qcd_a_cr'
	if options.outtag=='' :
		options.outtag='qcd_a_cr_sb'
	else :
		options.outtag=options.outtag+'_qcd_a_cr_sb'
	options.skimcut='qcd_A_CR_selection==1'
	options.MC_weights = qcda_weights
elif mode in ['qcd_b_cr','qcdbcr','qcd_b_cr'] :
	mode = 'qcd_b_cr'
	if options.outtag=='' :
		options.outtag='qcd_b_cr_sb'
	else :
		options.outtag=options.outtag+'_qcd_b_cr_sb'
	options.skimcut='qcd_B_CR_selection==1'
elif mode in ['qcd_c_cr','qcdccr','qcd_c_cr'] :
	mode = 'qcd_c_cr'
	if options.outtag=='' :
		options.outtag='qcd_c_cr_sb'
	else :
		options.outtag=options.outtag+'_qcd_c_cr_sb'
	options.skimcut='qcd_C_CR_selection==1'
	options.MC_weights = qcdc_weights

estimate_qcd = (options.skimcut=='def' and options.MC_weights=='def') or mode in ['wjets_cr','wjcr','wjetscr']

if options.skimcut=='def' :
	options.skimcut='fullselection==1'
	#options.skimcut='fullselection==1 && ((ak41_csvv2>0.8484 && (ak42_csvv2>0.8484 || ak43_csvv2>0.8484 || ak44_csvv2>0.8484)) || (ak42_csvv2>0.8484 && (ak43_csvv2>0.8484 || ak44_csvv2>0.8484 || ak41_csvv2>0.8484)) || (ak43_csvv2>0.8484 && (ak44_csvv2>0.8484 || ak41_csvv2>0.8484 || ak42_csvv2>0.8484)) || (ak44_csvv2>0.8484 && (ak41_csvv2>0.8484 || ak42_csvv2>0.8484 || ak43_csvv2>0.8484)))'
	#options.skimcut='fullselection==1 && tt_pt<100.'
	#options.skimcut='fullselection==1 && lep_Iso!=0.'
	#options.skimcut='fullselection==1 && ak41_pt>55. && ak42_pt>45.'
if options.MC_weights=='def' :
	options.MC_weights=def_weights

if options.usetopptrw=='no' :
	if options.outtag!='' :
		options.outtag+='_'
	options.outtag+='no_top_pt_rw'

lepstring = ''; shortlepstring = ''
if leptype=='muons' :
	lepstring = 'muon'; shortlepstring = '#mu'
elif leptype=='electrons' :
	lepstring = 'electron'; shortlepstring = 'e'
elif leptype=='all_leptons' :
	lepstring = 'lepton'; shortlepstring = 'lepton'

samplenames = []
shortnames = []
#WJets
samplenames.append('WJets_HT-70to100'); 		 shortnames.append('WJets')
samplenames.append('WJets_HT-100to200'); 		 shortnames.append('WJets')
samplenames.append('WJets_HT-200to400'); 		 shortnames.append('WJets')
samplenames.append('WJets_HT-400to600'); 		 shortnames.append('WJets')
samplenames.append('WJets_HT-600to800'); 		 shortnames.append('WJets')
samplenames.append('WJets_HT-800to1200'); 		 shortnames.append('WJets')
samplenames.append('WJets_HT-1200to2500'); 		 shortnames.append('WJets')
samplenames.append('WJets_HT-2500toInf'); 		 shortnames.append('WJets')
#DYJets
samplenames.append('DYJets_M-50_HT-70to100'); 	 shortnames.append('DYJets')
samplenames.append('DYJets_M-50_HT-100to200'); 	 shortnames.append('DYJets')
samplenames.append('DYJets_M-50_HT-200to400'); 	 shortnames.append('DYJets')
samplenames.append('DYJets_M-50_HT-400to600'); 	 shortnames.append('DYJets')
samplenames.append('DYJets_M-50_HT-600to800'); 	 shortnames.append('DYJets')
samplenames.append('DYJets_M-50_HT-800to1200');  shortnames.append('DYJets')
samplenames.append('DYJets_M-50_HT-1200to2500'); shortnames.append('DYJets')
samplenames.append('DYJets_M-50_HT-2500toInf');  shortnames.append('DYJets')
#Single top
samplenames.append('ST_s-c'); 					 shortnames.append('Single top')
samplenames.append('ST_t-c_top'); 				 shortnames.append('Single top')
samplenames.append('ST_t-c_antitop'); 			 shortnames.append('Single top')
samplenames.append('ST_tW-c_top'); 				 shortnames.append('Single top')
samplenames.append('ST_tW-c_antitop'); 			 shortnames.append('Single top')
#POWHEG TT
if generator == 'powheg' :
	samplenames.append('powheg_TT'); 			 shortnames.append('t#bar{t} other')
	samplenames.append('powheg_TT'); 			 shortnames.append('t#bar{t} semilep')
#MCATNLO TT
elif generator == 'mcatnlo' :
	samplenames.append('mcatnlo_TT'); 			 shortnames.append('t#bar{t} other')
	samplenames.append('mcatnlo_TT'); 			 shortnames.append('t#bar{t} semilep')
else :
	print 'generator type '+generator+' not recognized! Exiting.'
	exit()

#data
data_samplenames = []
if leptype=='muons' or leptype=='all_leptons' :
	data_samplenames.append('SingleMu_Run2016Bv2')
	data_samplenames.append('SingleMu_Run2016C')
	data_samplenames.append('SingleMu_Run2016D')
	data_samplenames.append('SingleMu_Run2016E')
	data_samplenames.append('SingleMu_Run2016F')
	data_samplenames.append('SingleMu_Run2016G')
	data_samplenames.append('SingleMu_Run2016Hv2')
	data_samplenames.append('SingleMu_Run2016Hv3')
if leptype=='electrons' or leptype=='all_leptons' :
	data_samplenames.append('SingleEl_Run2016Bv2')
	data_samplenames.append('SingleEl_Run2016C')
	data_samplenames.append('SingleEl_Run2016D')
	data_samplenames.append('SingleEl_Run2016E')
	data_samplenames.append('SingleEl_Run2016F')
	data_samplenames.append('SingleEl_Run2016G')
	data_samplenames.append('SingleEl_Run2016Hv2')
	data_samplenames.append('SingleEl_Run2016Hv3')

#QCD MC
qcd_mc_samplenames = []
qcd_mc_samplenames.append('QCD_HT-100to200')
qcd_mc_samplenames.append('QCD_HT-200to300')
qcd_mc_samplenames.append('QCD_HT-300to500')
qcd_mc_samplenames.append('QCD_HT-500to700')
qcd_mc_samplenames.append('QCD_HT-700to1000')
qcd_mc_samplenames.append('QCD_HT-1000to1500')
qcd_mc_samplenames.append('QCD_HT-1500to2000')
qcd_mc_samplenames.append('QCD_HT-2000toInf')

#set up output file
outname = 'qcd_crosscheck_plots_'+leptype+'_'+generator+'_'+str(date.today())
if options.outtag!='' :
	outname+='_'+options.outtag
garbageFileName = outname+'_garbage.root'
garbageFile = TFile(garbageFileName,'recreate') 

#Chain up the files
MC_chains = []
data_chain = TChain('tree')
qcd_mc_chain = TChain('tree')
MC_chains_qcd_sb = []
data_chain_qcd_sb = TChain('tree')
qcd_mc_chain_qcd_sb = TChain('tree')
shortnames_done = []
for i in range(len(samplenames)) :
	index = 0
	if shortnames[i] not in shortnames_done :
		MC_chains.append(TChain('tree'))
		MC_chains_qcd_sb.append(TChain('tree'))
		index = len(MC_chains)-1
		shortnames_done.append(shortnames[i])
	else :
		index = shortnames_done.index(shortnames[i])
	if options.usepreskim == 'no' :
		filenamelist = glob('../'+samplenames[i]+'/aggregated_'+samplenames[i]+'_*.root')
	elif options.usepreskim == 'yes' :
		filenamelist = glob('../total_ttree_files/'+samplenames[i]+'_skim_all.root')
	for filename in filenamelist :
		if filename.find('JES')==-1 and filename.find('JER')==-1 :
		#if filename.find('JES_up')!=-1 :
			MC_chains[index].Add(filename)
			MC_chains_qcd_sb[index].Add(filename)
			print 'adding %s to chain for shortname %s'%(filename,shortnames[i]) #DEBUG
#print 'MC_chains = %s'%(MC_chains) #DEBUG
for data_samplename in data_samplenames :
	if options.usepreskim == 'no' :
		filenamelist = glob('../'+data_samplename+'/aggregated_'+data_samplename+'_*.root')
	elif options.usepreskim == 'yes' :
		filenamelist = glob('../total_ttree_files/'+data_samplename+'_skim_all.root')
	for filename in filenamelist :
		data_chain.Add(filename)
		data_chain_qcd_sb.Add(filename)
		print 'adding %s to data chain'%(filename) #DEBUG
for qcd_mc_samplename in qcd_mc_samplenames :
	if options.usepreskim == 'no' :
		filenamelist = glob('../'+qcd_mc_samplename+'/aggregated_'+qcd_mc_samplename+'_*.root')
	elif options.usepreskim == 'yes' :
		filenamelist = glob('../total_ttree_files/'+qcd_mc_samplename+'_skim_all.root')
	for filename in filenamelist :
		qcd_mc_chain.Add(filename)
		qcd_mc_chain_qcd_sb.Add(filename)
		print 'adding %s to QCD MC chain'%(filename) #DEBUG

#skim the chains into reduced trees and save the reduced trees in the garbage file
manager = multiprocessing.Manager()
lock = multiprocessing.Lock()
com_cuts = '('+options.skimcut+')'
if leptype=='muons' :
	com_cuts+=' && lepflavor==1'
elif leptype=='electrons' :
	com_cuts+=' && lepflavor==2'
print 'Skim cuts: %s'%(com_cuts)
procs = []
for i in range(len(MC_chains)) :
	if len(procs)>=options.n_threads :
		for p in procs :
			p.join()
		procs = []
	realcuts = '('+com_cuts+')'
	if shortnames_done[i].find('semilep')!=-1 :
		realcuts = '(('+com_cuts+') && eventType<2)'
	elif shortnames_done[i].find('other')!=-1 :
		realcuts = '(('+com_cuts+') && (eventType==2 || eventType==3))'
	p = multiprocessing.Process(target=treeskimmer, args=(MC_chains[i],shortnames_done[i],realcuts,lock))
	p.start()
	procs.append(p)
newp = multiprocessing.Process(target=treeskimmer, args=(data_chain,'DATA',com_cuts,lock))
newp.start()
procs.append(newp)
newp2 = multiprocessing.Process(target=treeskimmer, args=(qcd_mc_chain,'QCD MC',com_cuts,lock))
newp2.start()
procs.append(newp2)
for p in procs :
	p.join()

#skim sideband chains and get the QCD transfer factors per topology
QCD_transfer_factors = {'t1':0.,'t2':0.,'t3':0.}
QCD_transfer_factors_up = {'t1':0.,'t2':0.,'t3':0.}
QCD_transfer_factors_MC_driven = {'t1':0.,'t2':0.,'t3':0.}
print 'SKIMMING CHAINS FOR QCD ESTIMATE'
if mode in ['wjets_cr','wjcr','wjetscr'] :
	com_cuts = '(qcd_A_CR_selection==1 || qcd_B_CR_selection==1 || qcd_C_CR_selection==1)'
	#com_cuts+=' && ((ak41_csvv2>0.8484 && (ak42_csvv2>0.8484 || ak43_csvv2>0.8484 || ak44_csvv2>0.8484)) || (ak42_csvv2>0.8484 && (ak43_csvv2>0.8484 || ak44_csvv2>0.8484 || ak41_csvv2>0.8484)) || (ak43_csvv2>0.8484 && (ak44_csvv2>0.8484 || ak41_csvv2>0.8484 || ak42_csvv2>0.8484)) || (ak44_csvv2>0.8484 && (ak41_csvv2>0.8484 || ak42_csvv2>0.8484 || ak43_csvv2>0.8484)))'
else :
	com_cuts = '(qcd_A_SR_selection==1 || qcd_B_SR_selection==1 || qcd_C_SR_selection==1)'
	#com_cuts+=' && ((ak41_csvv2>0.8484 && (ak42_csvv2>0.8484 || ak43_csvv2>0.8484 || ak44_csvv2>0.8484)) || (ak42_csvv2>0.8484 && (ak43_csvv2>0.8484 || ak44_csvv2>0.8484 || ak41_csvv2>0.8484)) || (ak43_csvv2>0.8484 && (ak44_csvv2>0.8484 || ak41_csvv2>0.8484 || ak42_csvv2>0.8484)) || (ak44_csvv2>0.8484 && (ak41_csvv2>0.8484 || ak42_csvv2>0.8484 || ak43_csvv2>0.8484)))'
if leptype=='muons' :
	com_cuts+=' && lepflavor==1'
elif leptype=='electrons' :
	com_cuts+=' && lepflavor==2'
print 'Skim cuts: %s'%(com_cuts)
procs = []
for i in range(len(MC_chains)) :
	if len(procs)>=options.n_threads :
		for p in procs :
			p.join()
		procs = []
	realcuts = '('+com_cuts+')'
	if shortnames_done[i].find('semilep')!=-1 :
		realcuts = '(('+com_cuts+') && eventType<2)'
	elif shortnames_done[i].find('other')!=-1 :
		realcuts = '(('+com_cuts+') && (eventType==2 || eventType==3))'
	p = multiprocessing.Process(target=treeskimmer, args=(MC_chains_qcd_sb[i],shortnames_done[i],realcuts,lock,True))
	p.start()
	procs.append(p)
newp = multiprocessing.Process(target=treeskimmer, args=(data_chain_qcd_sb,'DATA',com_cuts,lock,True))
newp.start()
procs.append(newp)
newp2 = multiprocessing.Process(target=treeskimmer, args=(qcd_mc_chain_qcd_sb,'QCD_MC',com_cuts,lock,True))
newp2.start()
procs.append(newp2)
for p in procs :
	p.join()
print 'DONE SKIMMING CHAINS FOR QCD ESTIMATE'
#draw chains into histograms
qcd_a_sb_t1_histo = TH1D('qcd_a_sb_t1_histo','',20,-1.,1.); qcd_a_sb_t1_histo_MC_driven = TH1D('qcd_a_sb_t1_histo_MC_driven','',20,-1.,1.)
qcd_b_sb_t1_histo = TH1D('qcd_b_sb_t1_histo','',20,-1.,1.); qcd_b_sb_t1_histo_MC_driven = TH1D('qcd_b_sb_t1_histo_MC_driven','',20,-1.,1.)
qcd_a_sb_t2_histo = TH1D('qcd_a_sb_t2_histo','',20,-1.,1.); qcd_a_sb_t2_histo_MC_driven = TH1D('qcd_a_sb_t2_histo_MC_driven','',20,-1.,1.)
qcd_b_sb_t2_histo = TH1D('qcd_b_sb_t2_histo','',20,-1.,1.); qcd_b_sb_t2_histo_MC_driven = TH1D('qcd_b_sb_t2_histo_MC_driven','',20,-1.,1.)
qcd_a_sb_t3_histo = TH1D('qcd_a_sb_t3_histo','',20,-1.,1.); qcd_a_sb_t3_histo_MC_driven = TH1D('qcd_a_sb_t3_histo_MC_driven','',20,-1.,1.)
qcd_b_sb_t3_histo = TH1D('qcd_b_sb_t3_histo','',20,-1.,1.); qcd_b_sb_t3_histo_MC_driven = TH1D('qcd_b_sb_t3_histo_MC_driven','',20,-1.,1.)
qcd_a_cut = 'qcd_A_SR_selection==1'# && ((ak41_csvv2>0.8484 && (ak42_csvv2>0.8484 || ak43_csvv2>0.8484 || ak44_csvv2>0.8484)) || (ak42_csvv2>0.8484 && (ak43_csvv2>0.8484 || ak44_csvv2>0.8484 || ak41_csvv2>0.8484)) || (ak43_csvv2>0.8484 && (ak44_csvv2>0.8484 || ak41_csvv2>0.8484 || ak42_csvv2>0.8484)) || (ak44_csvv2>0.8484 && (ak41_csvv2>0.8484 || ak42_csvv2>0.8484 || ak43_csvv2>0.8484)))'
qcd_b_cut = 'qcd_B_SR_selection==1'# && ((ak41_csvv2>0.8484 && (ak42_csvv2>0.8484 || ak43_csvv2>0.8484 || ak44_csvv2>0.8484)) || (ak42_csvv2>0.8484 && (ak43_csvv2>0.8484 || ak44_csvv2>0.8484 || ak41_csvv2>0.8484)) || (ak43_csvv2>0.8484 && (ak44_csvv2>0.8484 || ak41_csvv2>0.8484 || ak42_csvv2>0.8484)) || (ak44_csvv2>0.8484 && (ak41_csvv2>0.8484 || ak42_csvv2>0.8484 || ak43_csvv2>0.8484)))'
qcd_c_cut = 'qcd_C_SR_selection==1'# && ((ak41_csvv2>0.8484 && (ak42_csvv2>0.8484 || ak43_csvv2>0.8484 || ak44_csvv2>0.8484)) || (ak42_csvv2>0.8484 && (ak43_csvv2>0.8484 || ak44_csvv2>0.8484 || ak41_csvv2>0.8484)) || (ak43_csvv2>0.8484 && (ak44_csvv2>0.8484 || ak41_csvv2>0.8484 || ak42_csvv2>0.8484)) || (ak44_csvv2>0.8484 && (ak41_csvv2>0.8484 || ak42_csvv2>0.8484 || ak43_csvv2>0.8484)))'
if mode in ['wjets_cr','wjcr','wjetscr'] :
	qcd_a_cut = 'qcd_A_CR_selection==1'
	qcd_b_cut = 'qcd_B_CR_selection==1'
	qcd_c_cut = 'qcd_C_CR_selection==1'
if leptype=='muons' :
	qcd_a_cut+=' && lepflavor==1'
	qcd_b_cut+=' && lepflavor==1'
elif leptype=='electrons' :
	qcd_a_cut+=' && lepflavor==2'
	qcd_b_cut+=' && lepflavor==2'
data_chain_qcd_sb.Draw('cstar>>data_qcd_a_sb_t1(20,-1.,1.)','('+qcd_a_cut+' && eventTopology==1)')
qcd_a_sb_t1_histo.Add(gROOT.FindObject('data_qcd_a_sb_t1'))
data_chain_qcd_sb.Draw('cstar>>data_qcd_b_sb_t1(20,-1.,1.)','('+qcd_b_cut+' && eventTopology==1)')
qcd_b_sb_t1_histo.Add(gROOT.FindObject('data_qcd_b_sb_t1'))
data_chain_qcd_sb.Draw('cstar>>data_qcd_a_sb_t2(20,-1.,1.)','('+qcd_a_cut+' && eventTopology==2)')
qcd_a_sb_t2_histo.Add(gROOT.FindObject('data_qcd_a_sb_t2'))
data_chain_qcd_sb.Draw('cstar>>data_qcd_b_sb_t2(20,-1.,1.)','('+qcd_b_cut+' && eventTopology==2)')
qcd_b_sb_t2_histo.Add(gROOT.FindObject('data_qcd_b_sb_t2'))
data_chain_qcd_sb.Draw('cstar>>data_qcd_a_sb_t3(20,-1.,1.)','('+qcd_a_cut+' && eventTopology==3)')
qcd_a_sb_t3_histo.Add(gROOT.FindObject('data_qcd_a_sb_t3'))
data_chain_qcd_sb.Draw('cstar>>data_qcd_b_sb_t3(20,-1.,1.)','('+qcd_b_cut+' && eventTopology==3)')
qcd_b_sb_t3_histo.Add(gROOT.FindObject('data_qcd_b_sb_t3'))
print 'qcd sb numbers after adding data: type-1 A=%.2f, B=%.2f; type-2 A=%.2f, B=%.2f; type-3 A=%.2f, B=%.2f'%(qcd_a_sb_t1_histo.Integral(),qcd_b_sb_t1_histo.Integral(),qcd_a_sb_t2_histo.Integral(),qcd_b_sb_t2_histo.Integral(),qcd_a_sb_t3_histo.Integral(),qcd_b_sb_t3_histo.Integral()) #DEBUG
for i in range(len(MC_chains)) :
	mc_qcd_a_cut = qcd_a_cut
	mc_qcd_b_cut = qcd_b_cut
	if shortnames_done[i].find('semilep')!=-1 :
		mc_qcd_a_cut+=' && eventType<2'
		mc_qcd_b_cut+=' && eventType<2'
	elif shortnames_done[i].find('other')!=-1 :
		mc_qcd_a_cut+=' && (eventType==2 || eventType==3)'
		mc_qcd_b_cut+=' && (eventType==2 || eventType==3)'
	mc_qcda_weights = qcda_weights
	mc_qcdb_weights = qcdb_weights
	if shortnames_done[i].find('t#bar{t}')!=-1 and options.usetopptrw=='yes' :
		mc_qcda_weights+='*sf_top_pt_rw_v2'
		mc_qcdb_weights+='*sf_top_pt_rw_v2'
	MC_chains_qcd_sb[i].Draw('cstar>>mc_'+str(i)+'_a_sb_t1(20,-1.,1.)','('+mc_qcd_a_cut+' && eventTopology==1)*('+qcda_weights+')*(-1.)')
	qcd_a_sb_t1_histo.Add(gROOT.FindObject('mc_'+str(i)+'_a_sb_t1'))
	MC_chains_qcd_sb[i].Draw('cstar>>mc_'+str(i)+'_b_sb_t1(20,-1.,1.)','('+mc_qcd_b_cut+' && eventTopology==1)*('+qcdb_weights+')*(-1.)')
	qcd_b_sb_t1_histo.Add(gROOT.FindObject('mc_'+str(i)+'_b_sb_t1'))
	MC_chains_qcd_sb[i].Draw('cstar>>mc_'+str(i)+'_a_sb_t2(20,-1.,1.)','('+mc_qcd_a_cut+' && eventTopology==2)*('+qcda_weights+')*(-1.)')
	qcd_a_sb_t2_histo.Add(gROOT.FindObject('mc_'+str(i)+'_a_sb_t2'))
	MC_chains_qcd_sb[i].Draw('cstar>>mc_'+str(i)+'_b_sb_t2(20,-1.,1.)','('+mc_qcd_b_cut+' && eventTopology==2)*('+qcdb_weights+')*(-1.)')
	qcd_b_sb_t2_histo.Add(gROOT.FindObject('mc_'+str(i)+'_b_sb_t2'))
	MC_chains_qcd_sb[i].Draw('cstar>>mc_'+str(i)+'_a_sb_t3(20,-1.,1.)','('+mc_qcd_a_cut+' && eventTopology==3)*('+qcda_weights+')*(-1.)')
	qcd_a_sb_t3_histo.Add(gROOT.FindObject('mc_'+str(i)+'_a_sb_t3'))
	MC_chains_qcd_sb[i].Draw('cstar>>mc_'+str(i)+'_b_sb_t3(20,-1.,1.)','('+mc_qcd_b_cut+' && eventTopology==3)*('+qcdb_weights+')*(-1.)')
	qcd_b_sb_t3_histo.Add(gROOT.FindObject('mc_'+str(i)+'_b_sb_t3'))
	print 'qcd sb numbers after subtracting %s MC: type-1 A=%.2f, B=%.2f, type-2 A=%.2f, B=%.2f, type-3 A=%.2f, B=%.2f'%(shortnames_done[i],qcd_a_sb_t1_histo.Integral(),qcd_b_sb_t1_histo.Integral(),qcd_a_sb_t2_histo.Integral(),qcd_b_sb_t2_histo.Integral(),qcd_a_sb_t3_histo.Integral(),qcd_b_sb_t3_histo.Integral()) #DEBUG
for i in range(1,21) :
	if qcd_a_sb_t1_histo.GetBinContent(i)<=0. : qcd_a_sb_t1_histo.SetBinContent(i,0.)
	if qcd_b_sb_t1_histo.GetBinContent(i)<=0. : qcd_b_sb_t1_histo.SetBinContent(i,0.)
	if qcd_a_sb_t2_histo.GetBinContent(i)<=0. : qcd_a_sb_t2_histo.SetBinContent(i,0.)
	if qcd_b_sb_t2_histo.GetBinContent(i)<=0. : qcd_b_sb_t2_histo.SetBinContent(i,0.)
	if qcd_a_sb_t3_histo.GetBinContent(i)<=0. : qcd_a_sb_t3_histo.SetBinContent(i,0.)
	if qcd_b_sb_t3_histo.GetBinContent(i)<=0. : qcd_b_sb_t3_histo.SetBinContent(i,0.)
QCD_transfer_factors['t1']=qcd_a_sb_t1_histo.Integral()/qcd_b_sb_t1_histo.Integral() if qcd_a_sb_t1_histo.Integral()!=0. else 0.
QCD_transfer_factors['t2']=qcd_a_sb_t2_histo.Integral()/qcd_b_sb_t2_histo.Integral() if qcd_a_sb_t2_histo.Integral()!=0. else 0.
QCD_transfer_factors['t3']=qcd_a_sb_t3_histo.Integral()/qcd_b_sb_t3_histo.Integral() if qcd_a_sb_t3_histo.Integral()!=0. else 0.
QCD_transfer_factors_up['t1']=QCD_transfer_factors['t1']*(1.+sqrt(1./qcd_a_sb_t1_histo.Integral()+1./qcd_b_sb_t1_histo.Integral())) if qcd_a_sb_t1_histo.Integral()!=0. and qcd_b_sb_t1_histo.Integral()!=0. else 0.
QCD_transfer_factors_up['t2']=QCD_transfer_factors['t2']*(1.+sqrt(1./qcd_a_sb_t2_histo.Integral()+1./qcd_b_sb_t2_histo.Integral())) if qcd_a_sb_t2_histo.Integral()!=0. and qcd_b_sb_t2_histo.Integral()!=0. else 0.
QCD_transfer_factors_up['t3']=QCD_transfer_factors['t3']*(1.+sqrt(1./qcd_a_sb_t3_histo.Integral()+1./qcd_b_sb_t3_histo.Integral())) if qcd_a_sb_t3_histo.Integral()!=0. and qcd_b_sb_t3_histo.Integral()!=0. else 0.
print 'FINAL QCD TRANSFER FACTORS = %s'%(QCD_transfer_factors)
#and for the MC-driven factors:
qcd_mc_chain_qcd_sb.Draw('cstar>>qcd_mc_qcd_a_sb_t1(20,-1.,1.)','('+qcd_a_cut+' && eventTopology==1)')
qcd_a_sb_t1_histo_MC_driven.Add(gROOT.FindObject('qcd_mc_qcd_a_sb_t1'))
qcd_mc_chain_qcd_sb.Draw('cstar>>qcd_mc_qcd_b_sb_t1(20,-1.,1.)','('+qcd_b_cut+' && eventTopology==1)')
qcd_b_sb_t1_histo_MC_driven.Add(gROOT.FindObject('qcd_mc_qcd_b_sb_t1'))
qcd_mc_chain_qcd_sb.Draw('cstar>>qcd_mc_qcd_a_sb_t2(20,-1.,1.)','('+qcd_a_cut+' && eventTopology==2)')
qcd_a_sb_t2_histo_MC_driven.Add(gROOT.FindObject('qcd_mc_qcd_a_sb_t2'))
qcd_mc_chain_qcd_sb.Draw('cstar>>qcd_mc_qcd_b_sb_t2(20,-1.,1.)','('+qcd_b_cut+' && eventTopology==2)')
qcd_b_sb_t2_histo_MC_driven.Add(gROOT.FindObject('qcd_mc_qcd_b_sb_t2'))
qcd_mc_chain_qcd_sb.Draw('cstar>>qcd_mc_qcd_a_sb_t3(20,-1.,1.)','('+qcd_a_cut+' && eventTopology==3)')
qcd_a_sb_t3_histo_MC_driven.Add(gROOT.FindObject('qcd_mc_qcd_a_sb_t3'))
qcd_mc_chain_qcd_sb.Draw('cstar>>qcd_mc_qcd_b_sb_t3(20,-1.,1.)','('+qcd_b_cut+' && eventTopology==3)')
qcd_b_sb_t3_histo_MC_driven.Add(gROOT.FindObject('qcd_mc_qcd_b_sb_t3'))
for i in range(1,21) :
	if qcd_a_sb_t1_histo_MC_driven.GetBinContent(i)<=0. : qcd_a_sb_t1_histo_MC_driven.SetBinContent(i,0.)
	if qcd_b_sb_t1_histo_MC_driven.GetBinContent(i)<=0. : qcd_b_sb_t1_histo_MC_driven.SetBinContent(i,0.)
	if qcd_a_sb_t2_histo_MC_driven.GetBinContent(i)<=0. : qcd_a_sb_t2_histo_MC_driven.SetBinContent(i,0.)
	if qcd_b_sb_t2_histo_MC_driven.GetBinContent(i)<=0. : qcd_b_sb_t2_histo_MC_driven.SetBinContent(i,0.)
	if qcd_a_sb_t3_histo_MC_driven.GetBinContent(i)<=0. : qcd_a_sb_t3_histo_MC_driven.SetBinContent(i,0.)
	if qcd_b_sb_t3_histo_MC_driven.GetBinContent(i)<=0. : qcd_b_sb_t3_histo_MC_driven.SetBinContent(i,0.)
QCD_transfer_factors_MC_driven['t1']=qcd_a_sb_t1_histo_MC_driven.Integral()/qcd_b_sb_t1_histo_MC_driven.Integral() if qcd_a_sb_t1_histo_MC_driven.Integral()!=0. else 0.
QCD_transfer_factors_MC_driven['t2']=qcd_a_sb_t2_histo_MC_driven.Integral()/qcd_b_sb_t2_histo_MC_driven.Integral() if qcd_a_sb_t2_histo_MC_driven.Integral()!=0. else 0.
QCD_transfer_factors_MC_driven['t3']=qcd_a_sb_t3_histo_MC_driven.Integral()/qcd_b_sb_t3_histo_MC_driven.Integral() if qcd_a_sb_t3_histo_MC_driven.Integral()!=0. else 0.
print 'FINAL QCD TRANSFER FACTORS FROM QCD MC = %s'%(QCD_transfer_factors_MC_driven)

#plot class
class Plot(object) :

	def __init__(self,name,varlist,title,nBins,low,hi,MC_weights=options.MC_weights,addl_cuts='',iPos=11,lPos=3,logy=False,resid_zoom=1) :
		#save init info
		self._name = name
		self._varlist = varlist
		self._title = title
		self._nBins = nBins
		self._low = low
		self._hi = hi
		self._oneline = TLine(self._low,1.,self._hi,1.)
		self._oneline.SetLineWidth(3)
		self._oneline.SetLineStyle(2)
		self._MC_weights = MC_weights
		self._addl_cuts = addl_cuts
		self._cutstring = 'weight!=0' #vanilla (charge-summed)
		#self._cutstring = 'weight!=0 && lep_Q>0.' #positively charged leptons only
		#self._cutstring = 'weight!=0 && lep_Q<0.' #negatively charged leptons only
		for var in varlist :
			self._cutstring+=' && '+var+'!=-900.'
		if addl_cuts!='' :
			self._cutstring+=' && '+addl_cuts
		#print 'cuts for '+name+' = '+self._cutstring #DEBUG
		#luminosity object position
		self._iPos = iPos #iPos = 10*(alignment) + position (1/2/3 = left/center/right)
		#legend position
		self._lPos = lPos #1=left, 2=middle, 3=right (default), 4=above everything in 3 columns
		#whether to put the y-axis on a log scale
		self._logy = logy
		#how far to zoom in on the residual plot axes
		self._resid_zoom = resid_zoom

	def plot(self,treedict,outfilep) :
		outfilep.cd()
		#declare histograms
		self._qcd_mc_histo = TH1D(self._name+'_qcd_mc',self._title,self._nBins,self._low,self._hi)
		self._qcd_dd_histo = TH1D(self._name+'_qcd_dd',self._title,self._nBins,self._low,self._hi)
		#declare the data-driven QCD uncertainty histogram
		self._qcd_dd_err_histo = TH1D(self._name+'_qcd_dd_err_histo',self._title,self._nBins,self._low,self._hi)
		#declare the residual plot
		self._resid = TH1D(self._name+'_resid',self._title,self._nBins,self._low,self._hi)
		#declare the data-driven QCD uncertainty residuals histogram
		self._qcd_dd_err_resid = TH1D(self._name+'_qcd_dd_err_resid',self._title,self._nBins,self._low,self._hi)
		#put histograms in memory
		self._qcd_mc_histo.SetDirectory(0)
		self._qcd_dd_histo.SetDirectory(0)
		#set histogram attributes 
		self._qcd_mc_histo.SetMarkerStyle(20)
		self._resid.SetStats(0); self._qcd_dd_err_resid.SetStats(0)
		self._resid.SetMarkerStyle(20)
		self._qcd_dd_err_histo.SetMarkerStyle(1); self._qcd_dd_err_resid.SetMarkerStyle(1)
		self._qcd_dd_err_histo.SetFillStyle(3013); self._qcd_dd_err_resid.SetFillStyle(3013)
		self._qcd_dd_err_histo.SetFillColor(kBlack); self._qcd_dd_err_resid.SetFillColor(kBlack)
		self._qcd_dd_err_resid.SetTitle(';%s;MC/data pred.'%(self._qcd_dd_err_resid.GetXaxis().GetTitle()))
		self._qcd_dd_histo.SetMarkerStyle(21)
		self._qcd_dd_histo.SetMarkerColor(kYellow)
		self._qcd_dd_histo.SetLineColor(kYellow)
		self._qcd_dd_histo.SetFillColor(kYellow)
		self._qcd_dd_histo.GetXaxis().SetLabelSize(0.)
		outfilep.cd()
		#declare the canvas
		self._canv = TCanvas(self._name+'_canv',self._name+'_canv',1100,900)
		self._canv.cd()
		#plot and get the MC-driven QCD histograms
		for var in self._varlist :
			interactivename = self._name+'_qcd_mc_'+var.replace('(','').replace(')','')
			#print '	drawing %s...'%(interactivename) #DEBUG
			treedict['QCD_MC_QCD_SB'].Draw('%s>>%s(%d,%f,%f)'%(var,interactivename,self._nBins,self._low,self._hi),'('+self._MC_weights+')*('+qcd_c_cut+' && '+self._cutstring+')*('+str(QCD_transfer_factors[self._name.split('_')[-1]])+')')
			newinthisto=gROOT.FindObject(interactivename)
			newinthisto.SetTitle(self._title)
			newinthisto.SetBinContent(1,newinthisto.GetBinContent(1)+newinthisto.GetBinContent(0))
			newinthisto.SetBinContent(newinthisto.GetNbinsX(),newinthisto.GetBinContent(newinthisto.GetNbinsX())+newinthisto.GetBinContent(newinthisto.GetNbinsX()+1))
			self._qcd_mc_histo.Add(newinthisto)
		self._qcd_mc_histo.SetTitle(self._title)
		#plot and get the data-driven QCD histograms
		dummy_qcd_up_histo=TH1D(self._name+'_QCD_up',self._title,self._nBins,self._low,self._hi)
		for var in self._varlist :
			interactivename = self._name+'_data_qcd_c_sb_'+var.replace('(','').replace(')','')
			#print '	drawing %s...'%(interactivename) #DEBUG
			treedict['DATA_QCD_SB'].Draw('%s>>%s(%d,%f,%f)'%(var,interactivename,self._nBins,self._low,self._hi),'('+qcd_c_cut+' && '+self._cutstring+')*('+str(QCD_transfer_factors[self._name.split('_')[-1]])+')')
			treedict['DATA_QCD_SB'].Draw('%s>>%s(%d,%f,%f)'%(var,interactivename+'_up',self._nBins,self._low,self._hi),'('+qcd_c_cut+' && '+self._cutstring+')*('+str(QCD_transfer_factors_up[self._name.split('_')[-1]])+')')
			newinthisto=gROOT.FindObject(interactivename)
			newinthisto.SetTitle(self._title)
			newinthisto.SetBinContent(1,newinthisto.GetBinContent(1)+newinthisto.GetBinContent(0))
			newinthisto.SetBinContent(newinthisto.GetNbinsX(),newinthisto.GetBinContent(newinthisto.GetNbinsX())+newinthisto.GetBinContent(newinthisto.GetNbinsX()+1))
			self._qcd_dd_histo.Add(newinthisto)
			newinthistoup=gROOT.FindObject(interactivename+'_up')
			newinthistoup.SetBinContent(1,newinthistoup.GetBinContent(1)+newinthistoup.GetBinContent(0))
			newinthistoup.SetBinContent(newinthistoup.GetNbinsX(),newinthistoup.GetBinContent(newinthistoup.GetNbinsX())+newinthistoup.GetBinContent(newinthistoup.GetNbinsX()+1))
			dummy_qcd_up_histo.Add(newinthistoup)
		for i in range(len(shortnames_done)) :
			if not treedict[shortnames_done[i].replace(' ','_').replace('#','')+'_QCD_SB'].GetEntries()>0 :
				continue
			for var in self._varlist :
				interactivename = self._name+'_mc_qcd_c_sb_'+str(i)+'_'+var.replace('(','').replace(')','')
				#print '	drawing %s for %s samples...'%(interactivename,shortnames_done[i]) #DEBUG
				if shortnames_done[i].find('t#bar{t}')!=-1 and options.usetopptrw=='yes' :
					treedict[shortnames_done[i].replace(' ','_').replace('#','')+'_QCD_SB'].Draw('%s>>%s(%d,%f,%f)'%(var,interactivename,self._nBins,self._low,self._hi),'('+self._MC_weights+')*('+qcd_c_cut+' && '+self._cutstring+')*(-1.*'+str(QCD_transfer_factors[self._name.split('_')[-1]])+')*sf_top_pt_rw_v2')
					treedict[shortnames_done[i].replace(' ','_').replace('#','')+'_QCD_SB'].Draw('%s>>%s(%d,%f,%f)'%(var,interactivename+'_up',self._nBins,self._low,self._hi),'('+self._MC_weights+')*('+qcd_c_cut+' && '+self._cutstring+')*(-1.*'+str(QCD_transfer_factors_up[self._name.split('_')[-1]])+')*sf_top_pt_rw_v2')
				else :
					treedict[shortnames_done[i].replace(' ','_').replace('#','')+'_QCD_SB'].Draw('%s>>%s(%d,%f,%f)'%(var,interactivename,self._nBins,self._low,self._hi),'('+self._MC_weights+')*('+qcd_c_cut+' && '+self._cutstring+')*(-1.*'+str(QCD_transfer_factors[self._name.split('_')[-1]])+')')
					treedict[shortnames_done[i].replace(' ','_').replace('#','')+'_QCD_SB'].Draw('%s>>%s(%d,%f,%f)'%(var,interactivename+'_up',self._nBins,self._low,self._hi),'('+self._MC_weights+')*('+qcd_c_cut+' && '+self._cutstring+')*(-1.*'+str(QCD_transfer_factors_up[self._name.split('_')[-1]])+')')
				#print '	findobject for name %s returns %s'%(interactivename,gROOT.FindObject(interactivename)) #DEBUG
				newinthisto=gROOT.FindObject(interactivename)
				newinthisto.SetTitle(self._title)
				newinthisto.SetBinContent(1,newinthisto.GetBinContent(1)+newinthisto.GetBinContent(0))
				newinthisto.SetBinContent(newinthisto.GetNbinsX(),newinthisto.GetBinContent(newinthisto.GetNbinsX())+newinthisto.GetBinContent(newinthisto.GetNbinsX()+1))
				self._qcd_dd_histo.Add(newinthisto)
				newinthistoup=gROOT.FindObject(interactivename+'_up')
				newinthistoup.SetBinContent(1,newinthistoup.GetBinContent(1)+newinthistoup.GetBinContent(0))
				newinthistoup.SetBinContent(newinthistoup.GetNbinsX(),newinthistoup.GetBinContent(newinthistoup.GetNbinsX())+newinthistoup.GetBinContent(newinthistoup.GetNbinsX()+1))
				dummy_qcd_up_histo.Add(newinthistoup)
			self._qcd_dd_histo.SetTitle(self._title)
		for i in range(1,self._nBins+1) :
			if self._qcd_dd_histo.GetBinContent(i)<=0. :
				self._qcd_dd_histo.SetBinContent(i,0.)
			#binerr = sqrt((dummy_qcd_up_histo.GetBinContent(i)-self._qcd_dd_histo.GetBinContent(i))**2+self._qcd_dd_histo.GetBinContent(i))
			binerr = sqrt((dummy_qcd_up_histo.GetBinContent(i)-self._qcd_dd_histo.GetBinContent(i))**2+(0.3*self._qcd_dd_histo.GetBinContent(i))**2)
			self._qcd_dd_histo.SetBinError(i,binerr)
		#set the contents and errors for the data-driven QCD uncertainty histograms
		for i in range(1,self._nBins+1) :
			content = self._qcd_dd_histo.GetBinContent(i)
			err2 = 0.
			entries = self._qcd_dd_histo.GetEffectiveEntries()
			integral = self._qcd_dd_histo.Integral(0,self._nBins+1)
			thishisterr = max(self._qcd_dd_histo.GetBinError(i),sqrt(abs(self._qcd_dd_histo.GetBinContent(i))))
			if entries!=0. :
				err2+=(thishisterr**2)#*(integral/abs(entries))
			err = sqrt(abs(err2)) if err2!=0. else 0.
			self._qcd_dd_err_histo.SetBinContent(i,content)
			self._qcd_dd_err_histo.SetBinError(i,err)
			if content!=0. :
				self._qcd_dd_err_resid.SetBinContent(i,content/content)
				self._qcd_dd_err_resid.SetBinError(i,err/content)
		#set the contents/errors of the residual histogram
		for i in range(1,self._nBins+1) :
			MCdrivencont = self._qcd_mc_histo.GetBinContent(i)
			dd_cont = self._qcd_dd_histo.GetBinContent(i)
			resid_cont = 0.; resid_err = 0;
			if dd_cont!=0. :
				resid_cont=MCdrivencont/abs(dd_cont)
				resid_err = sqrt(MCdrivencont)/abs(dd_cont)
			self._resid.SetBinContent(i,resid_cont)
			self._resid.SetBinError(i,resid_err)
		#set the maximum of the data-driven histogram to show the whole thing
		dd_max=max(self._qcd_dd_histo.GetMaximum()+sqrt(self._qcd_dd_histo.GetMaximum()),self._qcd_mc_histo.GetMaximum()+sqrt(self._qcd_mc_histo.GetMaximum()))
		if self._lPos==4 :
			self._qcd_dd_histo.SetMaximum(1.25*dd_max)
		else :
			self._qcd_dd_histo.SetMaximum(1.05*dd_max) 
		#make a legend
		legwidth = 0.2 if self._lPos!=4 else 0.4
		legheight = (2)*0.05 if self._lPos!=4 else ((2)*0.05)/3.5
		x2 = 0.9; y2 = 0.83
		if self._lPos==1 :
			x2 = 0.44
		elif self._lPos==2 :
			x2 = 0.7
		self._leg = TLegend(x2-legwidth,y2-legheight,x2,y2)
		if self._lPos==4 :
			self._leg.SetNColumns(3)
		#self._leg.AddEntry(self._data_histo,'DATA ('+str(self._data_histo.GetEntries())+' entries)','PE')
		self._leg.AddEntry(self._qcd_mc_histo,'MC QCD','PE')
		self._leg.AddEntry(self._qcd_dd_histo,'data-driven QCD','F')
		self._leg.AddEntry(self._qcd_dd_err_histo,'data-driven QCD uncertainty','F')
		#build the final plot
		self._canv.cd()
		gStyle.SetTitleFontSize(0.00001)
		self._histo_pad = TPad(self._name+'_histo_pad',self._name+'_histo_pad',0,0.25,1,1)
		self._resid_pad = TPad(self._name+'_resid_pad',self._name+'_resid_pad',0,0,1.,0.25)
		self._histo_pad.SetLeftMargin(0.16); self._histo_pad.SetRightMargin(0.05) 
		self._histo_pad.SetTopMargin(0.11);	 self._histo_pad.SetBottomMargin(0.02)
		self._histo_pad.SetBorderMode(0)
		self._resid_pad.SetLeftMargin(0.16); self._resid_pad.SetRightMargin(0.05)
		self._resid_pad.SetTopMargin(0.0);   self._resid_pad.SetBottomMargin(0.42)
		self._resid_pad.SetBorderMode(0)
		self._histo_pad.Draw(); self._resid_pad.Draw()
		self._histo_pad.cd()
		if self._logy :
			self._histo_pad.SetLogy(1)
		self._qcd_dd_histo.Draw('HIST'); self._qcd_dd_err_histo.Draw("SAME E2"); self._qcd_mc_histo.Draw("SAME PE")
		#print 'MC stack title = %s, MC stack histo title = %s, MC_err title = %s, _data_histo title = %s'%(self._MC_stack.GetTitle(),self._MC_stack.GetHistogram().GetTitle(),self._MC_err_histo.GetTitle(), self._data_histo.GetTitle()) #DEBUG
		if self._lPos!=0 :
			self._leg.Draw("SAME")
		#make the channel identifier text
		chanTxt = shortlepstring+' + jets'
		if mode=='t' :
			chanTxt+=' (trigger sideband)'
		elif mode=='wjetscr' :
			chanTxt+=' (W+Jets CR)'
		elif mode=='qcd_a_sr' :
			chanTxt+=' (QCD transfer numerator sb)'
		elif mode=='qcd_b_sr' :
			chanTxt+=' (QCD transfer denominator sb)'
		elif mode=='qcd_c_sr' :
			chanTxt+=' (QCD shape sb)'
		elif mode=='qcd_a_cr' :
			chanTxt+=' (QCD transfer numerator sb, W+Jets CR)'
		elif mode=='qcd_b_cr' :
			chanTxt+=' (QCD transfer denominator sb, W+Jets CR)'
		elif mode=='qcd_c_cr' :
			chanTxt+=' (QCD shape sb, W+Jets CR)'
		chantxt = TLatex()
		chantxt.SetNDC()
		chantxt.SetTextAngle(0)
		chantxt.SetTextColor(kBlack)
		chantxt.SetTextFont(42)
		chantxt.SetTextAlign(11) 
		chantxt.SetTextSize(0.6*0.11)
		chantxt.DrawLatex(0.16,1-0.11+0.2*0.11,chanTxt)
		self._qcd_dd_histo.GetXaxis().SetLabelSize(0.)
		#self._MC_stack.GetHistogram().GetXaxis().SetTitle('')
		self._resid_pad.cd()
		self._qcd_dd_err_resid.Draw("E2"); self._oneline.Draw("SAME"); self._resid.Draw("SAME PE")
		self._qcd_dd_err_resid.GetXaxis().SetLabelSize(0.15)
		self._qcd_dd_err_resid.GetYaxis().SetLabelSize(0.15)
		self._qcd_dd_err_resid.GetYaxis().SetTitleOffset(0.25)
		self._qcd_dd_err_resid.GetXaxis().SetTitleSize((0.75/0.25)*self._qcd_dd_err_resid.GetXaxis().GetTitleSize())
		self._qcd_dd_err_resid.GetYaxis().SetTitleSize((0.75/0.25)*self._qcd_dd_err_resid.GetYaxis().GetTitleSize())
		self._qcd_dd_err_resid.GetYaxis().SetRangeUser(0.2,1.8)
		self._qcd_dd_err_resid.GetYaxis().SetNdivisions(504)
		self._canv.Update()
		#plot the CMS_Lumi lines on the canvases
		all_lumi_objs.append(CMS_lumi.CMS_lumi(self._histo_pad, iPeriod, self._iPos))
		outfilep.cd()
		self._canv.Write()
		self._canv.Print(self._name+'_canv.pdf','pdf')
		#put it on the total canvas
		#total_canv.cd(all_plots.index(self)+1)
		#self._MC_stack.Draw(); self._MC_err_histo.Draw("SAME E2"); self._data_histo.Draw("SAME PE"); self._leg.Draw("SAME")
		#total_canv.Update()

	def getCanv(self) :
		return self._canv
	def getName(self) :
		return self._name

#list of CMS_lumi objects so they don't get deleted
all_lumi_objs = []
#define all the plots
all_plots = []
#total_canv = TCanvas('total_canv','total_canv',1320,900)
#total_canv.Divide(6,5)
#control plots
if mode=='' or mode=='wjetscr' or mode.startswith('qcd') :
	#if not (mode.startswith('qcd') and leptype=='muons') :
	#type-1
	if mode=='wjetscr' :
		all_plots.append(Plot('npv_t1',['npv'],'; # primary vertices; Events/bin',61,-0.5,60.5,addl_cuts='eventTopology==1'))
		all_plots.append(Plot('ngv_t1',['ngoodvtx'],'; # good vertices; Events/bin',41,-0.5,40.5,addl_cuts='eventTopology==1',lPos=4))
		all_plots.append(Plot('lbsf_t1',['par_2'],'; #lambda_{2}; Events/0.05',20,0.2,1.2,addl_cuts='eventTopology==1'))
		all_plots.append(Plot('hsfs_t1',['par_3','par_4'],'; #lambda_{3/4}; Events/0.05',20,0.2,1.2,addl_cuts='eventTopology==1',iPos=33,lPos=1))
		all_plots.append(Plot('chi2_t1',['chi2'],'; #chi^{2}; Events/5',50,-65.,185.,addl_cuts='eventTopology==1',iPos=33,lPos=0))
	else :
		all_plots.append(Plot('npv_t1',['npv'],'; # primary vertices; Events/bin',61,-0.5,60.5,addl_cuts='eventTopology==1',iPos=33,lPos=2))
		all_plots.append(Plot('ngv_t1',['ngoodvtx'],'; # good vertices; Events/bin',41,-0.5,40.5,addl_cuts='eventTopology==1',lPos=4))
		all_plots.append(Plot('lbsf_t1',['par_2'],'; #lambda_{2}; Events/0.01',40,0.75,1.15,addl_cuts='eventTopology==1',lPos=4))
		all_plots.append(Plot('hsfs_t1',['par_3','par_4'],'; #lambda_{3/4}; Events/0.01',40,0.8,1.2,addl_cuts='eventTopology==1',iPos=33,lPos=1))
		all_plots.append(Plot('chi2_t1',['chi2'],'; #chi^{2}; Events/1',28,-43.,-15.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('vpz_t1',['par_0'],'; p^{#nu}_{Z} [GeV]; Events/25 GeV',40,-500.,500.,addl_cuts='eventTopology==1'))
	all_plots.append(Plot('lsf_t1',['par_1'],'; #lambda_{1}; Events/0.0005',40,0.985,1.005,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('leppT_t1',['lep_pt'],';'+lepstring+' p_{T} [GeV]; Events/10 GeV',40,50.,450.,addl_cuts='eventTopology==1',iPos=33,lPos=2))
	all_plots.append(Plot('lepeta_t1',['lep_eta'],';'+lepstring+' #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==1',lPos=4))
	if leptype=='muons' or leptype=='all_leptons' :
		all_plots.append(Plot('ak41pT_t1',['ak41_pt'],'; AK4 jet1 p_{T} [GeV]; Events/20 GeV',40,150.,950.,addl_cuts='eventTopology==1',lPos=4))
		all_plots.append(Plot('ak42pT_t1',['ak42_pt'],'; AK4 jet2 p_{T} [GeV]; Events/10 GeV',55,50.,600.,addl_cuts='eventTopology==1',lPos=4))
	elif leptype=='electrons' :
		all_plots.append(Plot('ak41pT_t1',['ak41_pt'],'; AK4 jet1 p_{T} [GeV]; Events/20 GeV',35,250.,950.,addl_cuts='eventTopology==1',lPos=4))
		all_plots.append(Plot('ak42pT_t1',['ak42_pt'],'; AK4 jet2 p_{T} [GeV]; Events/10 GeV',53,70.,600.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak41eta_t1',['ak41_eta'],'; AK4 jet1 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak41csv_t1',['ak41_csvv2'],'; AK4 jet1 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==1',lPos=2))
	all_plots.append(Plot('ak42eta_t1',['ak42_eta'],'; AK4 jet2 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak42csv_t1',['ak42_csvv2'],'; AK4 jet2 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==1',lPos=2))
	all_plots.append(Plot('ak8pT_t1',['ak81_pt'],'; top-tagged AK8 jet p_{T} [GeV]; Events/10 GeV',50,400.,900.,addl_cuts='eventTopology==1',iPos=33,lPos=2))
	all_plots.append(Plot('ak8eta_t1',['ak81_eta'],'; top-tagged AK8 jet #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak8M_t1',['ak81_M'],'; top-tagged AK8 jet mass [GeV]; Events/5 GeV',30,100.,250.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak8tau32_t1',['ak81_tau32'],'; top-tagged AK8 jet #tau_{32}; Events/0.05',16,0.,0.8,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak8SDM_t1',['ak81_SDM'],'; top-tagged AK8 jet softdrop mass [GeV]; Events/5 GeV',23,105.,220.,addl_cuts='eventTopology==1',lPos=4))
	if leptype=='muons' or leptype=='all_leptons' :
		all_plots.append(Plot('MET_t1',['met_E'],'; MET [GeV]; Events/20 GeV',40,50.,850.,addl_cuts='eventTopology==1',lPos=4))
	elif leptype=='electrons' :
		all_plots.append(Plot('MET_t1',['met_E'],'; MET [GeV]; Events/20 GeV',40,100.,900.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('METphi_t1',['met_phi'],'; MET #phi; Events/0.2',31,-3.2,3.2,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('scaledMETphi_t1',['scaled_met_phi'],'; MET #phi (postfit); Events/0.2',31,-3.2,3.2,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('lepWHT_t1',['met_E+lep_pt'],'; p_{T}^{lep}+MET [GeV]; Events/20',50,0.,1000.,addl_cuts='eventTopology==1',lPos=0))
	all_plots.append(Plot('nak4_t1',['nak4jets'],'; # AK4 jets; Events/bin',9,1.5,10.5,addl_cuts='eventTopology==1',iPos=33,lPos=4))
	all_plots.append(Plot('nak8_t1',['nak8jets'],'; # AK8 jets; Events/bin',11,-0.5,10.5,addl_cuts='eventTopology==1',iPos=33,lPos=2))
	all_plots.append(Plot('nttags_t1',['ntTags'],'; # top-tagged AK8 jets; Events/bin',11,-0.5,10.5,addl_cuts='eventTopology==1',iPos=33,lPos=2))
	all_plots.append(Plot('nbtags_t1',['nLbTags'],'; # b-tagged AK4 jets; Events/bin',5,0.5,5.5,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('lepWpT_t1',['scaled_lepW_pt'],'; reconstructed W_{lep} p_{T} [GeV]; Events/20 GeV',40,0.,800.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('lepWM_t1',['scaled_lepW_M'],'; reconstructed W_{lep} mass [GeV]; Events/0.01 GeV',40,80.25,80.65,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('leptpT_t1',['scaled_lept_pt'],'; reconstructed t_{lep} p_{T} [GeV]; Events/40 GeV',25,0.,1000.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('leptM_t1',['scaled_lept_M'],'; reconstructed t_{lep} mass [GeV]; Events/5 GeV',20,110.,210.,addl_cuts='eventTopology==1',iPos=33,lPos=1))
	all_plots.append(Plot('hadtpT_t1',['scaled_hadt_pt'],'; reconstructed t_{had} p_{T} [GeV]; Events/20 GeV',35,350.,1050.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('hadtM_t1',['scaled_hadt_M'],'; reconstructed t_{had} mass [GeV]; Events/5 GeV',24,110.,230.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('cstar_t1',['cstar'],'; c*; Events/0.1',20,-1.,1.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('x_F_t1',['abs(x_F)'],'; |x_{F}|; Events/0.02',30,0.,0.6,addl_cuts='eventTopology==1',iPos=33,lPos=2))
	all_plots.append(Plot('M_t1',['M'],'; M; Events/100 GeV',25,500.,3000.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ttpt_t1',['tt_pt'],'; top pair p_{T} [GeV]; Events/20 GeV',30,0.,600.,addl_cuts='eventTopology==1',iPos=33,lPos=2))
	all_plots.append(Plot('lep_hadt_dR_t1',['lep_hadt_deltaR'],'; #Delta R ('+shortlepstring+', t_{had}); Events/0.1',30,1.5,4.5,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('lep_v_dPhi_t1',['lep_v_deltaPhi'],'; #Delta #phi ('+shortlepstring+', #nu); Events/0.1',30,-1.5,1.5,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('lept_hadt_dR_t1',['lept_hadt_deltaR'],'; #Delta R (t_{lep}, t_{had}); Events/0.1',20,2.,4.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak43pT_t1',['ak43_pt'],'; AK4 jet3 p_{T} [GeV]; Events/10 GeV',30,30.,330.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak43eta_t1',['ak43_eta'],'; AK4 jet3 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak43csv_t1',['ak43_csvv2'],'; AK4 jet3 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==1',lPos=2))
	all_plots.append(Plot('ak44pT_t1',['ak44_pt'],'; AK4 jet4 p_{T} [GeV]; Events/10 GeV',20,30.,230.,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak44eta_t1',['ak44_eta'],'; AK4 jet4 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==1',lPos=4))
	all_plots.append(Plot('ak44csv_t1',['ak44_csvv2'],'; AK4 jet4 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==1',lPos=2))
	all_plots.append(Plot('ak8tau21_t1',['ak81_tau21'],'; top-tagged AK8 jet #tau_{21}; Events/0.05',20,0.,1.,addl_cuts='eventTopology==1',lPos=4))
	#type-2
	if mode=='wjetscr' :
		all_plots.append(Plot('npv_t2',['npv'],'; # primary vertices; Events/bin',61,-0.5,60.5,addl_cuts='eventTopology==2'))
		all_plots.append(Plot('ngv_t2',['ngoodvtx'],'; # good vertices; Events/bin',41,-0.5,40.5,addl_cuts='eventTopology==2',lPos=4))
		all_plots.append(Plot('lbsf_t2',['par_2'],'; #lambda_{2}; Events/0.05',20,0.2,1.2,addl_cuts='eventTopology==2',iPos=33,lPos=1))
		all_plots.append(Plot('hsfs_t2',['par_3','par_4','par_5'],'; #lambda_{3/4/5}; Events/0.05',20,0.2,1.2,addl_cuts='eventTopology==2',iPos=33,lPos=1))
		all_plots.append(Plot('chi2_t2',['chi2'],'; #chi^{2}; Events/5',50,-65.,185.,addl_cuts='eventTopology==2'))
	else :
		all_plots.append(Plot('npv_t2',['npv'],'; # primary vertices; Events/bin',61,-0.5,60.5,addl_cuts='eventTopology==2',iPos=33,lPos=2))
		all_plots.append(Plot('ngv_t2',['ngoodvtx'],'; # good vertices; Events/bin',41,-0.5,40.5,addl_cuts='eventTopology==2',lPos=4))
		all_plots.append(Plot('lbsf_t2',['par_2'],'; #lambda_{2}; Events/0.02',25,0.7,1.2,addl_cuts='eventTopology==2',iPos=33,lPos=1))
		all_plots.append(Plot('hsfs_t2',['par_3','par_4','par_5'],'; #lambda_{3/4/5}; Events/0.02',25,0.7,1.2,addl_cuts='eventTopology==2',iPos=33,lPos=1))
		all_plots.append(Plot('chi2_t2',['chi2'],'; #chi^{2}; Events/1',28,-43.,-15.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('vpz_t2',['par_0'],'; p^{#nu}_{Z} [GeV]; Events/25 GeV',40,-500.,500.,addl_cuts='eventTopology==2'))
	all_plots.append(Plot('lsf_t2',['par_1'],'; #lambda_{1}; Events/0.0005',40,0.985,1.005,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('leppT_t2',['lep_pt'],';'+lepstring+' p_{T} [GeV]; Events/10 GeV',30,50.,350.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('lepeta_t2',['lep_eta'],';'+lepstring+' #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==2',lPos=4))
	if leptype=='muons' or leptype=='all_leptons' :
		all_plots.append(Plot('ak41pT_t2',['ak41_pt'],'; AK4 jet1 p_{T} [GeV]; Events/10 GeV',35,150.,500.,addl_cuts='eventTopology==2',lPos=4))
		all_plots.append(Plot('ak42pT_t2',['ak42_pt'],'; AK4 jet2 p_{T} [GeV]; Events/10 GeV',30,50.,350.,addl_cuts='eventTopology==2',lPos=4))
	elif leptype=='electrons' :
		all_plots.append(Plot('ak41pT_t2',['ak41_pt'],'; AK4 jet1 p_{T} [GeV]; Events/10 GeV',35,250.,600.,addl_cuts='eventTopology==2',lPos=4))
		all_plots.append(Plot('ak42pT_t2',['ak42_pt'],'; AK4 jet2 p_{T} [GeV]; Events/10 GeV',33,70.,400.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak41eta_t2',['ak41_eta'],'; AK4 jet1 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak41csv_t2',['ak41_csvv2'],'; AK4 jet1 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==2',lPos=2))
	all_plots.append(Plot('ak42eta_t2',['ak42_eta'],'; AK4 jet2 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak42csv_t2',['ak42_csvv2'],'; AK4 jet2 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==2',lPos=2))
	all_plots.append(Plot('ak8pT_t2',['ak81_pt'],'; AK8 jet1 p_{T} [GeV]; Events/10 GeV',50,200.,700.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak8eta_t2',['ak81_eta'],'; AK8 jet1 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak8M_t2',['ak81_M'],'; AK8 jet1 mass [GeV]; Events/5 GeV',40,40.,240.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak8tau32_t2',['ak81_tau32'],'; AK8 jet1 #tau_{32}; Events/0.05',20,0.,1.,addl_cuts='eventTopology==2',lPos=2))
	all_plots.append(Plot('ak8SDM_t2',['ak81_SDM'],'; AK8 jet1 softdrop mass [GeV]; Events/5 GeV',40,40.,240.,addl_cuts='eventTopology==2',lPos=4))
	if leptype=='muons' or leptype=='all_leptons' :
		all_plots.append(Plot('MET_t2',['met_E'],'; MET [GeV]; Events/10 GeV',45,50.,500.,addl_cuts='eventTopology==2',lPos=4))
	elif leptype=='electrons' :
		all_plots.append(Plot('MET_t2',['met_E'],'; MET [GeV]; Events/10 GeV',45,100.,550.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('METphi_t2',['met_phi'],'; MET #phi; Events/0.2',31,-3.2,3.2,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('scaledMETphi_t2',['scaled_met_phi'],'; MET #phi (postfit); Events/0.2',31,-3.2,3.2,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('lepWHT_t2',['met_E+lep_pt'],'; p_{T}^{lep}+MET [GeV]; Events/20',50,0.,1000.,addl_cuts='eventTopology==2',lPos=0))
	all_plots.append(Plot('nak4_t2',['nak4jets'],'; # AK4 jets; Events/bin',6,3.5,9.5,addl_cuts='eventTopology==2',iPos=33,lPos=2))
	all_plots.append(Plot('nak8_t2',['nak8jets'],'; # AK8 jets; Events/bin',11,-0.5,10.5,addl_cuts='eventTopology==2',iPos=33,lPos=2))
	all_plots.append(Plot('nttags_t2',['ntTags'],'; # top-tagged AK8 jets; Events/bin',10,0.5,10.5,addl_cuts='eventTopology==2',iPos=33,lPos=2))
	all_plots.append(Plot('nbtags_t2',['nLbTags'],'; # b-tagged AK4 jets; Events/bin',5,0.5,5.5,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('lepWpT_t2',['scaled_lepW_pt'],'; reconstructed W_{lep} p_{T} [GeV]; Events/20 GeV',30,0.,600.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('lepWM_t2',['scaled_lepW_M'],'; reconstructed W_{lep} mass [GeV]; Events/0.01 GeV',40,80.25,80.65,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('leptpT_t2',['scaled_lept_pt'],'; reconstructed t_{lep} p_{T} [GeV]; Events/20 GeV',40,0.,800.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('leptM_t2',['scaled_lept_M'],'; reconstructed t_{lep} mass [GeV]; Events/5 GeV',20,110.,210.,addl_cuts='eventTopology==2',iPos=33,lPos=1))
	all_plots.append(Plot('hadtpT_t2',['scaled_hadt_pt'],'; reconstructed t_{had} p_{T} [GeV]; Events/20 GeV',30,0.,600.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('hadtM_t2',['scaled_hadt_M'],'; reconstructed t_{had} mass [GeV]; Events/5 GeV',30,100.,250.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('cstar_t2',['cstar'],'; c*; Events/0.1',20,-1.,1.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('x_F_t2',['abs(x_F)'],'; |x_{F}|; Events/0.05',20,0.,0.4,addl_cuts='eventTopology==2',iPos=33,lPos=2))
	all_plots.append(Plot('M_t2',['M'],'; M; Events/100 GeV',20,300.,2300.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ttpt_t2',['tt_pt'],'; top pair p_{T} [GeV]; Events/20 GeV',30,0.,600.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('lep_hadt_dR_t2',['lep_hadt_deltaR'],'; #Delta R ('+shortlepstring+', t_{had}); Events/0.2',25,0.,5.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('lep_v_dPhi_t2',['lep_v_deltaPhi'],'; #Delta #phi ('+shortlepstring+', #nu); Events/0.2',20,-2.,2.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('lept_hadt_dR_t2',['lept_hadt_deltaR'],'; #Delta R (t_{lep}, t_{had}); Events/0.1',25,0.,5.,addl_cuts='eventTopology==2',iPos=33,lPos=1))
	all_plots.append(Plot('ak43pT_t2',['ak43_pt'],'; AK4 jet3 p_{T} [GeV]; Events/5 GeV',34,30.,200.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak43eta_t2',['ak43_eta'],'; AK4 jet3 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak43csv_t2',['ak43_csvv2'],'; AK4 jet3 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==2',lPos=2))
	all_plots.append(Plot('ak44pT_t2',['ak44_pt'],'; AK4 jet4 p_{T} [GeV]; Events/5 GeV',20,30.,130.,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak44eta_t2',['ak44_eta'],'; AK4 jet4 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==2',lPos=4))
	all_plots.append(Plot('ak44csv_t2',['ak44_csvv2'],'; AK4 jet4 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==2',lPos=2))
	all_plots.append(Plot('ak8tau21_t2',['ak81_tau21'],'; AK8 jet1 #tau_{21}; Events/0.05',20,0.,1.,addl_cuts='eventTopology==2',lPos=0))
	if not mode=='wjetscr' :
		#type-3
		all_plots.append(Plot('npv_t3',['npv'],'; # primary vertices; Events/bin',61,-0.5,60.5,addl_cuts='eventTopology==3',iPos=33,lPos=2))
		all_plots.append(Plot('ngv_t3',['ngoodvtx'],'; # good vertices; Events/bin',41,-0.5,40.5,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('vpz_t3',['par_0'],'; p^{#nu}_{Z} [GeV]; Events/25 GeV',40,-500.,500.,addl_cuts='eventTopology==3'))
		all_plots.append(Plot('lsf_t3',['par_1'],'; #lambda_{1}; Events/0.0005',40,0.985,1.005,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('lbsf_t3',['par_2'],'; #lambda_{2}; Events/0.02',30,0.6,1.2,addl_cuts='eventTopology==3',iPos=33,lPos=1))
		all_plots.append(Plot('hsfs_t3',['par_3','par_4','par_5'],'; #lambda_{3/4/5}; Events/0.02',30,0.6,1.2,addl_cuts='eventTopology==3',iPos=33,lPos=1,resid_zoom=2))
		all_plots.append(Plot('chi2_t3',['chi2'],'; #chi^{2}; Events/5',34,-45.,125.,addl_cuts='eventTopology==3',lPos=4,logy=True))
		all_plots.append(Plot('leppT_t3',['lep_pt'],';'+lepstring+' p_{T} [GeV]; Events/10 GeV',34,30.,200.,addl_cuts='eventTopology==3',lPos=4,resid_zoom=2))
		all_plots.append(Plot('lepeta_t3',['lep_eta'],';'+lepstring+' #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak41pT_t3',['ak41_pt'],'; AK4 jet1 p_{T} [GeV]; Events/5 GeV',40,30.,230.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak41eta_t3',['ak41_eta'],'; AK4 jet1 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak41csv_t3',['ak41_csvv2'],'; AK4 jet1 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==3',lPos=4,logy=True))
		all_plots.append(Plot('ak42pT_t3',['ak42_pt'],'; AK4 jet2 p_{T} [GeV]; Events/5 GeV',30,30.,180.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak42eta_t3',['ak42_eta'],'; AK4 jet2 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak42csv_t3',['ak42_csvv2'],'; AK4 jet2 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==3',lPos=4,logy=True))
		#all_plots.append(Plot('ak8pT_t3',['ak81_pt'],'; AK8 jet1 p_{T} [GeV]; Events/20 GeV',40,100.,900.,addl_cuts='eventTopology==3'))
		#all_plots.append(Plot('ak8eta_t3',['ak81_eta'],'; AK8 jet1 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==3',lPos=4))
		#all_plots.append(Plot('ak8M_t3',['ak81_M'],'; AK8 jet1 mass [GeV]; Events/10 GeV',40,0.0,400.0,addl_cuts='eventTopology==3'))
		#all_plots.append(Plot('ak8tau32_t3',['ak81_tau32'],'; AK8 jet1 #tau_{32}; Events/0.05 GeV',20,0.,1.,addl_cuts='eventTopology==3',lPos=0))
		#all_plots.append(Plot('ak8SDM_t3',['ak81_SDM'],'; AK8 jet1 softdrop mass [GeV]; Events/10 GeV',40,0.,400.,addl_cuts='eventTopology==3'))
		all_plots.append(Plot('MET_t3',['met_E'],'; MET [GeV]; Events/5 GeV',60,40.,340.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('METphi_t3',['met_phi'],'; MET #phi; Events/0.2',31,-3.2,3.2,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('scaledMETphi_t3',['scaled_met_phi'],'; MET #phi (postfit); Events/0.2',31,-3.2,3.2,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('lepWHT_t3',['met_E+lep_pt'],'; p_{T}^{lep}+MET [GeV]; Events/20',50,0.,1000.,addl_cuts='eventTopology==3',lPos=0))
		all_plots.append(Plot('nak4_t3',['nak4jets'],'; # AK4 jets; Events/bin',5,3.5,8.5,addl_cuts='eventTopology==3',iPos=33,lPos=2,resid_zoom=2))
		#all_plots.append(Plot('nak8_t3',['nak8jets'],'; # AK8 jets; Events/bin',11,-0.5,10.5,addl_cuts='eventTopology==3',iPos=33,lPos=2))
		all_plots.append(Plot('nttags_t3',['ntTags'],'; # top-tagged AK8 jets; Events/bin',1,-0.5,0.5,addl_cuts='eventTopology==3',iPos=33,lPos=0,logy=True))
		all_plots.append(Plot('nbtags_t3',['nMbTags'],'; # b-tagged AK4 jets; Events/bin',3,1.5,4.5,addl_cuts='eventTopology==3',lPos=4,logy=True))
		all_plots.append(Plot('lepWpT_t3',['scaled_lepW_pt'],'; reconstructed W_{lep} p_{T} [GeV]; Events/5 GeV',60,0.,300.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('lepWM_t3',['scaled_lepW_M'],'; reconstructed W_{lep} mass [GeV]; Events/0.02 GeV',50,80.,81.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('leptpT_t3',['scaled_lept_pt'],'; reconstructed t_{lep} p_{T} [GeV]; Events/10 GeV',40,0.,400.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('leptM_t3',['scaled_lept_M'],'; reconstructed t_{lep} mass [GeV]; Events/5 GeV',30,100.,250.,addl_cuts='eventTopology==3'))
		all_plots.append(Plot('hadtpT_t3',['scaled_hadt_pt'],'; reconstructed t_{had} p_{T} [GeV]; Events/5 GeV',60,0.,300.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('hadtM_t3',['scaled_hadt_M'],'; reconstructed t_{had} mass [GeV]; Events/5 GeV',40,100.,300.,addl_cuts='eventTopology==3'))
		all_plots.append(Plot('cstar_t3',['cstar'],'; c*; Events/0.1',20,-1.,1.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('x_F_t3',['abs(x_F)'],'; |x_{F}|; Events/0.05',15,0.,0.3,addl_cuts='eventTopology==3',iPos=33,lPos=2))
		all_plots.append(Plot('M_t3',['M'],'; M; Events/100 GeV',12,300.,1500.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ttpt_t3',['tt_pt'],'; top pair p_{T} [GeV]; Events/10 GeV',30,0.,300.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('lep_hadt_dR_t3',['lep_hadt_deltaR'],'; #Delta R ('+shortlepstring+', t_{had}); Events/0.2',30,0.,6.,addl_cuts='eventTopology==3'))
		all_plots.append(Plot('lep_v_dPhi_t3',['lep_v_deltaPhi'],'; #Delta #phi ('+shortlepstring+', #nu); Events/0.2',32,-3.2,3.2,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('lept_hadt_dR_t3',['lept_hadt_deltaR'],'; #Delta R (t_{lep}, t_{had}); Events/0.1',30,0.,6.,addl_cuts='eventTopology==3'))
		all_plots.append(Plot('ak43pT_t3',['ak43_pt'],'; AK4 jet3 p_{T} [GeV]; Events/5 GeV',20,30.,130.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak43eta_t3',['ak43_eta'],'; AK4 jet3 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak43csv_t3',['ak43_csvv2'],'; AK4 jet3 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==3',lPos=4,logy=True))
		all_plots.append(Plot('ak44pT_t3',['ak44_pt'],'; AK4 jet4 p_{T} [GeV]; Events/2 GeV',30,30.,90.,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak44eta_t3',['ak44_eta'],'; AK4 jet4 #eta; Events/0.1',48,-2.4,2.4,addl_cuts='eventTopology==3',lPos=4))
		all_plots.append(Plot('ak44csv_t3',['ak44_csvv2'],'; AK4 jet4 CSVv2; Events/0.025',40,0.0,1.0,addl_cuts='eventTopology==3',lPos=4,logy=True))
		#all_plots.append(Plot('ak8tau21_t3',['ak81_tau21'],'; AK8 jet1 #tau_{21}; Events/0.05 GeV',20,0.,1.,addl_cuts='eventTopology==3'))

#plot plots on canvases in parallel 
def plotparallel(thisplotlist,outfilenames,i) :
	#declare the tree dictionary
	thistreedict = {}
	#open tree files and add to dictionary
	fileps = []
	for sn in shortnames_done :
		thisfilep_SR = TFile(outname+'_'+sn.replace(' ','_').replace('#','')+'_sr_skimmed_tree.root')
		thistreedict[sn.replace(' ','_').replace('#','')+'_SR'] = thisfilep_SR.Get('tree')
		fileps.append(thisfilep_SR)
	datafilep_SR = TFile(outname+'_DATA_sr_skimmed_tree.root')
	thistreedict['DATA_SR'] = datafilep_SR.Get('tree')
	fileps.append(datafilep_SR)
	qcdmcfilep_SR = TFile(outname+'_QCD_MC_sr_skimmed_tree.root')
	thistreedict['QCD_MC_SR'] = qcdmcfilep_SR.Get('tree')
	fileps.append(qcdmcfilep_SR)
	for sn in shortnames_done :
		thisfilep_QCD_SB = TFile(outname+'_'+sn.replace(' ','_').replace('#','')+'_qcd_sb_skimmed_tree.root')
		thistreedict[sn.replace(' ','_').replace('#','')+'_QCD_SB'] = thisfilep_QCD_SB.Get('tree')
		fileps.append(thisfilep_QCD_SB)
	datafilep_QCD_SB = TFile(outname+'_DATA_qcd_sb_skimmed_tree.root')
	thistreedict['DATA_QCD_SB'] = datafilep_QCD_SB.Get('tree')
	fileps.append(datafilep_QCD_SB)
	qcdmcfilep_QCD_SB = TFile(outname+'_QCD_MC_qcd_sb_skimmed_tree.root')
	thistreedict['QCD_MC_QCD_SB'] = qcdmcfilep_QCD_SB.Get('tree')
	fileps.append(qcdmcfilep_QCD_SB)
	#start a new file for this subset of the plots
	thisoutfilename = outname+'_'+str(i)+'.root'
	outfilenames.append(thisoutfilename)
	thisoutfile = TFile(thisoutfilename,'recreate')
	#plot plots
	for thisplot in thisplotlist :
		print 'Plotting for plot '+thisplot.getName()+' ('+str(all_plots.index(thisplot)+1)+' out of '+str(len(all_plots))+')...'
		thisplot.plot(thistreedict,thisoutfile)
	#close tree files
	for filep in fileps :
		filep.Close()
	#close the file of plots
	thisoutfile.Close()
n_procs = options.n_threads
procs = []
outfilenames = manager.list()
for i in range(n_procs) :
	n_plots = len(all_plots)/n_procs if i<n_procs-1 else len(all_plots)-((n_procs-1)*(len(all_plots)/n_procs))
	firstplot = i*(len(all_plots)/n_procs)
	thisplotlist = all_plots[firstplot:(firstplot+n_plots)]
	#print 'plots for batch %d (firstplot=%d, n_plots=%d):'%(i,firstplot,n_plots) #DEBUG
	#for p in thisplotlist : #DEBUG
	#	print p.getName() #DEBUG
	p = multiprocessing.Process(target=plotparallel, args=(thisplotlist,outfilenames,i))
	p.start()
	procs.append(p)
for p in procs :
	p.join()
cmd = 'hadd -f '+outname+'.root '
for ofn in outfilenames :
	cmd+=ofn+' '
print cmd
os.system(cmd)
if os.path.isfile(outname+'.root') :
	for ofn in outfilenames :
		os.system('rm -rf '+ofn)

##plot plots on canvases
#for plot in all_plots :
#	plot.plot()

if not os.path.isdir(outname+'_pdfs') :
	os.system('mkdir '+outname+'_pdfs')
os.system('mv *_canv.pdf '+outname+'_pdfs')
for sn in shortnames_done :
	os.system('rm -rf '+outname+'_'+sn.replace(' ','_').replace('#','')+'_*_skimmed_tree.root')
os.system('rm -rf '+outname+'_DATA_*_skimmed_tree.root')
os.system('rm -rf '+outname+'_QCD_MC_*_skimmed_tree.root')
os.system('rm -rf '+garbageFileName)
