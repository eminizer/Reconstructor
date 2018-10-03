#Global variables
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0
#Trigger paths
#MU_TRIG_PATH = 'HLT_Mu30_eta2p1_PFJet150_PFJet50'
#MU_TRIG_PATH = 'HLT_Mu45_eta2p1'
MU_TRIG_PATHS_BOOSTED = ['HLT_Mu50','HLT_TkMu50']
MU_TRIG_PATHS_RESOLVED = ['HLT_IsoMu24','HLT_IsoTkMu24']
#EL_TRIG_PATH = 'HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50'
#EL_TRIG_PATHS = ['HLT_Ele45_WPLoose_Gsf']
EL_TRIG_PATHS_BOOSTED = ['HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165','HLT_Ele115_CaloIdVT_GsfTrkIdT']
EL_TRIG_PATHS_RESOLVED = ['HLT_Ele27_WPTight_Gsf']

##########								   Imports  								##########

from math import pi, log
import multiprocessing
from ROOT import TFile, TTree, TLorentzVector
from branch import Branch
from eventTypeHelper import getEventType, findInitialPartons, findMCParticles
from jet import AK4Jet, AK8Jet
from lepton import Muon, Electron
from metHelper import setupMET
from kinfit import reconstruct
from angleReconstructor import getObservables, getMCRWs
from corrector import Corrector
import sys
import gc

################################   addBranch function  #################################
def AddBranch(readname=None,writename=None,ttreetype='F',inival=-900.,size='1',dictlist=None) :
	newBranch = Branch(readname=readname,writename=writename,ttreetype=ttreetype,inival=inival,size=size)
	for d in dictlist :
		if writename!=None :
			d[writename] = newBranch
		else :
			d[readname] = newBranch
	return newBranch

##########							   Treemaker Class 								##########

class Reconstructor(object) :


	##################################  Branches  ##################################
	allBranches = {}
	ak4JetBranches = {}
	ak8JetBranches = {}
	genUncBranches = {}
	triggerBranches = {}
	mcGenEventBranches = {}
	vtxBranches = {}
	filterBranches = {}
	metBranches = {}
	muonBranches = {}
	electronBranches = {}
	#GenWeight info
	genWeight = AddBranch('evt_Gen_Weight','genWeight','F',1.0,'1',[allBranches])
	rho 	  = AddBranch('evt_rho','rho','D',1.0,'1',[allBranches,ak4JetBranches,ak8JetBranches])
	#pileup information
	ntrueint = AddBranch('pu_NtrueInt','ntrueint','I',-1,'1',[allBranches,ak4JetBranches,ak8JetBranches])
	npv 	 = AddBranch('vtx_npv','npv','I',-1,'1',[allBranches,ak4JetBranches,ak8JetBranches])
	ngoodvtx = AddBranch('evt_NGoodVtx','ngoodvtx','I',-1,'1',[allBranches,ak4JetBranches,ak8JetBranches])
	#top pT reweighting
	toppTreweight = AddBranch('evt_top_pt_rw','toppTrw','F',1.0,'1',[allBranches])
	#vertex information
	thisdictlist=[allBranches,vtxBranches]
	vtx_size = AddBranch(readname='vtx_size',ttreetype='i',dictlist=thisdictlist)
	vtx_ndof = AddBranch(readname='vtx_ndof',ttreetype='I',size='vtx_size',dictlist=thisdictlist)
	vtx_z = AddBranch(readname='vtx_z',size='vtx_size',dictlist=thisdictlist)
	vtx_rho = AddBranch(readname='vtx_rho',size='vtx_size',dictlist=thisdictlist)
	#Scale, PDF, and alpha_s weights
	thisdictlist=[allBranches,genUncBranches]
	scale_size = AddBranch(readname='scale_size',ttreetype='i',dictlist=thisdictlist)
	scale_weights = AddBranch(readname='scale_Weights',ttreetype='F',size='scale_size',dictlist=thisdictlist)
	pdf_size = AddBranch(readname='pdf_size',ttreetype='i',dictlist=thisdictlist)
	pdf_weights = AddBranch(readname='pdf_Weights',ttreetype='F',size='pdf_size',dictlist=thisdictlist)
	alphas_size = AddBranch(readname='alphas_size',ttreetype='i',dictlist=thisdictlist)
	alphas_weights = AddBranch(readname='alphas_Weights',ttreetype='F',size='alphas_size',dictlist=thisdictlist)
	#MC GenEvent info
	thisdictlist = [allBranches,mcGenEventBranches]
	MC_part1_factor = AddBranch(readname='MC_part1_factor',dictlist=thisdictlist)
	MC_part1_ID 	= AddBranch(readname='MC_part1_ID',dictlist=thisdictlist)
	MC_part2_factor = AddBranch(readname='MC_part2_factor',dictlist=thisdictlist)
	MC_part2_ID 	= AddBranch(readname='MC_part2_ID',dictlist=thisdictlist)
	MC_lep_ID 		= AddBranch(readname='MC_lep_ID',dictlist=thisdictlist)
	MC_cstar 		= AddBranch(readname='MC_cstar',dictlist=thisdictlist)
	MC_x_F 			= AddBranch(readname='MC_x_F',dictlist=thisdictlist)
	MC_Mtt 			= AddBranch(readname='MC_Mtt',dictlist=thisdictlist)
	fourvectornames = ['t','tbar','lep','nu','lepb','hadW','hadb']
	for n in fourvectornames :
		AddBranch(readname='MC_'+n+'_pt',dictlist=thisdictlist)
		AddBranch(readname='MC_'+n+'_eta',dictlist=thisdictlist)
		AddBranch(readname='MC_'+n+'_phi',dictlist=thisdictlist)
		AddBranch(readname='MC_'+n+'_E',dictlist=thisdictlist)
	#Trigger Information
	thisdictlist = [allBranches,triggerBranches]
	for trigName in MU_TRIG_PATHS_BOOSTED :
		AddBranch(trigName,trigName,'I',-1,'1',thisdictlist)
	for trigName in MU_TRIG_PATHS_RESOLVED :
		AddBranch(trigName,trigName,'I',-1,'1',thisdictlist)
	for trigName in EL_TRIG_PATHS_BOOSTED :
		AddBranch(trigName,trigName,'I',-1,'1',thisdictlist)
	for trigName in EL_TRIG_PATHS_RESOLVED :
		AddBranch(trigName,trigName,'I',-1,'1',thisdictlist)
	#MET Filter information
	thisdictlist = [allBranches,filterBranches]
	noise 				= AddBranch('Flag_HBHENoiseFilter','noisefilter','I',-1,'1',thisdictlist)
	noiseIso 			= AddBranch('Flag_HBHENoiseIsoFilter','noiseIsofilter','I',-1,'1',thisdictlist)
	ecaldeadcell 		= AddBranch('Flag_EcalDeadCellTriggerPrimitiveFilter','ecaldeadcellfilter','I',-1,'1',thisdictlist)
	goodVertices 		= AddBranch('Flag_goodVertices','goodVerticesfilter','I',-1,'1',thisdictlist)
	eebadsc 			= AddBranch('Flag_eeBadScFilter','eebadscfilter','I',-1,'1',thisdictlist)
	globaltighthalo2016 = AddBranch('Flag_globalTightHalo2016Filter','globaltighthalo2016filter','I',-1,'1',thisdictlist)
	badpfmuon 			= AddBranch('Flag_BadPFMuonFilter','badpfmuonfilter','I',-1,'1',thisdictlist)
	badchargedcandidate = AddBranch('Flag_BadChargedCandidateFilter','badchargedcandidatefilter','I',-1,'1',thisdictlist)
	#MET
	thisdictlist = [allBranches,metBranches]
	met_size 	= AddBranch(readname='met_size',ttreetype='i',dictlist=thisdictlist)
	met_pts 	= AddBranch(readname='met_Pt',size='met_size',dictlist=thisdictlist)
	met_phis 	= AddBranch(readname='met_Phi',size='met_size',dictlist=thisdictlist)
	#muons
	thisdictlist = [allBranches,muonBranches]
	mu_size 	  = AddBranch(readname='mu_size',ttreetype='i',dictlist=thisdictlist)
	mu_pts 		  = AddBranch(readname='mu_Pt',size='mu_size',dictlist=thisdictlist)
	mu_etas 	  = AddBranch(readname='mu_Eta',size='mu_size',dictlist=thisdictlist)
	mu_phis 	  = AddBranch(readname='mu_Phi',size='mu_size',dictlist=thisdictlist)
	mu_es 		  = AddBranch(readname='mu_E',size='mu_size',dictlist=thisdictlist)
	mu_charges 	  = AddBranch(readname='mu_Charge',size='mu_size',dictlist=thisdictlist)
	mu_isos 	  = AddBranch(readname='mu_Iso04',size='mu_size',dictlist=thisdictlist)
	mu_miniIsos   = AddBranch(readname='mu_MiniIso',size='mu_size',dictlist=thisdictlist)
	mus_isMed     = AddBranch(readname='mu_IsMediumMuon',size='mu_size',dictlist=thisdictlist)
	mus_isMed2016 = AddBranch(readname='mu_IsMediumMuon2016',size='mu_size',dictlist=thisdictlist)
	mu_Keys 	  = AddBranch(readname='mu_Key',size='mu_size',dictlist=thisdictlist)
	#electrons
	thisdictlist = [allBranches,electronBranches]
	el_size 	  = AddBranch(readname='el_size',ttreetype='i',dictlist=thisdictlist)
	el_pts 		  = AddBranch(readname='el_Pt',size='el_size',dictlist=thisdictlist)
	el_etas 	  = AddBranch(readname='el_Eta',size='el_size',dictlist=thisdictlist)
	el_scetas 	  = AddBranch(readname='el_SCEta',size='el_size',dictlist=thisdictlist)
	el_phis 	  = AddBranch(readname='el_Phi',size='el_size',dictlist=thisdictlist)
	el_es 		  = AddBranch(readname='el_E',size='el_size',dictlist=thisdictlist)
	el_charges 	  = AddBranch(readname='el_Charge',size='el_size',dictlist=thisdictlist)
	el_isos 	  = AddBranch(readname='el_Iso03',size='el_size',dictlist=thisdictlist)
	el_miniIsos   = AddBranch(readname='el_MiniIso',size='el_size',dictlist=thisdictlist)
	el_Dxys 	  = AddBranch(readname='el_Dxy',size='el_size',dictlist=thisdictlist)
	el_Dzs 		  = AddBranch(readname='el_Dz',size='el_size',dictlist=thisdictlist)
	el_id 		  = AddBranch(readname='el_IDMedium_NoIso',ttreetype='I',size='el_size',dictlist=thisdictlist)
	el_id_tight   = AddBranch(readname='el_IDTight_NoIso',ttreetype='I',size='el_size',dictlist=thisdictlist)
	el_Keys 	  = AddBranch(readname='el_Key',size='el_size',dictlist=thisdictlist)
	#AK4 Jets
	thisdictlist = [allBranches,ak4JetBranches]
	ak4_size 					= AddBranch(readname='jetAK4CHS_size',ttreetype='i',dictlist=thisdictlist)
	ak4_pts 					= AddBranch(readname='jetAK4CHS_Pt',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_etas 					= AddBranch(readname='jetAK4CHS_Eta',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_phis 					= AddBranch(readname='jetAK4CHS_Phi',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_es 						= AddBranch(readname='jetAK4CHS_E',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_genpts 					= AddBranch(readname='jetAK4CHS_GenJetPt',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_genetas 				= AddBranch(readname='jetAK4CHS_GenJetEta',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_genphis 				= AddBranch(readname='jetAK4CHS_GenJetPhi',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_genEs 					= AddBranch(readname='jetAK4CHS_GenJetE',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_csvv2s 					= AddBranch(readname='jetAK4CHS_CSVv2',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_jec0s 					= AddBranch(readname='jetAK4CHS_jecFactor0',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_jecuncs 				= AddBranch(readname='jetAK4CHS_jecUncertainty',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_jetAs 					= AddBranch(readname='jetAK4CHS_jetArea',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_jetPtRess 				= AddBranch(readname='jetAK4CHS_PtResolution',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_keylists 				= AddBranch(readname='jetAK4CHS_Keys',ttreetype='vi',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_neutralMultiplicity 	= AddBranch(readname='jetAK4CHS_neutralMultiplicity',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_neutralHadronEnergyFrac = AddBranch(readname='jetAK4CHS_neutralHadronEnergyFrac',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_neutralEmEnergyFrac 	= AddBranch(readname='jetAK4CHS_neutralEmEnergyFrac',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_chargedHadronEnergyFrac = AddBranch(readname='jetAK4CHS_chargedHadronEnergyFrac',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_chargedEmEnergyFrac 	= AddBranch(readname='jetAK4CHS_chargedEmEnergyFrac',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_chargedMultiplicity 	= AddBranch(readname='jetAK4CHS_chargedMultiplicity',size='jetAK4CHS_size',dictlist=thisdictlist)
	ak4_flavors 				= AddBranch(readname='jetAK4CHS_HadronFlavour',size='jetAK4CHS_size',dictlist=thisdictlist)
	#AK8 Jets
	thisdictlist = [allBranches,ak8JetBranches]
	ak8_size 					= AddBranch(readname='jetAK8CHS_size',ttreetype='i',dictlist=thisdictlist)
	ak8_pts 					= AddBranch(readname='jetAK8CHS_Pt',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_etas 					= AddBranch(readname='jetAK8CHS_Eta',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_phis 					= AddBranch(readname='jetAK8CHS_Phi',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_es 						= AddBranch(readname='jetAK8CHS_E',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_genpts 					= AddBranch(readname='jetAK8CHS_GenJetPt',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_genetas 				= AddBranch(readname='jetAK8CHS_GenJetEta',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_genphis 				= AddBranch(readname='jetAK8CHS_GenJetPhi',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_genEs 					= AddBranch(readname='jetAK8CHS_GenJetE',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_csvv2s 					= AddBranch(readname='jetAK8CHS_CSVv2',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_jec0s 					= AddBranch(readname='jetAK8CHS_jecFactor0',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_jecuncs 				= AddBranch(readname='jetAK8CHS_jecUncertainty',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_jetAs 					= AddBranch(readname='jetAK8CHS_jetArea',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_jetPtRess 				= AddBranch(readname='jetAK8CHS_PtResolution',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_tau1s 					= AddBranch(readname='jetAK8CHS_tau1CHS',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_tau2s 					= AddBranch(readname='jetAK8CHS_tau2CHS',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_tau3s 					= AddBranch(readname='jetAK8CHS_tau3CHS',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_sdms 					= AddBranch(readname='jetAK8CHS_softDropMassCHS',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_keylists 				= AddBranch(readname='jetAK8CHS_Keys',ttreetype='vi',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_neutralMultiplicity 	= AddBranch(readname='jetAK8CHS_neutralMultiplicity',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_neutralHadronEnergyFrac = AddBranch(readname='jetAK8CHS_neutralHadronEnergyFrac',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_neutralEmEnergyFrac 	= AddBranch(readname='jetAK8CHS_neutralEmEnergyFrac',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_chargedHadronEnergyFrac = AddBranch(readname='jetAK8CHS_chargedHadronEnergyFrac',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_chargedEmEnergyFrac 	= AddBranch(readname='jetAK8CHS_chargedEmEnergyFrac',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_chargedMultiplicity 	= AddBranch(readname='jetAK8CHS_chargedMultiplicity',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_subjetIndex0 			= AddBranch(readname='jetAK8CHS_vSubjetIndex0',size='jetAK8CHS_size',dictlist=thisdictlist)
	ak8_subjetIndex1 			= AddBranch(readname='jetAK8CHS_vSubjetIndex1',size='jetAK8CHS_size',dictlist=thisdictlist)
	#AK8 Subjets
	thisdictlist = [allBranches,ak8JetBranches]
	subjak8_size 					= AddBranch(readname='subjetAK8CHS_size',ttreetype='i',dictlist=thisdictlist)
	subjak8_pts 					= AddBranch(readname='subjetAK8CHS_Pt',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_etas 					= AddBranch(readname='subjetAK8CHS_Eta',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_phis 					= AddBranch(readname='subjetAK8CHS_Phi',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_es 						= AddBranch(readname='subjetAK8CHS_E',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_genpts 					= AddBranch(readname='subjetAK8CHS_GenJetPt',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_genetas 				= AddBranch(readname='subjetAK8CHS_GenJetEta',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_genphis 				= AddBranch(readname='subjetAK8CHS_GenJetPhi',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_genEs 					= AddBranch(readname='subjetAK8CHS_GenJetE',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_csvv2s 					= AddBranch(readname='subjetAK8CHS_CSVv2',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_jec0s 					= AddBranch(readname='subjetAK8CHS_jecFactor0',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_keylists 				= AddBranch(readname='subjetAK8CHS_Keys',ttreetype='vi',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_neutralMultiplicity 	= AddBranch(readname='subjetAK8CHS_neutralMultiplicity',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_neutralHadronEnergyFrac = AddBranch(readname='subjetAK8CHS_neutralHadronEnergyFrac',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_neutralEmEnergyFrac 	= AddBranch(readname='subjetAK8CHS_neutralEmEnergyFrac',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_chargedHadronEnergyFrac = AddBranch(readname='subjetAK8CHS_chargedHadronEnergyFrac',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_chargedEmEnergyFrac 	= AddBranch(readname='subjetAK8CHS_chargedEmEnergyFrac',size='subjetAK8CHS_size',dictlist=thisdictlist)
	subjak8_chargedMultiplicity 	= AddBranch(readname='subjetAK8CHS_chargedMultiplicity',size='subjetAK8CHS_size',dictlist=thisdictlist)
	#PRODUCES
	weightBranches = {}
	scalefactorBranches = {}
	physobjectBranches = {}
	observableBranches = {}
	mctruthBranches = {}
	cut_branches = {}
	cr_sb_cut_branches = {}
	kinfit_branches = {}
	#Weights
	thisdictlist = [allBranches,weightBranches]
	weight 	 = AddBranch(writename='weight',inival=1.,dictlist=thisdictlist)
	wg1 	 = AddBranch(writename='wg1',inival=1.,dictlist=thisdictlist)
	wg2 	 = AddBranch(writename='wg2',inival=1.,dictlist=thisdictlist)
	wg3 	 = AddBranch(writename='wg3',inival=1.,dictlist=thisdictlist)
	wg4 	 = AddBranch(writename='wg4',inival=1.,dictlist=thisdictlist)
	wqs1 	 = AddBranch(writename='wqs1',inival=1.,dictlist=thisdictlist)
	wqs2 	 = AddBranch(writename='wqs2',inival=1.,dictlist=thisdictlist)
	wqa0 	 = AddBranch(writename='wqa0',inival=1.,dictlist=thisdictlist)
	wqa1 	 = AddBranch(writename='wqa1',inival=1.,dictlist=thisdictlist)
	wqa2 	 = AddBranch(writename='wqa2',inival=1.,dictlist=thisdictlist)
	wg1_opp  = AddBranch(writename='wg1_opp',inival=1.,dictlist=thisdictlist)
	wg2_opp  = AddBranch(writename='wg2_opp',inival=1.,dictlist=thisdictlist)
	wg3_opp  = AddBranch(writename='wg3_opp',inival=1.,dictlist=thisdictlist)
	wg4_opp  = AddBranch(writename='wg4_opp',inival=1.,dictlist=thisdictlist)
	wqs1_opp = AddBranch(writename='wqs1_opp',inival=1.,dictlist=thisdictlist)
	wqs2_opp = AddBranch(writename='wqs2_opp',inival=1.,dictlist=thisdictlist)
	wqa0_opp = AddBranch(writename='wqa0_opp',inival=1.,dictlist=thisdictlist)
	wqa1_opp = AddBranch(writename='wqa1_opp',inival=1.,dictlist=thisdictlist)
	wqa2_opp = AddBranch(writename='wqa2_opp',inival=1.,dictlist=thisdictlist)
	wega 	 = AddBranch(writename='wega',inival=1.,dictlist=thisdictlist)
	wegc 	 = AddBranch(writename='wegc',inival=1.,dictlist=thisdictlist)
	#scalefactors
	thisdictlist = [allBranches,scalefactorBranches]
	sf_pileup 			  = AddBranch(writename='sf_pileup',inival=1.,dictlist=thisdictlist)
	sf_pileup_low 		  = AddBranch(writename='sf_pileup_low',inival=1.,dictlist=thisdictlist)
	sf_pileup_hi 		  = AddBranch(writename='sf_pileup_hi',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_BtoF 	  = AddBranch(writename='sf_trig_eff_BtoF',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_BtoF_low  = AddBranch(writename='sf_trig_eff_BtoF_low',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_BtoF_hi   = AddBranch(writename='sf_trig_eff_BtoF_hi',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_GH 		  = AddBranch(writename='sf_trig_eff_GH',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_GH_low 	  = AddBranch(writename='sf_trig_eff_GH_low',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_GH_hi 	  = AddBranch(writename='sf_trig_eff_GH_hi',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_BtoF 		  = AddBranch(writename='sf_lep_ID_BtoF',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_BtoF_low 	  = AddBranch(writename='sf_lep_ID_BtoF_low',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_BtoF_hi 	  = AddBranch(writename='sf_lep_ID_BtoF_hi',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_GH 		  = AddBranch(writename='sf_lep_ID_GH',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_GH_low 	  = AddBranch(writename='sf_lep_ID_GH_low',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_GH_hi 	  = AddBranch(writename='sf_lep_ID_GH_hi',inival=1.,dictlist=thisdictlist)
	sf_lep_iso_BtoF 	  = AddBranch(writename='sf_lep_iso_BtoF',inival=1.,dictlist=thisdictlist)
	sf_lep_iso_BtoF_low   = AddBranch(writename='sf_lep_iso_BtoF_low',inival=1.,dictlist=thisdictlist)
	sf_lep_iso_BtoF_hi 	  = AddBranch(writename='sf_lep_iso_BtoF_hi',inival=1.,dictlist=thisdictlist)
	sf_lep_iso_GH 		  = AddBranch(writename='sf_lep_iso_GH',inival=1.,dictlist=thisdictlist)
	sf_lep_iso_GH_low 	  = AddBranch(writename='sf_lep_iso_GH_low',inival=1.,dictlist=thisdictlist)
	sf_lep_iso_GH_hi 	  = AddBranch(writename='sf_lep_iso_GH_hi',inival=1.,dictlist=thisdictlist)
	sf_lep_mini_iso 	  = AddBranch(writename='sf_lep_mini_iso',inival=1.,dictlist=thisdictlist)
	sf_lep_mini_iso_low   = AddBranch(writename='sf_lep_mini_iso_low',inival=1.,dictlist=thisdictlist)
	sf_lep_mini_iso_hi 	  = AddBranch(writename='sf_lep_mini_iso_hi',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_flavb 	  = AddBranch(writename='sf_btag_eff_flavb',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_flavb_low = AddBranch(writename='sf_btag_eff_flavb_low',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_flavb_hi  = AddBranch(writename='sf_btag_eff_flavb_hi',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_flavc 	  = AddBranch(writename='sf_btag_eff_flavc',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_flavc_low = AddBranch(writename='sf_btag_eff_flavc_low',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_flavc_hi  = AddBranch(writename='sf_btag_eff_flavc_hi',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_heavy 	  = AddBranch(writename='sf_btag_eff_heavy',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_heavy_low = AddBranch(writename='sf_btag_eff_heavy_low',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_heavy_hi  = AddBranch(writename='sf_btag_eff_heavy_hi',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_light 	  = AddBranch(writename='sf_btag_eff_light',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_light_low = AddBranch(writename='sf_btag_eff_light_low',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_light_hi  = AddBranch(writename='sf_btag_eff_light_hi',inival=1.,dictlist=thisdictlist)
	sf_ttag_eff 		  = AddBranch(writename='sf_ttag_eff',inival=1.,dictlist=thisdictlist)
	sf_ttag_eff_low 	  = AddBranch(writename='sf_ttag_eff_low',inival=1.,dictlist=thisdictlist)
	sf_ttag_eff_hi 		  = AddBranch(writename='sf_ttag_eff_hi',inival=1.,dictlist=thisdictlist)
	sf_mu_R 			  = AddBranch(writename='sf_mu_R',inival=1.,dictlist=thisdictlist)
	sf_mu_R_low 		  = AddBranch(writename='sf_mu_R_low',inival=1.,dictlist=thisdictlist)
	sf_mu_R_hi 			  = AddBranch(writename='sf_mu_R_hi',inival=1.,dictlist=thisdictlist)
	sf_mu_F 			  = AddBranch(writename='sf_mu_F',inival=1.,dictlist=thisdictlist)
	sf_mu_F_low 		  = AddBranch(writename='sf_mu_F_low',inival=1.,dictlist=thisdictlist)
	sf_mu_F_hi 			  = AddBranch(writename='sf_mu_F_hi',inival=1.,dictlist=thisdictlist)
	sf_scale_comb 		  = AddBranch(writename='sf_scale_comb',inival=1.,dictlist=thisdictlist)
	sf_scale_comb_low 	  = AddBranch(writename='sf_scale_comb_low',inival=1.,dictlist=thisdictlist)
	sf_scale_comb_hi 	  = AddBranch(writename='sf_scale_comb_hi',inival=1.,dictlist=thisdictlist)
	sf_pdf_alphas 		  = AddBranch(writename='sf_pdf_alphas',inival=1.,dictlist=thisdictlist)
	sf_pdf_alphas_low 	  = AddBranch(writename='sf_pdf_alphas_low',inival=1.,dictlist=thisdictlist)
	sf_pdf_alphas_hi 	  = AddBranch(writename='sf_pdf_alphas_hi',inival=1.,dictlist=thisdictlist)
	sf_top_pt_rw 		  = AddBranch(writename='sf_top_pt_rw',inival=1.,dictlist=thisdictlist)
	sf_top_pt_rw_low 	  = AddBranch(writename='sf_top_pt_rw_low',inival=1.,dictlist=thisdictlist)
	sf_top_pt_rw_hi 	  = AddBranch(writename='sf_top_pt_rw_hi',inival=1.,dictlist=thisdictlist)
	#physics objects
	thisdictlist = [allBranches,physobjectBranches]
	#fourvectors
	fourvectornames = ['muon1','muon2','ele1','ele2','lep','met','ak41','ak42','ak43','ak44','ak81','ak82']
	for fourvecname in fourvectornames :
		AddBranch(writename=fourvecname+'_pt',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_eta',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_phi',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_M',dictlist=thisdictlist)
	metE = AddBranch(writename='met_E',dictlist=thisdictlist)
	ak41_csvv2 = AddBranch(writename='ak41_csvv2',dictlist=thisdictlist)
	ak42_csvv2 = AddBranch(writename='ak42_csvv2',dictlist=thisdictlist)
	ak43_csvv2 = AddBranch(writename='ak43_csvv2',dictlist=thisdictlist)
	ak44_csvv2 = AddBranch(writename='ak44_csvv2',dictlist=thisdictlist)
	#rescaled fourvectors
	fourvectornames = ['lep','met','lepW','lepb','lept','hadt','had1','had2','had3']
	for fourvecname in fourvectornames :
		AddBranch(writename='scaled_'+fourvecname+'_pt',dictlist=thisdictlist)
		AddBranch(writename='scaled_'+fourvecname+'_eta',dictlist=thisdictlist)
		AddBranch(writename='scaled_'+fourvecname+'_phi',dictlist=thisdictlist)
		AddBranch(writename='scaled_'+fourvecname+'_M',dictlist=thisdictlist)
	#others
	leptonnames = ['lep','muon1','muon2','ele1','ele2']
	for lepname in leptonnames :
		AddBranch(writename=lepname+'_Q',ttreetype='I',inival=0,dictlist=thisdictlist)
		AddBranch(writename=lepname+'_relPt',dictlist=thisdictlist)
		AddBranch(writename=lepname+'_dR',dictlist=thisdictlist)
		AddBranch(writename=lepname+'_Iso',dictlist=thisdictlist)
		AddBranch(writename=lepname+'_MiniIso',dictlist=thisdictlist)
	ak8names = ['ak81','ak82','scaled_hadt']
	for ak8name in ak8names :
		AddBranch(writename=ak8name+'_tau32',dictlist=thisdictlist)
		AddBranch(writename=ak8name+'_tau21',dictlist=thisdictlist)
		AddBranch(writename=ak8name+'_SDM',dictlist=thisdictlist)
		AddBranch(writename=ak8name+'_isttagged',ttreetype='i',inival=2,dictlist=thisdictlist)
		AddBranch(writename=ak8name+'_isWtagged',ttreetype='i',inival=2,dictlist=thisdictlist)
	#miscellaneous stuff
	thisdictlist = [allBranches]
	#lepton type in the event (1 for muon, 2 for electron)
	lepflavor = AddBranch(writename='lepflavor',ttreetype='i',inival=0,dictlist=thisdictlist)
	#how mamy unique MET solutions were there
	nMETs = AddBranch(writename='nMETs',ttreetype='i',inival=0,dictlist=thisdictlist)
	#whether or not this event should be added twice and have its weight halved based on its initial state
	addTwice = AddBranch(writename='addTwice',ttreetype='i',inival=0,dictlist=thisdictlist)
	#event type (0=qqbar semilep TT, 1=qg/gg semilep TT, 2=dileptonic TT, 3=hadronic TT, 4=other background, 5=UNKNOWN)
	event_type = AddBranch(writename='eventType',ttreetype='I',inival=5,dictlist=thisdictlist)
	#event topology (1=fully merged, 2=partially merged, 3=fully resolved)
	event_topology = AddBranch(writename='eventTopology',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of t-tagged AK8 jets in the event
	ntTags = AddBranch(writename='ntTags',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of W-tagged AK8 jets in the event
	nWTags = AddBranch(writename='nWTags',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of loose/medium b-tagged AK4 jets in the event
	nLbTags = AddBranch(writename='nLbTags',ttreetype='I',inival=0,dictlist=thisdictlist)
	nMbTags = AddBranch(writename='nMbTags',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of b-tagged AK4 jets used in the chosen jet assignment hypothesis
	nbTagsUsed = AddBranch(writename='nbTagsUsed',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of AK4 jets
	nak4jets = AddBranch(writename='nak4jets',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of AK8 jets
	nak8jets = AddBranch(writename='nak8jets',ttreetype='I',inival=0,dictlist=thisdictlist)
	#deltaR between postfit lepton/hadronic top
	lep_hadt_deltaR = AddBranch(writename='lep_hadt_deltaR',dictlist=thisdictlist)
	#deltaPhi between postfit lepton/neutrino
	lep_v_deltaPhi = AddBranch(writename='lep_v_deltaPhi',dictlist=thisdictlist)
	#deltaR between postfit leptonic/hadronic top
	lept_hadt_deltaR = AddBranch(writename='lept_hadt_deltaR',dictlist=thisdictlist)
	#top pair pT
	tt_pt = AddBranch(writename='tt_pt',dictlist=thisdictlist)
	#Obervables
	thisdictlist = [allBranches,observableBranches]
	#cosine(theta)
	cstar 			= AddBranch(writename='cstar',dictlist=thisdictlist)
	cstar_prefit 	= AddBranch(writename='cstar_prefit',dictlist=thisdictlist)
	cstar_corprefit = AddBranch(writename='cstar_corprefit',dictlist=thisdictlist)
	#Feynman x
	x_F 		  = AddBranch(writename='x_F',dictlist=thisdictlist)
	x_F_prefit 	  = AddBranch(writename='x_F_prefit',dictlist=thisdictlist)
	x_F_corprefit = AddBranch(writename='x_F_corprefit',dictlist=thisdictlist)
	#ttbar invariant mass
	M 			= AddBranch(writename='M',dictlist=thisdictlist)
	M_prefit 	= AddBranch(writename='M_prefit',dictlist=thisdictlist)
	M_corprefit = AddBranch(writename='M_corprefit',dictlist=thisdictlist)
	#initial quark vector
	thisdictlist = [allBranches,mctruthBranches]
	mctruthfourvectornames = ['q','qbar','MCt','MCtbar','MClep','MCv','MClepb','MChadb','MChadW']	
	for fourvecname in mctruthfourvectornames :
		AddBranch(writename=fourvecname+'_pt',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_eta',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_phi',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_M',dictlist=thisdictlist)
	#MC truth observables
	cstar_MC = AddBranch(writename='cstar_MC',dictlist=thisdictlist)
	x_F_MC 	 = AddBranch(writename='x_F_MC',dictlist=thisdictlist)
	M_MC 	 = AddBranch(writename='M_MC',dictlist=thisdictlist)
	#cut variables for full selection
	thisdictlist = [allBranches,cut_branches]
	cutnames = ['goodpv','metfilters','trigger','onelepton','isolepton','btags','ak4jetmult','jetcuts','METcuts','lepcuts','kinfitchi2','recoleptM','validminimization','fullselection']
	for cutname in cutnames :
		AddBranch(writename=cutname,ttreetype='i',inival=2,dictlist=thisdictlist)
	#cut variables for various control regions and sideband selections
	thisdictlist = [allBranches,cr_sb_cut_branches]
	cutnames = ['wjets_cr_selection','qcd_A_SR_selection','qcd_B_SR_selection','qcd_C_SR_selection','qcd_A_CR_selection','qcd_B_CR_selection','qcd_C_CR_selection']
	for cutname in cutnames :
		AddBranch(writename=cutname,ttreetype='i',inival=2,dictlist=thisdictlist)
	#cut variables etc. for electron trigger efficiency measurement
	#eltrig_cut_branches = {}
	#thisdictlist = [allBranches,eltrig_cut_branches]
	#cutnames = ['mutrigger','eltrigger','isomu','isoel','twoleptons','opplepcharge','btags','lepWpT','fullselection']
	#for cutname in cutnames :
	#	AddBranch(writename='eltrig_'+cutname,ttreetype='i',inival=2,dictlist=thisdictlist)
	#eltrig_themuon_pt = AddBranch(writename='eltrig_themuon_pt',dictlist=[allBranches])
	#eltrig_themuon_eta = AddBranch(writename='eltrig_themuon_eta',dictlist=[allBranches])
	#eltrig_theelectron_pt = AddBranch(writename='eltrig_theelectron_pt',dictlist=[allBranches])
	#eltrig_theelectron_eta = AddBranch(writename='eltrig_theelectron_eta',dictlist=[allBranches])
	#cut variables for alternate lepton isolation algorithm selections
	alt_lep_iso_cut_branches = {}
	thisdictlist = [allBranches,alt_lep_iso_cut_branches]
	cutnames = ['miniisolepton','oneminiisolepton','fullminiisoselection']
	for cutname in cutnames :
		AddBranch(writename=cutname,ttreetype='i',inival=2,dictlist=thisdictlist)
	#kinfit variables
	thisdictlist=[kinfit_branches,allBranches]
	chi2 = AddBranch(writename='chi2',inival=900.,dictlist=thisdictlist)
	par_0 = AddBranch(writename='par_0',dictlist=thisdictlist)
	par_1 = AddBranch(writename='par_1',dictlist=thisdictlist)
	par_2 = AddBranch(writename='par_2',dictlist=thisdictlist)
	par_3 = AddBranch(writename='par_3',dictlist=thisdictlist)
	par_4 = AddBranch(writename='par_4',dictlist=thisdictlist)
	par_5 = AddBranch(writename='par_5',dictlist=thisdictlist)
	nhypotheses = AddBranch(writename='nhypotheses',ttreetype='i',inival=0,dictlist=thisdictlist)
	ismatchable = AddBranch(writename='ismatchable',ttreetype='i',inival=0,dictlist=thisdictlist)
	iscorrect = AddBranch(writename='iscorrect',ttreetype='i',inival=0,dictlist=thisdictlist)
	ismatchedpostfit = AddBranch(writename='ismatchedpostfit',ttreetype='i',inival=0,dictlist=thisdictlist)
	lepWcorprefitM = AddBranch(writename='lepWcorprefitM',dictlist=thisdictlist)
	leptcorprefitM = AddBranch(writename='leptcorprefitM',dictlist=thisdictlist)
	hadtcorprefitM = AddBranch(writename='hadtcorprefitM',dictlist=thisdictlist)
	corlep_sf 	   = AddBranch(writename='corlep_sf',dictlist=thisdictlist)
	corlepb_sf 	   = AddBranch(writename='corlepb_sf',dictlist=thisdictlist)
	nearestJetPt = AddBranch(writename='nearestJetPt',dictlist=thisdictlist)
	wasCleanedFromNearestJet = AddBranch(writename='wasCleanedFromNearestJet',dictlist=thisdictlist)


	##################################  ANALYZE FUNCTION  ##################################
	def analyze(self,eventnumber) :

		#get the event in the tree
		self.inputTree.GetEntry(eventnumber)

		#-----------------------------------------------------Below here is a bunch of preselection and object assignment-----------------------------------------------------#
		#light preskim for requisite physics objects
		if not self.__hasRequisitePhysicsObjects__() :
			#print 'EVENT %d NOT VALID; MISSING REQUISITE PHYSICS OBJECTS!'%(eventnumber) #DEBUG
			return
		#MC stuff (event type, addTwice, eventweight, MC truth fourvectors for ttbar events)
		if not self.is_data :
			mctruthfourvecs = self.__getMCTruthProperties__()
		#if not self.event_type.getWriteValue()==2 :
		#	return
		#print '----------------------- event number %d -------------------------'%(eventnumber) #DEBUG
		#For the record, trigger information is handled automatically
		#initial raw MET
		#print '	Handling MET...' #DEBUG
		met = TLorentzVector(); met.SetPtEtaPhiM(self.met_pts.getReadValue(), 0., self.met_phis.getReadValue(), 0.)
		#initial muons and electrons
		muons, electrons = self.__makeLeptonLists__()
		allleps = muons+electrons
		allleps.sort(key=lambda x: x.getPt(), reverse=True)
		if len(allleps)<1 :
			#print 'EVENT NUMBER %d NOT VALID; MISSING LEPTONS (# muons = %d, # electrons = %d)'%(eventnumber,len(muons),len(electrons)) #DEBUG
			return
		#jets (adjusting met and counting btagged AK4 jets also)
		ak4jets, ak4jetsforisocalc, ak8jets, newmet, nLbtags, nMbtags, nttags = self.__makeJetLists__(met,allleps)
		met.SetPtEtaPhiM(newmet.Pt(),newmet.Eta(),newmet.Phi(),newmet.M())
		#if the lepton cleaning got rid of too many jets toss the event
		if not len(ak4jets)>0 :
			#print 'EVENT NUMBER %d NOT VALID; MISSING JETS (# AK4jets = %d, # AK8Jets = %d)'%(eventnumber,len(ak4jets),len(ak8jets)) #DEBUG
			return
		#calculate lepton isolation variables now that we know the list of valid jets
		#print '	Calculating lepton isolation values...' #DEBUG
		for lep in (allleps) :
			lep.calculateIsolation(ak4jetsforisocalc)
		#print '		Done.'#DEBUG
		#Set the event toplogy
		#print '	Setting event topology...' #DEBUG
		topology, ttags = self.__setEventTopology__(ak8jets)
		#print '		There are %d top-tagged jets; the event topology is type %d'%(len(ttags),topology) #DEBUG
		#is it possibble to reconstruct this event using a lepton+jets topology?
		canreconstruct = topology==1 or len(ak4jets)>=4
		#figure out whether the event is muonic or electronic, assign lepton
		lep = self.__assignLepton__(allleps,topology)
		#print 'LEP CHARGE = %d'%(int(lep.getQ())) #DEBUG
		#if not lep.getType()=='mu' : return #DEBUG
		#if not lep.getType()=='el' : return #DEBUG
		#Set physics object fourvectors
		self.__writePhysObjFourvecs__(muons,electrons,ak4jets,ttags,ak8jets,topology)
		#neutrino handling and setup for fit
		met1_vec, met2_vec = setupMET(lep.getFourVector(),met)
		self.nMETs.setWriteValue(2) if met1_vec.Pz() != met2_vec.Pz() else self.nMETs.setWriteValue(1)
		self.metE.setWriteValue(met1_vec.E())
		#-----------------------------------------------------------------Below here is event reconstruction-----------------------------------------------------------------#
		if canreconstruct :
			#print '	Reconstructing event...' #DEBUG
			#build the list of jet assignment hypotheses
			#checkmem('','before_hypotheses') #DEBUG
			hypotheses = getHypothesisList(topology,lep,met1_vec,met2_vec,self.nMETs.getWriteValue(),ak4jets[:5],ttags)
			#checkmem('',' after_hypotheses') #DEBUG
			if len(hypotheses)==0 :
				#if there's no valid type-1 hypotheses, try again as type-2
				if topology==1 and len(ak4jets)>3 :
					topology=2
					self.event_topology.setWriteValue(topology)
					hypotheses = getHypothesisList(topology,lep,met1_vec,met2_vec,self.nMETs.getWriteValue(),ak4jets[:5],ttags)
				#otherwise it's not reconstructable
				else :
					#print 'event %d invalid; no separated top tagged jets and not enough AK4 jets'%(eventnumber) #DEBUG
					self.cut_branches['fullselection'].setWriteValue(0)
					return
			#print 'size of %d hypotheses = %d, in event %d (topology %d) with %d AK4 jets '%(len(hypotheses),sys.getsizeof(hypotheses),eventnumber,topology,len(ak4jets)) #DEBUG
			#print '		Will try %d jet assignment hypotheses...'%(len(hypotheses)) #DEBUG
			self.nhypotheses.setWriteValue(len(hypotheses))
			#do the monte carlo matching
			corrhypindex=-1
			if self.event_type.getWriteValue()<2 and not self.is_data :
				corrhypindex = self.__matchJetAssignmentHypotheses__(topology,hypotheses,mctruthfourvecs)
			#send the hypotheses to the kinematic fit 
			hypindex, scaledlep, scaledmet, scaledlepb, scaledhadt, scaledhad1, scaledhad2, scaledhad3, fitchi2, finalpars = self.__writeKinfitResults__(hypotheses,topology)
			#print '		Done.' #DEBUG 
		else :
			self.cut_branches['validminimization'].setWriteValue(0)
			self.cut_branches['fullselection'].setWriteValue(0)
		#--------------------------------------------------------------------Below here is event selection--------------------------------------------------------------------#
		#print '	Calculating cut variables...' #DEBUG
		#FULL SELECTION CRITERIA
		if canreconstruct :
			self.__assignFullSelectionCutVars__(canreconstruct,topology,nLbtags,nMbtags,fitchi2,scaledlep+scaledmet+scaledlepb,lep,electrons,muons,ak4jets,ak8jets,scaledmet)
		else :
			self.__assignFullSelectionCutVars__(canreconstruct,topology,nLbtags,nMbtags,1000000.,None,lep,electrons,muons,ak4jets,ak8jets,met)
		#print 'FULL SELECTION = %d'%self.cut_branches['fullselection'].getWriteValue() #DEBUG
		#W+Jets Control Region Selections
		self.__assignWJetsCRSelectionCutVars__()
		#QCD sidebands for ABCD method background estimation (pass all cuts but MET and lepton isolation)
		self.__assignQCDSBSelectionCutVars__()
		#ELECTRON TRIGGER SAMPLE SELECTIONS
		#self.__assignElectronTriggerCRSelectionCutVars__(muons,electrons,topology,met1_vec,met2_vec)
		#ALTERNATE LEPTON ISOLATION (MiniIsolation) ALGORITHM SELECTIONS
		self.__assignMiniIsolationFullSelectionCutVars__(lep,muons,electrons)
		#print '		Done. (fullselection=%d)'%(self.cut_branches['fullselection'].getWriteValue()) #DEBUG
		#-----------------------------------------------------Below here is a bunch of variable and weight calculation-----------------------------------------------------# 
		if canreconstruct :
			if not self.is_data : 
				#see if the kinematic fit returned the matched hypothesis
				self.__checkForCorrectJetAssignment__(topology,hypindex,corrhypindex,hypotheses)
				#try the MC matching again with the postfit quantities 
				if self.event_type.getWriteValue()<2 : 
					hypothesis=hypotheses[hypindex] 
					MClep_vec = mctruthfourvecs['MClep']; MCv_vec = mctruthfourvecs['MCv']; MClepb_vec = mctruthfourvecs['MClepb']
					MChadt_vec = mctruthfourvecs['MCtbar'] if mctruthfourvecs['MClep_charge']>0 else mctruthfourvecs['MCt']
					self.ismatchedpostfit.setWriteValue(1) if (scaledlep.DeltaR(MClep_vec)<0.1 and scaledmet.DeltaPhi(MCv_vec)<0.3 and scaledlepb.DeltaR(MClepb_vec)<0.4 and scaledhadt.DeltaR(MChadt_vec)<0.8) else self.ismatchedpostfit.setWriteValue(0) 
				#	print 'ismatchedpostfit = %d'%(self.ismatchedpostfit.getWriteValue()) #DEBUG 
			#Write Kinematic fit debugging variables 
			self.chi2.setWriteValue(fitchi2) 
			locs = [self.par_0,self.par_1,self.par_2,self.par_3,self.par_4,self.par_5] 
			for i in range(len(finalpars)) : 
				locs[i].setWriteValue(finalpars[i]) 
			#write the scaled fourvectors
			self.__writeScaledFourvectors__(topology,scaledlep,scaledmet,scaledlepb,scaledhadt,scaledhad1,scaledhad2,scaledhad3,hypotheses[hypindex])
			#reconstruct the observables in multiple ways
			if not self.is_data :
				self.__reconstructObservables__(scaledlep,scaledmet,scaledlepb,scaledhadt,lep,hypotheses,hypindex,corrhypindex,mctruthfourvecs['MClep'],mctruthfourvecs['MClepb'])
			else :
				self.__reconstructObservables__(scaledlep,scaledmet,scaledlepb,scaledhadt,lep,hypotheses,hypindex,corrhypindex)
		#MC Truth observable and reweighting calculation 
		if not self.is_data :
			self.__calculateReweights__(mctruthfourvecs,topology,lep,ak4jets,nttags)
		#Finally write the event to the tree and closeout
		self.__closeout__() #yay! A complete event!

	##################################  Event Analyzer helper functions  ##################################

	def __hasRequisitePhysicsObjects__(self) :
		#Make sure the event has some met, at least one lepton, and at least one AK4 jet
		metsize = self.met_size.getReadValue()
		musize = self.mu_size.getReadValue()
		elsize = self.el_size.getReadValue()
		ak4jetsize = self.ak4_size.getReadValue()
		ak8jetsize = self.ak8_size.getReadValue()
		if not (metsize>0 and (musize>0 or elsize>0) and ak4jetsize>0) :
			#print 'EVENT NUMBER %d NOT VALID; MISSING REQUISITE PHYSICS OBJECTS (metsize = %d, musize = %d, elsize = %d, ak4jetsize = %d, ak8jetsize = %d)'%(eventnumber,metsize,musize,elsize,ak4jetsize,ak8jetsize)
			return False
		return True

	def __getMCTruthProperties__(self) :
		#event type split
		self.event_type.setWriteValue(getEventType(self.mcGenEventBranches))
		#set the addTwice value (true for qqbar and symmetric gg events)
		self.addTwice.setWriteValue(self.event_type.getWriteValue()==0 or (self.event_type.getWriteValue()==1 and self.mcGenEventBranches['MC_part1_ID'].getReadValue()==self.mcGenEventBranches['MC_part2_ID'].getReadValue()))
		#set the eventweight
		eventweight = self.genWeight.getReadValue()*(self.kfac*self.xsec/self.totalweight)
		self.weight.setWriteValue(eventweight)
		#Mother particle and MC truth top assignment
		vecdict = {}
		if self.event_type!=4 :
			q_vec, qbar_vec = findInitialPartons(self.mcGenEventBranches)
			MCt_vec, MCtbar_vec, MClep_vec, MCv_vec, MClepb_vec, MChadW_vec, MChadb_vec, MClep_charge = findMCParticles(self.mcGenEventBranches)				
			vecdict = {'q':q_vec,'qbar':qbar_vec,'MCt':MCt_vec,'MCtbar':MCtbar_vec,'MClep':MClep_vec,'MCv':MCv_vec,
					   'MClepb':MClepb_vec,'MChadW':MChadW_vec,'MChadb':MChadb_vec,'MClep_charge':MClep_charge}
		#write out fourvectors of MC particles
		for name,vec in vecdict.iteritems() :
			if vec!=None and name!='MClep_charge' :
				self.__setFourVectorBranchValues__(name,vec)
		return vecdict

	def __makeLeptonLists__(self) :
		#muons
		#print '	Handling Muons (%d total)...'%(self.mu_size.getReadValue()) #DEBUG
		#allmuons = [] #DEBUG
		muons = []
		for i in range(self.mu_size.getReadValue()) :
			newmuon=Muon(self.muonBranches,i,self.run_era)
		#	allmuons.append(newmuon) #DEBUG
			if newmuon.isValid() :
				muons.append(newmuon)
		#print '		Added %d Muons.'%(len(muons)) #DEBUG
		#electrons
		#print '	Handling Electrons (%d total)...'%(self.el_size.getReadValue()) #DEBUG
		#allelectrons = [] #DEBUG
		electrons = []
		for i in range(self.el_size.getReadValue()) :
			newele = Electron(self.electronBranches,i)
		#	allelectrons.append(newele) #DEBUG
			if newele.isValid() :
				electrons.append(newele)
		#print '		Added %d Electrons.'%(len(electrons)) #DEBUG
		#print '-----------------------------------------------------' #DEBUG
		#es = 'all eles: ' #DEBUG
		#for e in allelectrons : #DEBUG
		#	es+='(%.2f,%d) '%(e.getPt(),e.getQ()) #DEBUG
		#print es #DEBUG
		#es = 'selected eles: ' #DEBUG
		#for e in electrons : #DEBUG
		#	es+='(%.2f,%d) '%(e.getPt(),e.getQ()) #DEBUG
		#print es #DEBUG
		#ms = 'all mus: ' #DEBUG
		#for m in allmuons : #DEBUG
		#	ms+='(%.2f,%d) '%(m.getPt(),m.getQ()) #DEBUG
		#print ms #DEBUG
		#ms = 'selected mus: ' #DEBUG
		#for m in muons : #DEBUG
		#	ms+='(%.2f,%d) '%(m.getPt(),m.getQ()) #DEBUG
		#print ms #DEBUG
		return muons,electrons

	def __makeJetLists__(self,met,leplist) :
		#print '	Adding AK4 jets (%d total)...'%(self.ak4_size.getReadValue())#DEBUG
		ak4jets = []; ak4jetsforisocalc = []; ak8jets = []; metcorrvecs = []
		for i in range(self.ak4_size.getReadValue()) :
			newJet = AK4Jet(self.ak4JetBranches,i,self.JES,self.JER,leplist,self.corrector,self.is_data)
			metcorrvecs.append(newJet.getMETCorrectionVec())
			if newJet.isValidForIsoCalc() :
				ak4jetsforisocalc.append(newJet)
			if newJet.isValid() :
				ak4jets.append(newJet)
		#print '		Added %d AK4 Jets (%d for the isolation calcuations).'%(len(ak4jets),len(ak4jetsforisocalc)) #DEBUG
		#count the number of valid loose/medium b-tagged AK4 jets
		nLbtags = 0; nMbtags = 0
		for ak4jet in ak4jets :
			if ak4jet.isLbTagged() :
				nLbtags+=1
			if ak4jet.isMbTagged() :
				nMbtags+=1
		#print '		%d of the AK4 jets are loosely b-tagged.'%(nLbtags) #DEBUG
		self.nLbTags.setWriteValue(nLbtags)
		self.nMbTags.setWriteValue(nMbtags)
		#print '	Adding AK8 jets (%d total)...'%(self.ak8_size.getReadValue())#DEBUG
		for i in range(self.ak8_size.getReadValue()) :
			newJet = AK8Jet(self.ak8JetBranches,i,self.JES,self.JER,leplist,self.corrector,self.is_data)
			metcorrvecs.append(newJet.getMETCorrectionVec())
			if newJet.isValid() :
				ak8jets.append(newJet)
		#count the number of top tagged jets
		n_top_tags = 0
		for ak8jet in ak8jets :
			if ak8jet.isTopTagged() :
				n_top_tags+=1
		#sort the lists of jets by pT
		ak4jets.sort(key=lambda x: x.getPt(), reverse=True)
		ak4jetsforisocalc.sort(key=lambda x: x.getPt(), reverse=True)
		ak8jets.sort(key=lambda x: x.getPt(), reverse=True)
		self.nak4jets.setWriteValue(len(ak4jets)); self.nak8jets.setWriteValue(len(ak8jets));
		#print '		Added %d AK8 Jets.'%(len(ak8jets)) #DEBUG
		#recorrect the met
		#print '		met before = (%.1f,%.1f,%.1f,%.1f)'%(met.Pt(),met.Eta(),met.Phi(),met.M()) #DEBUG
		for metcorrvec in metcorrvecs :
			met+=metcorrvec
		met.SetPz(0.); met.SetE(met.Vect().Mag()) #reset the eta and mass of the MET
		#print '		met after = (%.1f,%.1f,%.1f,%.1f)'%(met.Pt(),met.Eta(),met.Phi(),met.M()) #DEBUG
		#now that the MET is adjusted set its branch values
		self.__setFourVectorBranchValues__('met',met)
		self.metE.setWriteValue(met.E())
		#return the stuff I'll need later
		return ak4jets,ak4jetsforisocalc,ak8jets,met,nLbtags,nMbtags,n_top_tags

	def __setEventTopology__(self,ak8jets) :
		#find top- and W-tags
		ttags = []; wtags = []
		for ak8jet in ak8jets :
			if ak8jet.isTopTagged() : 
				ttags.append(ak8jet)
			if ak8jet.isWTagged() : 
				wtags.append(ak8jet)
		n_ttags = len(ttags); n_Wtags = len(wtags)
		self.ntTags.setWriteValue(n_ttags); self.nWTags.setWriteValue(n_Wtags)
		if n_ttags>=1 : #If there is a t-tagged AK8 jet than this event has a type-1 (fully merged) topology
			topology=1
		elif len(ak8jets)>=1 : #If there are no t-tagged AK8 jets, but there IS at least one AK8 jet and four AK4 jets, the event is "boosted untagged" (type-2 topology)
			topology=2
		else : #If there are no AK8 jets at all, but there are at least 4 AK4 jets, the event is "resolved" (type-3 topology)
			topology=3
		#print 'n_ttags = %d, nAK8Jets = %d, topology = %d'%(n_ttags,len(ak8jets),topology)
		self.event_topology.setWriteValue(topology)
		return topology, ttags

	def __assignLepton__(self,allleps,topology) :
		lep = allleps[0]
		lepisiso = lep.is2DIso(topology) and (topology<3 or topology==3 and ((lep.getType=='mu' and lep.isTightIso()) or (lep.getType=='el' and lep.isIso())))
		if not lepisiso :
			for lepcand in allleps[1:] :
				lepcandisiso = lepcand.is2DIso(topology) and (topology<3 or topology==3 and ((lepcand.getType=='mu' and lepcand.isTightIso()) or (lepcand.getType=='el' and lepcand.isIso())))
				if lepcandisiso :
					lep = lepcand
					break
		if lep.getType()=='mu' : self.lepflavor.setWriteValue(1)
		elif lep.getType()=='el' : self.lepflavor.setWriteValue(2)
		else :
			#print 'EVENT NUMBER %d NOT VALID; LEPTON FLAVOR %s UNIDENTIFIED'%(eventnumber, lep.getType()) #DEBUG
			return
		self.nearestJetPt.setWriteValue(lep.getNearestJetPt())
		self.wasCleanedFromNearestJet.setWriteValue(lep.wasCleanedFromNearestJet())
		#write the analysis lepton variables
		self.__setLeptonBranchValues__('lep',lep)
		#print '	Lepton assigned (leading isolated %s).'%('muon' if self.lepflavor.getWriteValue()==1 else 'electron') #DEBUG
		return lep

	def __writePhysObjFourvecs__(self,muons,electrons,ak4jets,ttags,ak8jets,topology) :
		mulen = len(muons)
		if mulen>0 :
			self.__setLeptonBranchValues__('muon1',muons[0])
			if mulen>1 :
				self.__setLeptonBranchValues__('muon2',muons[1])
		ellen = len(electrons)
		if ellen>0 :
			self.__setLeptonBranchValues__('ele1',electrons[0])
			if ellen>1 :
				self.__setLeptonBranchValues__('ele2',electrons[1])
		self.__setFourVectorBranchValues__('ak41',ak4jets[0].getFourVector())
		self.ak41_csvv2.setWriteValue(ak4jets[0].getCSVv2())
		if len(ak4jets)>1 :
			self.__setFourVectorBranchValues__('ak42',ak4jets[1].getFourVector())
			self.ak42_csvv2.setWriteValue(ak4jets[1].getCSVv2())
			if len(ak4jets)>2 :
				self.__setFourVectorBranchValues__('ak43',ak4jets[2].getFourVector())
				self.ak43_csvv2.setWriteValue(ak4jets[2].getCSVv2())
				if len(ak4jets)>3 :
					self.__setFourVectorBranchValues__('ak44',ak4jets[3].getFourVector())
					self.ak44_csvv2.setWriteValue(ak4jets[3].getCSVv2())
		ak8jetwritelist = ttags if topology==1 else ak8jets
		if len(ak8jetwritelist)>0 : 
			self.__setFourVectorBranchValues__('ak81',ak8jetwritelist[0].getFourVector())
			self.allBranches['ak81_tau32'].setWriteValue(ak8jetwritelist[0].getTau32())
			self.allBranches['ak81_tau21'].setWriteValue(ak8jetwritelist[0].getTau21())
			self.allBranches['ak81_SDM'].setWriteValue(ak8jetwritelist[0].getSDM())
			self.allBranches['ak81_isttagged'].setWriteValue(1) if ak8jetwritelist[0].isTopTagged() else self.allBranches['ak81_isttagged'].setWriteValue(0)
			self.allBranches['ak81_isWtagged'].setWriteValue(1) if ak8jetwritelist[0].isWTagged() else self.allBranches['ak81_isWtagged'].setWriteValue(0)
			if len(ak8jetwritelist)>1 :
				self.__setFourVectorBranchValues__('ak82',ak8jetwritelist[1].getFourVector())
				self.allBranches['ak82_tau32'].setWriteValue(ak8jetwritelist[1].getTau32())
				self.allBranches['ak82_tau21'].setWriteValue(ak8jetwritelist[1].getTau21())
				self.allBranches['ak82_SDM'].setWriteValue(ak8jetwritelist[1].getSDM())
				self.allBranches['ak82_isttagged'].setWriteValue(1) if ak8jetwritelist[1].isTopTagged() else self.allBranches['ak82_isttagged'].setWriteValue(0)
				self.allBranches['ak82_isWtagged'].setWriteValue(1) if ak8jetwritelist[1].isWTagged() else self.allBranches['ak82_isWtagged'].setWriteValue(0)

	def __matchJetAssignmentHypotheses__(self,topology,hypotheses,mcvecs) :
		matchedhypindex=-1
		#print '----------------event %d--------------------'%(eventnumber) #DEBUG
		manager = multiprocessing.Manager()
		hypdict = manager.dict() #keys: hypothesis list indices. Values: sumdRdM
		all_parallel_matching_groups = []; j=-1
		batchsize = 1 if self.onGrid=='yes' else 5
		for i in range(len(hypotheses)) :
			if i%batchsize==0 :
				j+=1
				all_parallel_matching_groups.append([])
			newmatchinggroup = (hypotheses[i],i)
			all_parallel_matching_groups[j].append(newmatchinggroup)
		procs = []
		for pmg in all_parallel_matching_groups :
			for thisjob in pmg :
				p = multiprocessing.Process(target=self.__checkHypothesisMatchingParallel__,args=(topology,thisjob[0],mcvecs,thisjob[1],hypdict))
				p.start(); procs.append(p)
			for proc in procs :
				proc.join()
		#find the best-matched hypothesis
		mindRdM = 1000000.
		for index in hypdict.keys() :
			sumdRdM = hypdict[index]
			if sumdRdM<mindRdM :
				mindRdM = sumdRdM
				matchedhypindex = index
		#print '------------------------------------------------------------------------' #DEBUG
		if matchedhypindex==-1 :
			#print 'EVENT IS NOT MATCHABLE' #DEBUG
			self.ismatchable.setWriteValue(0)
		else : 
			#print 'Event has a correct assignment hypothesis at index %d'%(matchedhypindex) #DEBUG
			self.ismatchable.setWriteValue(1)
		#return the index of the best-matched hypothesis
		return matchedhypindex
	def __checkHypothesisMatchingParallel__(self,topology,hypothesis,mcvecs,i,hypdict) :
		#check the lepton 
		hyplvec = hypothesis[0].getFourVector()
		ldR = hyplvec.DeltaR(mcvecs['MClep'])
		ldM = hyplvec.M()-mcvecs['MClep'].M()
		if not ldR<0.1 : 
			hypdict[i]=1000000.; return
		#neutrino
		vdP = abs(hypothesis[1].DeltaPhi(mcvecs['MCv']))
		if not vdP<0.3 : 
			hypdict[i]=1000000.; return
		#leptonic b
		hyplbvec = hypothesis[2].getFourVector()
		lbdR = hyplbvec.DeltaR(mcvecs['MClepb'])
		lbdM = hyplbvec.M()-mcvecs['MClepb'].M()
		if not lbdR<0.4 : 
			hypdict[i]=1000000.; return
		#get the MC truth top quark vectors
		MClept_vec = mcvecs['MCt'] if mcvecs['MClep_charge']>0 else mcvecs['MCtbar']
		MChadt_vec = mcvecs['MCtbar'] if mcvecs['MClep_charge']>0 else mcvecs['MCt']
		#leptonic top
		hypltvec = hyplvec+hypothesis[1]+hyplbvec
		ltdR = hypltvec.DeltaR(MClept_vec)
		ltdM = hypltvec.M()-MClept_vec.M()
		#the rest are topology-dependent
		hyphtvec = hypothesis[3][0].getFourVector()
		if topology>1 : 
			hyphtvec = hyphtvec+hypothesis[3][1].getFourVector()+hypothesis[3][2].getFourVector()
		htdR = hyphtvec.DeltaR(MChadt_vec)
		htdM = hyphtvec.M()-MChadt_vec.M()
		sumdRdM = ltdR*abs(ltdM/MClept_vec.M())+htdR*abs(htdM/MChadt_vec.M())
		#print 'hypothesis number %d: ldR/dM=%.2f/%.2f, vdP=%.2f, lbdR/dM=%.2f/%.2f, ltdR/dM=%.2f/%.2f, htdR/dM=%.2f/%.2f, sumdRdM=%.2f, mindRdM=%.2f'%(i,ldR,ldM,vdP,lbdR,lbdM,ltdR,ltdM,htdR,htdM,sumdRdM,mindRdM) #DEBUG					
		#if the matching deltaR requirements are satisfied put this hypothesis index in the list with its sumdRdM
		if htdR<0.8 :
			hypdict[i]=sumdRdM
		else :
			hypdict[i]=1000000.
		return

	def __writeKinfitResults__(self,hypotheses,topology) :
		hypindex, scaledlep, scaledmet, scaledlepb, scaledhadt, scaledhad1, scaledhad2, scaledhad3, fitchi2, finalpars = reconstruct(hypotheses,topology,self.onGrid) 
		if scaledlep==None : 
			#print 'EVENT NUMBER '+str(eventnumber)+' NOT VALID; NO KINEMATIC FITS CONVERGED' #DEBUG 
			self.cut_branches['validminimization'].setWriteValue(0) 
			#print 'hypotheses = %s'%(hypotheses) #DEBUG
			hypindex = 0; scaledlep=hypotheses[0][0].getFourVector(); scaledmet=hypotheses[0][1]; scaledlepb=hypotheses[0][2].getFourVector()
			scaledhadt=hypotheses[0][3][0].getFourVector()
			if topology>1 :
				scaledhad1=hypotheses[0][3][0].getFourVector(); scaledhad2=hypotheses[0][3][1].getFourVector(); scaledhad3=hypotheses[0][3][2].getFourVector()
				scaledhadt=scaledhad1+scaledhad2+scaledhad3
		else :
			self.cut_branches['validminimization'].setWriteValue(1) 
		#set the number of bTags used 
		nusedbtags=0 
		if (topology<3 and hypotheses[hypindex][2].isLbTagged()) or hypotheses[hypindex][2].isMbTagged() : nusedbtags+=1 #leptonic side b-jet 
		for thisjet in hypotheses[hypindex][3] : #hadronic side jets 
			if (topology<3 and thisjet.isLbTagged()) or thisjet.isMbTagged() : nusedbtags+=1 
		self.nbTagsUsed.setWriteValue(nusedbtags) 
		#Set the scaled MET energy
		self.metE.setWriteValue(scaledmet.E())
		#set the postfit topology variables
		self.lep_hadt_deltaR.setWriteValue(scaledlep.DeltaR(scaledhadt))
		self.lep_v_deltaPhi.setWriteValue(scaledlep.DeltaPhi(scaledmet))
		self.lept_hadt_deltaR.setWriteValue((scaledlep+scaledmet+scaledlepb).DeltaR(scaledhadt))
		self.tt_pt.setWriteValue((scaledlep+scaledmet+scaledlepb+scaledhadt).Pt())
		return hypindex, scaledlep, scaledmet, scaledlepb, scaledhadt, scaledhad1, scaledhad2, scaledhad3, fitchi2, finalpars

	def __checkForCorrectJetAssignment__(self,topology,hypindex,corrhypindex,hypotheses) :
		#if the fit returned the correct hypothesis with either MET solution, it was correct 
		self.iscorrect.setWriteValue(0)
		if hypindex==corrhypindex :
			self.iscorrect.setWriteValue(1) 
		elif corrhypindex!=-1 :
			thishyp = hypotheses[hypindex]
			corrhyp = hypotheses[corrhypindex]
			if thishyp[0]==corrhyp[0] and thishyp[2]==corrhyp[2] and thishyp[3][0]==corrhyp[3][0] and (topology==1 or (thishyp[3][1]==corrhyp[3][1] and thishyp[3][2]==corrhyp[3][2])) :  
				self.iscorrect.setWriteValue(1) 

	def __writeScaledFourvectors__(self,topology,slep,smet,slepb,shadt,shad1,shad2,shad3,hyp) :
		#fill the TTree with the scaled fourvector variables 
		self.__setFourVectorBranchValues__('scaled_lep',slep) 
		self.__setFourVectorBranchValues__('scaled_met',smet) 
		self.__setFourVectorBranchValues__('scaled_lepW',slep+smet) 
		self.__setFourVectorBranchValues__('scaled_lepb',slepb) 
		self.__setFourVectorBranchValues__('scaled_lept',slep+smet+slepb) 
		self.__setFourVectorBranchValues__('scaled_hadt',shadt) 
		if topology==1 :
			self.allBranches['scaled_hadt_tau32'].setWriteValue(hyp[3][0].getTau32()) 
			self.allBranches['scaled_hadt_tau21'].setWriteValue(hyp[3][0].getTau21()) 
			self.allBranches['scaled_hadt_SDM'].setWriteValue(hyp[3][0].getSDM()) 
			self.docut('scaled_hadt_isttagged',hyp[3][0].isTopTagged(),self.allBranches)
			self.docut('scaled_hadt_isWtagged',hyp[3][0].isWTagged(),self.allBranches)
		if shad1!=None: self.__setFourVectorBranchValues__('scaled_had1',shad1) 
		if shad2!=None: self.__setFourVectorBranchValues__('scaled_had2',shad2) 
		if shad3!=None: self.__setFourVectorBranchValues__('scaled_had3',shad3) 

	def __reconstructObservables__(self,slep,smet,slepb,shadt,lep,hyps,hypi,chypi,mclepv=None,mclepbv=None) :	
		#first using the chosen postfit 
		cstar_s,xF_s,M_s = getObservables(slep+smet+slepb,shadt,lep.getQ())  
		self.cstar.setWriteValue(cstar_s); self.x_F.setWriteValue(xF_s); self.M.setWriteValue(M_s) 
		#again using the chosen prefit instead
		hyp = hyps[hypi] 
		prefitlept = hyp[0].getFourVector()+hyp[1]+hyp[2].getFourVector() 
		prefithadt = hyp[3][0].getFourVector() 
		for otherjet in hyp[3][1:] : 
			prefithadt = prefithadt+otherjet.getFourVector() 
		cstar,xF,M = getObservables(prefitlept,prefithadt,hyp[0].getQ())  
		self.cstar_prefit.setWriteValue(cstar); self.x_F_prefit.setWriteValue(xF); self.M_prefit.setWriteValue(M) 
		#finally using the correct assignment prefit
		#print 'isData = %s, corrhypindex=%d'%(self.is_data,corrhypindex) #DEBUG
		if (not self.is_data) and chypi!=-1 : 
			hyp = hyps[chypi]	
			prefitlepW = hyp[0].getFourVector()+hyp[1]
			prefitlept = prefitlepW+hyp[2].getFourVector()
			prefithadt = hyp[3][0].getFourVector()
			for otherjet in hyp[3][1:] :
				prefithadt = prefithadt+otherjet.getFourVector()
			#print 'correct hyp leptonic W mass = %.2f, leptonic top mass = %.2f, hadronic top mass = %.2f'%(prefitlepW.M(),prefitlept.M(),prefithadt.M()) #DEBUG
			self.lepWcorprefitM.setWriteValue(prefitlepW.M())
			self.leptcorprefitM.setWriteValue(prefitlept.M())
			self.hadtcorprefitM.setWriteValue(prefithadt.M())
			self.corlep_sf.setWriteValue(mclepv.P()/slep.P())
			self.corlepb_sf.setWriteValue(mclepbv.P()/slepb.P())
			cstar,xF,M = getObservables(prefitlept,prefithadt,hyp[0].getQ()) 
			#print 'setting correct prefit observable values: cstar = %.2f, x_F = %.2f, M = %.2f'%(cstar,x_F,M) #DEBUG
			self.cstar_corprefit.setWriteValue(cstar); self.x_F_corprefit.setWriteValue(xF); self.M_corprefit.setWriteValue(M)

	def __calculateReweights__(self,mctvecs,topology,lep,ak4jets,ntoptags) :
		#template reweights
		if self.event_type.getWriteValue()<2 : 
			self.cstar_MC.setWriteValue(self.mcGenEventBranches['MC_cstar'].getReadValue()) 
			self.x_F_MC.setWriteValue(self.mcGenEventBranches['MC_x_F'].getReadValue()) 
			self.M_MC.setWriteValue(self.mcGenEventBranches['MC_Mtt'].getReadValue()) 
			if self.event_type.getWriteValue()!=4 : 
				( wg1,wg2,wg3,wg4,wqs1,wqs2,wqa0,wqa1,wqa2, 
					wg1_opp,wg2_opp,wg3_opp,wg4_opp,wqs1_opp,wqs2_opp,wqa0_opp,wqa1_opp,wqa2_opp, 
					wega, wegc ) = getMCRWs(self.cstar_MC.getWriteValue(),mctvecs['MCt'],mctvecs['MCtbar'],self.corrector)  
			self.wg1.setWriteValue(wg1) 
			self.wg2.setWriteValue(wg2) 
			self.wg3.setWriteValue(wg3) 
			self.wg4.setWriteValue(wg4) 
			self.wqs1.setWriteValue(wqs1) 
			self.wqs2.setWriteValue(wqs2) 
			self.wqa0.setWriteValue(wqa0) 
			self.wqa1.setWriteValue(wqa1) 
			self.wqa2.setWriteValue(wqa2) 
			self.wg1_opp.setWriteValue(wg1_opp) 
			self.wg2_opp.setWriteValue(wg2_opp) 
			self.wg3_opp.setWriteValue(wg3_opp) 
			self.wg4_opp.setWriteValue(wg4_opp) 
			self.wqs1_opp.setWriteValue(wqs1_opp) 
			self.wqs2_opp.setWriteValue(wqs2_opp) 
			self.wqa0_opp.setWriteValue(wqa0_opp) 
			self.wqa1_opp.setWriteValue(wqa1_opp) 
			self.wqa2_opp.setWriteValue(wqa2_opp) 
			self.wega.setWriteValue(wega) 
			self.wegc.setWriteValue(wegc) 
		#scale factor and reweighting calculations
		#Pileup reweighting
		pu_sf, pu_sf_up, pu_sf_down = self.corrector.getPileupReweight(self.allBranches['ntrueint'].getReadValue())
		self.sf_pileup.setWriteValue(pu_sf); self.sf_pileup_hi.setWriteValue(pu_sf_up); self.sf_pileup_low.setWriteValue(pu_sf_down)
		#Trigger efficiency reweighting
		( trig_eff_BtoF_sf, trig_eff_BtoF_sf_up, trig_eff_BtoF_sf_down, 
			trig_eff_GH_sf, trig_eff_GH_sf_up, trig_eff_GH_sf_down ) = self.corrector.getTrigEff(topology,lep)
		self.sf_trig_eff_BtoF.setWriteValue(trig_eff_BtoF_sf); self.sf_trig_eff_BtoF_hi.setWriteValue(trig_eff_BtoF_sf_up); self.sf_trig_eff_BtoF_low.setWriteValue(trig_eff_BtoF_sf_down)
		self.sf_trig_eff_GH.setWriteValue(trig_eff_GH_sf); 	   self.sf_trig_eff_GH_hi.setWriteValue(trig_eff_GH_sf_up); 	self.sf_trig_eff_GH_low.setWriteValue(trig_eff_GH_sf_down)
		#Lepton ID efficiency reweighting
		( lep_ID_BtoF_sf, lep_ID_BtoF_sf_up, lep_ID_BtoF_sf_down, 
			lep_ID_GH_sf, lep_ID_GH_sf_up, lep_ID_GH_sf_down ) = self.corrector.getIDEff(self.allBranches['npv'].getReadValue(),lep)
		self.sf_lep_ID_BtoF.setWriteValue(lep_ID_BtoF_sf); self.sf_lep_ID_BtoF_hi.setWriteValue(lep_ID_BtoF_sf_up); self.sf_lep_ID_BtoF_low.setWriteValue(lep_ID_BtoF_sf_down)
		self.sf_lep_ID_GH.setWriteValue(lep_ID_GH_sf); 	   self.sf_lep_ID_GH_hi.setWriteValue(lep_ID_GH_sf_up); 	self.sf_lep_ID_GH_low.setWriteValue(lep_ID_GH_sf_down)
		#Lepton isolation efficiency reweighting
		( lep_iso_BtoF_sf, lep_iso_BtoF_sf_up, lep_iso_BtoF_sf_down, 
			lep_iso_GH_sf, lep_iso_GH_sf_up, lep_iso_GH_sf_down ) = self.corrector.getIsoEff(self.allBranches['npv'].getReadValue(),topology,lep)
		self.sf_lep_iso_BtoF.setWriteValue(lep_iso_BtoF_sf); self.sf_lep_iso_BtoF_hi.setWriteValue(lep_iso_BtoF_sf_up); self.sf_lep_iso_BtoF_low.setWriteValue(lep_iso_BtoF_sf_down)
		self.sf_lep_iso_GH.setWriteValue(lep_iso_GH_sf); 	 self.sf_lep_iso_GH_hi.setWriteValue(lep_iso_GH_sf_up); 	self.sf_lep_iso_GH_low.setWriteValue(lep_iso_GH_sf_down)
		#Lepton mini isolation efficiency reweighting
		lep_mini_iso_sf, lep_mini_iso_sf_up, lep_mini_iso_sf_down = self.corrector.getMiniIsoEff(self.allBranches['npv'].getReadValue(),lep)
		self.sf_lep_mini_iso.setWriteValue(lep_mini_iso_sf); self.sf_lep_mini_iso_hi.setWriteValue(lep_mini_iso_sf_up); self.sf_lep_mini_iso_low.setWriteValue(lep_mini_iso_sf_down)
		#b-tagging efficiency reweighting
		( btag_eff_flavb_sf, btag_eff_flavb_sf_up, btag_eff_flavb_sf_down, 
			btag_eff_flavc_sf, btag_eff_flavc_sf_up, btag_eff_flavc_sf_down, 
			btag_eff_heavy_sf, btag_eff_heavy_sf_up, btag_eff_heavy_sf_down, 
			btag_eff_light_sf, btag_eff_light_sf_up, btag_eff_light_sf_down ) = self.corrector.getBTagEff(topology,ak4jets)
		self.sf_btag_eff_flavb.setWriteValue(btag_eff_flavb_sf); self.sf_btag_eff_flavb_hi.setWriteValue(btag_eff_flavb_sf_up); self.sf_btag_eff_flavb_low.setWriteValue(btag_eff_flavb_sf_down)
		self.sf_btag_eff_flavc.setWriteValue(btag_eff_flavc_sf); self.sf_btag_eff_flavc_hi.setWriteValue(btag_eff_flavc_sf_up); self.sf_btag_eff_flavc_low.setWriteValue(btag_eff_flavc_sf_down)
		self.sf_btag_eff_heavy.setWriteValue(btag_eff_heavy_sf); self.sf_btag_eff_heavy_hi.setWriteValue(btag_eff_heavy_sf_up); self.sf_btag_eff_heavy_low.setWriteValue(btag_eff_heavy_sf_down)
		self.sf_btag_eff_light.setWriteValue(btag_eff_light_sf); self.sf_btag_eff_light_hi.setWriteValue(btag_eff_light_sf_up); self.sf_btag_eff_light_low.setWriteValue(btag_eff_light_sf_down)
		#top-tagging efficiency reweighting (only for type-1/2 events since type-3 events have no AK8 jets)
		if topology!=3 and ntoptags!=0 :
			#SF and uncertainties come from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetTopTagging
			ttag_eff_sf, ttag_eff_sf_up, ttag_eff_sf_down = ntoptags*(1.06), ntoptags*(1.06+0.09), ntoptags*(1.06-0.04)
			self.sf_ttag_eff.setWriteValue(ttag_eff_sf); self.sf_ttag_eff_hi.setWriteValue(ttag_eff_sf_up); self.sf_ttag_eff_low.setWriteValue(ttag_eff_sf_down)
		#Scale, pdf/alpha_s, and top pT reweights (ttbar only)
		if self.event_type.getWriteValue()<4 :
			( mu_R_sf, mu_R_sf_up, mu_R_sf_down, mu_F_sf, mu_F_sf_up, mu_F_sf_down,
			 scale_comb_sf, scale_comb_sf_up, scale_comb_sf_down, pdf_alphas_sf, pdf_alphas_sf_up, pdf_alphas_sf_down ) = self.corrector.getGenReweights(self.genUncBranches)
			self.sf_mu_R.setWriteValue(mu_R_sf); self.sf_mu_R_hi.setWriteValue(mu_R_sf_up); self.sf_mu_R_low.setWriteValue(mu_R_sf_down)
			self.sf_mu_F.setWriteValue(mu_F_sf); self.sf_mu_F_hi.setWriteValue(mu_F_sf_up); self.sf_mu_F_low.setWriteValue(mu_F_sf_down)
			self.sf_scale_comb.setWriteValue(scale_comb_sf); self.sf_scale_comb_hi.setWriteValue(scale_comb_sf_up); self.sf_scale_comb_low.setWriteValue(scale_comb_sf_down)
			self.sf_pdf_alphas.setWriteValue(pdf_alphas_sf); self.sf_pdf_alphas_hi.setWriteValue(pdf_alphas_sf_up); self.sf_pdf_alphas_low.setWriteValue(pdf_alphas_sf_down)
			nominal_top_pt_rw_value = self.toppTreweight.getReadValue()/self.toppTreweightmean
			self.toppTreweight.setWriteValue(nominal_top_pt_rw_value)
			self.sf_top_pt_rw.setWriteValue(1.)
			self.sf_top_pt_rw_hi.setWriteValue(nominal_top_pt_rw_value)
			self.sf_top_pt_rw_low.setWriteValue(1.+(1.-nominal_top_pt_rw_value))

	##################################  selection cut helper functions  ##################################

	def docut(self,cutname,cutbool,cutdict=None) :
		if cutdict==None :
			cutdict=self.cut_branches
		cutdict[cutname].setWriteValue(1) if cutbool else cutdict[cutname].setWriteValue(0)

	def __assignFullSelectionCutVars__(self,canreconstruct,topology,nLbtags,nMbtags,fitchi2,scaledlept,lep,electrons,muons,ak4jets,ak8jets,met) :
		#number of AK4 jets
		self.docut('ak4jetmult',canreconstruct)
		#good primary vertex
		self.docut('goodpv',(self.vtxBranches['vtx_ndof'].getReadValue()>=4 and abs(self.vtxBranches['vtx_z'].getReadValue())<24. and abs(self.vtxBranches['vtx_rho'].getReadValue())<2.))
		#met filtering
		metfiltercuts = []
		for branch in self.filterBranches.values() :
			if self.is_data==0 and branch.getReadName()=='Flag_eeBadScFilter' : #the one filter that shouldn't be applied to Monte Carlo
				continue
			if branch.getWriteValue()!=1 :
				metfiltercuts.append(False)
		self.docut('metfilters',metfiltercuts.count(False)==0)
		#number of btags
		self.docut('btags',((topology<3 and nLbtags>0) or (topology==3 and nMbtags>1)))
		#kinematic fit chi2
		self.docut('kinfitchi2',(canreconstruct and (topology==3 or fitchi2<-15)))
		#reconstructed leptonic top mass for boosted events
		self.docut('recoleptM',(canreconstruct and (topology==3 or scaledlept.M()<210.)))
		#other cuts are lepton flavor specific
		other_leps = []; allcuts = []
		if self.lepflavor.getWriteValue()==1 :
			#trigger
			trigvals = []
			if topology==1 or topology==2 :
				for trigName in MU_TRIG_PATHS_BOOSTED :
					trigvals.append(self.triggerBranches[trigName].getWriteValue())
			elif topology==3 :
				for trigName in MU_TRIG_PATHS_RESOLVED :
					trigvals.append(self.triggerBranches[trigName].getWriteValue())
			self.docut('trigger',trigvals.count(1)>0)
			#isolated lepton
			self.docut('isolepton',(lep.is2DIso(topology) and (topology<3 or lep.isTightIso())))
			#exactly one isolated lepton
			other_leps+=electrons; other_leps+=muons[1:]; noOtherLeps = True
			for olep in other_leps :
				if olep.is2DIso(topology) and (topology<3 or olep.isLooseIso()) :
					noOtherLeps = False
					break
			self.docut('onelepton',noOtherLeps)
			#leading ak4 jets
			self.docut('jetcuts',(topology==3 or ((len(ak4jets)>1 and ak4jets[0].getPt()>150. and ak4jets[1].getPt()>50.) and (topology==1 or ak8jets[0].getSDM()>40.))))
			#MET cuts
			self.docut('METcuts',((topology==3 and met.E()>40.) or met.E()>50.))
			#boosted lepton pT cuts
			self.docut('lepcuts',(topology==3 or lep.getPt()>50.))
		elif self.lepflavor.getWriteValue()==2 :
			#trigger
			trigvals = []
			if topology==1 or topology==2 :
				for trigName in EL_TRIG_PATHS_BOOSTED :
					trigvals.append(self.triggerBranches[trigName].getWriteValue())
			elif topology==3 :
				for trigName in EL_TRIG_PATHS_RESOLVED :
					trigvals.append(self.triggerBranches[trigName].getWriteValue())
			self.docut('trigger',trigvals.count(1)>0)
			#isolated lepton
			self.docut('isolepton',(lep.is2DIso(topology) and (topology<3 or lep.isIso())))
			#exactly one isolated lepton
			other_leps+=muons; other_leps+=electrons[1:]; noOtherLeps = True
			for lep in other_leps :
				if (topology<3 and lep.is2DIso(topology)) or (topology==3 and lep.is2DIso(topology) and lep.isLooseIso()) :
					noOtherLeps = False
					break
			self.docut('onelepton',noOtherLeps)
			#leading ak4 jets
			self.docut('jetcuts',(topology==3 or ((len(ak4jets)>1 and ak4jets[0].getPt()>250. and ak4jets[1].getPt()>70.) and (topology==1 or ak8jets[0].getSDM()>40.))))
			#MET cuts
			self.docut('METcuts',((topology==3 and met.E()>40.) or met.E()>100.))
			#boosted lepton pT cuts
			self.docut('lepcuts',(topology==3 or lep.getPt()>50.))
		#full selection
		for cutbranch in self.cut_branches.values() :
			if cutbranch.getWriteValue()==0 :
				allcuts.append(False)
		self.docut('fullselection',allcuts.count(False)==0)

	def __assignWJetsCRSelectionCutVars__(self) :
		#W+Jets control region (fails kinematic fit chi2 cuts OR reconstructed leptonic top mass cuts) for boosted events
		wjets_cr_pass_cutlist = ['goodpv','metfilters','trigger','onelepton','isolepton','btags','ak4jetmult','jetcuts','METcuts','lepcuts','validminimization']
		wjets_cr_fail_cutlist = ['kinfitchi2','recoleptM']
		#npassedcuts = sum([int(self.cut_branches[cutname].getWriteValue()==1) for cutname in wjets_cr_pass_cutlist]) #DEBUG
		#nfailedcuts = sum([int(self.cut_branches[cutname].getWriteValue()==0) for cutname in wjets_cr_fail_cutlist]) #DEBUG
		#print 'WJets CR: passed %d/%d pass cuts; failed %d/%d fail cuts'%(npassedcuts,len(wjets_cr_pass_cutlist),nfailedcuts,len(wjets_cr_fail_cutlist)) #DEBUG
		for cutname in wjets_cr_pass_cutlist :
			if not self.cut_branches[cutname].getWriteValue()==1 :
				self.cr_sb_cut_branches['wjets_cr_selection'].setWriteValue(0)
				break
		if self.cr_sb_cut_branches['wjets_cr_selection'].getWriteValue()!=0 :
			for cutname in wjets_cr_fail_cutlist :
				if self.cut_branches[cutname].getWriteValue()==0 :
					self.cr_sb_cut_branches['wjets_cr_selection'].setWriteValue(1)
					break
		if self.cr_sb_cut_branches['wjets_cr_selection'].getWriteValue()==2 :
			self.cr_sb_cut_branches['wjets_cr_selection'].setWriteValue(0)
		#if npassedcuts==len(wjets_cr_pass_cutlist) and nfailedcuts>0 : #DEBUG
		#	print 'WJetsCR: selection = %d'%(self.cr_sb_cut_branches['wjets_cr_selection'].getWriteValue()) #DEBUG

	def __assignQCDSBSelectionCutVars__(self) :
		qcd_sb_pass_cutlist = ['goodpv','metfilters','trigger','onelepton','btags','ak4jetmult','jetcuts','lepcuts','validminimization']
		for cutname in qcd_sb_pass_cutlist :
			if not self.cut_branches[cutname].getWriteValue()==1 :
				self.cr_sb_cut_branches['qcd_A_SR_selection'].setWriteValue(0)
				self.cr_sb_cut_branches['qcd_B_SR_selection'].setWriteValue(0)
				self.cr_sb_cut_branches['qcd_C_SR_selection'].setWriteValue(0)
				self.cr_sb_cut_branches['qcd_A_CR_selection'].setWriteValue(0)
				self.cr_sb_cut_branches['qcd_B_CR_selection'].setWriteValue(0)
				self.cr_sb_cut_branches['qcd_C_CR_selection'].setWriteValue(0)
				break
		#separate QCD sideband regions in 2D space
		if self.cr_sb_cut_branches['qcd_A_SR_selection'].getWriteValue()!=0 :
			issignalregion = self.cut_branches['kinfitchi2'].getWriteValue()==1 and self.cut_branches['recoleptM'].getWriteValue()==1
			isolepval = self.cut_branches['isolepton'].getWriteValue()
			metcutsval = self.cut_branches['METcuts'].getWriteValue()
			if issignalregion :
				self.docut('qcd_A_SR_selection',isolepval==1 and metcutsval==0,self.cr_sb_cut_branches)
				self.docut('qcd_B_SR_selection',isolepval==0 and metcutsval==0,self.cr_sb_cut_branches)
				self.docut('qcd_C_SR_selection',isolepval==0 and metcutsval==1,self.cr_sb_cut_branches)
			else :
				self.docut('qcd_A_CR_selection',isolepval==1 and metcutsval==0,self.cr_sb_cut_branches)
				self.docut('qcd_B_CR_selection',isolepval==0 and metcutsval==0,self.cr_sb_cut_branches)
				self.docut('qcd_C_CR_selection',isolepval==0 and metcutsval==1,self.cr_sb_cut_branches)

	def __assignElectronTriggerCRSelectionCutVars__(self,muons,electrons,topology,met1,met2) :
		return
		##muon trigger
		#trigvals = []
		#if topology==1 or topology==2 :
		#	for trigName in MU_TRIG_PATHS_BOOSTED :
		#		trigvals.append(self.triggerBranches[trigName].getWriteValue())
		#elif topology==3 :
		#	for trigName in MU_TRIG_PATHS_RESOLVED :
		#		trigvals.append(self.triggerBranches[trigName].getWriteValue())
		#self.docut('eltrig_mutrigger',trigvals.count(1)>0,self.eltrig_cut_branches)
		##electron trigger
		#trigvals = []
		#if topology==1 or topology==2 :
		#	for trigName in EL_TRIG_PATHS_BOOSTED :
		#		trigvals.append(self.triggerBranches[trigName].getWriteValue())
		#elif topology==3 :
		#	for trigName in EL_TRIG_PATHS_RESOLVED :
		#		trigvals.append(self.triggerBranches[trigName].getWriteValue())
		#self.docut('eltrig_eltrigger',trigvals.count(1)>0,self.eltrig_cut_branches)
		##find the hardest isolated (if possible) muon and electron and set the "two leptons" variable
		#allisoleps = []
		#themuon = None; theelectron = None
		#if len(muons)>0 :
		#	for m in muons :
		#		if m.is2DIso(topology) and (topology<3 or m.isLooseIso()) :
		#			allisoleps.append(m)
		#			if themuon==None :
		#				themuon = m
		#				break
		#	if themuon==None :
		#		themuon = muons[0]
		#if len(electrons)>0 :
		#	for e in electrons :
		#		if e.is2DIso(topology) and (topology<3 or e.isIso()) :
		#			allisoleps.append(e)
		#			if theelectron==None :
		#				theelectron = e
		#				break
		#	if theelectron==None :
		#		theelectron = electrons[0]
		#if themuon!=None :
		#	self.eltrig_themuon_pt.setWriteValue(themuon.getPt())
		#	self.eltrig_themuon_eta.setWriteValue(themuon.getEta())
		#if theelectron!=None :
		#	self.eltrig_theelectron_pt.setWriteValue(theelectron.getPt())
		#	self.eltrig_theelectron_eta.setWriteValue(theelectron.getEta())
		#self.docut('eltrig_twoleptons',themuon!=None and theelectron!=None and len(allisoleps)==2,self.eltrig_cut_branches)
		##isolated muon
		#self.docut('eltrig_isomu',themuon!=None and themuon.is2DIso(topology) and (topology<3 or themuon.isTightIso()),self.eltrig_cut_branches)
		##isolated electron
		#self.docut('eltrig_isoel',theelectron!=None and theelectron.is2DIso(topology) and (topology<3 or theelectron.isIso()),self.eltrig_cut_branches)
		##opposite lepton charge
		#self.docut('eltrig_opplepcharge',themuon!=None and theelectron!=None and themuon.getQ()*theelectron.getQ()<0,self.eltrig_cut_branches)
		##number of btags
		#self.eltrig_cut_branches['eltrig_btags'].setWriteValue(self.cut_branches['btags'].getWriteValue())
		##leptonic-side W pT
		#themuonvec = None if themuon==None else themuon.getFourVector()
		#theelectronvec = None if theelectron==None else theelectron.getFourVector()
		#self.docut('eltrig_lepWpT',(themuon!=None and ((themuonvec+met1).Pt()>50. or (themuonvec+met2).Pt()>50.)) or (theelectronvec!=None and ((theelectronvec+met1).Pt()>50. or (theelectronvec+met2).Pt()>50.)),self.eltrig_cut_branches)
		##full selection
		#self.eltrig_cut_branches['eltrig_fullselection'].setWriteValue(1)
		#for cutkey, cutbranch in self.eltrig_cut_branches.iteritems() :
		#	if cutbranch.getWriteValue()==0 :
		#		if cutkey!='eltrig_eltrigger' :
		#			self.eltrig_cut_branches['eltrig_fullselection'].setWriteValue(0)

	def __assignMiniIsolationFullSelectionCutVars__(self,lep,muons,electrons) :
		otherpasscuts = ['goodpv','metfilters','trigger','btags','ak4jetmult','jetcuts','METcuts','lepcuts','kinfitchi2','recoleptM','validminimization','fullselection']
		#start with lepton flavor-specific cuts
		other_leps = []
		if self.lepflavor.getWriteValue()==1 :
			#lepton miniisolation
			self.docut('miniisolepton',lep.isMedMiniIso(),self.alt_lep_iso_cut_branches)
			#other loose lepton veto
			other_leps+=electrons; other_leps+=muons[1:]; noOtherLeps = True
			for olep in other_leps :
				if olep.isLooseMiniIso() :
					noOtherLeps = False
					break
			self.docut('oneminiisolepton',noOtherLeps,self.alt_lep_iso_cut_branches)
		elif self.lepflavor.getWriteValue()==2 :
			#lepton miniisolation
			self.docut('miniisolepton',lep.isTightMiniIso(),self.alt_lep_iso_cut_branches)
			#other loose lepton veto
			other_leps+=electrons; other_leps+=muons[1:]; noOtherLeps = True
			for olep in other_leps :
				if olep.isLooseMiniIso() :
					noOtherLeps = False
					break
			self.docut('oneminiisolepton',noOtherLeps,self.alt_lep_iso_cut_branches)
		#full selection just uses what we already have
		for passcut in otherpasscuts :
			if self.cut_branches[passcut].getWriteValue()!=1 :
				self.alt_lep_iso_cut_branches['fullminiisoselection'].setWriteValue(0)
				break
		if self.alt_lep_iso_cut_branches['fullminiisoselection'].getWriteValue()!=0 and self.alt_lep_iso_cut_branches['miniisolepton'].getWriteValue()==1 and self.alt_lep_iso_cut_branches['oneminiisolepton'].getWriteValue()==1 :
			self.alt_lep_iso_cut_branches['fullminiisoselection'].setWriteValue(1)
		else :
			self.alt_lep_iso_cut_branches['fullminiisoselection'].setWriteValue(0)

	##################################  #__init__ function  ##################################
	def __init__(self,fileName,tree,isData,xsec,kfac,jes,jer,onGrid,pu_histo,totweight,renormdict) :
		#output file
		self.outfile_name = fileName
		self.outfile = TFile(self.outfile_name,'recreate')
		#output tree
		self.tree = TTree('tree','recreate')
		#Run era
		self.run_era = None
		fnsplit = fileName.split('Run2016')
		if len(fnsplit)>1 :
			self.run_era=fnsplit[1][0]
		#input tree
		self.inputTree = tree
		#Set input and output branch locations
		for branch in sorted(self.allBranches.values()) :
			branch.initialize(self.inputTree,self.tree)
		#data or MC?
		self.is_data = isData
		#cross section
		self.xsec = xsec
		#k factor
		self.kfac = kfac
		#JEC systematics?
		self.JES = jes
		self.JER = jer
		#are we running on the grid?
		self.onGrid = onGrid
		#Set the total weight
		self.totalweight = totweight
		#Set the corrector that does event weights and JEC calculations
		self.corrector = Corrector(self.is_data,onGrid,pu_histo,self.run_era,renormdict)
		#set the mean value of the top pT reweighting factor for this sample
		self.toppTreweightmean = renormdict['topptrwmean']

	##################################   reset function   ##################################
	#########  sets all relevant values back to initial values to get ready for next event  ##########
	def reset(self) :
		for branch in self.allBranches.values() :
			branch.reset()
		gc.collect()

	##############################   tree filling functions  ###############################		

	def __setFourVectorBranchValues__(self,name,vec) :
		self.allBranches[name+'_pt'].setWriteValue(vec.Pt())
		self.allBranches[name+'_eta'].setWriteValue(vec.Eta())
		self.allBranches[name+'_phi'].setWriteValue(vec.Phi())
		self.allBranches[name+'_M'].setWriteValue(vec.M())

	def __setLeptonBranchValues__(self,name,lep) :
		self.__setFourVectorBranchValues__(name,lep.getFourVector())
		self.allBranches[name+'_Q'].setWriteValue(int(lep.getQ()))
		self.allBranches[name+'_relPt'].setWriteValue(lep.getRelPt())
		self.allBranches[name+'_dR'].setWriteValue(lep.getDR())
		self.allBranches[name+'_Iso'].setWriteValue(lep.getIso())
		self.allBranches[name+'_MiniIso'].setWriteValue(lep.getMiniIso())

	########## function to close out the event, called before kicking back to runner #########
	def __closeout__(self) :
		#fill ttree
		#print 'written M = %.4f, written M_prefit = %.4f'%(self.M.getWriteValue(),self.M_prefit.getWriteValue()) #DEBUG
		self.outfile.cd()
		self.tree.Fill()

	################################## __del__ function  ###################################
	def __del__(self) :
		self.outfile.cd()
		self.outfile.Write()
		self.outfile.Close()

def getHypothesisList(topology,lep,met1_vec,met2_vec,nMETs,ak4jets,ttags) :
	#first figure out how many btags there are in this event
	nbtags = sum([x.isLbTagged() for x in ak4jets]) if topology<3 else sum([x.isMbTagged() for x in ak4jets])
	hypotheses = []
	met_options = [met1_vec,met2_vec]
	for i in range(nMETs) :
		#The leptonic side is the same for everything
		#there's only one lepton, so first choose the MET hypothesis
		thismet = met_options[i]
		#for any choice of leptonic-side b-jet
		for j in range(len(ak4jets)) :
			lepbCandJet = ak4jets[j]
			lepbistagged = lepbCandJet.isLbTagged() if topology<3 else lepbCandJet.isMbTagged()
			#the rest is topology-dependent
			#FULLY-MERGED EVENTS: 
			#the hadronic top candidate is the AK8 jet
			#just choose the leptonic b candidate from all AK4 jets
			#hypotheses are [lepton, neutrino, leptonic b-jet, [hadronic top jet]]
			if topology==1 : 
				#if there are btags the leptonic b candidate must be one of them
				if nbtags>0 and not lepbistagged :
					continue
				#check the reco leptonic side top
				thisleptop = lep.getFourVector()+thismet+lepbCandJet.getFourVector()
				#for all the choices of a merged top
				for ttag in ttags :
					#check the separation
					if thisleptop.DeltaR(ttag.getFourVector())>2. :
						#append this hypothesis
						hypotheses.append([lep,thismet,lepbCandJet,[ttag]])
				continue
			#other events need more than one AK4 jet
			#loop over the other ak4 jets to find the first hadronic-side AK4 jet
			for k in range(len(ak4jets)) :
				if j==k :
					continue
				had1CandJet = ak4jets[k]
				had1istagged = had1CandJet.isLbTagged() if topology<3 else had1CandJet.isMbTagged()
				hadsidehasbtag = had1istagged					
				#loop over the remaining AK4 jets for a second hadronic-side jet
				for m in range(k+1,len(ak4jets)) :
					if j==m :
						continue
					had2CandJet = ak4jets[m]
					had2istagged = had2CandJet.isLbTagged() if topology<3 else had2CandJet.isMbTagged()
					hadsidehasbtag = had1istagged or had2istagged
					#loop one last time over the remaining AK4 jets for a third and final hadronic-side jet
					for n in range(m+1,len(ak4jets)) :
						if j==n :
							continue
						had3CandJet = ak4jets[n]
						had3istagged = had3CandJet.isLbTagged() if topology<3 else had3CandJet.isMbTagged()
						hadsidehasbtag = had1istagged or had2istagged or had3istagged
						#BOOSTED UNTAGGED WITH THREE HADRONIC SIDE JETS, AND FULLY RESOLVED
						#the hadronic side of the decay has three jets 
						#hypotheses are [lepton, neutrino, leptonic b-jet, [list of three hadronic-side jets]]
						if topology==2 or topology==3 :
							#make sure the btags are in place
							if (nbtags==1 and not (hadsidehasbtag or lepbistagged)) or (nbtags==2 and not (hadsidehasbtag and lepbistagged)) :
								continue
							#append this hypothesis
							hypotheses.append([lep,thismet,lepbCandJet,[had1CandJet,had2CandJet,had3CandJet]])
	if hypotheses==[] : #DEBUG
		print 'found no hypotheses for a type-%d event with %d MET options, %d AK4 jets, %d btags!'%(topology,nMETs,len(ak4jets),nbtags) #DEBUG
	return hypotheses