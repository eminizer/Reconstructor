#Global variables
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0
#Trigger paths
MU_TRIG_PATH = 'HLT_Mu30_eta2p1_PFJet150_PFJet50'
EL_TRIG_PATH = 'HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50'

##########								   Imports  								##########

from math import pi, log
import copy
from ROOT import TFile, TTree, TLorentzVector
from branch import Branch
from eventTypeHelper import keepEventType, findInitialPartons, findMCParticles
from jet import AK4Jet, AK8Jet
from lepton import Muon, Electron
from metHelper import setupMET
from ttbarReconstructor import reconstruct
from angleReconstructor import getObservables, getMCObservables
from corrector import Corrector

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
	#GenWeight info
	genWeight = AddBranch('evt_Gen_Weight','genWeight','F',1.0,'1',[allBranches])
	rho 	  = AddBranch('evt_rho','rho','D',1.0,'1',[allBranches,ak4JetBranches,ak8JetBranches])
	npv 	  = AddBranch('evt_npv','npv','I',-1,'1',[allBranches,ak4JetBranches,ak8JetBranches])
	#MC GenEvent info
	mcGenEventBranches = {}
	thisdictlist = [allBranches,mcGenEventBranches]
	gen_size 	   = AddBranch(readname='gen_size',ttreetype='i',dictlist=thisdictlist)
	gen_ID 		   = AddBranch(readname='gen_ID',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Status 	   = AddBranch(readname='gen_Status',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Mom0ID 	   = AddBranch(readname='gen_Mom0ID',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Mom0Status = AddBranch(readname='gen_Mom0Status',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Mom1ID 	   = AddBranch(readname='gen_Mom1ID',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Mom1Status = AddBranch(readname='gen_Mom1Status',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Dau0ID 	   = AddBranch(readname='gen_Dau0ID',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Dau0Status = AddBranch(readname='gen_Dau0Status',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Dau1ID 	   = AddBranch(readname='gen_Dau1ID',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Dau1Status = AddBranch(readname='gen_Dau1Status',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Dau1Status = AddBranch(readname='gen_Dau1Status',ttreetype='I',size='gen_size',dictlist=thisdictlist)
	gen_Pt 		   = AddBranch(readname='gen_Pt',ttreetype='F',size='gen_size',dictlist=thisdictlist)
	gen_Eta 	   = AddBranch(readname='gen_Eta',ttreetype='F',size='gen_size',dictlist=thisdictlist)
	gen_Phi 	   = AddBranch(readname='gen_Phi',ttreetype='F',size='gen_size',dictlist=thisdictlist)
	gen_E 		   = AddBranch(readname='gen_E',ttreetype='F',size='gen_size',dictlist=thisdictlist)
	#Trigger Information
	triggerBranches = {}
	thisdictlist = [allBranches,triggerBranches]
	muTrig = AddBranch(MU_TRIG_PATH,'muTrig','I',-1,'1',thisdictlist)
	elTrig = AddBranch(EL_TRIG_PATH,'elTrig','I',-1,'1',thisdictlist)
	#MET Filter information
	filterBranches = {}
	thisdictlist = [allBranches,filterBranches]
	noise 				= AddBranch('Flag_HBHENoiseFilter','noisefilter','I',-1,'1',thisdictlist)
	noiseIso 			= AddBranch('Flag_HBHENoiseIsoFilter','noiseIsofilter','I',-1,'1',thisdictlist)
	ecaldeadcell 		= AddBranch('Flag_EcalDeadCellTriggerPrimitiveFilter','ecaldeadcellfilter','I',-1,'1',thisdictlist)
	goodVertices 		= AddBranch('Flag_goodVertices','goodVerticesfilter','I',-1,'1',thisdictlist)
	eebadsc 			= AddBranch('Flag_eeBadScFilter','eebadscfilter','I',-1,'1',thisdictlist)
	globaltighthalo2016 = AddBranch('Flag_globalTightHalo2016Filter','globaltighthalo2016filter','I',-1,'1',thisdictlist)
	#pileup
	pileupBranches = {}
	thisdictlist = [allBranches,pileupBranches]
	pu = AddBranch(readname='pu_NtrueInt',ttreetype='I',dictlist=thisdictlist)
	#MET
	metBranches = {}
	thisdictlist = [allBranches,metBranches]
	met_size 	= AddBranch(readname='met_size',ttreetype='i',dictlist=thisdictlist)
	met_pts 	= AddBranch(readname='met_Pt',size='met_size',dictlist=thisdictlist)
	met_phis 	= AddBranch(readname='met_Phi',size='met_size',dictlist=thisdictlist)
	#muons
	muonBranches = {}
	thisdictlist = [allBranches,muonBranches]
	mu_size 	  = AddBranch(readname='mu_size',ttreetype='i',dictlist=thisdictlist)
	mu_pts 		  = AddBranch(readname='mu_Pt',size='mu_size',dictlist=thisdictlist)
	mu_etas 	  = AddBranch(readname='mu_Eta',size='mu_size',dictlist=thisdictlist)
	mu_phis 	  = AddBranch(readname='mu_Phi',size='mu_size',dictlist=thisdictlist)
	mu_es 		  = AddBranch(readname='mu_E',size='mu_size',dictlist=thisdictlist)
	mu_charges 	  = AddBranch(readname='mu_Charge',size='mu_size',dictlist=thisdictlist)
	mus_isLoose   = AddBranch(readname='mu_IsLooseMuon',size='mu_size',dictlist=thisdictlist)
	mu_Keys 	  = AddBranch(readname='mu_Key',size='mu_size',dictlist=thisdictlist)
	#electrons
	electronBranches = {}
	thisdictlist = [allBranches,electronBranches]
	el_size 	  = AddBranch(readname='el_size',ttreetype='i',dictlist=thisdictlist)
	el_pts 		  = AddBranch(readname='el_Pt',size='el_size',dictlist=thisdictlist)
	el_etas 	  = AddBranch(readname='el_Eta',size='el_size',dictlist=thisdictlist)
	el_scetas 	  = AddBranch(readname='el_SCEta',size='el_size',dictlist=thisdictlist)
	el_phis 	  = AddBranch(readname='el_Phi',size='el_size',dictlist=thisdictlist)
	el_es 		  = AddBranch(readname='el_E',size='el_size',dictlist=thisdictlist)
	el_charges 	  = AddBranch(readname='el_Charge',size='el_size',dictlist=thisdictlist)
	el_idLoose 	  = AddBranch(readname='el_IDTight_NoIso',size='el_size',dictlist=thisdictlist)
	el_Keys 	  = AddBranch(readname='el_Key',size='el_size',dictlist=thisdictlist)
	#AK4 Jets
	thisdictlist = [allBranches,ak4JetBranches]
	ak4_size 					= AddBranch(readname='jetAK4_size',ttreetype='i',dictlist=thisdictlist)
	ak4_pts 					= AddBranch(readname='jetAK4_Pt',size='jetAK4_size',dictlist=thisdictlist)
	ak4_etas 					= AddBranch(readname='jetAK4_Eta',size='jetAK4_size',dictlist=thisdictlist)
	ak4_phis 					= AddBranch(readname='jetAK4_Phi',size='jetAK4_size',dictlist=thisdictlist)
	ak4_es 						= AddBranch(readname='jetAK4_E',size='jetAK4_size',dictlist=thisdictlist)
	ak4_genpts 					= AddBranch(readname='jetAK4_GenJetPt',size='jetAK4_size',dictlist=thisdictlist)
	ak4_genetas 				= AddBranch(readname='jetAK4_GenJetEta',size='jetAK4_size',dictlist=thisdictlist)
	ak4_genphis 				= AddBranch(readname='jetAK4_GenJetPhi',size='jetAK4_size',dictlist=thisdictlist)
	ak4_csvv2s 					= AddBranch(readname='jetAK4_CSVv2',size='jetAK4_size',dictlist=thisdictlist)
	ak4_jec0s 					= AddBranch(readname='jetAK4_jecFactor0',size='jetAK4_size',dictlist=thisdictlist)
	ak4_jecuncs 				= AddBranch(readname='jetAK4_jecUncertainty',size='jetAK4_size',dictlist=thisdictlist)
	ak4_jetAs 					= AddBranch(readname='jetAK4_jetArea',size='jetAK4_size',dictlist=thisdictlist)
	ak4_keylists 				= AddBranch(readname='jetAK4_Keys',ttreetype='vi',size='jetAK4_size',dictlist=thisdictlist)
	ak4_neutralMultiplicity 	= AddBranch(readname='jetAK4_neutralMultiplicity',size='jetAK4_size',dictlist=thisdictlist)
	ak4_neutralHadronEnergyFrac = AddBranch(readname='jetAK4_neutralHadronEnergyFrac',size='jetAK4_size',dictlist=thisdictlist)
	ak4_neutralEmEnergyFrac 	= AddBranch(readname='jetAK4_neutralEmEnergyFrac',size='jetAK4_size',dictlist=thisdictlist)
	ak4_chargedHadronEnergyFrac = AddBranch(readname='jetAK4_chargedHadronEnergyFrac',size='jetAK4_size',dictlist=thisdictlist)
	ak4_chargedEmEnergyFrac 	= AddBranch(readname='jetAK4_chargedEmEnergyFrac',size='jetAK4_size',dictlist=thisdictlist)
	ak4_chargedMultiplicity 	= AddBranch(readname='jetAK4_chargedMultiplicity',size='jetAK4_size',dictlist=thisdictlist)
	#AK8 Jets
	thisdictlist = [allBranches,ak8JetBranches]
	ak8_size 					= AddBranch(readname='jetAK8_size',ttreetype='i',dictlist=thisdictlist)
	ak8_pts 					= AddBranch(readname='jetAK8_Pt',size='jetAK8_size',dictlist=thisdictlist)
	ak8_etas 					= AddBranch(readname='jetAK8_Eta',size='jetAK8_size',dictlist=thisdictlist)
	ak8_phis 					= AddBranch(readname='jetAK8_Phi',size='jetAK8_size',dictlist=thisdictlist)
	ak8_es 						= AddBranch(readname='jetAK8_E',size='jetAK8_size',dictlist=thisdictlist)
	ak8_genpts 					= AddBranch(readname='jetAK8_GenJetPt',size='jetAK8_size',dictlist=thisdictlist)
	ak8_genetas 				= AddBranch(readname='jetAK8_GenJetEta',size='jetAK8_size',dictlist=thisdictlist)
	ak8_genphis 				= AddBranch(readname='jetAK8_GenJetPhi',size='jetAK8_size',dictlist=thisdictlist)
	ak8_csvv2s 					= AddBranch(readname='jetAK8_CSVv2',size='jetAK8_size',dictlist=thisdictlist)
	ak8_jec0s 					= AddBranch(readname='jetAK8_jecFactor0',size='jetAK8_size',dictlist=thisdictlist)
	ak8_jecuncs 				= AddBranch(readname='jetAK8_jecUncertainty',size='jetAK8_size',dictlist=thisdictlist)
	ak8_jetAs 					= AddBranch(readname='jetAK8_jetArea',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau1s 					= AddBranch(readname='jetAK8_tau1',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau2s 					= AddBranch(readname='jetAK8_tau2',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau3s 					= AddBranch(readname='jetAK8_tau3',size='jetAK8_size',dictlist=thisdictlist)
	ak8_sdms 					= AddBranch(readname='jetAK8_softDropMass',size='jetAK8_size',dictlist=thisdictlist)
	ak8_keylists 				= AddBranch(readname='jetAK8_Keys',ttreetype='vi',size='jetAK8_size',dictlist=thisdictlist)
	ak8_neutralMultiplicity 	= AddBranch(readname='jetAK8_neutralMultiplicity',size='jetAK8_size',dictlist=thisdictlist)
	ak8_neutralHadronEnergyFrac = AddBranch(readname='jetAK8_neutralHadronEnergyFrac',size='jetAK8_size',dictlist=thisdictlist)
	ak8_neutralEmEnergyFrac 	= AddBranch(readname='jetAK8_neutralEmEnergyFrac',size='jetAK8_size',dictlist=thisdictlist)
	ak8_chargedHadronEnergyFrac = AddBranch(readname='jetAK8_chargedHadronEnergyFrac',size='jetAK8_size',dictlist=thisdictlist)
	ak8_chargedEmEnergyFrac 	= AddBranch(readname='jetAK8_chargedEmEnergyFrac',size='jetAK8_size',dictlist=thisdictlist)
	ak8_chargedMultiplicity 	= AddBranch(readname='jetAK8_chargedMultiplicity',size='jetAK8_size',dictlist=thisdictlist)
	ak8_subjetIndex0 			= AddBranch(readname='jetAK8_vSubjetIndex0',size='jetAK8_size',dictlist=thisdictlist)
	ak8_subjetIndex1 			= AddBranch(readname='jetAK8_vSubjetIndex1',size='jetAK8_size',dictlist=thisdictlist)
	#AK8 Subjets
	thisdictlist = [allBranches,ak8JetBranches]
	subjak8_size 					= AddBranch(readname='subjetAK8_size',ttreetype='i',dictlist=thisdictlist)
	subjak8_pts 					= AddBranch(readname='subjetAK8_Pt',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_etas 					= AddBranch(readname='subjetAK8_Eta',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_phis 					= AddBranch(readname='subjetAK8_Phi',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_es 						= AddBranch(readname='subjetAK8_E',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_genpts 					= AddBranch(readname='subjetAK8_GenJetPt',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_genetas 				= AddBranch(readname='subjetAK8_GenJetEta',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_genphis 				= AddBranch(readname='subjetAK8_GenJetPhi',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_csvv2s 					= AddBranch(readname='subjetAK8_CSVv2',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_jec0s 					= AddBranch(readname='subjetAK8_jecFactor0',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_keylists 				= AddBranch(readname='subjetAK8_Keys',ttreetype='vi',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_neutralMultiplicity 	= AddBranch(readname='subjetAK8_neutralMultiplicity',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_neutralHadronEnergyFrac = AddBranch(readname='subjetAK8_neutralHadronEnergyFrac',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_neutralEmEnergyFrac 	= AddBranch(readname='subjetAK8_neutralEmEnergyFrac',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_chargedHadronEnergyFrac = AddBranch(readname='subjetAK8_chargedHadronEnergyFrac',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_chargedEmEnergyFrac 	= AddBranch(readname='subjetAK8_chargedEmEnergyFrac',size='subjetAK8_size',dictlist=thisdictlist)
	subjak8_chargedMultiplicity 	= AddBranch(readname='subjetAK8_chargedMultiplicity',size='subjetAK8_size',dictlist=thisdictlist)
	#PRODUCES
	#Weights
	weightBranches = {}
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
	scalefactorBranches = {}
	thisdictlist = [allBranches,scalefactorBranches]
	sf_top_pT 		  = AddBranch(writename='sf_top_pT',inival=1.,dictlist=thisdictlist)
	sf_top_pT_low 	  = AddBranch(writename='sf_top_pT_low',inival=1.,dictlist=thisdictlist)
	sf_top_pT_hi 	  = AddBranch(writename='sf_top_pT_hi',inival=1.,dictlist=thisdictlist)
	sf_btag_eff 	  = AddBranch(writename='sf_btag_eff',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_low   = AddBranch(writename='sf_btag_eff_low',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_hi 	  = AddBranch(writename='sf_btag_eff_hi',inival=1.,dictlist=thisdictlist)
	sf_pileup 		  = AddBranch(writename='sf_pileup',inival=1.,dictlist=thisdictlist)
	sf_pileup_low 	  = AddBranch(writename='sf_pileup_low',inival=1.,dictlist=thisdictlist)
	sf_pileup_hi 	  = AddBranch(writename='sf_pileup_hi',inival=1.,dictlist=thisdictlist)
	sf_lep_ID 		  = AddBranch(writename='sf_lep_ID',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_low 	  = AddBranch(writename='sf_lep_ID_low',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_hi 	  = AddBranch(writename='sf_lep_ID_hi',inival=1.,dictlist=thisdictlist)
	sf_trig_eff 	  = AddBranch(writename='sf_trig_eff',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_low   = AddBranch(writename='sf_trig_eff_low',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_hi 	  = AddBranch(writename='sf_trig_eff_hi',inival=1.,dictlist=thisdictlist)
	#physics objects
	physobjectBranches = {}
	thisdictlist = [allBranches,physobjectBranches]
	#fourvectors
	fourvectornames = ['muon1','muon2','ele1','ele2','met','ak41','ak42','ak8']
	for fourvecname in fourvectornames :
		AddBranch(writename=fourvecname+'_pt',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_eta',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_phi',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_M',dictlist=thisdictlist)
	#rescaled fourvectors
	fourvectornames = ['lep','met','lepW','lepb','lept','hadW','hadb','hadt']
	for fourvecname in fourvectornames :
		AddBranch(writename='scaled_'+fourvecname+'_pt',dictlist=thisdictlist)
		AddBranch(writename='scaled_'+fourvecname+'_eta',dictlist=thisdictlist)
		AddBranch(writename='scaled_'+fourvecname+'_phi',dictlist=thisdictlist)
		AddBranch(writename='scaled_'+fourvecname+'_M',dictlist=thisdictlist)
	#others
	Q_l 		  = AddBranch(writename='Q_l',ttreetype='I',inival=0,dictlist=thisdictlist)
	leptonnames = ['muon1','muon2','ele1','ele2']
	for lepname in leptonnames :
		AddBranch(writename=lepname+'_Q',ttreetype='I',inival=0,dictlist=thisdictlist)
		AddBranch(writename=lepname+'_ID',ttreetype='i',inival=2,dictlist=thisdictlist)
		AddBranch(writename=lepname+'_relPt',dictlist=thisdictlist)
		AddBranch(writename=lepname+'_dR',dictlist=thisdictlist)
	hadt_tau32 = AddBranch(writename='hadt_tau32',dictlist=thisdictlist)
	hadt_tau21 = AddBranch(writename='hadt_tau21',dictlist=thisdictlist)
	hadt_SDM   = AddBranch(writename='hadt_SDM',dictlist=thisdictlist)
	hadW_tau32 = AddBranch(writename='hadW_tau32',dictlist=thisdictlist)
	hadW_tau21 = AddBranch(writename='hadW_tau21',dictlist=thisdictlist)
	hadW_SDM   = AddBranch(writename='hadW_SDM',dictlist=thisdictlist)
	#kinematic fit stuff
	thisdictlist = [allBranches]
	ttype = AddBranch(writename='ttype',ttreetype='i',inival=0,dictlist=thisdictlist)
	chi2 = AddBranch(writename='chi2',inival=0.0,dictlist=thisdictlist)
	nFits = AddBranch(writename='nFits',ttreetype='i',inival=0,dictlist=thisdictlist)
	#whether or not this event should be added twice and have its weight halved based on its initial state
	addTwice = AddBranch(writename='addTwice',ttreetype='i',inival=0,dictlist=thisdictlist)
	#Obervables
	observableBranches = {}
	thisdictlist = [allBranches,observableBranches]
	#cosine(theta)
	cstar 		 = AddBranch(writename='cstar',dictlist=thisdictlist)
	cstar_scaled = AddBranch(writename='cstar_scaled',dictlist=thisdictlist)
	cstar_corprefit = AddBranch(writename='cstar_corprefit',dictlist=thisdictlist)
	#Feynman x
	x_F 	   = AddBranch(writename='x_F',dictlist=thisdictlist)
	x_F_scaled = AddBranch(writename='x_F_scaled',dictlist=thisdictlist)
	x_F_corprefit = AddBranch(writename='x_F_corprefit',dictlist=thisdictlist)
	#ttbar invariant mass
	M 		 = AddBranch(writename='M',dictlist=thisdictlist)
	M_scaled = AddBranch(writename='M_scaled',dictlist=thisdictlist)
	M_corprefit = AddBranch(writename='M_corprefit',dictlist=thisdictlist)
	#initial quark vector
	mctruthBranches = {}
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
	#cut variables
	cut_branches = {}
	thisdictlist = [allBranches,cut_branches]
	cutnames = ['metfilters','onelepton','isolepton','smalljets','fullselection']
	for cutname in cutnames :
		AddBranch(writename=cutname,ttreetype='i',inival=2,dictlist=thisdictlist)
	#debugging variables
	kinfit_debug_branches = {}
	thisdictlist=[kinfit_debug_branches,allBranches]
	par_0 = AddBranch(writename='par_0',dictlist=thisdictlist)
	par_1 = AddBranch(writename='par_1',dictlist=thisdictlist)
	par_2 = AddBranch(writename='par_2',dictlist=thisdictlist)
	par_3 = AddBranch(writename='par_3',dictlist=thisdictlist)
	par_4 = AddBranch(writename='par_4',dictlist=thisdictlist)
	par_5 = AddBranch(writename='par_5',dictlist=thisdictlist)
	lnLtl = AddBranch(writename='lnLtl',dictlist=thisdictlist)
	lnLth = AddBranch(writename='lnLth',dictlist=thisdictlist)
	lnLwl = AddBranch(writename='lnLwl',dictlist=thisdictlist)
	lnLwh = AddBranch(writename='lnLwh',dictlist=thisdictlist)
	lnLlsf = AddBranch(writename='lnLlsf',dictlist=thisdictlist)
	lnLblsf = AddBranch(writename='lnLblsf',dictlist=thisdictlist)
	lnLhsf1 = AddBranch(writename='lnLhsf1',dictlist=thisdictlist)
	lnLhsf2 = AddBranch(writename='lnLhsf2',dictlist=thisdictlist)
	lnLhsf3 = AddBranch(writename='lnLhsf3',dictlist=thisdictlist)
	nhypotheses = AddBranch(writename='nhypotheses',ttreetype='i',inival=0,dictlist=thisdictlist)
	ismatchable = AddBranch(writename='ismatchable',ttreetype='i',inival=2,dictlist=thisdictlist)
	iscorrect = AddBranch(writename='iscorrect',ttreetype='i',inival=2,dictlist=thisdictlist)
	ismatchedpostfit = AddBranch(writename='ismatchedpostfit',ttreetype='i',inival=2,dictlist=thisdictlist)
	lepWcorprefitM = AddBranch(writename='lepWcorprefitM',dictlist=thisdictlist)
	hadWcorprefitM = AddBranch(writename='hadWcorprefitM',dictlist=thisdictlist)
	leptcorprefitM = AddBranch(writename='leptcorprefitM',dictlist=thisdictlist)
	hadtcorprefitM = AddBranch(writename='hadtcorprefitM',dictlist=thisdictlist)

	##################################  ANALYZE FUNCTION  ##################################
	def analyze(self,eventnumber) :
		#get the event in the tree
		self.inputTree.GetEntry(eventnumber)

		#Make sure the event has some met, at least one lepton, and at least one of each type of jet
		metsize = self.met_size.getReadValue()
		musize = self.mu_size.getReadValue()
		elsize = self.el_size.getReadValue()
		ak4jetsize = self.ak4_size.getReadValue()
		ak8jetsize = self.ak8_size.getReadValue()
		if not (metsize>0 and (musize>0 or elsize>0) and ak4jetsize>0 and ak8jetsize>0) :
			#print 'EVENT NUMBER %d NOT VALID; MISSING REQUISITE PHYSICS OBJECTS (metsize = %d, musize = %d, elsize = %d, ak4jetsize = %d, ak8jetsize = %d)'%(eventnumber,
			#	metsize,musize,elsize,ak4jetsize,ak8jetsize)
			return

		#MC stuff
		if not self.is_data :
			#event type split
			keepevent,addtwice = keepEventType(self.mcGenEventBranches,self.event_type,self.MC_generator)
			if not keepevent :
				return
			#set the addTwice value
			self.addTwice.setWriteValue(addtwice)
			#set the eventweight
			eventweight = self.genWeight.getReadValue()*(self.xsec/self.totalweight)
			self.weight.setWriteValue(eventweight)
			#Mother particle and MC truth top assignment
			vecnames = []; vecobjs = []
			if self.event_type!=4 :
				q_vec, qbar_vec = findInitialPartons(self.mcGenEventBranches,self.MC_generator)
				vecnames += ['q','qbar']
				vecobjs  += [q_vec,qbar_vec]
			if self.event_type<2 :
				MCt_vec, MCtbar_vec, MClep_vec, MClep_charge, MCv_vec, MClepb_vec, MChadW_vec, MChadb_vec = findMCParticles(self.mcGenEventBranches,self.MC_generator)
				vecnames += ['MCt','MCtbar','MClep','MCv','MClepb','MChadW','MChadb']
				vecobjs  += [MCt_vec,MCtbar_vec,MClep_vec,MCv_vec,MClepb_vec,MChadW_vec,MChadb_vec]
			#write out fourvectors of MC particles
			for i in range(len(vecnames)) :
				self.__setFourVectorBranchValues__(vecnames[i],vecobjs[i])

		#For the record, trigger information is handled automatically

		#MET
		met = TLorentzVector(); met.SetPtEtaPhiM(self.met_pts.getReadValue(), 0., self.met_phis.getReadValue(), 0.)
		self.__setFourVectorBranchValues__('met',met)

		#muons
		muons = []
		for i in range(musize) :
			newmuon=Muon(self.muonBranches,i)
			if newmuon.getPt()>50. and abs(newmuon.getEta())<2.1 and newmuon.getID()==1 : muons.append(newmuon)
		muons.sort(key=lambda x: x.getPt(), reverse=True)

		#electrons
		electrons = []
		for i in range(elsize) :
			newele = Electron(self.electronBranches,i)
			if newele.getPt()>40 and abs(newele.getEtaSC())<2.5 and newele.getID()==1 : electrons.append(newele)
		electrons.sort(key=lambda x: x.getPt(), reverse=True)

		#figure out whether the event is muonic or electronic, assign lep
		lep = None
		if len(muons)>0 and (len(electrons)==0 or muons[0].getPt()>electrons[0].getPt()) :
			lep = muons[0]; leptype='mu'
		elif len(electrons)>0 and (len(muons)==0 or electrons[0].getPt()>muons[0].getPt()) :
			lep = electrons[0]; leptype='el'
		if lep==None :
			#print 'EVENT NUMBER %d NOT VALID; NO LEPTONS'%(eventnumber) #DEBUG
			return
		self.Q_l.setWriteValue(int(lep.getQ()))

		#jets
		ak4jets = []; ak8jets = [];
		for i in range(ak4jetsize) :
			newJet = AK4Jet(self.ak4JetBranches,i,self.JES,self.JER,lep,self.corrector,self.is_data)
			if newJet.getFourVector()!=None and newJet.getPt()>15. and abs(newJet.getEta())<3.0 and newJet.isIDed() :
				ak4jets.append(newJet)
		for i in range(ak8jetsize) :
			newJet = AK8Jet(self.ak8JetBranches,i,self.JES,self.JER,lep,self.corrector,self.is_data)
			if newJet.getFourVector()!=None and newJet.getPt()>500. and abs(newJet.getEta())<2.4 and newJet.isIDed() :
				ak8jets.append(newJet)
		ak4jets.sort(key=lambda x: x.getPt(), reverse=True)
		ak8jets.sort(key=lambda x: x.getPt(), reverse=True)
		#if the lepton cleaning got rid of too many jets, or if the AK8 jet doesn't have at least 2 subjets, toss the event
		if not (len(ak4jets)>0 and len(ak8jets)>0 and ak8jets[0].getNSubjets()>1) :
			#s = 'EVENT NUMBER %d NOT VALID; MISSING JETS (# AK4jets = %d, # AK8Jets = %d'%(eventnumber,len(ak4jets),len(ak8jets)) #DEBUG
			#if (len(ak8jets)>0) : #DEBUG
			#	s+= ', # AK8 subjets = %d'%(ak8jets[0].getNSubjets()) #DEBUG
			#print s+')' #DEBUG
			return

		#calculate lepton isolation variables
		for muon in muons :
			muon.calculateIsolation(ak4jets)
		for electron in electrons :
			electron.calculateIsolation(ak4jets)

		#remove the ak4 jets we needed just for isolation calculations
		i=0
		while len(ak4jets)>0 and i<len(ak4jets) :
			ak4jet = ak4jets[i]
			if not (ak4jet.getPt()>30. and abs(ak4jet.getEta())<2.4) :
				ak4jets.pop(i)
			else :
				i+=1
		if not len(ak4jets)>1 :
			#print 'EVENT NUMBER %d NOT VALID; MISSING JETS (# AK4jets = %d)'%(eventnumber,len(ak4jets)) #DEBUG
			return

		#Set physics object fourvectors
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
		ak4len = len(ak4jets)
		if ak4len>0 :
			self.__setFourVectorBranchValues__('ak41',ak4jets[0].getFourVector())
			if ak4len>1 :
				self.__setFourVectorBranchValues__('ak42',ak4jets[1].getFourVector())
		self.__setFourVectorBranchValues__('ak8',ak8jets[0].getFourVector())
		self.hadt_tau32.setWriteValue(ak8jets[0].getTau32())
		self.hadt_tau21.setWriteValue(ak8jets[0].getTau21())
		self.hadt_SDM.setWriteValue(ak8jets[0].getSDM())

		#neutrino handling and setup for fit
		met1_vec, met2_vec = setupMET(lep.getFourVector(),met)
		self.nFits.setWriteValue(2)
		if met1_vec.Pz() == met2_vec.Pz() :
			self.nFits.setWriteValue(1)

		#--------------------------------------------------------------------Below here is event selection--------------------------------------------------------------------#
		['metfilters','onelepton','isolepton','smalljets','fullselection']
		#met filtering
		metfiltercuts = []
		for branch in self.filterBranches.values() :
			if branch.getWriteValue()!=1 :
				metfiltercuts.append(False)
		if metfiltercuts.count(False)==0 : self.cut_branches['metfilters'].setWriteValue(1)
		else : self.cut_branches['metfilters'].setWriteValue(0)
		#exactly one lepton
		other_leps = []
		if leptype=='mu' :
			other_leps+=electrons
			for i in range(1,len(muons)) :
				other_leps.append(muons[i])
		elif leptype=='el' :
			other_leps+=muons
			for i in range(1,len(electrons)) :
				other_leps.append(electrons[i])
		if len(other_leps)==0 : self.cut_branches['onelepton'].setWriteValue(1)
		else : self.cut_branches['onelepton'].setWriteValue(0)
		#isolated lepton
		if lep.getRelPt()>20. or lep.getDR()>0.4 : self.cut_branches['isolepton'].setWriteValue(1)
		else : self.cut_branches['isolepton'].setWriteValue(0)
		#two leading ak4 jets
		if leptype=='mu' :
			if ak4jets[0].getPt()>150. and abs(ak4jets[0].getEta())<2.4 and ak4jets[1].getPt()>50. and abs(ak4jets[1].getEta())<2.4 : self.cut_branches['smalljets'].setWriteValue(1)
			else : self.cut_branches['smalljets'].setWriteValue(0)
		elif leptype=='el' :
			if ak4jets[0].getPt()>250. and abs(ak4jets[0].getEta())<2.4 and ak4jets[1].getPt()>70. and abs(ak4jets[1].getEta())<2.4 : self.cut_branches['smalljets'].setWriteValue(1)
			else : self.cut_branches['smalljets'].setWriteValue(0)
		#full selection
		allcuts = []
		for cutbranch in self.cut_branches.values() :
			if cutbranch.getWriteValue()==0 :
				allcuts.append(False)
		if leptype=='mu' :
			allcuts.append((lep.getPt()+met.E())>150.)
			allcuts.append(met.E()>50.)
		elif leptype=='el' :
			allcuts.append(met.E()>120.)
		if allcuts.count(False)==0 : self.cut_branches['fullselection'].setWriteValue(1)
		else : self.cut_branches['fullselection'].setWriteValue(0)
		#-----------------------------------------------------------------Below here is event reconstruction-----------------------------------------------------------------#
		
		#list of jet assignment hypotheses
		hypotheses = []
		met_options = [met1_vec,met2_vec]
		#figure out whether the event is a type1 (fully merged) or type2 (partially merged) event
		if ak8jets[0].isTopTagged() : #if the event has a top tag, it's fully merged
			self.ttype.setWriteValue(1)
			#top is the hardest ak8 jet
			hadtCandJet = ak8jets[0]			
			for i in range(self.nFits.getWriteValue()) :
				thismet = met_options[i]
				#leptonic b is an ak4 jet on the opposite hemisphere
				for j in range(len(ak4jets)) :
					ak4jet = ak4jets[j]
					dRcheck = (ak4jet.getFourVector()+lep.getFourVector()+thismet).DeltaR(hadtCandJet.getFourVector())
					if dRcheck > pi/2. : 
						hypotheses.append([lep,thismet,ak4jet,hadtCandJet])
		else : #otherwise it's partially merged
			self.ttype.setWriteValue(2)
			#hadronic W is the hardest ak8 jet
			hadWCandJet = ak8jets[0]
			#make the list of hypotheses
			hypotheses = []; met_options = [met1_vec,met2_vec]
			for i in range(self.nFits.getWriteValue()) :
				thismet = met_options[i]
				#hadronic/leptonic b's are the hardest ak4 jets on the appropriate hemisphere
				hadbcands = []; lepbcands = []
				for i in range(len(ak4jets)) : 
					lepbak4jet = ak4jets[i]
					for j in range(len(ak4jets)) :
						if j==i : continue
						hadbak4jet = ak4jets[j]
						#make sure that the hadronic ak4 jet is outside the hadronic W jet
						if hadbak4jet.getFourVector().DeltaR(hadWCandJet.getFourVector())<1.2 : continue
						#if this jet combination produces a viable hypothesis, assign the candidates
						leptopcandvec = lep.getFourVector()+thismet+ak4jets[i].getFourVector()
						hadtopcandvec = hadWCandJet.getFourVector()+ak4jets[j].getFourVector()
						if leptopcandvec.DeltaR(hadtopcandvec) > pi/2 :
							hypotheses.append([lep,thismet,ak4jets[i],hadWCandJet,ak4jets[j]])
		#if no hypotheses had a valid assignment, chuck it
		if len(hypotheses)==0 :
			print 'EVENT NUMBER %d NOT VALID; NO AK4 JET ASSIGNMENT CREATES HEMISPHERICALLY SEPARATED TOPS'%(eventnumber)
			return
		self.nhypotheses.setWriteValue(len(hypotheses))
		#do the monte carlo matching
		corrhypindex=-1
		if self.event_type<2 and not self.is_data :
			for i in range(len(hypotheses)) :
				hypothesis=hypotheses[i]
				#check the lepton, neutrino, and leptonic b
				if hypothesis[0].getFourVector().DeltaR(MClep_vec)<0.1 and hypothesis[1].DeltaPhi(MCv_vec)<0.3 and hypothesis[2].getFourVector().DeltaR(MClepb_vec)<0.4 :
					if len(hypothesis)==4 : #type 1
						#check the hadronic top
						if hypothesis[3].getFourVector().DeltaR(MChadW_vec)<0.8 and hypothesis[3].getFourVector().DeltaR(MChadb_vec)<0.8 :
							corrhypindex = i; self.ismatchable.setWriteValue(1); break
					else : #type 2
						#check the hadronic W and b
						if hypothesis[3].getFourVector().DeltaR(MChadW_vec)<0.8 and hypothesis[4].getFourVector().DeltaR(MChadb_vec)<0.4 :
							corrhypindex = i; self.ismatchable.setWriteValue(1); break
		#	if corrhypindex==-1 : #DEBUG
		#		print 'EVENT NUMBER %d IS NOT MATCHABLE'%(eventnumber) #DEBUG
		#	else : #DEBUG
		#		print 'Event number %d has a correct assignment hypothesis at index %d'%(eventnumber,corrhypindex) #DEBUG
		#send the hypotheses to the kinematic fit
		scaledlep = TLorentzVector(); scaledmet = TLorentzVector()
		scaledlepb = TLorentzVector(); scaledhadt = TLorentzVector()
		scaledhadW = TLorentzVector(); scaledhadb = TLorentzVector()
		hypindex, scaledlep, scaledmet, scaledlepb, scaledhadW, scaledhadb, scaledhadt, fitchi2, finalpars = reconstruct(hypotheses)
		if scaledlep==None :
			print 'EVENT NUMBER '+str(eventnumber)+' NOT VALID; NO KINEMATIC FITS CONVERGED'
			return
		if not self.is_data :
			if hypindex==corrhypindex : self.iscorrect.setWriteValue(1)
			elif corrhypindex!=-1 : self.iscorrect.setWriteValue(0)
		#try the MC matching again with the postfit quantities
		if self.event_type<2 and not self.is_data :
			hypothesis=hypotheses[hypindex]
			#check the lepton, neutrino, and leptonic b
			if scaledlep.DeltaR(MClep_vec)<0.1 and scaledmet.DeltaPhi(MCv_vec)<0.3 and scaledlepb.DeltaR(MClepb_vec)<0.4 :
				if len(hypothesis)==4 : #type 1
					#check the hadronic top
					if scaledhadt.DeltaR(MChadW_vec)<0.8 and scaledhadt.DeltaR(MChadb_vec)<0.8 :
						self.ismatchedpostfit.setWriteValue(1)
				else : #type 2
					#check the hadronic W and b
					if scaledhadW.DeltaR(MChadW_vec)<0.8 and scaledhadb.DeltaR(MChadb_vec)<0.4 :
						self.ismatchedpostfit.setWriteValue(1)
			if self.ismatchedpostfit.getWriteValue()!=1 : self.ismatchedpostfit.setWriteValue(0)
							
		#Kinematic fit debugging variables
		self.chi2.setWriteValue(fitchi2)
		locs = [self.par_0,self.par_1,self.par_2,self.par_3,self.par_4,self.par_5]
		for i in range(len(finalpars)) :
			locs[i].setWriteValue(finalpars[i])
		MW = 80.4; MW_ht2 = 80.4; MT = 172.5; MT_lt1 = 174.2; MT_lt2 = 173.4; MT_ht1 = 187.5; MT_ht2 = 173.8; 
		WW = 2.0; WW_ht2 = 2.0; WT = 1.4; WT_lt1 = 31.80; WT_lt2 = 26.03; WT_ht1 = 18.78; WT_ht2 = 16.39; 
		SIGMAJ  = 0.10; SIGMAL  = 0.03; 
		self.lnLlsf.setWriteValue((finalpars[1]-1.)*(finalpars[1]-1.)/(SIGMAL*SIGMAL))
		self.lnLblsf.setWriteValue((finalpars[2]-1.)*(finalpars[2]-1.)/(SIGMAJ*SIGMAJ))
		self.lnLhsf1.setWriteValue((finalpars[3]-1.)*(finalpars[3]-1.)/(SIGMAJ*SIGMAJ))
		self.lnLhsf2.setWriteValue((finalpars[4]-1.)*(finalpars[4]-1.)/(SIGMAJ*SIGMAJ))
		wl = scaledmet + scaledlep; tl = wl + scaledlepb
		if self.ttype.getWriteValue()==1 :
			th = scaledhadt
			mwl2 = wl.M2(); mtl2 = tl.M2(); mth2 = th.M2()
			self.lnLtl.setWriteValue(1./((mtl2-MT_lt1*MT_lt1)**2+MT_lt1*MT_lt1*WT_lt1*WT_lt1))
			self.lnLth.setWriteValue(1./((mth2-MT_ht1*MT_ht1)**2+MT_ht1*MT_ht1*WT_ht1*WT_ht1))
			self.lnLwl.setWriteValue(((MT_lt1*MT_lt1-mwl2)*(MT_lt1*MT_lt1-mwl2)*(2*MT_lt1*MT_lt1+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW))
		elif self.ttype.getWriteValue()==2 :
			wh = scaledhadW; th = wh+scaledhadb
			mwl2 = wl.M2(); mtl2 = tl.M2(); mwh2 = wh.M2(); mth2 = th.M2()
			self.lnLtl.setWriteValue(1./((mtl2-MT_lt2*MT_lt2)**2+MT_lt2*MT_lt2*WT_lt2*WT_lt2))
			self.lnLth.setWriteValue(1./((mth2-MT_ht2*MT_ht2)**2+MT_ht2*MT_ht2*WT_ht2*WT_ht2))
			self.lnLwl.setWriteValue(((MT_lt2*MT_lt2-mwl2)*(MT_lt2*MT_lt2-mwl2)*(2*MT_lt2*MT_lt2+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW))
			self.lnLwh.setWriteValue(((MT_ht2*MT_ht2-mwh2)*(MT_ht2*MT_ht2-mwh2)*(2*MT_ht2*MT_ht2+mwh2))/((mwh2-MW_ht2*MW_ht2)*(mwh2-MW_ht2*MW_ht2)+MW_ht2*MW_ht2*WW_ht2*WW_ht2))
			self.lnLhsf3.setWriteValue((finalpars[5]-1.)*(finalpars[5]-1.)/(SIGMAJ*SIGMAJ))

		#fill the TTree with the scaled fourvector variables
		self.__setFourVectorBranchValues__('scaled_lep',scaledlep)
		self.__setFourVectorBranchValues__('scaled_met',scaledmet)
		self.__setFourVectorBranchValues__('scaled_lepW',scaledlep+scaledmet)
		self.__setFourVectorBranchValues__('scaled_lepb',scaledlepb)
		self.__setFourVectorBranchValues__('scaled_lept',scaledlep+scaledmet+scaledlepb)
		if len(hypotheses[hypindex])==5 :
			self.__setFourVectorBranchValues__('scaled_hadW',scaledhadW)
			self.__setFourVectorBranchValues__('scaled_hadb',scaledhadb)
		elif len(hypotheses[hypindex])==4 :
			self.__setFourVectorBranchValues__('scaled_hadt',scaledhadt)

		#reconstruct the observables 
		#using the chosen postfit
		cstar_s,xF_s,M_s = getObservables(scaledlep+scaledmet+scaledlepb,scaledhadt,lep.getQ()) 
		self.cstar_scaled.setWriteValue(cstar_s); self.x_F_scaled.setWriteValue(xF_s); self.M_scaled.setWriteValue(M_s)
		#using the chosen prefit
		hypothesis = hypotheses[hypindex]
		prefitlept = hypothesis[0].getFourVector()+hypothesis[1]+hypothesis[2].getFourVector()
		prefithadt = hypothesis[3].getFourVector()
		if len(hypothesis)==5 : prefithadt+=hypothesis[4].getFourVector()
		cstar,xF,M = getObservables(prefitlept,prefithadt,hypothesis[0].getQ()) 
		self.cstar.setWriteValue(cstar); self.x_F.setWriteValue(xF); self.M.setWriteValue(M)
		#using the correct assignment prefit
		if (not self.is_data) and corrhypindex!=-1 : 
			hypothesis = hypotheses[corrhypindex]	
			prefitlepW = hypothesis[0].getFourVector()+hypothesis[1]
			prefitlept = prefitlepW+hypothesis[2].getFourVector()
			prefithadW = hypothesis[3].getFourVector()
			prefithadt = copy.deepcopy(prefithadW)
			if len(hypothesis)==5 : 
				prefithadt+=hypothesis[4].getFourVector()
				self.hadWcorprefitM.setWriteValue(prefithadW.M())
			self.lepWcorprefitM.setWriteValue(prefitlepW.M())
			self.hadWcorprefitM.setWriteValue(prefithadW.M())
			self.leptcorprefitM.setWriteValue(prefitlept.M())
			self.hadtcorprefitM.setWriteValue(prefithadt.M())
			cstar,xF,M = getObservables(prefitlept,prefithadt,hypothesis[0].getQ()) 
			self.cstar_corprefit.setWriteValue(cstar); self.x_F_corprefit.setWriteValue(xF); self.M_corprefit.setWriteValue(M)
		#MC Truth observable and reweighting calculation
		if self.event_type<2 and not self.is_data :
			if self.event_type!=4 :
				( cstar_MC,x_F_MC,M_MC,wg1,wg2,wg3,wg4,wqs1,wqs2,wqa0,wqa1,wqa2,
					wg1_opp,wg2_opp,wg3_opp,wg4_opp,wqs1_opp,wqs2_opp,wqa0_opp,wqa1_opp,wqa2_opp,
					wega, wegc ) = getMCObservables(q_vec,qbar_vec,MCt_vec,MCtbar_vec,self.event_type) 
			self.cstar_MC.setWriteValue(cstar_MC)
			self.x_F_MC.setWriteValue(x_F_MC)
			self.M_MC.setWriteValue(M_MC)
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
#			#scale factor and reweighting calculations
#			if self.lep_type==0 :
#				meas_lep_pt=muons[0].vec.Pt(); meas_lep_eta=muons[0].vec.Eta()
#			elif self.lep_type==1 :
#				meas_lep_pt=electrons[0].vec.Pt(); meas_lep_eta=electrons[0].vec.Eta()
#			#8TeV numbers
#			self.sf_top_pT[0], self.sf_top_pT_low[0], self.sf_top_pT_hi[0] = self.MCcorrector.getToppT_reweight(MCt_vec,MCtbar_vec,self.Q_l[0])
#			self.sf_top_pT_new[0], self.sf_top_pT_new_low[0], self.sf_top_pT_new_hi[0] = self.MCcorrector.getNewToppT_reweight(MCt_vec,MCtbar_vec)
#			self.sf_pileup[0], self.sf_pileup_low[0], self.sf_pileup_hi[0] = self.MCcorrector.getpileup_reweight(MCpileup)
#			( self.sf_lep_ID[0], self.sf_lep_ID_low[0], 
#				self.sf_lep_ID_hi[0] ) = self.MCcorrector.getID_eff(pileup,meas_lep_pt,meas_lep_eta,self.lep_type)
#			( self.sf_trig_eff[0], self.sf_trig_eff_low[0], 
#				self.sf_trig_eff_hi[0] ) = self.MCcorrector.gettrig_eff(pileup,meas_lep_pt,meas_lep_eta,self.lep_type)

		self.__closeout__() #yay! A complete event!

	##################################  #__init__ function  ##################################
	def __init__(self,fileName,tree,isData,xsec,generator,jes,jer,onGrid,totweight) :
		#output file
		self.outfile_name = fileName
		self.outfile = TFile(self.outfile_name,'recreate')
		#output tree
		self.tree = TTree('tree','recreate')
		#event type?
		if fileName.find('qq_semilep_TT')!=-1 :
			print 'only SEMILEPTONIC QQBAR EVENTS will be analyzed from this file'
			self.event_type = 0
		elif fileName.find('gg_semilep_TT')!=-1 :
			print 'only SEMILEPTONIC GG (qg,qiqbarj,etc.) EVENTS will be analyzed from this file'
			self.event_type = 1
		elif fileName.find('dilep_TT')!=-1 :
			print 'only DILEPTONIC EVENTS will be analyzed from this file'
			self.event_type = 2
		elif fileName.find('had_TT')!=-1 :
			print 'only HADRONIC EVENTS will be analyzed from this file'
			self.event_type = 3
		else :
			print 'ALL event types will be analyzed from this file'
			self.event_type = 4
		#input tree
		self.inputTree = tree
		#Set input and output branch locations
		for branch in sorted(self.allBranches.values()) :
			branch.initialize(self.inputTree,self.tree)
		#data or MC?
		self.is_data = isData
		#cross section
		self.xsec = xsec
		#MC generator?
		gen = generator.lower()
		if gen == 'powheg' or gen == 'madgraph' or gen == 'pythia8' or gen == 'mcatnlo' :
			self.MC_generator = gen
		elif gen == 'mg5' :
			self.MC_generator = 'madgraph'
		elif gen == 'mc@nlo' or gen == 'amcatnlo' :
			self.MC_generator = 'mcatnlo'
		elif gen == 'none' :
			self.MC_generator = 'none'
		#JEC systematics?
		self.JES = jes
		self.JER = jer
		#Set the total weight
		self.totalweight = totweight
		#Set the corrector that does event weights and JEC calculations
		self.corrector = Corrector(self.is_data,self.MC_generator,self.event_type,onGrid)

	##################################   reset function   ##################################
	#########  sets all relevant values back to initial values to get ready for next event  ##########
	def reset(self) :
		for branch in self.allBranches.values() :
			branch.reset()
		pass

	##############################   tree filling functions  ###############################		

	def __setFourVectorBranchValues__(self,name,vec) :
		self.allBranches[name+'_pt'].setWriteValue(vec.Pt())
		self.allBranches[name+'_eta'].setWriteValue(vec.Eta())
		self.allBranches[name+'_phi'].setWriteValue(vec.Phi())
		self.allBranches[name+'_M'].setWriteValue(vec.M())

	def __setLeptonBranchValues__(self,name,lep) :
		self.__setFourVectorBranchValues__(name,lep.getFourVector())
		self.allBranches[name+'_Q'].setWriteValue(int(lep.getQ()))
		self.allBranches[name+'_ID'].setWriteValue(int(lep.getID()))
		self.allBranches[name+'_relPt'].setWriteValue(lep.getRelPt())
		self.allBranches[name+'_dR'].setWriteValue(lep.getDR())

	########## function to close out the event, called before kicking back to runner #########
	def __closeout__(self) :
		#fill ttree
		#print 'written M = %.4f, written M_scaled = %.4f'%(self.M.getWriteValue(),self.M_scaled.getWriteValue()) #DEBUG
		self.tree.Fill()

	################################## __del__ function  ###################################
	def __del__(self) :
		self.outfile.cd()
		#self.tree.Write()
		self.outfile.Write()
		self.outfile.Close()
