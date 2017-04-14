#Global variables
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0
#Trigger paths
#MU_TRIG_PATH = 'HLT_Mu30_eta2p1_PFJet150_PFJet50'
#MU_TRIG_PATH = 'HLT_Mu45_eta2p1'
MU_TRIG_PATHS = ['HLT_Mu50','HLT_TkMu50']
#EL_TRIG_PATH = 'HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50'
EL_TRIG_PATHS = ['HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50']

##########								   Imports  								##########

from math import pi, log
from ROOT import TFile, TTree, TLorentzVector
from branch import Branch
from eventTypeHelper import getEventType, findInitialPartons, findMCParticles
from jet import AK4Jet, AK8Jet
from lepton import Muon, Electron
from metHelper import setupMET
from ttbarReconstructor import TTBarReconstructor
from angleReconstructor import getObservables, getMCRWs
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
	#pileup
	npv = AddBranch('pu_NtrueInt','npv','I',-1,'1',[allBranches,ak4JetBranches,ak8JetBranches])
	#Scale, PDF, and alpha_s weights
	genUncBranches = {}
	thisdictlist=[allBranches,genUncBranches]
	scale_size = AddBranch(readname='scale_size',ttreetype='i',dictlist=thisdictlist)
	scale_weights = AddBranch(readname='scale_Weights',ttreetype='F',size='scale_size',dictlist=thisdictlist)
	pdf_size = AddBranch(readname='pdf_size',ttreetype='i',dictlist=thisdictlist)
	pdf_weights = AddBranch(readname='pdf_Weights',ttreetype='F',size='pdf_size',dictlist=thisdictlist)
	alphas_size = AddBranch(readname='alphas_size',ttreetype='i',dictlist=thisdictlist)
	alphas_weights = AddBranch(readname='alphas_Weights',ttreetype='F',size='alphas_size',dictlist=thisdictlist)
	#MC GenEvent info
	mcGenEventBranches = {}
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
	triggerBranches = {}
	thisdictlist = [allBranches,triggerBranches]
	for trigName in MU_TRIG_PATHS :
		AddBranch(trigName,trigName,'I',-1,'1',thisdictlist)
	for trigName in EL_TRIG_PATHS :
		AddBranch(trigName,trigName,'I',-1,'1',thisdictlist)
	#MET Filter information
	filterBranches = {}
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
	mu_isos 	  = AddBranch(readname='mu_Iso04',size='mu_size',dictlist=thisdictlist)
	mus_isMed     = AddBranch(readname='mu_IsMediumMuon',size='mu_size',dictlist=thisdictlist)
	mus_isMed2016 = AddBranch(readname='mu_IsMediumMuon2016',size='mu_size',dictlist=thisdictlist)
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
	el_isos 	  = AddBranch(readname='el_Iso03',size='el_size',dictlist=thisdictlist)
	el_id 		  = AddBranch(readname='el_IDMedium_NoIso',ttreetype='I',size='el_size',dictlist=thisdictlist)
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
	sf_pileup 		  = AddBranch(writename='sf_pileup',inival=1.,dictlist=thisdictlist)
	sf_pileup_low 	  = AddBranch(writename='sf_pileup_low',inival=1.,dictlist=thisdictlist)
	sf_pileup_hi 	  = AddBranch(writename='sf_pileup_hi',inival=1.,dictlist=thisdictlist)
	sf_lep_ID 		  = AddBranch(writename='sf_lep_ID',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_low 	  = AddBranch(writename='sf_lep_ID_low',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_hi 	  = AddBranch(writename='sf_lep_ID_hi',inival=1.,dictlist=thisdictlist)
	sf_mu_R 		  = AddBranch(writename='sf_mu_R',inival=1.,dictlist=thisdictlist)
	sf_mu_R_low 	  = AddBranch(writename='sf_mu_R_low',inival=1.,dictlist=thisdictlist)
	sf_mu_R_hi 		  = AddBranch(writename='sf_mu_R_hi',inival=1.,dictlist=thisdictlist)
	sf_mu_F 		  = AddBranch(writename='sf_mu_F',inival=1.,dictlist=thisdictlist)
	sf_mu_F_low 	  = AddBranch(writename='sf_mu_F_low',inival=1.,dictlist=thisdictlist)
	sf_mu_F_hi 		  = AddBranch(writename='sf_mu_F_hi',inival=1.,dictlist=thisdictlist)
	sf_scale_comb 	  = AddBranch(writename='sf_scale_comb',inival=1.,dictlist=thisdictlist)
	sf_scale_comb_low = AddBranch(writename='sf_scale_comb_low',inival=1.,dictlist=thisdictlist)
	sf_scale_comb_hi  = AddBranch(writename='sf_scale_comb_hi',inival=1.,dictlist=thisdictlist)
	sf_pdf_alphas 	  = AddBranch(writename='sf_pdf_alphas',inival=1.,dictlist=thisdictlist)
	sf_pdf_alphas_low = AddBranch(writename='sf_pdf_alphas_low',inival=1.,dictlist=thisdictlist)
	sf_pdf_alphas_hi  = AddBranch(writename='sf_pdf_alphas_hi',inival=1.,dictlist=thisdictlist)
	#physics objects
	physobjectBranches = {}
	thisdictlist = [allBranches,physobjectBranches]
	#fourvectors
	fourvectornames = ['muon1','muon2','ele1','ele2','lep','met','ak41','ak42','ak43','ak44','ak81','ak82']
	for fourvecname in fourvectornames :
		AddBranch(writename=fourvecname+'_pt',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_eta',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_phi',dictlist=thisdictlist)
		AddBranch(writename=fourvecname+'_M',dictlist=thisdictlist)
	#rescaled fourvectors
	fourvectornames = ['lep','met','lepW','lepb','lept','bigjet','had1','had2','had3']
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
	ak8names = ['ak81','ak82','scaled_bigjet']
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
	#number of b-tagged AK4 jets in the event
	nbTags = AddBranch(writename='nbTags',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of b-tagged AK4 jets used in the chosen jet assignment hypothesis
	nbTagsUsed = AddBranch(writename='nbTagsUsed',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of hadronic side jets used in the 'correct' jet assignment hypothesis
	nHadSideJetsCorrect = AddBranch(writename='nHadSideJetsCorrect',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of hadronic side jets used in the chosen jet assignment hypothesis
	nHadSideJets = AddBranch(writename='nHadSideJets',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of AK4 jets
	nak4jets = AddBranch(writename='nak4jets',ttreetype='I',inival=0,dictlist=thisdictlist)
	#number of AK8 jets
	nak8jets = AddBranch(writename='nak8jets',ttreetype='I',inival=0,dictlist=thisdictlist)
	#Obervables
	observableBranches = {}
	thisdictlist = [allBranches,observableBranches]
	#cosine(theta)
	cstar 		 = AddBranch(writename='cstar',dictlist=thisdictlist)
	cstar_prefit = AddBranch(writename='cstar_prefit',dictlist=thisdictlist)
	cstar_corprefit = AddBranch(writename='cstar_corprefit',dictlist=thisdictlist)
	#Feynman x
	x_F 	   = AddBranch(writename='x_F',dictlist=thisdictlist)
	x_F_prefit = AddBranch(writename='x_F_prefit',dictlist=thisdictlist)
	x_F_corprefit = AddBranch(writename='x_F_corprefit',dictlist=thisdictlist)
	#ttbar invariant mass
	M 		 = AddBranch(writename='M',dictlist=thisdictlist)
	M_prefit = AddBranch(writename='M_prefit',dictlist=thisdictlist)
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
	cutnames = ['metfilters','trigger','onelepton','isolepton','jetcuts','fullselection','validminimization']
	for cutname in cutnames :
		AddBranch(writename=cutname,ttreetype='i',inival=2,dictlist=thisdictlist)
	#debugging variables
	kinfit_debug_branches = {}
	thisdictlist=[kinfit_debug_branches,allBranches]
	chi2 = AddBranch(writename='chi2',inival=0.0,dictlist=thisdictlist)
	par_0 = AddBranch(writename='par_0',dictlist=thisdictlist)
	par_1 = AddBranch(writename='par_1',dictlist=thisdictlist)
	par_2 = AddBranch(writename='par_2',dictlist=thisdictlist)
	par_3 = AddBranch(writename='par_3',dictlist=thisdictlist)
	par_4 = AddBranch(writename='par_4',dictlist=thisdictlist)
	par_5 = AddBranch(writename='par_5',dictlist=thisdictlist)
	nhypotheses = AddBranch(writename='nhypotheses',ttreetype='i',inival=0,dictlist=thisdictlist)
	ismatchable = AddBranch(writename='ismatchable',ttreetype='i',inival=2,dictlist=thisdictlist)
	iscorrect = AddBranch(writename='iscorrect',ttreetype='i',inival=2,dictlist=thisdictlist)
	ismatchedpostfit = AddBranch(writename='ismatchedpostfit',ttreetype='i',inival=2,dictlist=thisdictlist)
	lepWcorprefitM = AddBranch(writename='lepWcorprefitM',dictlist=thisdictlist)
	leptcorprefitM = AddBranch(writename='leptcorprefitM',dictlist=thisdictlist)
	hadtcorprefitM = AddBranch(writename='hadtcorprefitM',dictlist=thisdictlist)

	##################################  ANALYZE FUNCTION  ##################################
	def analyze(self,eventnumber) :
		#get the event in the tree
		self.inputTree.GetEntry(eventnumber)

		#-----------------------------------------------------Below here is a bunch of preselection and object assignment-----------------------------------------------------#

		#Make sure the event has some met, at least one lepton, and at least one AK4 jet
		metsize = self.met_size.getReadValue()
		musize = self.mu_size.getReadValue()
		elsize = self.el_size.getReadValue()
		ak4jetsize = self.ak4_size.getReadValue()
		ak8jetsize = self.ak8_size.getReadValue()
		if not (metsize>0 and (musize>0 or elsize>0) and ak4jetsize>0) :
			#print 'EVENT NUMBER %d NOT VALID; MISSING REQUISITE PHYSICS OBJECTS (metsize = %d, musize = %d, elsize = %d, ak4jetsize = %d, ak8jetsize = %d)'%(eventnumber,metsize,musize,elsize,ak4jetsize,ak8jetsize)
			return

		#MC stuff
		if not self.is_data :
			#event type split
			self.event_type.setWriteValue(getEventType(self.mcGenEventBranches))
			#set the addTwice value
			self.addTwice.setWriteValue(self.event_type.getWriteValue()==0 or (self.event_type.getWriteValue()<4 and self.mcGenEventBranches['MC_part1_ID'].getReadValue()==self.mcGenEventBranches['MC_part2_ID'].getReadValue()))
			#set the eventweight
			eventweight = self.genWeight.getReadValue()*(self.xsec/self.totalweight)
			self.weight.setWriteValue(eventweight)
			#Mother particle and MC truth top assignment
			vecnames = []; vecobjs = []
			if self.event_type!=4 :
				q_vec, qbar_vec = findInitialPartons(self.mcGenEventBranches)
				MCt_vec, MCtbar_vec, MClep_vec, MCv_vec, MClepb_vec, MChadW_vec, MChadb_vec, MClep_charge = findMCParticles(self.mcGenEventBranches)				
				vecnames += ['q','qbar','MCt','MCtbar','MClep','MCv','MClepb','MChadW','MChadb']
				vecobjs  += [q_vec,qbar_vec,MCt_vec,MCtbar_vec,MClep_vec,MCv_vec,MClepb_vec,MChadW_vec,MChadb_vec]
			#write out fourvectors of MC particles
			for i in range(len(vecnames)) :
				if vecobjs[i]!=None :
					self.__setFourVectorBranchValues__(vecnames[i],vecobjs[i])

		#print '------------------------------------------------' #DEBUG

		#For the record, trigger information is handled automatically

		#MET
		#print '	Handling MET...' #DEBUG
		met = TLorentzVector(); met.SetPtEtaPhiM(self.met_pts.getReadValue(), 0., self.met_phis.getReadValue(), 0.)
		self.__setFourVectorBranchValues__('met',met)

		#muons
		#print '	Handling Muons...' #DEBUG
		muons = []
		for i in range(musize) :
			newmuon=Muon(self.muonBranches,i,self.run_era)
			if newmuon.getPt()>55. and abs(newmuon.getEta())<2.5 and newmuon.getID()==1 : muons.append(newmuon)
		muons.sort(key=lambda x: x.getPt(), reverse=True)
		#print '		Added %d Muons.'%(len(muons)) #DEBUG

		#electrons
		#print '	Handling Electrons...' #DEBUG
		electrons = []
		for i in range(elsize) :
			newele = Electron(self.electronBranches,i)
			if newele.getPt()>55. and abs(newele.getEtaSC())<2.5 and newele.getID()==1 : electrons.append(newele)
		electrons.sort(key=lambda x: x.getPt(), reverse=True)
		#print '		Added %d Electrons.'%(len(electrons)) #DEBUG

		if len(muons)==0 and len(electrons)==0 :
			#print 'EVENT NUMBER %d NOT VALID; MISSING LEPTONS (# muons = %d, # electrons = %d)'%(eventnumber,len(ak4jets),len(ak8jets)) #DEBUG
			return

		#jets
		ak4jets = []; ak8jets = [];
		#print '	Adding AK4 Jets...' #DEBUG
		for i in range(ak4jetsize) :
			newJet = AK4Jet(self.ak4JetBranches,i,self.JES,self.JER,muons+electrons,self.corrector,self.is_data)
			if newJet.getFourVector()!=None and newJet.getPt()>15. and abs(newJet.getEta())<3.0 and newJet.isIDed() :
				ak4jets.append(newJet)
		#print '		Added %d AK4 Jets.'%(len(ak4jets)) #DEBUG
		#print '	Adding AK8 Jets...' #DEBUG
		for i in range(ak8jetsize) :
			newJet = AK8Jet(self.ak8JetBranches,i,self.JES,self.JER,muons+electrons,self.corrector,self.is_data)
			if newJet.getFourVector()!=None and newJet.getPt()>200. and abs(newJet.getEta())<2.4 and newJet.getNSubjets()>1 and newJet.isIDed() :
				ak8jets.append(newJet)
		ak4jets.sort(key=lambda x: x.getPt(), reverse=True)
		ak8jets.sort(key=lambda x: x.getPt(), reverse=True)
		#print '		Added %d AK8 Jets.'%(len(ak8jets)) #DEBUG
		#if the lepton cleaning got rid of too many jets toss the event
		if not len(ak4jets)>0 :
			#print 'EVENT NUMBER %d NOT VALID; MISSING JETS (# AK4jets = %d, # AK8Jets = %d)'%(eventnumber,len(ak4jets),len(ak8jets)) #DEBUG
			return
		#set the numbers of jets in the event
		self.nak4jets.setWriteValue(len(ak4jets)); self.nak8jets.setWriteValue(len(ak8jets));

		#calculate lepton isolation variables
		#print '	Calculating lepton isolation values...' #DEBUG
		for muon in muons :
			muon.calculateIsolation(ak4jets)
		for electron in electrons :
			electron.calculateIsolation(ak4jets)
		#print '		Done.'#DEBUG

		#remove the ak4 jets we needed just for isolation calculations
		i=0
		while len(ak4jets)>0 and i<len(ak4jets) :
			ak4jet = ak4jets[i]
			if not (ak4jet.getPt()>30. and abs(ak4jet.getEta())<2.4) :
				ak4jets.pop(i)
			else :
				i+=1
		#print '	Refined AK4 Jets (there are now %d in the event)'%(len(ak4jets)) #DEBUG
		if not len(ak4jets)>1 :
			#print 'EVENT NUMBER %d NOT VALID; MISSING JETS (# AK4jets = %d)'%(eventnumber,len(ak4jets)) #DEBUG
			return
		#count the number of b-tagged AK4 jets
		nbtags = 0
		for ak4jet in ak4jets :
			if ak4jet.isbTagged() :
				nbtags+=1
		self.nbTags.setWriteValue(nbtags)
		#print '		%d of the AK4 jets are b-tagged.'%(nbtags) #DEBUG

		#Set the event toplogy
		#print '	Setting event topology...' #DEBUG
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
		elif n_Wtags>=1 : #If there is instead a W-tagged AK8 jet than this event has a type-2 (partially merged) topology
			topology=2
		elif len(ak8jets)>=1 : #If there are no t- or W-tagged AK8 jets, but there IS at least one AK8 jet, the event is "boosted untagged" (type-3 topology)
			topology=3
		elif len(ak4jets)>=4 and nbtags>=2 : #If there are no AK8 jets at all, but there are at least 4 AK4 jets (at least two b-tagged), the event is "resolved" (type-4 topology)
			topology=4
		else : #otherwise this event is GARBAGE! GARBAGE I tell you!!
			#print 'EVENT NUMBER %d NOT VALID; EVENT TOPOLOGY CANNOT BE DETERMINED! (# top tags = %d, # W tags = %d, #AK4 jets=%d, #b-tags=%d)'%(eventnumber,ttags,wtags,len(ak4jets),nbtags) #DEBUG
			return
		self.event_topology.setWriteValue(topology)
		#print '		There are %d top-tagged and %d W-tagged jets; the event topology is type %d'%(ttags,wtags,topology) #DEBUG

		#figure out whether the event is muonic or electronic, assign lep
		allleps = muons+electrons
		allleps.sort(key=lambda x: x.getPt(), reverse=True)
		lep = allleps[0]
		for lepcand in allleps[1:] :
			if (topology<4 and (lepcand.getRelPt()>20. or lepcand.getDR()>0.4)) or (topology==4 and lepcand.isIso()) :
				lep = lepcand
				if lep.getType()=='mu' : self.lepflavor.setWriteValue(1)
				elif lep.getType()=='el' : self.lepflavor.setWriteValue(2)
				else :
					#print 'EVENT NUMBER %d NOT VALID; LEPTON FLAVOR %s UNIDENTIFIED'%(eventnumber, lep.getType()) #DEBUG
					return
				break
		#write the analysis lepton variables
		self.__setLeptonBranchValues__('lep',lep)
		#print '	Lepton assigned (leading isolated %s).'%('muon' if self.lepflavor.getWriteValue()==1 else 'electron') #DEBUG

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
		self.__setFourVectorBranchValues__('ak41',ak4jets[0].getFourVector())
		self.__setFourVectorBranchValues__('ak42',ak4jets[1].getFourVector())
		if len(ak4jets)>2 :
			self.__setFourVectorBranchValues__('ak43',ak4jets[2].getFourVector())
			if len(ak4jets)>3 :
				self.__setFourVectorBranchValues__('ak44',ak4jets[3].getFourVector())
		if len(ak8jets)>0 : 
			self.__setFourVectorBranchValues__('ak81',ak8jets[0].getFourVector())
			self.allBranches['ak81_tau32'].setWriteValue(ak8jets[0].getTau32())
			self.allBranches['ak81_tau21'].setWriteValue(ak8jets[0].getTau21())
			self.allBranches['ak81_SDM'].setWriteValue(ak8jets[0].getSDM())
			self.allBranches['ak81_isttagged'].setWriteValue(1) if ak8jets[0].isTopTagged() else self.allBranches['ak81_isttagged'].setWriteValue(0)
			self.allBranches['ak81_isWtagged'].setWriteValue(1) if ak8jets[0].isWTagged() else self.allBranches['ak81_isWtagged'].setWriteValue(0)
			if len(ak8jets)>1 :
				self.__setFourVectorBranchValues__('ak82',ak8jets[1].getFourVector())
				self.allBranches['ak82_tau32'].setWriteValue(ak8jets[1].getTau32())
				self.allBranches['ak82_tau21'].setWriteValue(ak8jets[1].getTau21())
				self.allBranches['ak82_SDM'].setWriteValue(ak8jets[1].getSDM())
				self.allBranches['ak82_isttagged'].setWriteValue(1) if ak8jets[1].isTopTagged() else self.allBranches['ak82_isttagged'].setWriteValue(0)
				self.allBranches['ak82_isWtagged'].setWriteValue(1) if ak8jets[1].isWTagged() else self.allBranches['ak82_isWtagged'].setWriteValue(0)

		#neutrino handling and setup for fit
		met1_vec, met2_vec = setupMET(lep.getFourVector(),met)
		self.nMETs.setWriteValue(2) if met1_vec.Pz() != met2_vec.Pz() else self.nMETs.setWriteValue(1)

		#--------------------------------------------------------------------Below here is event selection--------------------------------------------------------------------#
		
		#print '	Calculating cut variables...' #DEBUG

		#met filtering
		metfiltercuts = []
		for branch in self.filterBranches.values() :
			if self.is_data==0 and branch.getReadName()=='Flag_eeBadScFilter' : #the one filter that shouldn't be applied to Monte Carlo
				continue
			if branch.getWriteValue()!=1 :
				metfiltercuts.append(False)
		self.cut_branches['metfilters'].setWriteValue(1) if metfiltercuts.count(False)==0 else self.cut_branches['metfilters'].setWriteValue(0)
		#isolated lepton
		if topology<4 :
			self.cut_branches['isolepton'].setWriteValue(1) if (lep.getRelPt()>20. or lep.getDR()>0.4) else self.cut_branches['isolepton'].setWriteValue(0)
		elif topology==4 :
			self.cut_branches['isolepton'].setWriteValue(1) if lep.isIso() else self.cut_branches['isolepton'].setWriteValue(0)
		#other cuts are lepton flavor specific
		other_leps = []; allcuts = []
		if self.lepflavor.getWriteValue()==1 :
			#trigger
			for trigName in MU_TRIG_PATHS :
				if self.triggerBranches[trigName].getWriteValue()==1 :
					self.cut_branches['trigger'].setWriteValue(1)
					break
				self.cut_branches['trigger'].setWriteValue(0)
			#exactly one isolated lepton
			other_leps+=electrons
			other_leps+=muons[1:]
			i=0
			while len(other_leps)>0 and i<len(other_leps) :
				lepcand = other_leps[i]
				if (topology<4 and not (lep.getRelPt()>20. or lep.getDR()>0.4)) or (topology==4 and not lep.isIso()) :
					other_leps.pop(i)
				else :
					i+=1
			self.cut_branches['onelepton'].setWriteValue(1) if len(other_leps)==0 else self.cut_branches['onelepton'].setWriteValue(0)
			#leading ak4 jets
			self.cut_branches['jetcuts'].setWriteValue(1) if (topology==4 or (ak4jets[0].getPt()>150. and ak4jets[1].getPt()>50.)) else self.cut_branches['jetcuts'].setWriteValue(0)
			#add'l cuts
			allcuts.append(topology==4 or (lep.getPt()+met.E())>150.)
			allcuts.append(topology==4 or met.E()>50.)
		elif self.lepflavor.getWriteValue()==2 :
			#trigger
			for trigName in EL_TRIG_PATHS :
				if self.triggerBranches[trigName].getWriteValue()==1 :
					self.cut_branches['trigger'].setWriteValue(1)
					break
				self.cut_branches['trigger'].setWriteValue(0)
			other_leps+=muons
			other_leps+=electrons[1:]
			i=0
			while len(other_leps)>0 and i<len(other_leps) :
				lepcand = other_leps[i]
				if (topology<4 and not (lep.getRelPt()>20. or lep.getDR()>0.4)) or (topology==4 and not lep.isIso()) :
					other_leps.pop(i)
				else :
					i+=1
			self.cut_branches['onelepton'].setWriteValue(1) if len(other_leps)==0 else self.cut_branches['onelepton'].setWriteValue(0)
			#leading ak4 jets
			self.cut_branches['jetcuts'].setWriteValue(1) if (topology==4 or (ak4jets[0].getPt()>250. and ak4jets[1].getPt()>70.)) else self.cut_branches['jetcuts'].setWriteValue(0)
			#add'l cuts
			allcuts.append(topology==4 or met.E()>120.)
		#full selection
		for cutbranch in self.cut_branches.values() :
			if cutbranch.getWriteValue()==0 :
				allcuts.append(False)
		self.cut_branches['fullselection'].setWriteValue(1) if allcuts.count(False)==0 else self.cut_branches['fullselection'].setWriteValue(0)

		#print '		Done. (fullselection=%d)'%(self.cut_branches['fullselection'].getWriteValue()) #DEBUG

		#-----------------------------------------------------------------Below here is event reconstruction-----------------------------------------------------------------#

		#print '	Reconstructing event...' #DEBUG

		#build the list of jet assignment hypotheses
		hypotheses = []
		met_options = [met1_vec,met2_vec]
		for i in range(self.nMETs.getWriteValue()) :
			#The leptonic side is the same for everything
			#there's only one lepton, so first choose the MET hypothesis
			thismet = met_options[i]
			#for any choice of leptonic-side b-jet
			for j in range(len(ak4jets)) :
				lepbCandJet = ak4jets[j]
				lepbistagged = lepbCandJet.isbTagged()
				#the rest is topology-dependent
				#FULLY-MERGED EVENTS: 
				#the hadronic top candidate is the AK8 jet
				#just choose the leptonic b candidate from all AK4 jets
				#hypotheses are [lepton, neutrino, leptonic b-jet, [hadronic top jet]]
				if topology==1 : 
					#if there are btags the leptonic b candidate must be one of them
					if nbtags>0 and not lepbistagged :
						continue
					#for all the choices of a merged top
					for ttag in ttags :
						#append this hypothesis
						hypotheses.append([lep,thismet,lepbCandJet,[ttag]])
					continue
				#other events need more than one AK4 jet
				#loop over the other ak4 jets to find the first hadronic-side AK4 jet
				for k in range(len(ak4jets)) :
					if j==k :
						continue
					had1CandJet = ak4jets[k]
					had1istagged = had1CandJet.isbTagged()
					hadsidehasbtag = had1istagged
					#PARTIALLY-MERGED EVENTS: 
					#the hadronic top candidate is the W-tagged AK8 jet with an additional AK4 jet
					#hypotheses are [lepton, neutrino, leptonic b-jet, hadronic W jet, hadronic b jet]
					if topology==2 :
						#if the btags are in place
						if (nbtags==1 and not (hadsidehasbtag or lepbistagged)) or (nbtags>1 and not (hadsidehasbtag and lepbistagged)) :
							continue
						#for any choice of W-tag
						for wtag in wtags :
							#append this hypothesis
							hypotheses.append([lep,thismet,lepbCandJet,[wtag,had1CandJet]])
						continue
					#BOOSTED UNTAGGED AND RESOLVED EVENTS:
					#the hadronic side of the decay can have up to three jets 
					#hypotheses are [lepton, neutrino, leptonic b-jet, [list of up to three hadronic-side jets (exactly three for resolved events)]]
					#BOOSTED UNTAGGED WITH ONE HADRONIC_SIDE JET
					if topology==3 and (nbtags==0 or lepbistagged) :
						#append this hypothesis
						hypotheses.append([lep,thismet,lepbCandJet,[had1CandJet]])
					#loop over the remaining AK4 jets for a second hadronic-side jet
					for m in range(k+1,len(ak4jets)) :
						if j==m :
							continue
						had2CandJet = ak4jets[m]
						had2istagged = had2CandJet.isbTagged()
						hadsidehasbtag = had1istagged or had2istagged
						#BOOSTED UNTAGGED WITH TWO HADRONIC SIDE JETS
						if topology==3 and (nbtags==0 or (nbtags==1 and (hadsidehasbtag or lepbistagged)) or (nbtags>1 and (hadsidehasbtag and lepbistagged))) :
							#append this hypothesis
							hypotheses.append([lep,thismet,lepbCandJet,[had1CandJet,had2CandJet]])
						#loop one last time over the remaining AK4 jets for a third and final hadronic-side jet
						for n in range(m+1,len(ak4jets)) :
							if j==n :
								continue
							had3CandJet = ak4jets[n]
							had3istagged = had3CandJet.isbTagged()
							hadsidehasbtag = had1istagged or had2istagged or had3istagged
							#BOOSTED UNTAGGED WITH THREE HADRONIC SIDE JETS, AND FULLY RESOLVED
							if topology==3 or topology==4 :
								#make sure the btags are in place
								if (nbtags==1 and not (hadsidehasbtag or lepbistagged)) or (nbtags==2 and not (hadsidehasbtag and lepbistagged)) :
									continue
								#append this hypothesis
								hypotheses.append([lep,thismet,lepbCandJet,[had1CandJet,had2CandJet,had3CandJet]])
		#print '		Will try %d jet assignment hypotheses...'%(len(hypotheses)) #DEBUG
		self.nhypotheses.setWriteValue(len(hypotheses))
		#do the monte carlo matching
		corrhypindex=-1; mindMdR = 10000.
		if self.event_type.getWriteValue()<2 and not self.is_data :
			#print '----------------event %d--------------------'%(eventnumber) #DEBUG
			for i in range(len(hypotheses)) :
				hypothesis=hypotheses[i]
				#get the MC truth top quark vectors
				MClept_vec = MCt_vec if MClep_charge>0 else MCtbar_vec
				MChadt_vec = MCtbar_vec if MClep_charge>0 else MCt_vec
				#check the lepton 
				hyplvec = hypothesis[0].getFourVector()
				ldR = hyplvec.DeltaR(MClep_vec)
				ldM = hyplvec.M()-MClep_vec.M()
				if not ldR<0.1 : continue
				#neutrino
				vdP = abs(hypothesis[1].DeltaPhi(MCv_vec))
				if not vdP<0.3 : continue
				#leptonic b
				hyplbvec = hypothesis[2].getFourVector()
				lbdR = hyplbvec.DeltaR(MClepb_vec)
				lbdM = hyplbvec.M()-MClepb_vec.M()
				if not lbdR<0.4 : continue
				#leptonic top
				hypltvec = hyplvec+hypothesis[1]+hyplbvec
				ltdR = hypltvec.DeltaR(MClept_vec)
				ltdM = hypltvec.M()-MClept_vec.M()
				#the rest are topology-dependent
				hyphtvec = None; hyphWvec = None; hyphbvec = None
				if topology==1 : hyphtvec = hypothesis[3][0].getFourVector()
				elif topology==2 : 
					hyphtvec = hypothesis[3][0].getFourVector()+hypothesis[3][1].getFourVector()
					hyphWvec = hypothesis[3][0].getFourVector()
					hyphbvec = hypothesis[3][1].getFourVector()
				elif topology==3 :
					hyphtvec = hypothesis[3][0].getFourVector()
					for otherjet in hypothesis[3][1:] :
						hyphtvec = hyphtvec+otherjet.getFourVector()
				elif topology==4 : hyphtvec = hypothesis[3][0].getFourVector()+hypothesis[3][1].getFourVector()+hypothesis[3][2].getFourVector()
				if hyphtvec!=None :
					htdR = hyphtvec.DeltaR(MChadt_vec)
					htdM = hyphtvec.M()-MChadt_vec.M()
				if hyphWvec!=None and hyphbvec!=None :
					hWdR = hyphWvec.DeltaR(MChadW_vec)
					hWdM = hyphWvec.M()-MChadW_vec.M()
					hbdR = hyphbvec.DeltaR(MChadb_vec)
					hbdM = hyphbvec.M()-MChadb_vec.M()
				sumdMdR = ltdR*abs(ltdM/MClept_vec.M())+htdR*abs(htdM/MChadt_vec.M())
				#print 'hypothesis number %d: ldR/dM=%.2f/%.2f, vdP=%.2f, lbdR/dM=%.2f/%.2f, ltdR/dM=%.2f/%.2f, htdR/dM=%.2f/%.2f, sumdMdR=%.2f, mindMdR=%.2f'%(i,ldR,ldM,vdP,lbdR,lbdM,ltdR,ltdM,htdR,htdM,sumdMdR,mindMdR) #DEBUG					
				#if the matching deltaR requirements are satisfied and this hypothesis has the most accurate mass thus far it's the new matched hypothesis.
				if htdR<0.8 and (topology!=2 or (topology==2 and hWdR<0.8 and hbdR<0.4)) and sumdMdR<mindMdR :
					mindMdR = sumdMdR; corrhypindex = i; self.nHadSideJetsCorrect.setWriteValue(len(hypothesis[3])); self.ismatchable.setWriteValue(1)
			#print '------------------------------------------------------------------------' #DEBUG
			if corrhypindex==-1 :
				#print 'EVENT NUMBER %d IS NOT MATCHABLE'%(eventnumber) #DEBUG
				self.ismatchable.setWriteValue(0)
			#else : #DEBUG
				#print 'Event number %d has a correct assignment hypothesis at index %d'%(eventnumber,corrhypindex) #DEBUG
		#send the hypotheses to the kinematic fit
		scaledlep = TLorentzVector(); scaledmet = TLorentzVector(); scaledlepb = TLorentzVector(); 
		scaledhadWs1 = TLorentzVector(); scaledhadWs2 = TLorentzVector(); scaledhadW = TLorentzVector();
		scaledhadb = TLorentzVector(); scaledhadt = TLorentzVector()
		hypindex, scaledlep, scaledmet, scaledlepb, scaledbigjet, scaledhad1, scaledhad2, scaledhad3, fitchi2, finalpars = self.ttbarreconstructor.reconstruct(hypotheses,topology)
		if scaledlep==None :
			#print 'EVENT NUMBER '+str(eventnumber)+' NOT VALID; NO KINEMATIC FITS CONVERGED' #DEBUG
			self.cut_branches['validminimization'].setWriteValue(0)
			self.cut_branches['fullselection'].setWriteValue(0)
			self.__closeout__()
			return
		self.cut_branches['validminimization'].setWriteValue(1)
		#build the post-fit hadronic top
		scaledhadt = scaledbigjet if scaledbigjet!=None else scaledhad1
		if topology==2 :
			scaledhadt = scaledhadt+scaledhad3
		elif topology>2 :
			if scaledhad2!=None :
				scaledhadt = scaledhadt+scaledhad2
				if scaledhad3!=None :
					scaledhadt = scaledhadt+scaledhad3
		#set the number of bTags used
		nusedbtags=0
		if hypotheses[hypindex][2].isbTagged() : nusedbtags+=1 #leptonic side b-jet
		for thisjet in hypotheses[hypindex][3] : #hadronic side jets
			if thisjet.isbTagged() : nusedbtags+=1
		self.nbTagsUsed.setWriteValue(nusedbtags)
		self.nHadSideJets.setWriteValue(3-[scaledhad1,scaledhad2,scaledhad3].count(None))

		#print '		Done.' #DEBUG

		#-----------------------------------------------------Below here is a bunch of variable and weight calculation-----------------------------------------------------#

		if not self.is_data :
			#if the fit returned the correct hypothesis with either MET solution, it was correct
			if hypindex==corrhypindex or (corrhypindex!=-1 and hypotheses[hypindex][2]==hypotheses[corrhypindex][2]) : 
				self.iscorrect.setWriteValue(1)
			elif corrhypindex!=-1 : 
				self.iscorrect.setWriteValue(0)
		#try the MC matching again with the postfit quantities
		if self.event_type.getWriteValue()<2 and not self.is_data :
			hypothesis=hypotheses[hypindex]
			self.ismatchedpostfit.setWriteValue(1) if (scaledlep.DeltaR(MClep_vec)<0.1 and scaledmet.DeltaPhi(MCv_vec)<0.3 and scaledlepb.DeltaR(MClepb_vec)<0.6 and scaledhadt.DeltaR(MChadt_vec)<1.2) else self.ismatchedpostfit.setWriteValue(0)
		#	print 'ismatchedpostfit = %d'%(self.ismatchedpostfit.getWriteValue()) #DEBUG
							
		#Kinematic fit debugging variables
		self.chi2.setWriteValue(fitchi2)
		locs = [self.par_0,self.par_1,self.par_2,self.par_3,self.par_4,self.par_5]
		for i in range(len(finalpars)) :
			locs[i].setWriteValue(finalpars[i])

		#fill the TTree with the scaled fourvector variables
		self.__setFourVectorBranchValues__('scaled_lep',scaledlep)
		self.__setFourVectorBranchValues__('scaled_met',scaledmet)
		self.__setFourVectorBranchValues__('scaled_lepW',scaledlep+scaledmet)
		self.__setFourVectorBranchValues__('scaled_lepb',scaledlepb)
		self.__setFourVectorBranchValues__('scaled_lept',scaledlep+scaledmet+scaledlepb)
		if scaledbigjet!=None : 
			self.__setFourVectorBranchValues__('scaled_bigjet',scaledbigjet)
			self.allBranches['scaled_bigjet_tau32'].setWriteValue(hypotheses[hypindex][3][0].getTau32())
			self.allBranches['scaled_bigjet_tau21'].setWriteValue(hypotheses[hypindex][3][0].getTau21())
			self.allBranches['scaled_bigjet_SDM'].setWriteValue(hypotheses[hypindex][3][0].getSDM())
			self.allBranches['scaled_bigjet_isttagged'].setWriteValue(1) if hypotheses[hypindex][3][0].isTopTagged() else self.allBranches['scaled_bigjet_isttagged'].setWriteValue(0)
			self.allBranches['scaled_bigjet_isWtagged'].setWriteValue(1) if hypotheses[hypindex][3][0].isWTagged() else self.allBranches['scaled_bigjet_isWtagged'].setWriteValue(0)
		if scaledhad1!=None : self.__setFourVectorBranchValues__('scaled_had1',scaledhad1)
		if scaledhad2!=None : self.__setFourVectorBranchValues__('scaled_had2',scaledhad2)
		if scaledhad3!=None : self.__setFourVectorBranchValues__('scaled_had3',scaledhad3)

		#reconstruct the observables 
		#using the chosen postfit
		cstar_s,xF_s,M_s = getObservables(scaledlep+scaledmet+scaledlepb,scaledhadt,lep.getQ()) 
		self.cstar.setWriteValue(cstar_s); self.x_F.setWriteValue(xF_s); self.M.setWriteValue(M_s)
		#using the chosen prefit
		hypothesis = hypotheses[hypindex]
		prefitlept = hypothesis[0].getFourVector()+hypothesis[1]+hypothesis[2].getFourVector()
		prefithadt = hypothesis[3][0].getFourVector()
		for otherjet in hypothesis[3][1:] :
			prefithadt = prefithadt+otherjet.getFourVector()
		cstar,xF,M = getObservables(prefitlept,prefithadt,hypothesis[0].getQ()) 
		self.cstar_prefit.setWriteValue(cstar); self.x_F_prefit.setWriteValue(xF); self.M_prefit.setWriteValue(M)
		#using the correct assignment prefit
		#print 'isData = %s, corrhypindex=%d'%(self.is_data,corrhypindex) #DEBUG
		if (not self.is_data) and corrhypindex!=-1 : 
			hypothesis = hypotheses[corrhypindex]	
			prefitlepW = hypothesis[0].getFourVector()+hypothesis[1]
			prefitlept = prefitlepW+hypothesis[2].getFourVector()
			prefithadt = hypothesis[3][0].getFourVector()
			for otherjet in hypothesis[3][1:] :
				prefithadt = prefithadt+otherjet.getFourVector()
			#print 'correct hypothesis leptonic W mass = %.2f, leptonic top mass = %.2f, hadronic top mass = %.2f'%(prefitlepW.M(),prefitlept.M(),prefithadt.M()) #DEBUG
			self.lepWcorprefitM.setWriteValue(prefitlepW.M())
			self.leptcorprefitM.setWriteValue(prefitlept.M())
			self.hadtcorprefitM.setWriteValue(prefithadt.M())
			cstar,xF,M = getObservables(prefitlept,prefithadt,hypothesis[0].getQ()) 
			#print 'setting correct prefit observable values: cstar = %.2f, x_F = %.2f, M = %.2f'%(cstar,x_F,M) #DEBUG
			self.cstar_corprefit.setWriteValue(cstar); self.x_F_corprefit.setWriteValue(xF); self.M_corprefit.setWriteValue(M)
		#MC Truth observable and reweighting calculation
		if self.event_type.getWriteValue()<2 and not self.is_data :
			self.cstar_MC.setWriteValue(self.mcGenEventBranches['MC_cstar'].getReadValue())
			self.x_F_MC.setWriteValue(self.mcGenEventBranches['MC_x_F'].getReadValue())
			self.M_MC.setWriteValue(self.mcGenEventBranches['MC_Mtt'].getReadValue())
			if self.event_type.getWriteValue()!=4 :
				( wg1,wg2,wg3,wg4,wqs1,wqs2,wqa0,wqa1,wqa2,
					wg1_opp,wg2_opp,wg3_opp,wg4_opp,wqs1_opp,wqs2_opp,wqa0_opp,wqa1_opp,wqa2_opp,
					wega, wegc ) = getMCRWs(self.cstar_MC.getWriteValue(),MCt_vec,MCtbar_vec,self.alpha,self.epsilon) 
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
			pu_sf, pu_sf_up, pu_sf_down = self.corrector.getPileupReweight(self.allBranches['npv'].getReadValue())
			self.sf_pileup.setWriteValue(pu_sf); self.sf_pileup_hi.setWriteValue(pu_sf_up); self.sf_pileup_low.setWriteValue(pu_sf_down)
			#Lepton ID reweighting
			lep_ID_sf, lep_ID_sf_up, lep_ID_sf_down = self.corrector.getIDEff(self.allBranches['npv'].getReadValue(),lep)
			self.sf_lep_ID.setWriteValue(lep_ID_sf); self.sf_lep_ID_hi.setWriteValue(lep_ID_sf_up); self.sf_lep_ID_low.setWriteValue(lep_ID_sf_down)
			#Scale, and pdf/alpha_s reweights
			( mu_R_sf, mu_R_sf_up, mu_R_sf_down, mu_F_sf, mu_F_sf_up, mu_F_sf_down,
			 scale_comb_sf, scale_comb_sf_up, scale_comb_sf_down, pdf_alphas_sf, pdf_alphas_sf_up, pdf_alphas_sf_down ) = self.corrector.getGenReweights(self.genUncBranches)
			self.sf_mu_R.setWriteValue(mu_R_sf); self.sf_mu_R_hi.setWriteValue(mu_R_sf_up); self.sf_mu_R_low.setWriteValue(mu_R_sf_down)
			self.sf_mu_F.setWriteValue(mu_F_sf); self.sf_mu_F_hi.setWriteValue(mu_F_sf_up); self.sf_mu_F_low.setWriteValue(mu_F_sf_down)
			self.sf_scale_comb.setWriteValue(scale_comb_sf); self.sf_scale_comb_hi.setWriteValue(scale_comb_sf_up); self.sf_scale_comb_low.setWriteValue(scale_comb_sf_down)
			self.sf_pdf_alphas.setWriteValue(pdf_alphas_sf); self.sf_pdf_alphas_hi.setWriteValue(pdf_alphas_sf_up); self.sf_pdf_alphas_low.setWriteValue(pdf_alphas_sf_down)
#			( self.sf_trig_eff[0], self.sf_trig_eff_low[0], 
#				self.sf_trig_eff_hi[0] ) = self.MCcorrector.gettrig_eff(pileup,meas_lep_pt,meas_lep_eta,self.lep_type)

		self.__closeout__() #yay! A complete event!

	##################################  #__init__ function  ##################################
	def __init__(self,fileName,tree,isData,xsec,jes,jer,onGrid,pu_histo,totweight,renormdict) :
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
		#JEC systematics?
		self.JES = jes
		self.JER = jer
		#Set the total weight
		self.totalweight = totweight
		#Set the alpha and epsilon values used to calculate event reweights
		self.alpha = renormdict['alpha']
		self.epsilon = renormdict['epsilon']
		#Set the ttbar reconstructor object
		self.ttbarreconstructor = TTBarReconstructor()
		#Set the corrector that does event weights and JEC calculations
		self.corrector = Corrector(self.is_data,onGrid,pu_histo,self.run_era,renormdict)

	##################################   reset function   ##################################
	#########  sets all relevant values back to initial values to get ready for next event  ##########
	def reset(self) :
		for branch in self.allBranches.values() :
			branch.reset()

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

	########## function to close out the event, called before kicking back to runner #########
	def __closeout__(self) :
		#fill ttree
		#print 'written M = %.4f, written M_prefit = %.4f'%(self.M.getWriteValue(),self.M_prefit.getWriteValue()) #DEBUG
		self.tree.Fill()

	################################## __del__ function  ###################################
	def __del__(self) :
		self.outfile.cd()
		self.outfile.Write()
		self.outfile.Close()

 