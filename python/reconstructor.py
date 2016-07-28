#Global variables
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0
#Trigger paths
MU_TRIG_PATH = 'HLT_Mu30_eta2p1_PFJet150_PFJet50'
EL_TRIG_PATH = 'HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50'

##########								   Imports  								##########

from math import pi
from ROOT import TFile, TTree, TLorentzVector
from branch import Branch
from eventTypeHelper import keepEventType, findInitialPartons, findMCTops
from jet import AK4Jet, AK8Jet
from lepton import Muon, Electron
from metHelper import setupMET
from ttbarReconstructor import reconstruct
from mc_corrector import MC_corrector

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
	#GenWeight info
	genWeight = AddBranch('evt_Gen_Weight','genWeight','F',1.0,'1',[allBranches])
	xsec 	  = AddBranch('evt_XSec','xsec','F',1.0,'1',[allBranches])
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
	#pileup
	pileupBranches = {}
	thisdictlist = [allBranches,pileupBranches]
	pu = AddBranch(readname='pu_NtrueInt',ttreetype='I',dictlist=thisdictlist)
	#MET
	metBranches = {}
	thisdictlist = [allBranches,metBranches]
	met_size 	= AddBranch(readname='met_size',ttreetype='i',dictlist=thisdictlist)
	met_pts 	= AddBranch(readname='met_Pt',size='met_size',dictlist=thisdictlist)
	met_pxs 	= AddBranch(readname='met_Px',size='met_size',dictlist=thisdictlist)
	met_pys 	= AddBranch(readname='met_Py',size='met_size',dictlist=thisdictlist)
	met_phis 	= AddBranch(readname='met_Phi',size='met_size',dictlist=thisdictlist)
	#muons
	muonBranches = {}
	thisdictlist = [allBranches,muonBranches]
	mu_size 	= AddBranch(readname='mu_size',ttreetype='i',dictlist=thisdictlist)
	mu_pts 		= AddBranch(readname='mu_Pt',size='mu_size',dictlist=thisdictlist)
	mu_etas 	= AddBranch(readname='mu_Eta',size='mu_size',dictlist=thisdictlist)
	mu_phis 	= AddBranch(readname='mu_Phi',size='mu_size',dictlist=thisdictlist)
	mu_es 		= AddBranch(readname='mu_E',size='mu_size',dictlist=thisdictlist)
	mu_charges 	= AddBranch(readname='mu_Charge',size='mu_size',dictlist=thisdictlist)
	mus_isLoose = AddBranch(readname='mu_IsLooseMuon',size='mu_size',dictlist=thisdictlist)
	#electrons
	electronBranches = {}
	thisdictlist = [allBranches,electronBranches]
	el_size 	= AddBranch(readname='el_size',ttreetype='i',dictlist=thisdictlist)
	el_pts 		= AddBranch(readname='el_Pt',size='el_size',dictlist=thisdictlist)
	el_etas 	= AddBranch(readname='el_Eta',size='el_size',dictlist=thisdictlist)
	el_phis 	= AddBranch(readname='el_Phi',size='el_size',dictlist=thisdictlist)
	el_es 		= AddBranch(readname='el_E',size='el_size',dictlist=thisdictlist)
	el_charges 	= AddBranch(readname='el_Charge',size='el_size',dictlist=thisdictlist)
	el_detas_in = AddBranch(readname='el_dEtaIn',size='el_size',dictlist=thisdictlist)
	el_dphis_in = AddBranch(readname='el_dPhiIn',size='el_size',dictlist=thisdictlist)
	el_5x5siees = AddBranch(readname='el_full5x5siee',size='el_size',dictlist=thisdictlist)
	el_HoEs 	= AddBranch(readname='el_HoE',size='el_size',dictlist=thisdictlist)
	el_D0s 		= AddBranch(readname='el_D0',size='el_size',dictlist=thisdictlist)
	el_Dzs 		= AddBranch(readname='el_Dz',size='el_size',dictlist=thisdictlist)
	el_ooEmooPs = AddBranch(readname='el_ooEmooP',size='el_size',dictlist=thisdictlist)
	el_hasMCVs 	= AddBranch(readname='el_hasMatchedConVeto',size='el_size',dictlist=thisdictlist)
	el_misshits = AddBranch(readname='el_missHits',size='el_size',dictlist=thisdictlist)
	el_vidLoose = AddBranch(readname='el_vidLoose',size='el_size',dictlist=thisdictlist)
	#AK4 Jets
	ak4JetBranches = {}
	thisdictlist = [allBranches,ak4JetBranches]
	ak4_size 	   = AddBranch(readname='jetAK4_size',ttreetype='i',dictlist=thisdictlist)
	ak4_pts 	   = AddBranch(readname='jetAK4_Pt',size='jetAK4_size',dictlist=thisdictlist)
	ak4_etas 	   = AddBranch(readname='jetAK4_Eta',size='jetAK4_size',dictlist=thisdictlist)
	ak4_phis 	   = AddBranch(readname='jetAK4_Phi',size='jetAK4_size',dictlist=thisdictlist)
	ak4_es 		   = AddBranch(readname='jetAK4_E',size='jetAK4_size',dictlist=thisdictlist)
	ak4_smpts 	   = AddBranch(readname='jetAK4_SmearedPt',size='jetAK4_size',dictlist=thisdictlist)
	ak4_smetas 	   = AddBranch(readname='jetAK4_SmearedPEta',size='jetAK4_size',dictlist=thisdictlist)
	ak4_smphis 	   = AddBranch(readname='jetAK4_SmearedPhi',size='jetAK4_size',dictlist=thisdictlist)
	ak4_smes 	   = AddBranch(readname='jetAK4_SmearedE',size='jetAK4_size',dictlist=thisdictlist)
	ak4_csvv2s 	   = AddBranch(readname='jetAK4_CSVv2',size='jetAK4_size',dictlist=thisdictlist)
	ak4_jecuncs    = AddBranch(readname='jetAK4_jecUncertainty',size='jetAK4_size',dictlist=thisdictlist)
	ak4_JERSFs 	   = AddBranch(readname='jetAK4_JERSF',size='jetAK4_size',dictlist=thisdictlist)
	ak4_JERSFUps   = AddBranch(readname='jetAK4_JERSFUp',size='jetAK4_size',dictlist=thisdictlist)
	ak4_JERSFDowns = AddBranch(readname='jetAK4_JERSFDown',size='jetAK4_size',dictlist=thisdictlist)
	#AK8 Jets
	ak8JetBranches = {}
	thisdictlist = [allBranches,ak8JetBranches]
	ak8_size 	   = AddBranch(readname='jetAK8_size',ttreetype='i',dictlist=thisdictlist)
	ak8_pts 	   = AddBranch(readname='jetAK8_Pt',size='jetAK8_size',dictlist=thisdictlist)
	ak8_etas 	   = AddBranch(readname='jetAK8_Eta',size='jetAK8_size',dictlist=thisdictlist)
	ak8_phis 	   = AddBranch(readname='jetAK8_Phi',size='jetAK8_size',dictlist=thisdictlist)
	ak8_es 		   = AddBranch(readname='jetAK8_E',size='jetAK8_size',dictlist=thisdictlist)
	ak8_smpts 	   = AddBranch(readname='jetAK8_SmearedPt',size='jetAK8_size',dictlist=thisdictlist)
	ak8_smetas 	   = AddBranch(readname='jetAK8_SmearedPEta',size='jetAK8_size',dictlist=thisdictlist)
	ak8_smphis 	   = AddBranch(readname='jetAK8_SmearedPhi',size='jetAK8_size',dictlist=thisdictlist)
	ak8_smes 	   = AddBranch(readname='jetAK8_SmearedE',size='jetAK8_size',dictlist=thisdictlist)
	ak8_csvv2s 	   = AddBranch(readname='jetAK8_CSVv2',size='jetAK8_size',dictlist=thisdictlist)
	ak8_jecuncs    = AddBranch(readname='jetAK8_jecUncertainty',size='jetAK8_size',dictlist=thisdictlist)
	ak8_JERSFs 	   = AddBranch(readname='jetAK8_JERSF',size='jetAK8_size',dictlist=thisdictlist)
	ak8_JERSFUps   = AddBranch(readname='jetAK8_JERSFUp',size='jetAK8_size',dictlist=thisdictlist)
	ak8_JERSFDowns = AddBranch(readname='jetAK8_JERSFDown',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau1s 	   = AddBranch(readname='jetAK8_tau1',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau2s 	   = AddBranch(readname='jetAK8_tau2',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau3s 	   = AddBranch(readname='jetAK8_tau3',size='jetAK8_size',dictlist=thisdictlist)
	ak8_sdms 	   = AddBranch(readname='jetAK8_softDropMass',size='jetAK8_size',dictlist=thisdictlist)
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
	sf_top_pT_new 	  = AddBranch(writename='sf_top_pT_new',inival=1.,dictlist=thisdictlist)
	sf_top_pT_new_low = AddBranch(writename='sf_top_pT_new_low',inival=1.,dictlist=thisdictlist)
	sf_top_pT_new_hi  = AddBranch(writename='sf_top_pT_new_hi',inival=1.,dictlist=thisdictlist)
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
	#leptons
	physobjectBranches = {}
	thisdictlist = [allBranches,physobjectBranches]
	Q_l 		  = AddBranch(writename='Q_l',ttreetype='I',inival=0,dictlist=thisdictlist)
	muon1_pt 	  = AddBranch(writename='muon1_pt',dictlist=thisdictlist)
	muon1_eta 	  = AddBranch(writename='muon1_eta',dictlist=thisdictlist)
	muon1_phi 	  = AddBranch(writename='muon1_phi',dictlist=thisdictlist)
	muon1_M 	  = AddBranch(writename='muon1_M',dictlist=thisdictlist) 
	muon1_Q 	  = AddBranch(writename='muon1_Q',ttreetype='I',inival=0,dictlist=thisdictlist)
	muon1_ID 	  = AddBranch(writename='muon1_ID',ttreetype='i',inival=2,dictlist=thisdictlist)
	muon1_relPt   = AddBranch(writename='muon1_relPt',dictlist=thisdictlist)
	muon1_dR 	  = AddBranch(writename='muon1_dR',dictlist=thisdictlist)
	muon2_pt 	  = AddBranch(writename='muon2_pt', dictlist=thisdictlist)
	muon2_eta 	  = AddBranch(writename='muon2_eta',dictlist=thisdictlist)
	muon2_phi 	  = AddBranch(writename='muon2_phi',dictlist=thisdictlist)
	muon2_M 	  = AddBranch(writename='muon2_M',  dictlist=thisdictlist)
	muon2_Q 	  = AddBranch(writename='muon2_Q',ttreetype='I',inival=0,dictlist=thisdictlist)
	muon2_ID 	  = AddBranch(writename='muon2_ID',ttreetype='i',inival=2,dictlist=thisdictlist)
	muon2_relPt   = AddBranch(writename='muon2_relPt',dictlist=thisdictlist)
	muon2_dR 	  = AddBranch(writename='muon2_dR',dictlist=thisdictlist)
	ele1_pt 	  = AddBranch(writename='ele1_pt',dictlist=thisdictlist)
	ele1_eta 	  = AddBranch(writename='ele1_eta',dictlist=thisdictlist)
	ele1_phi 	  = AddBranch(writename='ele1_phi',dictlist=thisdictlist)
	ele1_M 		  = AddBranch(writename='ele1_M',dictlist=thisdictlist) 
	ele1_Q 		  = AddBranch(writename='ele1_Q',ttreetype='I',inival=0,dictlist=thisdictlist)
	ele1_ID 	  = AddBranch(writename='ele1_ID',ttreetype='i',inival=2,dictlist=thisdictlist)
	ele1_relPt    = AddBranch(writename='ele1_relPt',dictlist=thisdictlist)
	ele1_dR 	  = AddBranch(writename='ele1_dR',dictlist=thisdictlist)
	ele2_pt 	  = AddBranch(writename='ele2_pt', dictlist=thisdictlist)
	ele2_eta 	  = AddBranch(writename='ele2_eta',dictlist=thisdictlist)
	ele2_phi 	  = AddBranch(writename='ele2_phi',dictlist=thisdictlist)
	ele2_M 		  = AddBranch(writename='ele2_M',  dictlist=thisdictlist)
	ele2_Q 		  = AddBranch(writename='ele2_Q',ttreetype='I',inival=0,dictlist=thisdictlist)
	ele2_ID 	  = AddBranch(writename='ele2_ID',ttreetype='i',inival=2,dictlist=thisdictlist)
	ele2_relPt    = AddBranch(writename='ele2_relPt',dictlist=thisdictlist)
	ele2_dR 	  = AddBranch(writename='ele2_dR',dictlist=thisdictlist)
	#neutrino
	met_pt 	= AddBranch(writename='met_pt',dictlist=thisdictlist)
	met_eta = AddBranch(writename='met_eta',dictlist=thisdictlist)
	met_phi = AddBranch(writename='met_phi',dictlist=thisdictlist)
	met_M   = AddBranch(writename='met_M',dictlist=thisdictlist)
	#leptonic b
	lepb_pt  = AddBranch(writename='lepb_pt',dictlist=thisdictlist)
	lepb_eta = AddBranch(writename='lepb_eta',dictlist=thisdictlist)
	lepb_phi = AddBranch(writename='lepb_phi',dictlist=thisdictlist)
	lepb_M   = AddBranch(writename='lepb_M',dictlist=thisdictlist)
	#hadronic top
	hadt_pt    = AddBranch(writename='hadt_pt',dictlist=thisdictlist)
	hadt_eta   = AddBranch(writename='hadt_eta',dictlist=thisdictlist)
	hadt_phi   = AddBranch(writename='hadt_phi',dictlist=thisdictlist)
	hadt_M 	   = AddBranch(writename='hadt_M',dictlist=thisdictlist)
	hadt_tau32 = AddBranch(writename='hadt_tau32',dictlist=thisdictlist)
	hadt_tau21 = AddBranch(writename='hadt_tau21',dictlist=thisdictlist)
	hadt_SDM   = AddBranch(writename='hadt_SDM',dictlist=thisdictlist)
	#rescaled fourvectors
	scaled_lep_pt 	  = AddBranch(writename='scaled_lep_pt',dictlist=thisdictlist)
	scaled_lep_eta 	  = AddBranch(writename='scaled_lep_eta',dictlist=thisdictlist)
	scaled_lep_phi 	  = AddBranch(writename='scaled_lep_phi',dictlist=thisdictlist)
	scaled_lep_M 	  = AddBranch(writename='scaled_lep_M',dictlist=thisdictlist) 
	scaled_met_pt 	= AddBranch(writename='scaled_met_pt',dictlist=thisdictlist)
	scaled_met_eta = AddBranch(writename='scaled_met_eta',dictlist=thisdictlist)
	scaled_met_phi = AddBranch(writename='scaled_met_phi',dictlist=thisdictlist)
	scaled_met_M   = AddBranch(writename='scaled_met_M',dictlist=thisdictlist)
	scaled_lepW_pt  = AddBranch(writename='scaled_lepW_pt',dictlist=thisdictlist)
	scaled_lepW_eta = AddBranch(writename='scaled_lepW_eta',dictlist=thisdictlist)
	scaled_lepW_phi = AddBranch(writename='scaled_lepW_phi',dictlist=thisdictlist)
	scaled_lepW_M   = AddBranch(writename='scaled_lepW_M',dictlist=thisdictlist)
	scaled_lepb_pt  = AddBranch(writename='scaled_lepb_pt',dictlist=thisdictlist)
	scaled_lepb_eta = AddBranch(writename='scaled_lepb_eta',dictlist=thisdictlist)
	scaled_lepb_phi = AddBranch(writename='scaled_lepb_phi',dictlist=thisdictlist)
	scaled_lepb_M   = AddBranch(writename='scaled_lepb_M',dictlist=thisdictlist)
	scaled_lept_pt  = AddBranch(writename='scaled_lept_pt',dictlist=thisdictlist)
	scaled_lept_eta = AddBranch(writename='scaled_lept_eta',dictlist=thisdictlist)
	scaled_lept_phi = AddBranch(writename='scaled_lept_phi',dictlist=thisdictlist)
	scaled_lept_M   = AddBranch(writename='scaled_lept_M',dictlist=thisdictlist)
	scaled_hadt_pt  = AddBranch(writename='scaled_hadt_pt',dictlist=thisdictlist)
	scaled_hadt_eta = AddBranch(writename='scaled_hadt_eta',dictlist=thisdictlist)
	scaled_hadt_phi = AddBranch(writename='scaled_hadt_phi',dictlist=thisdictlist)
	scaled_hadt_M   = AddBranch(writename='scaled_hadt_M',dictlist=thisdictlist)
	#kinematic fit stuff
	thisdictlist = [allBranches]
	chi2 = AddBranch(writename='chi2',inival=0.0,dictlist=thisdictlist)
	nFits = AddBranch(writename='nFits',ttreetype='i',inival=0,dictlist=thisdictlist)
	#whether or not this event should be added twice and have its weight halved based on whether its initial state
	#was symmetric (this will only be nonzero for qqbar and some gg events)
	addTwice = AddBranch(writename='addTwice',ttreetype='i',inival=0,dictlist=thisdictlist)
	#Obervables
	observableBranches = {}
	thisdictlist = [allBranches,observableBranches]
	#cosine(theta)
	cstar 		  = AddBranch(writename='cstar',dictlist=thisdictlist)
	cstar_scaled = AddBranch(writename='cstar_scaled',dictlist=thisdictlist)
	#Feynman x
	x_F 		= AddBranch(writename='x_F',dictlist=thisdictlist)
	x_F_scaled = AddBranch(writename='x_F_scaled',dictlist=thisdictlist)
	#ttbar invariant mass
	M 		  = AddBranch(writename='M',dictlist=thisdictlist)
	M_scaled = AddBranch(writename='M_scaled',dictlist=thisdictlist)
	#initial quark vector
	mctruthBranches = {}
	thisdictlist = [allBranches,mctruthBranches]
	q_pt  = AddBranch(writename='q_pt',dictlist=thisdictlist)
	q_eta = AddBranch(writename='q_eta',dictlist=thisdictlist)
	q_phi = AddBranch(writename='q_phi',dictlist=thisdictlist)
	q_M   = AddBranch(writename='q_M',dictlist=thisdictlist)
	#initial antiquark vector
	qbar_pt  = AddBranch(writename='qbar_pt',dictlist=thisdictlist)
	qbar_eta = AddBranch(writename='qbar_eta',dictlist=thisdictlist)
	qbar_phi = AddBranch(writename='qbar_phi',dictlist=thisdictlist)
	qbar_M   = AddBranch(writename='qbar_M',dictlist=thisdictlist)
	#MC top vector
	MCt_pt  = AddBranch(writename='MCt_pt',dictlist=thisdictlist)
	MCt_eta = AddBranch(writename='MCt_eta',dictlist=thisdictlist)
	MCt_phi = AddBranch(writename='MCt_phi',dictlist=thisdictlist)
	MCt_M   = AddBranch(writename='MCt_M',dictlist=thisdictlist)
	#MC antitop vector
	MCtbar_pt  = AddBranch(writename='MCtbar_pt',dictlist=thisdictlist)
	MCtbar_eta = AddBranch(writename='MCtbar_eta',dictlist=thisdictlist)
	MCtbar_phi = AddBranch(writename='MCtbar_phi',dictlist=thisdictlist)
	MCtbar_M   = AddBranch(writename='MCtbar_M',dictlist=thisdictlist)
	#MC truth observables
	cstar_MC = AddBranch(writename='cstar_MC',dictlist=thisdictlist)
	x_F_MC 	 = AddBranch(writename='x_F_MC',dictlist=thisdictlist)
	M_MC 	 = AddBranch(writename='M_MC',dictlist=thisdictlist)

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
			eventweight = self.genWeight.getReadValue()*(self.xsec.getReadValue()/self.totalweight)
			self.weight.setWriteValue(eventweight)
			#Mother particle and MC truth top assignment
			q_vec, qbar_vec = findInitialPartons(self.mcGenEventBranches,self.MC_generator)
			MCt_vec, MCtbar_vec = findMCTops(self.mcGenEventBranches,self.MC_generator)
			self.__setFourVectorBranchValues__('q',q_vec)
			self.__setFourVectorBranchValues__('qbar',qbar_vec)
			self.__setFourVectorBranchValues__('MCt',MCt_vec)
			self.__setFourVectorBranchValues__('MCtbar',MCtbar_vec)

		#For the record, trigger information is handled automatically

		#Copy over and populate objects
		#MET
		met = TLorentzVector(); met.SetPtEtaPhiM(self.met_pts.getReadValue(), 0., self.met_phis.getReadValue(), 0.)
		self.__setFourVectorBranchValues__('met',met)
		#jets
		ak4jets = []; ak8jets = [];
		for i in range(ak4jetsize) :
			ak4jets.append(AK4Jet(self.ak4JetBranches,i,self.JEC))
		for i in range(ak8jetsize) :
			ak8jets.append(AK8Jet(self.ak8JetBranches,i,self.JEC))
		ak4jets.sort(key=lambda x: x.getPt(), reverse=True)
		ak8jets.sort(key=lambda x: x.getPt(), reverse=True)
		#muons
		muons = []
		for i in range(musize) :
			muons.append(Muon(self.muonBranches,i,ak4jets))
		muons.sort(key=lambda x: x.getPt(), reverse=True)
		mulen = len(muons)
		if mulen>0 :
			self.__setLeptonBranchValues__('muon1',muons[0])
			if mulen>1 :
				self.__setLeptonBranchValues__('muon2',muons[1])
		#electrons
		electrons = []
		for i in range(elsize) :
			electrons.append(Electron(self.electronBranches,i,ak4jets))
		electrons.sort(key=lambda x: x.getPt(), reverse=True)
		ellen = len(electrons)
		if ellen>0 :
			self.__setLeptonBranchValues__('ele1',electrons[0])
			if ellen>1 :
				self.__setLeptonBranchValues__('ele2',electrons[1])

		#figure out whether the event is muonic or electronic, assign lep
		lep = None
		if len(muons)>0 and (len(electrons)==0 or muons[0].getPt()>electrons[0].getPt()) :
			lep = muons[0]
		elif len(electrons)>0 and (len(muons)==0 or electrons[0].getPt()>muons[0].getPt()) :
			lep = electrons[0]
		self.Q_l.setWriteValue(int(lep.getQ()))

		#neutrino handling and setup for fit
		met1_vec, met2_vec = setupMET(lep.getFourVector(),met)
		self.nFits.setWriteValue(2)
		if met1_vec.Pz() == met2_vec.Pz() :
			self.nFits.setWriteValue(1)

		#assign the top and b candidates
		#top is the hardest ak8 jet
		hadtCandJet = ak8jets[0]
		#leptonic b is the hardest ak4 jet on the opposite hemisphere
		lepbCandJet = None
		furthest_distance = 0.
		furthest_index = 0
		for i in range(len(ak4jets)) :
			ak4jet = ak4jets[i]
			dRcheck = ak4jet.getFourVector().DeltaR(hadtCandJet.getFourVector())
			if dRcheck>furthest_distance :
				furthest_distance=dRcheck
				furthest_index=i
			if dRcheck < pi :
				continue
			else :
				lepbCandJet = ak4jet
				break
		if lepbCandJet == None :
			print 'EVENT NUMBER %d NOT VALID; NO AK4 JET IN LEPTONIC HEMISPHERE (furthest jet is %.2f away, index %d of %d)'%(eventnumber,
				furthest_distance,furthest_index,len(ak4jets))
			return
		#and fill the uncorrected fourvectors and other attributes
		self.__setFourVectorBranchValues__('hadt',hadtCandJet.getFourVector())
		self.hadt_tau32.setWriteValue(hadtCandJet.getTau32())
		self.hadt_tau21.setWriteValue(hadtCandJet.getTau21())
		self.hadt_SDM.setWriteValue(hadtCandJet.getSDM())
		self.__setFourVectorBranchValues__('lepb',lepbCandJet.getFourVector())

		#event reconstruction with kinematic fit
		scaledlep = TLorentzVector(); scaledmet = TLorentzVector() 
		scaledlepb = TLorentzVector(); scaledhadt = TLorentzVector()
		scaledlep, scaledmet, scaledlepb, scaledhadt, fitchi2 = reconstruct(lep.getFourVector(),met1_vec,met2_vec,lepbCandJet.getFourVector(),hadtCandJet.getFourVector())
		if scaledlep==None :
			print 'EVENT NUMBER '+str(eventnumber)+' NOT VALID; NEITHER KINEMATIC FIT CONVERGED'
			return
		self.chi2.setWriteValue(fitchi2)

		#fill the TTree with the scaled fourvector variables
		self.__setFourVectorBranchValues__('scaled_lep',scaledlep)
		self.__setFourVectorBranchValues__('scaled_met',scaledmet)
		self.__setFourVectorBranchValues__('scaled_lepW',scaledlep+scaledmet)
		self.__setFourVectorBranchValues__('scaled_lepb',scaledlepb)
		self.__setFourVectorBranchValues__('scaled_lept',scaledlep+scaledmet+scaledlepb)
		self.__setFourVectorBranchValues__('scaled_hadt',scaledhadt)

#		#reconstruct the observables using both the scaled and unscaled vectors
#		if self.lep_type == 0 :
#			self.cstar[0], self.x_F[0], self.M[0] = getObservables(muons[0].vec+met1_vec+jets[1].vec,jets[0].vec,self.Q_l[0]) 
#		elif self.lep_type == 1 :
#			self.cstar[0], self.x_F[0], self.M[0] = getObservables(electrons[0].vec+met1_vec+jets[1].vec,jets[0].vec,self.Q_l[0]) 
#		( self.cstar_scaled[0], self.x_F_scaled[0], 
#			self.M_scaled[0] ) = getObservables(scaledlep+scaledmet+scaledlepb,scaledhadt,self.Q_l[0]) 
#		#above function in angleReconstructor.py

#		#MC Truth observable and reweighting calculation
#		if self.is_data==0 :
#			if self.event_type!=4 :
#				( self.cstar_MC[0],self.x_F_MC[0],self.M_MC[0],
#					self.wg1[0],self.wg2[0],self.wg3[0],self.wg4[0],
#					self.wqs1[0],self.wqs2[0],self.wqa0[0],self.wqa1[0],self.wqa2[0],
#					self.wg1_opp[0],self.wg2_opp[0],self.wg3_opp[0],self.wg4_opp[0],
#					self.wqs1_opp[0],self.wqs2_opp[0],
#					self.wqa0_opp[0],self.wqa1_opp[0],self.wqa2_opp[0],
#					self.wega[0], self.wegc[0] ) = getMCObservables(q_vec,qbar_vec,MCt_vec,MCtbar_vec,self.event_type) 
#			#scale factor and reweighting calculations
#			if self.lep_type==0 :
#				meas_lep_pt=muons[0].vec.Pt(); meas_lep_eta=muons[0].vec.Eta()
#			elif self.lep_type==1 :
#				meas_lep_pt=electrons[0].vec.Pt(); meas_lep_eta=electrons[0].vec.Eta()
#			#8TeV numbers
#			self.sf_top_pT[0], self.sf_top_pT_low[0], self.sf_top_pT_hi[0] = self.corrector.getToppT_reweight(MCt_vec,MCtbar_vec,self.Q_l[0])
#			self.sf_top_pT_new[0], self.sf_top_pT_new_low[0], self.sf_top_pT_new_hi[0] = self.corrector.getNewToppT_reweight(MCt_vec,MCtbar_vec)
#			self.sf_pileup[0], self.sf_pileup_low[0], self.sf_pileup_hi[0] = self.corrector.getpileup_reweight(MCpileup)
#			( self.sf_lep_ID[0], self.sf_lep_ID_low[0], 
#				self.sf_lep_ID_hi[0] ) = self.corrector.getID_eff(pileup,meas_lep_pt,meas_lep_eta,self.lep_type)
#			( self.sf_trig_eff[0], self.sf_trig_eff_low[0], 
#				self.sf_trig_eff_hi[0] ) = self.corrector.gettrig_eff(pileup,meas_lep_pt,meas_lep_eta,self.lep_type)

		self.__closeout__() #yay! A complete event!

	##################################  #__init__ function  ##################################
	def __init__(self,fileName,tree,isData,generator,jec,onGrid,totweight) :
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
		self.JEC = jec
		#Set the total weight
		self.totalweight = totweight
		#Set Monte Carlo reweighter
		self.corrector = MC_corrector(self.MC_generator,self.event_type,onGrid)

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
		self.tree.Fill()

	################################## __del__ function  ###################################
	def __del__(self) :
		self.outfile.cd()
		self.tree.Write()
		self.outfile.Write()
		self.outfile.Close()