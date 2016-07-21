#Global variables
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0
#Trigger paths
MU_TRIG_PATH = 'HLT_Mu30_eta2p1_PFJet150_PFJet50'
EL_TRIG_PATH = 'HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50'

##########								   Imports  								##########

from ROOT import TFile
from branch import Branch

##########							   Treemaker Class 								##########

class Reconstructor(object) :

	##################################  Branches  ##################################
	allBranches = {}
	#GenWeight info
	genWeight = self.__AddBranch__('evt_Gen_Weight','genWeight','F',1.0,'1',[allBranches])
	#MC GenEvent info
	mcGenEventBranches = {}
	thisdictlist = [allBranches,mcGenEventBranches]
	gen_size 	   = self.__AddBranch__(readname='gen_size',ttreetype='i',dictlist=thisdictlist)
	gen_ID 		   = self.__AddBranch__(readname='gen_ID',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Status 	   = self.__AddBranch__(readname='gen_Status',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Mom0ID 	   = self.__AddBranch__(readname='gen_Mom0ID',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Mom0Status = self.__AddBranch__(readname='gen_Mom0Status',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Mom1ID 	   = self.__AddBranch__(readname='gen_Mom1ID',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Mom1Status = self.__AddBranch__(readname='gen_Mom1Status',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Dau0ID 	   = self.__AddBranch__(readname='gen_Dau0ID',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Dau0Status = self.__AddBranch__(readname='gen_Dau0Status',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Dau1ID 	   = self.__AddBranch__(readname='gen_Dau1ID',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	gen_Dau1Status = self.__AddBranch__(readname='gen_Dau1Status',ttreetype='I',,size='gen_size',dictlist=thisdictlist)
	#Trigger Information
	triggerBranches = {}
	thisdictlist = [allBranches,triggerBranches]
	muTrig = self.__AddBranch__(MU_TRIG_PATH,'muTrig','I',-1,'1',thisdictlist)
	elTrig = self.__AddBranch__(EL_TRIG_PATH,'elTrig','I',-1,'1',thisdictlist)
	#pileup
	pileupBranches = {}
	thisdictlist = [allBranches,pileupBranches]
	pu = self.__AddBranch__(readname='pu_NtrueInt',ttreetype='I',dictlist=thisdictlist)
	#MET
	metBranches = {}
	thisdictlist = [allBranches,metBranches]
	met_size 	= self.__AddBranch__(readname='met_size',ttreetype='i',dictlist=thisdictlist)
	met_pts 	= self.__AddBranch__(readname='met_Pt',size='met_size',dictlist=thisdictlist)
	met_pxs 	= self.__AddBranch__(readname='met_Px',size='met_size',dictlist=thisdictlist)
	met_pys 	= self.__AddBranch__(readname='met_Py',size='met_size',dictlist=thisdictlist)
	met_phis 	= self.__AddBranch__(readname='met_Phi',size='met_size',dictlist=thisdictlist)
	#muons
	muonBranches = {}
	thisdictlist = [allBranches,muonBranches]
	mu_size 	= self.__AddBranch__(readname='mu_size',ttreetype='i',dictlist=thisdictlist)
	mu_pts 		= self.__AddBranch__(readname='mu_Pt',size='mu_size',dictlist=thisdictlist)
	mu_etas 	= self.__AddBranch__(readname='mu_Eta',size='mu_size',dictlist=thisdictlist)
	mu_phis 	= self.__AddBranch__(readname='mu_Phi',size='mu_size',dictlist=thisdictlist)
	mu_es 		= self.__AddBranch__(readname='mu_E',size='mu_size',dictlist=thisdictlist)
	mu_charges 	= self.__AddBranch__(readname='mu_Charge',size='mu_size',dictlist=thisdictlist)
	mus_isLoose = self.__AddBranch__(readname='mu_IsLooseMuon',size='mu_size',dictlist=thisdictlist)
	#electrons
	electronBranches = {}
	thisdictlist = [allBranches,electronBranches]
	el_size 	= self.__AddBranch__(readname='el_size',ttreetype='i',dictlist=thisdictlist)
	el_pts 		= self.__AddBranch__(readname='el_Pt',size='el_size',dictlist=thisdictlist)
	el_etas 	= self.__AddBranch__(readname='el_Eta',size='el_size',dictlist=thisdictlist)
	el_phis 	= self.__AddBranch__(readname='el_Phi',size='el_size',dictlist=thisdictlist)
	el_es 		= self.__AddBranch__(readname='el_E',size='el_size',dictlist=thisdictlist)
	el_charges 	= self.__AddBranch__(readname='el_Charge',size='el_size',dictlist=thisdictlist)
	el_detas_in = self.__AddBranch__(readname='el_dEtaIn',size='el_size',dictlist=thisdictlist)
	el_dphis_in = self.__AddBranch__(readname='el_dPhiIn',size='el_size',dictlist=thisdictlist)
	el_5x5siees = self.__AddBranch__(readname='el_full5x5siee',size='el_size',dictlist=thisdictlist)
	el_HoEs 	= self.__AddBranch__(readname='el_HoE',size='el_size',dictlist=thisdictlist)
	el_isEBs 	= self.__AddBranch__(readname='el_isEB',size='el_size',dictlist=thisdictlist)
	els_isLoose = self.__AddBranch__(readname='el_isLoose',size='el_size',dictlist=thisdictlist)
	el_D0s 		= self.__AddBranch__(readname='el_D0',size='el_size',dictlist=thisdictlist)
	el_Dzs 		= self.__AddBranch__(readname='el_Dz',size='el_size',dictlist=thisdictlist)
	el_ooEmooPs = self.__AddBranch__(readname='el_ooEmooP',size='el_size',dictlist=thisdictlist)
	el_hasMCVs 	= self.__AddBranch__(readname='el_hasMatchedConVeto',size='el_size',dictlist=thisdictlist)
	el_misshits = self.__AddBranch__(readname='el_missHits',size='el_size',dictlist=thisdictlist)
	#AK4 Jets
	ak4JetBranches = {}
	thisdictlist = [allBranches,ak4JetBranches]
	ak4_size 	   = self.__AddBranch__(readname='jetAK4_size',ttreetype='i',dictlist=thisdictlist)
	ak4_pts 	   = self.__AddBranch__(readname='jetAK4_Pt',size='jetAK4_size',dictlist=thisdictlist)
	ak4_etas 	   = self.__AddBranch__(readname='jetAK4_Eta',size='jetAK4_size',dictlist=thisdictlist)
	ak4_phis 	   = self.__AddBranch__(readname='jetAK4_Phi',size='jetAK4_size',dictlist=thisdictlist)
	ak4_es 		   = self.__AddBranch__(readname='jetAK4_E',size='jetAK4_size',dictlist=thisdictlist)
	ak4_smpts 	   = self.__AddBranch__(readname='jetAK4_SmearedPt',size='jetAK4_size',dictlist=thisdictlist)
	ak4_smetas 	   = self.__AddBranch__(readname='jetAK4_SmearedPEta',size='jetAK4_size',dictlist=thisdictlist)
	ak4_smphis 	   = self.__AddBranch__(readname='jetAK4_SmearedPhi',size='jetAK4_size',dictlist=thisdictlist)
	ak4_smes 	   = self.__AddBranch__(readname='jetAK4_SmearedE',size='jetAK4_size',dictlist=thisdictlist)
	ak4_csvv2s 	   = self.__AddBranch__(readname='jetAK4_CSVv2',size='jetAK4_size',dictlist=thisdictlist)
	ak4_jecuncs    = self.__AddBranch__(readname='jetAK4_jecUncertainty',size='jetAK4_size',dictlist=thisdictlist)
	ak4_JERSFs 	   = self.__AddBranch__(readname='jetAK4_JERSF',size='jetAK4_size',dictlist=thisdictlist)
	ak4_JERSFUps   = self.__AddBranch__(readname='jetAK4_JERSFUp',size='jetAK4_size',dictlist=thisdictlist)
	ak4_JERSFDowns = self.__AddBranch__(readname='jetAK4_JERSFDown',size='jetAK4_size',dictlist=thisdictlist)
	#AK8 Jets
	ak8JetBranches = {}
	thisdictlist = [allBranches,ak8JetBranches]
	ak8_size 	   = self.__AddBranch__(readname='jetAK8_size',ttreetype='i',dictlist=thisdictlist)
	ak8_pts 	   = self.__AddBranch__(readname='jetAK8_Pt',size='jetAK8_size',dictlist=thisdictlist)
	ak8_etas 	   = self.__AddBranch__(readname='jetAK8_Eta',size='jetAK8_size',dictlist=thisdictlist)
	ak8_phis 	   = self.__AddBranch__(readname='jetAK8_Phi',size='jetAK8_size',dictlist=thisdictlist)
	ak8_es 		   = self.__AddBranch__(readname='jetAK8_E',size='jetAK8_size',dictlist=thisdictlist)
	ak8_smpts 	   = self.__AddBranch__(readname='jetAK8_SmearedPt',size='jetAK8_size',dictlist=thisdictlist)
	ak8_smetas 	   = self.__AddBranch__(readname='jetAK8_SmearedPEta',size='jetAK8_size',dictlist=thisdictlist)
	ak8_smphis 	   = self.__AddBranch__(readname='jetAK8_SmearedPhi',size='jetAK8_size',dictlist=thisdictlist)
	ak8_smes 	   = self.__AddBranch__(readname='jetAK8_SmearedE',size='jetAK8_size',dictlist=thisdictlist)
	ak8_csvv2s 	   = self.__AddBranch__(readname='jetAK8_CSVv2',size='jetAK8_size',dictlist=thisdictlist)
	ak8_jecuncs    = self.__AddBranch__(readname='jetAK8_jecUncertainty',size='jetAK8_size',dictlist=thisdictlist)
	ak8_JERSFs 	   = self.__AddBranch__(readname='jetAK8_JERSF',size='jetAK8_size',dictlist=thisdictlist)
	ak8_JERSFUps   = self.__AddBranch__(readname='jetAK8_JERSFUp',size='jetAK8_size',dictlist=thisdictlist)
	ak8_JERSFDowns = self.__AddBranch__(readname='jetAK8_JERSFDown',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau1s 	   = self.__AddBranch__(readname='jetAK8_tau1',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau2s 	   = self.__AddBranch__(readname='jetAK8_tau2',size='jetAK8_size',dictlist=thisdictlist)
	ak8_tau3s 	   = self.__AddBranch__(readname='jetAK8_tau3',size='jetAK8_size',dictlist=thisdictlist)
	ak8_sdms 	   = self.__AddBranch__(readname='jetAK8_softDropMass',size='jetAK8_size',dictlist=thisdictlist)
	#PRODUCES
	#Weights
	weightBranches = {}
	thisdictlist = [allBranches,weightBranches]
	wg1 	 = self.__AddBranch__(writename='wg1',inival=1.,dictlist=thisdictlist)
	wg2 	 = self.__AddBranch__(writename='wg2',inival=1.,dictlist=thisdictlist)
	wg3 	 = self.__AddBranch__(writename='wg3',inival=1.,dictlist=thisdictlist)
	wg4 	 = self.__AddBranch__(writename='wg4',inival=1.,dictlist=thisdictlist)
	wqs1 	 = self.__AddBranch__(writename='wqs1',inival=1.,dictlist=thisdictlist)
	wqs2 	 = self.__AddBranch__(writename='wqs2',inival=1.,dictlist=thisdictlist)
	wqa0 	 = self.__AddBranch__(writename='wqa0',inival=1.,dictlist=thisdictlist)
	wqa1 	 = self.__AddBranch__(writename='wqa1',inival=1.,dictlist=thisdictlist)
	wqa2 	 = self.__AddBranch__(writename='wqa2',inival=1.,dictlist=thisdictlist)
	wg1_opp  = self.__AddBranch__(writename='wg1_opp',inival=1.,dictlist=thisdictlist)
	wg2_opp  = self.__AddBranch__(writename='wg2_opp',inival=1.,dictlist=thisdictlist)
	wg3_opp  = self.__AddBranch__(writename='wg3_opp',inival=1.,dictlist=thisdictlist)
	wg4_opp  = self.__AddBranch__(writename='wg4_opp',inival=1.,dictlist=thisdictlist)
	wqs1_opp = self.__AddBranch__(writename='wqs1_opp',inival=1.,dictlist=thisdictlist)
	wqs2_opp = self.__AddBranch__(writename='wqs2_opp',inival=1.,dictlist=thisdictlist)
	wqa0_opp = self.__AddBranch__(writename='wqa0_opp',inival=1.,dictlist=thisdictlist)
	wqa1_opp = self.__AddBranch__(writename='wqa1_opp',inival=1.,dictlist=thisdictlist)
	wqa2_opp = self.__AddBranch__(writename='wqa2_opp',inival=1.,dictlist=thisdictlist)
	wega 	 = self.__AddBranch__(writename='wega',inival=1.,dictlist=thisdictlist)
	wegc 	 = self.__AddBranch__(writename='wegc',inival=1.,dictlist=thisdictlist)
	#scalefactors
	scalefactorBranches = {}
	thisdictlist = [allBranches,scalefactorBranches]
	sf_top_pT 		  = self.__AddBranch__(writename='sf_top_pT',inival=1.,dictlist=thisdictlist)
	sf_top_pT_low 	  = self.__AddBranch__(writename='sf_top_pT_low',inival=1.,dictlist=thisdictlist)
	sf_top_pT_hi 	  = self.__AddBranch__(writename='sf_top_pT_hi',inival=1.,dictlist=thisdictlist)
	sf_top_pT_new 	  = self.__AddBranch__(writename='sf_top_pT_new',inival=1.,dictlist=thisdictlist)
	sf_top_pT_new_low = self.__AddBranch__(writename='sf_top_pT_new_low',inival=1.,dictlist=thisdictlist)
	sf_top_pT_new_hi  = self.__AddBranch__(writename='sf_top_pT_new_hi',inival=1.,dictlist=thisdictlist)
	sf_btag_eff 	  = self.__AddBranch__(writename='sf_btag_eff',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_low   = self.__AddBranch__(writename='sf_btag_eff_low',inival=1.,dictlist=thisdictlist)
	sf_btag_eff_hi 	  = self.__AddBranch__(writename='sf_btag_eff_hi',inival=1.,dictlist=thisdictlist)
	sf_pileup 		  = self.__AddBranch__(writename='sf_pileup',inival=1.,dictlist=thisdictlist)
	sf_pileup_low 	  = self.__AddBranch__(writename='sf_pileup_low',inival=1.,dictlist=thisdictlist)
	sf_pileup_hi 	  = self.__AddBranch__(writename='sf_pileup_hi',inival=1.,dictlist=thisdictlist)
	sf_lep_ID 		  = self.__AddBranch__(writename='sf_lep_ID',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_low 	  = self.__AddBranch__(writename='sf_lep_ID_low',inival=1.,dictlist=thisdictlist)
	sf_lep_ID_hi 	  = self.__AddBranch__(writename='sf_lep_ID_hi',inival=1.,dictlist=thisdictlist)
	sf_trig_eff 	  = self.__AddBranch__(writename='sf_trig_eff',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_low   = self.__AddBranch__(writename='sf_trig_eff_low',inival=1.,dictlist=thisdictlist)
	sf_trig_eff_hi 	  = self.__AddBranch__(writename='sf_trig_eff_hi',inival=1.,dictlist=thisdictlist)
	#leptons
	physobjectBranches = {}
	thisdictlist = [allBranches,physobjectBranches]
	self.Q_l 	  = self.__AddBranch__(writename='Q_l',ttreetype='I',inival=0,dictlist=thisdictlist)
	muon1_pt 	  = self.__AddBranch__(writename='muon1_pt',dictlist=thisdictlist)
	muon1_eta 	  = self.__AddBranch__(writename='muon1_eta',dictlist=thisdictlist)
	muon1_phi 	  = self.__AddBranch__(writename='muon1_phi',dictlist=thisdictlist)
	muon1_M 	  = self.__AddBranch__(writename='muon1_M',dictlist=thisdictlist) 
	muon1_Q 	  = self.__AddBranch__(writename='muon1_Q',ttreetpye='I',inival=0,dictlist=thisdictlist)
	muon1_ID 	  = self.__AddBranch__(writename='muon1_ID',ttreetype='i',inival=2,dictlist=thisdictlist)
	muon1_isLoose = self.__AddBranch__(writename='muon1_isLoose',ttreetype='i',inival=2,dictlist=thisdictlist)
	muon1_relPt   = self.__AddBranch__(writename='muon1_relPt',dictlist=thisdictlist)
	muon1_dR 	  = self.__AddBranch__(writename='muon1_dR',dictlist=thisdictlist)
	muon2_pt 	  = self.__AddBranch__(writename='muon2_pt', dictlist=thisdictlist)
	muon2_eta 	  = self.__AddBranch__(writename='muon2_eta',dictlist=thisdictlist)
	muon2_phi 	  = self.__AddBranch__(writename='muon2_phi',dictlist=thisdictlist)
	muon2_M 	  = self.__AddBranch__(writename='muon2_M',  dictlist=thisdictlist)
	muon2_Q 	  = self.__AddBranch__(writename='muon2_Q',ttreetype='I',inival=0,dictlist=thisdictlist)
	muon2_isTight = self.__AddBranch__(writename='muon2_isTight',ttreetype='i',inival=2,dictlist=thisdictlist)
	muon2_isLoose = self.__AddBranch__(writename='muon2_isLoose',ttreetype='i',inival=2,dictlist=thisdictlist)
	muon2_relPt   = self.__AddBranch__(writename='muon2_relPt',dictlist=thisdictlist)
	muon2_dR 	  = self.__AddBranch__(writename='muon2_dR',dictlist=thisdictlist)
	ele1_pt 	  = self.__AddBranch__(writename='ele1_pt',dictlist=thisdictlist)
	ele1_eta 	  = self.__AddBranch__(writename='ele1_eta',dictlist=thisdictlist)
	ele1_phi 	  = self.__AddBranch__(writename='ele1_phi',dictlist=thisdictlist)
	ele1_M 		  = self.__AddBranch__(writename='ele1_M',dictlist=thisdictlist) 
	ele1_Q 		  = self.__AddBranch__(writename='ele1_Q',ttreetpye='I',inival=0,dictlist=thisdictlist)
	ele1_ID 	  = self.__AddBranch__(writename='ele1_ID',ttreetype='i',inival=2,dictlist=thisdictlist)
	ele1_isLoose  = self.__AddBranch__(writename='ele1_isLoose',ttreetype='i',inival=2,dictlist=thisdictlist)
	ele1_relPt    = self.__AddBranch__(writename='ele1_relPt',dictlist=thisdictlist)
	ele1_dR 	  = self.__AddBranch__(writename='ele1_dR',dictlist=thisdictlist)
	ele2_pt 	  = self.__AddBranch__(writename='ele2_pt', dictlist=thisdictlist)
	ele2_eta 	  = self.__AddBranch__(writename='ele2_eta',dictlist=thisdictlist)
	ele2_phi 	  = self.__AddBranch__(writename='ele2_phi',dictlist=thisdictlist)
	ele2_M 		  = self.__AddBranch__(writename='ele2_M',  dictlist=thisdictlist)
	ele2_Q 		  = self.__AddBranch__(writename='ele2_Q',ttreetype='I',inival=0,dictlist=thisdictlist)
	ele2_isTight  = self.__AddBranch__(writename='ele2_isTight',ttreetype='i',inival=2,dictlist=thisdictlist)
	ele2_isLoose  = self.__AddBranch__(writename='ele2_isLoose',ttreetype='i',inival=2,dictlist=thisdictlist)
	ele2_relPt    = self.__AddBranch__(writename='ele2_relPt',dictlist=thisdictlist)
	ele2_dR 	  = self.__AddBranch__(writename='ele2_dR',dictlist=thisdictlist)
	#neutrino
	met_pt 	= self.__AddBranch__(writename='met_pt',dictlist=thisdictlist)
	met_eta = self.__AddBranch__(writename='met_eta',dictlist=thisdictlist)
	met_phi = self.__AddBranch__(writename='met_phi',dictlist=thisdictlist)
	met_M   = self.__AddBranch__(writename='met_M',dictlist=thisdictlist)
	#leptonic W
	lepW_pt  = self.__AddBranch__(writename='lepW_pt',dictlist=thisdictlist)
	lepW_eta = self.__AddBranch__(writename='lepW_eta',dictlist=thisdictlist)
	lepW_phi = self.__AddBranch__(writename='lepW_phi',dictlist=thisdictlist)
	lepW_M   = self.__AddBranch__(writename='lepW_M',dictlist=thisdictlist)
	#leptonic b
	lepb_pt  = self.__AddBranch__(writename='lepb_pt',dictlist=thisdictlist)
	lepb_eta = self.__AddBranch__(writename='lepb_eta',dictlist=thisdictlist)
	lepb_phi = self.__AddBranch__(writename='lepb_phi',dictlist=thisdictlist)
	lepb_M   = self.__AddBranch__(writename='lepb_M',dictlist=thisdictlist)
	#leptonic top
	lept_pt  = self.__AddBranch__(writename='lept_pt',dictlist=thisdictlist)
	lept_eta = self.__AddBranch__(writename='lept_eta',dictlist=thisdictlist)
	lept_phi = self.__AddBranch__(writename='lept_phi',dictlist=thisdictlist)
	lept_M   = self.__AddBranch__(writename='lept_M',dictlist=thisdictlist)
	#hadronic top
	hadt_pt  = self.__AddBranch__(writename='hadt_pt',dictlist=thisdictlist)
	hadt_eta = self.__AddBranch__(writename='hadt_eta',dictlist=thisdictlist)
	hadt_phi = self.__AddBranch__(writename='hadt_phi',dictlist=thisdictlist)
	hadt_M   = self.__AddBranch__(writename='hadt_M',dictlist=thisdictlist)
	hadt_tau32  = self.__AddBranch__(writename='hadt_tau32',dictlist=thisdictlist)
	hadt_tau21  = self.__AddBranch__(writename='hadt_tau21',dictlist=thisdictlist)
	#rescaled fourvectors
	scaled_lep_pt 	  = self.__AddBranch__(writename='scaled_lep_pt',dictlist=thisdictlist)
	scaled_lep_eta 	  = self.__AddBranch__(writename='scaled_lep_eta',dictlist=thisdictlist)
	scaled_lep_phi 	  = self.__AddBranch__(writename='scaled_lep_phi',dictlist=thisdictlist)
	scaled_lep_M 	  = self.__AddBranch__(writename='scaled_lep_M',dictlist=thisdictlist) 
	scaled_met_pt 	= self.__AddBranch__(writename='scaled_met_pt',dictlist=thisdictlist)
	scaled_met_eta = self.__AddBranch__(writename='scaled_met_eta',dictlist=thisdictlist)
	scaled_met_phi = self.__AddBranch__(writename='scaled_met_phi',dictlist=thisdictlist)
	scaled_met_M   = self.__AddBranch__(writename='scaled_met_M',dictlist=thisdictlist)
	scaled_lepW_pt  = self.__AddBranch__(writename='scaled_lepW_pt',dictlist=thisdictlist)
	scaled_lepW_eta = self.__AddBranch__(writename='scaled_lepW_eta',dictlist=thisdictlist)
	scaled_lepW_phi = self.__AddBranch__(writename='scaled_lepW_phi',dictlist=thisdictlist)
	scaled_lepW_M   = self.__AddBranch__(writename='scaled_lepW_M',dictlist=thisdictlist)
	scaled_lepb_pt  = self.__AddBranch__(writename='scaled_lepb_pt',dictlist=thisdictlist)
	scaled_lepb_eta = self.__AddBranch__(writename='scaled_lepb_eta',dictlist=thisdictlist)
	scaled_lepb_phi = self.__AddBranch__(writename='scaled_lepb_phi',dictlist=thisdictlist)
	scaled_lepb_M   = self.__AddBranch__(writename='scaled_lepb_M',dictlist=thisdictlist)
	scaled_lept_pt  = self.__AddBranch__(writename='scaled_lept_pt',dictlist=thisdictlist)
	scaled_lept_eta = self.__AddBranch__(writename='scaled_lept_eta',dictlist=thisdictlist)
	scaled_lept_phi = self.__AddBranch__(writename='scaled_lept_phi',dictlist=thisdictlist)
	scaled_lept_M   = self.__AddBranch__(writename='scaled_lept_M',dictlist=thisdictlist)
	scaled_hadt_pt  = self.__AddBranch__(writename='scaled_hadt_pt',dictlist=thisdictlist)
	scaled_hadt_eta = self.__AddBranch__(writename='scaled_hadt_eta',dictlist=thisdictlist)
	scaled_hadt_phi = self.__AddBranch__(writename='scaled_hadt_phi',dictlist=thisdictlist)
	scaled_hadt_M   = self.__AddBranch__(writename='scaled_hadt_M',dictlist=thisdictlist)
	#kinematic fit stuff
	thisdictlist = [allBranches]
	chi2 = self.__AddBranch__(writename='chi2',inival=0.0,dictlist=thisdictlist)
	nFits = self.__AddBranch__('nFits',self.nFits,'/i',0)
	#whether or not this event should be added twice and have its weight halved based on whether its initial state
	#was symmetric (this will only be nonzero for qqbar and some gg events)
	addTwice = self.__AddBranch__(writename='addTwice',ttreetype='i',inival=0,dictlist=thisdictlist)
	#Obervables
	observableBranches = {}
	thisdictlist = [allBranches,observableBranches]
	#cosine(theta)
	cstar 		  = self.__AddBranch__(writename='cstar',dictlist=thisdictlist)
	cstar_scaled = self.__AddBranch__(writename='cstar_scaled',dictlist=thisdictlist)
	#Feynman x
	x_F 		= self.__AddBranch__(writename='x_F',dictlist=thisdictlist)
	x_F_scaled = self.__AddBranch__(writename='x_F_scaled',dictlist=thisdictlist)
	#ttbar invariant mass
	M 		  = self.__AddBranch__(writename='M',dictlist=thisdictlist)
	M_scaled = self.__AddBranch__(writename='M_scaled',dictlist=thisdictlist)
	#initial quark vector
	mctruthBranches = {}
	thisdictlist = [allBranches,mctruthBranches]
	q_pt  = self.__AddBranch__(writename='q_pt',dictlist=thisdictlist)
	q_eta = self.__AddBranch__(writename='q_eta',dictlist=thisdictlist)
	q_phi = self.__AddBranch__(writename='q_phi',dictlist=thisdictlist)
	q_M   = self.__AddBranch__(writename='q_M',dictlist=thisdictlist)
	#initial antiquark vector
	qbar_pt  = self.__AddBranch__(writename='qbar_pt',dictlist=thisdictlist)
	qbar_eta = self.__AddBranch__(writename='qbar_eta',dictlist=thisdictlist)
	qbar_phi = self.__AddBranch__(writename='qbar_phi',dictlist=thisdictlist)
	qbar_M   = self.__AddBranch__(writename='qbar_M',dictlist=thisdictlist)
	#MC top vector
	MCt_pt  = self.__AddBranch__(writename='MCt_pt',dictlist=thisdictlist)
	MCt_eta = self.__AddBranch__(writename='MCt_eta',dictlist=thisdictlist)
	MCt_phi = self.__AddBranch__(writename='MCt_phi',dictlist=thisdictlist)
	MCt_M   = self.__AddBranch__(writename='MCt_M',dictlist=thisdictlist)
	#MC antitop vector
	MCtbar_pt  = self.__AddBranch__(writename='MCtbar_pt',dictlist=thisdictlist)
	MCtbar_eta = self.__AddBranch__(writename='MCtbar_eta',dictlist=thisdictlist)
	MCtbar_phi = self.__AddBranch__(writename='MCtbar_phi',dictlist=thisdictlist)
	MCtbar_M   = self.__AddBranch__(writename='MCtbar_M',dictlist=thisdictlist)
	#MC truth observables
	cstar_MC 	  = self.__AddBranch__(writename='cstar_MC',dictlist=thisdictlist)
	x_F_MC 	= self.__AddBranch__(writename='x_F_MC',dictlist=thisdictlist)
	M_MC 	  = self.__AddBranch__(writename='M_MC',dictlist=thisdictlist)

	##################################  ANALYZE FUNCTION  ##################################
	def analyze(self,eventnumber) :
		#keep track of whether event has been cut
		keepEvent = True
		#event type split
		if self.is_data == 0 :
			#GenParticles
			event.getByLabel(self.genLabel,self.genHandle)
			if not self.genHandle.isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			GenParticles = self.genHandle.product()
			#pythia8 nTuple genParticles
			genPartVars = []
			for i in range(len(self.genPartHandles)) :
				event.getByLabel(self.genPartLabels[i],self.genPartHandles[i])
				if not self.genPartHandles[i].isValid() :
					self.ERR_CODE = ERR_INVALID_HANDLE
					return self.ERR_CODE
				genPartVars.append(self.genPartHandles[i].product())
			if self.event_type != 4 :
				keepEvent,add_twice = eventTypeCheck(self.MC_generator,GenParticles,genPartVars,self.event_type) 
							#abovefunction in eventTypeHelper.py
				if add_twice :
					self.addTwice[0] = 1
		if not keepEvent :
			return self.ERR_CODE

		#if it's mcatnlo, check the sign of the event weight
		if self.MC_generator == 'mcatnlo' :
			event.getByLabel(self.GenEventLabel,self.GenEventHandle) 
			if not self.GenEventHandle.isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			GenEvent = self.GenEventHandle.product()
			if GenEvent.weight()<0 :
				self.weight[0]*=-1

		#Trigger information
		event.getByLabel(self.trigLabel,self.trigHandle)
		if not self.trigHandle.isValid() :
			self.ERR_CODE = ERR_INVALID_HANDLE
			return self.ERR_CODE
		trigResults = self.trigHandle.product()
		trigNames = event.object().triggerNames(trigResults)
		for i in range(trigResults.size()) :
			s = str(trigNames.triggerName(i))
			if s.startswith(MU_TRIG_PATH) :
				if trigResults.accept(i) :
					self.mu_trigger[0] = 1
				else :
					self.mu_trigger[0] = 0
				if self.el_trigger[0] != 2 :
					break
			elif s.startswith(EL_TRIG_PATH) :
				if trigResults.accept(i) :
					self.el_trigger[0] = 1
				else :
					self.el_trigger[0] = 0
				if self.mu_trigger[0] != 2 :
					break

		#PDF information
		if self.is_data == 0 :
			event.getByLabel(self.CT10Label,self.CT10Handle)
			event.getByLabel(self.cteqLabel,self.cteqHandle)
			event.getByLabel(self.GJRLabel,self.GJRHandle)
			CT10ws = self.CT10Handle.product()
			cteqws = self.cteqHandle.product()
			GJRws  = self.GJRHandle.product()
			for i in range(len(CT10ws)) :
				self.CT10_weights[i] = CT10ws[i]/CT10ws[0]
			for i in range(len(cteqws)) :
				self.cteq_weights[i] = cteqws[i]/cteqws[0]
			for i in range(len(GJRws)) :
				self.GJR_weights[i] = GJRws[i]/CT10ws[0]
		else :
			for i in range(len(self.CT10_weights)) :
				self.CT10_weights[i] = 1.0
			for i in range(len(self.cteq_weights)) :
				self.cteq_weights[i] = 1.0
			for i in range(len(self.GJR_weights)) :
				self.GJR_weights[i] = 1.0

		#Mother particle and MC truth top assignment
		if self.is_data == 0 : #MC truth values only relevant for semileptonic qqbar->ttbar
			q_vec 	 = findInitialQuark(self.MC_generator,GenParticles,genPartVars) #function in eventTypeHelper.py
			qbar_vec = ROOT.TLorentzVector(q_vec.X(),q_vec.Y(),-1.0*q_vec.Z(),q_vec.E())
			MCt_vec, MCtbar_vec = findMCTops(self.MC_generator,GenParticles) #function in eventTypeHelper.py
		else : #if we don't have the MC truth information, we have to assign which is which later when we do the boost
			q_vec 	 = ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
			qbar_vec = ROOT.TLorentzVector(1.0,0.0,-1.0*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
			MCt_vec    = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
			MCtbar_vec = ROOT.TLorentzVector(-1.0,0.0,0.0,1.0)
		self.__fillMC__(q_vec,qbar_vec,MCt_vec,MCtbar_vec)

		#get all the info from the event
		#MET
		metVars = []
		for i in range(len(self.metHandles)) :
			event.getByLabel(self.metLabels[i],self.metHandles[i])
			if not self.metHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			metVars.append(self.metHandles[i].product())
		met = ROOT.TLorentzVector(); met.SetPtEtaPhiM(metVars[0][0], 0., metVars[0][1], 0.)
		self.__fillMET__(met)
		#jets
		alljets = []; jets = []; jetVars = []
		for i in range(len(self.jetHandles)) :
			event.getByLabel(self.jetLabels[i],self.jetHandles[i])
			if not self.jetHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			jetVars.append(self.jetHandles[i].product())
		#adjust the jets as per the JEC
		if self.JES != 'nominal' or self.JER != 'nominal' :
			for i in range(len(jetVars[0])) :
				newJet = adjustJEC(jetVars[0][i],jetVars[1][i],jetVars[2][i],jetVars[3][i],jetVars[4][i],jetVars[5][i],jetVars[6][i],jetVars[7][i],self.JES,self.JER)
				jetVars[0][i].SetPt(newJet.Pt())
				jetVars[0][i].SetEta(newJet.Eta())
				jetVars[0][i].SetPhi(newJet.Phi())
				jetVars[0][i].SetM(newJet.M())
			for i in range(len(jetVars[8])) :
				newJet = adjustJEC(jetVars[8][i],jetVars[9][i],jetVars[10][i],jetVars[11][i],jetVars[12][i],jetVars[13][i],jetVars[14][i],jetVars[15][i],self.JES,self.JER)
				jetVars[8][i].SetPt(newJet.Pt())
				jetVars[8][i].SetEta(newJet.Eta())
				jetVars[8][i].SetPhi(newJet.Phi())
				jetVars[8][i].SetM(newJet.M())
		#build the list of analysis jets
		for i in range(len(jetVars[0])) :
			flavor = -1
			if len(jetVars)>20 :
				flavor = jetVars[20][i]
			newJet = jet(jetVars[0][i],jetVars[8],jetVars[16],jetVars[17],jetVars[18],jetVars[19][i],flavor)
			jets.append(newJet); alljets.append(newJet)
		#separate the jets into the top and b candidates
		jets = selectJets(jets)
		self.__fillJets__(jets)
		if len(jets)<2 :
			return self.ERR_CODE
		self.sf_btag_eff[0] = jets[1].btagSF 
		self.sf_btag_eff_low[0] = jets[1].btagSFlow 
		self.sf_btag_eff_hi[0] = jets[1].btagSFhigh
		#muons
		muons = []; muVars = []
		for i in range(len(self.muHandles)) :
			event.getByLabel(self.muLabels[i],self.muHandles[i])
			if not self.muHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			muVars.append(self.muHandles[i].product())
		for i in range(len(muVars[0])) :
			newMuon = muon(muVars[0][i],muVars[1][i],muVars[2][i],muVars[3][i],alljets)
			muons.append(newMuon)
		muons.sort(key = lambda x: x.vec.Pt(),reverse=True)
		#electrons
		electrons = []; elVars = []
		for i in range(len(self.elHandles)) :
			event.getByLabel(self.elLabels[i],self.elHandles[i])
			if not self.elHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			elVars.append(self.elHandles[i].product())
		for i in range(len(elVars[0])) :
			newElectron = electron(elVars[0][i],elVars[1][i],elVars[2][i],elVars[3][i],met,alljets)
			electrons.append(newElectron)
		electrons.sort(key = lambda x: x.vec.Pt(),reverse=True)
		#assign the lepton type and fill the tree
		if len(muons)<1 and len(electrons)<1 :
			return self.ERR_CODE
		self.lep_type = 0
		if len(muons)<1 or (len(electrons)>0 and len(muons)>0 and electrons[0].vec.Pt()>muons[0].vec.Pt()) :
			self.lep_type = 1
		self.__fillMuons__(muons)
		self.__fillElectrons__(electrons)
		#pileup
		event.getByLabel(self.pileupLabel,self.pileupHandle)
		if not self.pileupHandle.isValid() :
			self.ERR_CODE = ERR_INVALID_HANDLE
			return self.ERR_CODE
		pileup = self.pileupHandle.product()[0]
		MCpileup = 0
		if self.is_data == 0 :
			event.getByLabel(self.MCpileupLabel,self.MCpileupHandle)
			if not self.MCpileupHandle.isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			MCpileup = self.MCpileupHandle.product()[0]
		self.pileup[0] = int(pileup); self.MC_pileup[0] = int(MCpileup)
		
		#neutrino handling and setup for fit
		if self.lep_type==0 :
			met1_vec, met2_vec = setupMET(muons[0].vec,metVars) #function in metHelper.py
		elif self.lep_type==1 :
			met1_vec, met2_vec = setupMET(electrons[0].vec,metVars) #function in metHelper.py
		self.met_pt[0], self.met_eta[0] = met1_vec.Pt(),  met1_vec.Eta() 
		self.met_phi[0], self.met_M[0]  = met1_vec.Phi(), met1_vec.M()
		self.nFits[0] = 2
		if met1_vec.Pz() == met2_vec.Pz() :
			self.nFits[0] = 1
		
		#fill the rest of the leptonic fourvectors
		if self.lep_type==0 :
			self.__fillLepSide__(muons[0].vec,met1_vec,jets[1].vec)
		elif self.lep_type==1 :
			self.__fillLepSide__(electrons[0].vec,met1_vec,jets[1].vec)

		#event reconstruction with kinematic fit
		scaledlep = ROOT.TLorentzVector(); scaledmet = ROOT.TLorentzVector() 
		scaledlepb = ROOT.TLorentzVector(); scaledhadt = ROOT.TLorentzVector()
		if self.lep_type == 0 :
			scaledlep, scaledmet, scaledlepb, scaledhadt, self.chi2[0] = reconstruct(muons[0].vec,met1_vec,met2_vec,jets) 
		elif self.lep_type == 1 :
			scaledlep, scaledmet, scaledlepb, scaledhadt, self.chi2[0] = reconstruct(electrons[0].vec,met1_vec,met2_vec,jets) 
		#above function in ttbarReconstructor.py

		#fill the TTree with the scaled fourvector variables
		self.__fillScaledFourVecs__(scaledlep,scaledmet,scaledlepb,scaledhadt)

		#reconstruct the observables using both the scaled and unscaled vectors
		if self.lep_type == 0 :
			self.cstar[0], self.x_F[0], self.M[0] = getObservables(muons[0].vec+met1_vec+jets[1].vec,jets[0].vec,self.Q_l[0]) 
		elif self.lep_type == 1 :
			self.cstar[0], self.x_F[0], self.M[0] = getObservables(electrons[0].vec+met1_vec+jets[1].vec,jets[0].vec,self.Q_l[0]) 
		( self.cstar_scaled[0], self.x_F_scaled[0], 
			self.M_scaled[0] ) = getObservables(scaledlep+scaledmet+scaledlepb,scaledhadt,self.Q_l[0]) 
		#above function in angleReconstructor.py

		#MC Truth observable and reweighting calculation
		if self.is_data==0 :
			if self.event_type!=4 :
				( self.cstar_MC[0],self.x_F_MC[0],self.M_MC[0],
					self.wg1[0],self.wg2[0],self.wg3[0],self.wg4[0],
					self.wqs1[0],self.wqs2[0],self.wqa0[0],self.wqa1[0],self.wqa2[0],
					self.wg1_opp[0],self.wg2_opp[0],self.wg3_opp[0],self.wg4_opp[0],
					self.wqs1_opp[0],self.wqs2_opp[0],
					self.wqa0_opp[0],self.wqa1_opp[0],self.wqa2_opp[0],
					self.wega[0], self.wegc[0] ) = getMCObservables(q_vec,qbar_vec,MCt_vec,MCtbar_vec,self.event_type) 
			#scale factor and reweighting calculations
			if self.lep_type==0 :
				meas_lep_pt=muons[0].vec.Pt(); meas_lep_eta=muons[0].vec.Eta()
			elif self.lep_type==1 :
				meas_lep_pt=electrons[0].vec.Pt(); meas_lep_eta=electrons[0].vec.Eta()
			#8TeV numbers
			self.sf_top_pT[0], self.sf_top_pT_low[0], self.sf_top_pT_hi[0] = self.corrector.getToppT_reweight(MCt_vec,MCtbar_vec,self.Q_l[0])
			self.sf_top_pT_new[0], self.sf_top_pT_new_low[0], self.sf_top_pT_new_hi[0] = self.corrector.getNewToppT_reweight(MCt_vec,MCtbar_vec)
			self.sf_pileup[0], self.sf_pileup_low[0], self.sf_pileup_hi[0] = self.corrector.getpileup_reweight(MCpileup)
			( self.sf_lep_ID[0], self.sf_lep_ID_low[0], 
				self.sf_lep_ID_hi[0] ) = self.corrector.getID_eff(pileup,meas_lep_pt,meas_lep_eta,self.lep_type)
			( self.sf_trig_eff[0], self.sf_trig_eff_low[0], 
				self.sf_trig_eff_hi[0] ) = self.corrector.gettrig_eff(pileup,meas_lep_pt,meas_lep_eta,self.lep_type)
		self.__closeout__() #yay! A successful event!

	##################################  #__init__ function  ##################################
	def __init__(self,fileName,tree,isData,generator,jec,onGrid) :
		#output file
		self.outfile_name = fileName
		self.outfile = TFile(self.outfile_name,'recreate')
		#output tree
		self.tree = TTree('tree','recreate')
		#total weight histogram
		self.tot_weight_histo = TH1F('totweight','Total Sum of Weights; ; total weight',1,0.,1.)
		self.tot_weight_histo.SetDirectory(0)
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
		elif eventType == 'had_TT' :
			print 'only HADRONIC EVENTS will be analyzed from this file'
			self.event_type = 3
		else :
			print 'ALL event types will be analyzed from this file'
			self.event_type = 4
		#input tree
		self.inputTree = tree
		#Set input and output branch locations
		for branch in self.allBranches.values() :
			branch.initialize(self.inputTree,self.tree)
		#data or MC?
		if isData == 'no' :
			self.is_data = False
		elif isData == 'yes' :
			self.is_data = True
		else :
			print 'ERROR: cannot determine if inputted file is a data or MC file!'
			print '	options.data = '+isData+''
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
		#Set Monte Carlo reweighter
		self.corrector = MC_corrector(self.MC_generator,self.event_type,onGrid)

	##################################   reset function   ##################################
	#########  sets all relevant values back to initial values to get ready for next event  ##########
	def reset(self) :
		for branch in self.allBranches.values() :
			branch.reset()
		pass

	################################   addBranch function  #################################
	def __AddBranch__(self,readname=None,writename=None,ttreetype='F',inival=-900.,size='1',dictlist=None) :
		newBranch = Branch(readname=readname,writename=writename,ttreetype=ttreetype,inival=inival,size=size)
		for d in dictlist :
			if writename==None :
				d[readname] = newBranch
			else :
				d[writename] = newBranch
		return newBranch

	##############################   tree filling functions  ###############################		
	def __fillMC__(self,qvec,qbarvec,MCtvec,MCtbarvec) :
		self.q_pt[0], 		self.q_eta[0] 	   = qvec.Pt(), 	  qvec.Eta()
		self.q_phi[0], 		self.q_M[0] 	   = qvec.Phi(), 	  qvec.M()
		self.qbar_pt[0], 	self.qbar_eta[0]   = qbarvec.Pt(),    qbarvec.Eta()
		self.qbar_phi[0], 	self.qbar_M[0] 	   = qbarvec.Phi(),   qbarvec.M()
		self.MCt_pt[0], 	self.MCt_eta[0]    = MCtvec.Pt(), 	  MCtvec.Eta()
		self.MCt_phi[0], 	self.MCt_M[0] 	   = MCtvec.Phi(),    MCtvec.M()
		self.MCtbar_pt[0], 	self.MCtbar_eta[0] = MCtbarvec.Pt(),  MCtbarvec.Eta()
		self.MCtbar_phi[0], self.MCtbar_M[0]   = MCtbarvec.Phi(), MCtbarvec.M()

	def __fillMET__(self,metvec) :
		self.met_pt[0]  = metvec.Pt(); self.met_eta[0] = metvec.Eta()
		self.met_phi[0] = metvec.Phi(); self.met_M[0] 	= metvec.M()

	def __fillJets__(self,jetList) :
		if len(jetList) > 0 :
			self.hadt_pt[0]  = jetList[0].vec.Pt(); self.hadt_eta[0] = jetList[0].vec.Eta()
			self.hadt_phi[0] = jetList[0].vec.Phi(); self.hadt_M[0]   = jetList[0].vec.M()
			self.hadt_tau32[0]  = jetList[0].tau32; self.hadt_tau21[0]  = jetList[0].tau21
			self.hadt_csv[0] 	 = jetList[0].csv; self.hadt_flavor[0] = jetList[0].flavor
		if len(jetList) > 1 :
			self.lepb_pt[0]  = jetList[1].vec.Pt(); self.lepb_eta[0] = jetList[1].vec.Eta()
			self.lepb_phi[0] = jetList[1].vec.Phi(); self.lepb_M[0]   = jetList[1].vec.M()
			self.lepb_tau32[0]  = jetList[1].tau32; self.lepb_tau21[0]  = jetList[1].tau21
			self.lepb_csv[0] 	 = jetList[1].csv; self.lepb_flavor[0] = jetList[1].flavor

	def __fillMuons__(self,mulist) :
		if len(mulist) > 0 :
			self.muon1_pt[0]  = mulist[0].vec.Pt(); self.muon1_eta[0] = mulist[0].vec.Eta()
			self.muon1_phi[0] = mulist[0].vec.Phi(); self.muon1_M[0]   = mulist[0].vec.M()
			self.muon1_Q[0] 	   = mulist[0].charge; self.muon1_isTight[0] = mulist[0].isTight
			self.muon1_isLoose[0] = mulist[0].isLoose; self.muon1_relPt[0]   = mulist[0].relPt
			self.muon1_dR[0] 	   = mulist[0].dR
			if self.lep_type == 0 :
				self.Q_l[0] = mulist[0].charge
		if len(mulist) > 1 :
			self.muon2_pt[0]  = mulist[1].vec.Pt(); self.muon2_eta[0] = mulist[1].vec.Eta()
			self.muon2_phi[0] = mulist[1].vec.Phi(); self.muon2_M[0]   = mulist[1].vec.M()
			self.muon2_Q[0] 	   = mulist[1].charge; self.muon2_isTight[0] = mulist[1].isTight
			self.muon2_isLoose[0] = mulist[1].isLoose; self.muon2_relPt[0]   = mulist[1].relPt
			self.muon2_dR[0] 	   = mulist[1].dR

	def __fillElectrons__(self,ellist) :
		if len(ellist) > 0 :
			self.ele1_pt[0]  = ellist[0].vec.Pt(); self.ele1_eta[0] = ellist[0].vec.Eta()
			self.ele1_phi[0] = ellist[0].vec.Phi(); self.ele1_M[0]   = ellist[0].vec.M()
			self.ele1_Q[0] 	   = ellist[0].charge; self.ele1_isTight[0] = ellist[0].isTight
			self.ele1_isLoose[0] = ellist[0].isLoose; self.ele1_relPt[0]   = ellist[0].relPt
			self.ele1_dR[0] 	   = ellist[0].dR
			self.ele1_tri_el_val[0]  = ellist[0].triangle_el_val
			self.ele1_tri_jet_val[0] = ellist[0].triangle_jet_val
			self.ele1_tri_cut_val[0] = ellist[0].triangle_cut_val
			if self.lep_type == 1 :
				self.Q_l[0] = ellist[0].charge
		if len(ellist) > 1 :
			self.ele2_pt[0]  = ellist[1].vec.Pt(); self.ele2_eta[0] = ellist[1].vec.Eta()
			self.ele2_phi[0] = ellist[1].vec.Phi(); self.ele2_M[0]   = ellist[1].vec.M()
			self.ele2_Q[0] 	   = ellist[1].charge; self.ele2_isTight[0] = ellist[1].isTight
			self.ele2_isLoose[0] = ellist[1].isLoose; self.ele2_relPt[0]   = ellist[1].relPt
			self.ele2_dR[0] 	   = ellist[1].dR
			self.ele2_tri_el_val[0]  = ellist[1].triangle_el_val
			self.ele2_tri_jet_val[0] = ellist[1].triangle_jet_val
			self.ele2_tri_cut_val[0] = ellist[1].triangle_cut_val

	def __fillLepSide__(self,lepton,met,lepb) :
		lepW = lepton+met
		self.lepW_pt[0] 	= lepW.Pt(); 	self.lepW_eta[0] = lepW.Eta()
		self.lepW_phi[0] = lepW.Phi(); 	self.lepW_M[0] 	= lepW.M()
		lept = lepW+lepb
		self.lept_pt[0] 	= lept.Pt(); 	self.lept_eta[0] = lept.Eta()
		self.lept_phi[0] = lept.Phi(); 	self.lept_M[0] 	= lept.M()

	def __fillScaledFourVecs__(self,lepton,met,lepb,hadt) :
		self.scaled_lep_pt[0] 	= lepton.Pt(); 	self.scaled_lep_eta[0] = lepton.Eta()
		self.scaled_lep_phi[0] = lepton.Phi(); self.scaled_lep_M[0] 	= lepton.M()
		self.scaled_met_pt[0] 	= met.Pt(); 	self.scaled_met_eta[0] = met.Eta()
		self.scaled_met_phi[0] = met.Phi(); 	self.scaled_met_M[0] 	= met.M()
		lepW = lepton+met
		self.scaled_lepW_pt[0] 	= lepW.Pt(); 	self.scaled_lepW_eta[0] = lepW.Eta()
		self.scaled_lepW_phi[0] = lepW.Phi(); 	self.scaled_lepW_M[0] 	= lepW.M()
		self.scaled_lepb_pt[0]  = lepb.Pt();  self.scaled_lepb_eta[0] = lepb.Eta()
		self.scaled_lepb_phi[0] = lepb.Phi(); self.scaled_lepb_M[0] 	= lepb.M()
		lept = lepW+lepb
		self.scaled_lept_pt[0] 	= lept.Pt(); 	self.scaled_lept_eta[0] = lept.Eta()
		self.scaled_lept_phi[0] = lept.Phi(); 	self.scaled_lept_M[0] 	= lept.M()
		self.scaled_hadt_pt[0]  = hadt.Pt();  self.scaled_hadt_eta[0] = hadt.Eta()
		self.scaled_hadt_phi[0] = hadt.Phi(); self.scaled_hadt_M[0] 	= hadt.M()

	########## function to close out the event, called before kicking back to runner #########
	def __closeout__(self) :
		#fill ttree
		self.tree.Fill()

	################################## __del__ function  ###################################
	def __del__(self) :
		self.f.cd()
		self.f.Write()
		self.f.Close()