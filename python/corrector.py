#Global variables
#Pileup
DATA_PU_HISTO_NOMINAL_FILENAME = 'DataPileupHistogram_69200.root'
DATA_PU_HISTO_XS_UP_FILENAME   = 'DataPileupHistogram_72383.root'
DATA_PU_HISTO_XS_DOWN_FILENAME = 'DataPileupHistogram_66017.root'
#Muon trigger efficiency
MUON_TRIG_EFF_BTOF_ROOT_FILENAME = 'EfficienciesAndSF_RunBtoF_muon_trigger.root'
MUON_TRIG_EFF_GH_ROOT_FILENAME 	 = 'EfficienciesAndSF_Period4_muon_trigger.root'
MUON_TRIG_EFF_BOOSTED_PT_HISTONAME 	 = 'Mu50_OR_TkMu50_PtBins/pt_ratio'
MUON_TRIG_EFF_BOOSTED_ETA_HISTONAME  = 'Mu50_OR_TkMu50_EtaBins/eta_ratio'
MUON_TRIG_EFF_RESOLVED_PT_HISTONAME  = 'IsoMu24_OR_IsoTkMu24_PtBins/pt_ratio'
MUON_TRIG_EFF_RESOLVED_ETA_HISTONAME = 'IsoMu24_OR_IsoTkMu24_EtaBins/eta_ratio'
#Electron trigger efficiency
ELE_TRIG_EFF_BOOSTED_ROOT_FILENAME = 'SFs_elec_eta_pt_DATA_Run2016_Ele50PFJet165_OR_Ele115.root'
ELE_TRIG_EFF_BOOSTED_HISTONAME 	   = 'EGamma_SF2D'
ELE_TRIG_EFF_RESOLVED_BTOF_ROOT_FILENAME = 'TriggerSF_Run2016BCDEF_v1.root'
ELE_TRIG_EFF_RESOLVED_GH_ROOT_FILENAME 	 = 'TriggerSF_Run2016GH_v1.root'
ELE_TRIG_EFF_RESOLVED_HISTONAME 		 = 'Ele27_WPTight_Gsf'
#Muon tracking efficiency
MUON_TRK_EFF_ROOT_FILENAME = 'Tracking_EfficienciesAndSF_BCDEFGH_muon_tracking.root'
MUON_TRK_EFF_ETA_GRAPHNAME = 'ratio_eff_eta3_dr030e030_corr'
MUON_TRK_EFF_PU_GRAPHNAME  = 'ratio_eff_vtx_dr030e030_corr'
#Muon ID Efficiency
MUON_ID_EFF_BTOF_ROOT_FILENAME 	  = 'EfficienciesAndSF_BCDEF_muon_ID.root'
MUON_ID_EFF_BTOF_PT_ETA_HISTONAME = 'MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio'
MUON_ID_EFF_BTOF_PU_HISTONAME 	  = 'MC_NUM_MediumID2016_DEN_genTracks_PAR_vtx/tag_nVertices_ratio'
MUON_ID_EFF_GH_ROOT_FILENAME 	  = 'EfficienciesAndSF_GH_muon_ID.root'
MUON_ID_EFF_GH_PT_ETA_HISTONAME   = 'MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio'
MUON_ID_EFF_GH_PU_HISTONAME 	  = 'MC_NUM_MediumID_DEN_genTracks_PAR_vtx/tag_nVertices_ratio'
#Electron ID Efficiency
ELE_ID_EFF_ROOT_FILENAME   = 'egamma_ID_ISO_eff_SFs.root'
ELE_ID_EFF_2D_HISTONAME    = 'GsfElectronToCutBasedSpring15T'
ELE_RECO_EFF_ROOT_FILENAME = 'egamma_reco_eff_SFs.root'
ELE_RECO_EFF_2D_HISTONAME  = 'EGamma_SF2D'
#Muon isolation efficiency
MUON_ISO_EFF_BTOF_ROOT_FILENAME = 'EfficienciesAndSF_BCDEF_muon_iso.root'
MUON_ISO_EFF_GH_ROOT_FILENAME   = 'EfficienciesAndSF_GH_muon_iso.root'
MUON_ISO_EFF_PT_ETA_HISTONAME  = 'TightISO_MediumID_pt_eta/pt_abseta_ratio'
MUON_ISO_EFF_PU_HISTONAME  = 'TightISO_MediumID_vtx/tag_nVertices_ratio'
#Electron MiniIsolation Efficiency
ELE_MINIISO_EFF_ROOT_FILENAME   = 'egamma_ID_ISO_eff_SFs.root'
ELE_MINIISO_EFF_2D_HISTONAME    = 'MVAVLooseElectronToMini'
#b-Tagging efficiency
BTAGGING_MC_ROOT_FILENAME 		 = 'MC_btagging_efficiency.root'
BTAGGING_MC_B_RATIO_HISTONAME 	 = 'bjet_ratio'
BTAGGING_MC_C_RATIO_HISTONAME 	 = 'cjet_ratio'
BTAGGING_MC_UDSG_RATIO_HISTONAME = 'udsg_ratio'
BTAGGING_SF_CSV_FILENAME 		 = 'CSVv2_Moriond17_B_H.csv'

#Imports
from ROOT import vector, FactorizedJetCorrector, JetCorrectorParameters, JetCorrectionUncertainty 
from ROOT import TFile, Double, gSystem, BTagCalibration, BTagCalibrationReader, std
import ROOT
from random import gauss
from math import sqrt

class Corrector(object) :

	##################################  #__init__ function  ##################################
	def __init__(self,isdata,onGrid,pu_histo,runera,renormdict) :
		#JEC setup
		self.__ak4JetCorrector, self.__ak4JecUncertainty = setupJECCorrector(onGrid,isdata,'AK4PFchs',runera)
		self.__ak8JetCorrector, self.__ak8JecUncertainty = setupJECCorrector(onGrid,isdata,'AK8PFchs',runera)
		#pileup
		self.__MC_pu_histo, self.__data_pu_histo_nom, self.__data_pu_histo_up, self.__data_pu_histo_down = self.setupPileupHistos(onGrid,pu_histo)
		#trigger efficiency
		( self.__muon_trig_eff_boosted_BtoF_vs_pt, self.__muon_trig_eff_boosted_BtoF_vs_eta, 
			self.__muon_trig_eff_boosted_GH_vs_pt, self.__muon_trig_eff_boosted_GH_vs_eta, 
			self.__muon_trig_eff_resolved_BtoF_vs_pt, self.__muon_trig_eff_resolved_BtoF_vs_eta, 
			self.__muon_trig_eff_resolved_GH_vs_pt, self.__muon_trig_eff_resolved_GH_vs_eta ) = self.setupMuonTriggerHistos(onGrid)
		self.__muon_trig_b_pt_low  = self.__muon_trig_eff_boosted_BtoF_vs_pt.GetXaxis().GetXmin()
		self.__muon_trig_b_pt_hi   = self.__muon_trig_eff_boosted_BtoF_vs_pt.GetXaxis().GetXmax()
		self.__muon_trig_b_eta_low = self.__muon_trig_eff_boosted_BtoF_vs_eta.GetXaxis().GetXmin()
		self.__muon_trig_b_eta_hi  = self.__muon_trig_eff_boosted_BtoF_vs_eta.GetXaxis().GetXmax()
		self.__muon_trig_r_pt_low  = self.__muon_trig_eff_resolved_BtoF_vs_pt.GetXaxis().GetXmin()
		self.__muon_trig_r_pt_hi   = self.__muon_trig_eff_resolved_BtoF_vs_pt.GetXaxis().GetXmax()
		self.__muon_trig_r_eta_low = self.__muon_trig_eff_resolved_BtoF_vs_eta.GetXaxis().GetXmin()
		self.__muon_trig_r_eta_hi  = self.__muon_trig_eff_resolved_BtoF_vs_eta.GetXaxis().GetXmax()
		self.__ele_trig_eff_boosted, self.__ele_trig_eff_resolved_BtoF, self.__ele_trig_eff_resolved_GH = self.setupEleTriggerHistos(onGrid)
		etxa_b = self.__ele_trig_eff_boosted.GetXaxis(); etya_b = self.__ele_trig_eff_boosted.GetYaxis()
		self.__ele_trig_eff_b_pt_low  = etxa_b.GetXmin(); self.__ele_trig_eff_b_pt_hi  = etxa_b.GetXmax()
		self.__ele_trig_eff_b_eta_low = etya_b.GetXmin(); self.__ele_trig_eff_eta_hi = etya_b.GetXmax()
		etxa_r = self.__ele_trig_eff_resolved_BtoF.GetXaxis(); etya_r = self.__ele_trig_eff_resolved_BtoF.GetYaxis()
		self.__ele_trig_eff_r_pt_low  = etxa_r.GetXmin(); self.__ele_trig_eff_r_pt_hi  = etxa_r.GetXmax()
		self.__ele_trig_eff_r_eta_low = etya_r.GetXmin(); self.__ele_trig_eff_r_eta_hi = etya_r.GetXmax()
		#Tracking efficiency
		self.__muon_trk_eff_vs_eta, self.__muon_trk_eff_vs_pu = self.setupMuonTrackingGraphs(onGrid)
		tempx = Double(0.); tempy = Double(0.)
		self.__muon_trk_eff_vs_eta_nbins = self.__muon_trk_eff_vs_eta.GetN()
		self.__muon_trk_eff_vs_eta_lowxs = []; self.__muon_trk_eff_vs_eta_hixs = []
		self.__muon_trk_eff_vs_eta_ys = []; self.__muon_trk_eff_vs_eta_yerrups = []; self.__muon_trk_eff_vs_eta_yerrdowns = []
		for i in range(self.__muon_trk_eff_vs_eta_nbins) :
			self.__muon_trk_eff_vs_eta.GetPoint(i,tempx,tempy)
			self.__muon_trk_eff_vs_eta_lowxs.append(float(tempx-self.__muon_trk_eff_vs_eta.GetErrorXlow(i)))
			self.__muon_trk_eff_vs_eta_hixs.append(float(tempx+self.__muon_trk_eff_vs_eta.GetErrorXhigh(i)))
			self.__muon_trk_eff_vs_eta_ys.append(float(tempy))
			self.__muon_trk_eff_vs_eta_yerrups.append(float(self.__muon_trk_eff_vs_eta.GetErrorYhigh(i)))
			self.__muon_trk_eff_vs_eta_yerrdowns.append(float(self.__muon_trk_eff_vs_eta.GetErrorYlow(i)))
		self.__muon_trk_eff_vs_pu_nbins = self.__muon_trk_eff_vs_pu.GetN()
		self.__muon_trk_eff_vs_pu_lowxs = []; self.__muon_trk_eff_vs_pu_hixs = []
		self.__muon_trk_eff_vs_pu_ys = []; self.__muon_trk_eff_vs_pu_yerrups = []; self.__muon_trk_eff_vs_pu_yerrdowns = []
		for i in range(self.__muon_trk_eff_vs_pu_nbins) :
			self.__muon_trk_eff_vs_pu.GetPoint(i,tempx,tempy)
			self.__muon_trk_eff_vs_pu_lowxs.append(float(tempx-self.__muon_trk_eff_vs_pu.GetErrorXlow(i)))
			self.__muon_trk_eff_vs_pu_hixs.append(float(tempx+self.__muon_trk_eff_vs_pu.GetErrorXhigh(i)))
			self.__muon_trk_eff_vs_pu_ys.append(float(tempy))
			self.__muon_trk_eff_vs_pu_yerrups.append(float(self.__muon_trk_eff_vs_pu.GetErrorYhigh(i)))
			self.__muon_trk_eff_vs_pu_yerrdowns.append(float(self.__muon_trk_eff_vs_pu.GetErrorYlow(i)))
		#ID efficiency
		( self.__muon_id_eff_BtoF_abseta_vs_pt, self.__muon_id_eff_BtoF_vs_pu, 
			self.__muon_id_eff_GH_abseta_vs_pt, self.__muon_id_eff_GH_vs_pu ) = self.setupMuonIDHistos(onGrid)
		mixa = self.__muon_id_eff_BtoF_abseta_vs_pt.GetXaxis(); miya = self.__muon_id_eff_BtoF_abseta_vs_pt.GetYaxis()
		self.__muon_id_eff_pt_low  = mixa.GetXmin(); 									 self.__muon_id_eff_pt_hi  = mixa.GetXmax()
		self.__muon_id_eff_eta_low = miya.GetXmin(); 									 self.__muon_id_eff_eta_hi = miya.GetXmax()
		self.__muon_id_eff_pu_low  = self.__muon_id_eff_BtoF_vs_pu.GetXaxis().GetXmin(); self.__muon_id_eff_pu_hi  = self.__muon_id_eff_BtoF_vs_pu.GetXaxis().GetXmax()
		self.__ele_id_eff_eta_vs_pt, self.__ele_reco_eff_pt_vs_eta = self.setupEleIDHistos(onGrid)
		idxa = self.__ele_id_eff_eta_vs_pt.GetXaxis(); idya = self.__ele_id_eff_eta_vs_pt.GetYaxis()
		self.__ele_id_eff_pt_low  = idxa.GetXmin(); self.__ele_id_eff_pt_hi  = idxa.GetXmax()
		self.__ele_id_eff_eta_low = idya.GetXmin(); self.__ele_id_eff_eta_hi = idya.GetXmax()
		recoxa = self.__ele_reco_eff_pt_vs_eta.GetXaxis(); recoya = self.__ele_reco_eff_pt_vs_eta.GetYaxis()
		self.__ele_reco_eff_eta_low = recoxa.GetXmin(); self.__ele_reco_eff_eta_hi = recoxa.GetXmax()
		self.__ele_reco_eff_pt_low  = recoya.GetXmin(); self.__ele_reco_eff_pt_hi  = recoya.GetXmax()
		#Isolation efficiency
		( self.__muon_iso_eff_BtoF_abseta_vs_pt, self.__muon_iso_eff_BtoF_vs_pu, 
			self.__muon_iso_eff_GH_abseta_vs_pt, self.__muon_iso_eff_GH_vs_pu ) = self.setupMuonIsoHistos(onGrid)
		misoxa = self.__muon_iso_eff_BtoF_abseta_vs_pt.GetXaxis(); misoya = self.__muon_iso_eff_BtoF_abseta_vs_pt.GetYaxis()
		self.__muon_iso_eff_pt_low  = misoxa.GetXmin(); self.__muon_iso_eff_pt_hi  = misoxa.GetXmax()
		self.__muon_iso_eff_eta_low = misoya.GetXmin(); self.__muon_iso_eff_eta_hi = misoya.GetXmax()
		self.__muon_iso_eff_pu_low  = self.__muon_iso_eff_BtoF_vs_pu.GetXaxis().GetXmin();  self.__muon_iso_eff_pu_hi  = self.__muon_iso_eff_BtoF_vs_pu.GetXaxis().GetXmax()
		#MiniIsolation efficiency
		self.__ele_miniiso_eff_abseta_vs_pt = self.setupEleMiniIsoHistos(onGrid)
		miniisoxa = self.__ele_miniiso_eff_abseta_vs_pt.GetXaxis(); miniisoya = self.__ele_miniiso_eff_abseta_vs_pt.GetYaxis()
		self.__ele_miniiso_eff_pt_low  = miniisoxa.GetXmin(); self.__ele_miniiso_eff_pt_hi  = miniisoxa.GetXmax()
		self.__ele_miniiso_eff_eta_low = miniisoya.GetXmin(); self.__ele_miniiso_eff_eta_hi = miniisoya.GetXmax()
		#b-tagging efficiency
		self.__setupBTaggingEff__(onGrid)
		#PDF/alpha_s stuff
		self.__renormdict = renormdict
		self.__alphalist = getList(renormdict,'alpha')
		self.__epsilonlist = getList(renormdict,'epsilon')

	def __setupBTaggingEff__(self,onGrid) :
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		#Get the MC eff file histograms and limits
		btag_MC_eff_file = TFile.Open(pp+BTAGGING_MC_ROOT_FILENAME)
		self.__btag_MC_eff_bjet_histo_l  = btag_MC_eff_file.Get(BTAGGING_MC_B_RATIO_HISTONAME+'_loose')
		self.__btag_MC_eff_cjet_histo_l  = btag_MC_eff_file.Get(BTAGGING_MC_C_RATIO_HISTONAME+'_loose')
		self.__btag_MC_eff_udsg_histo_l  = btag_MC_eff_file.Get(BTAGGING_MC_UDSG_RATIO_HISTONAME+'_loose')
		self.__btag_MC_eff_bjet_histo_m  = btag_MC_eff_file.Get(BTAGGING_MC_B_RATIO_HISTONAME+'_medium')
		self.__btag_MC_eff_cjet_histo_m  = btag_MC_eff_file.Get(BTAGGING_MC_C_RATIO_HISTONAME+'_medium')
		self.__btag_MC_eff_udsg_histo_m  = btag_MC_eff_file.Get(BTAGGING_MC_UDSG_RATIO_HISTONAME+'_medium')
		xaxis = self.__btag_MC_eff_bjet_histo_l.GetXaxis(); yaxis = self.__btag_MC_eff_bjet_histo_l.GetYaxis()
		self.__btag_MC_eff_eta_low = xaxis.GetXmax()
		self.__btag_MC_eff_eta_hi  = xaxis.GetXmax()
		self.__btag_MC_eff_pt_low = yaxis.GetXmax()
		self.__btag_MC_eff_pt_hi  = yaxis.GetXmax()
		self.__btag_MC_eff_bjet_histo.SetDirectory(0); self.__btag_MC_eff_cjet_histo.SetDirectory(0)
		self.__btag_MC_eff_udsg_histo.SetDirectory(0)
		btag_MC_eff_file.Close()
		#Load the stuff we need for interpreting the SF file
		#This code from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration
		gSystem.Load('libCondFormatsBTauObjects') 
		gSystem.Load('libCondToolsBTau') 
		# get the sf data loaded
		calib = BTagCalibration('csvv2',pp+BTAGGING_SF_CSV_FILENAME)
		# making a std::vector<std::string>> in python is a bit awkward, 
		# but works with root (needed to load other sys types):
		#v_sys = std.vector('string') #changed this line but I think this works
		v_sys = getattr(ROOT, 'vector<string>')()
		v_sys.push_back('up')
		v_sys.push_back('down')
		# make a reader instance and load the sf data
		self.__reader_l = BTagCalibrationReader(0, 		 # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
		    								  "central", # central systematic type
		    								  v_sys) 	 # vector of other sys. types    
		self.__reader_l.load(calib, 
						   0, 	   # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
						   "comb") # measurement type
		self.__reader_l.load(calib,1,"comb")
		self.__reader_l.load(calib,2,"incl")
		self.__reader_m = BTagCalibrationReader(1,"central",v_sys) #1: medium
		self.__reader_m.load(calib,0,"comb")
		self.__reader_m.load(calib,1,"comb")
		self.__reader_m.load(calib,2,"incl")

	def __getTrkEff__(self,pileup,lepflav_or_lep,pt=-111.,eta=-111.) :
		nomfac = 1.; uperr = 1.; downerr = 1.
		doubleetaerr = False; doublepuerr = False
		lepflav = lepflav_or_lep
		if pt==-111. :
			lepflav = lepflav_or_lep.getType()
			pt = lepflav_or_lep.getPt()
			eta = lepflav_or_lep.getEta()
		if lepflav=='mu' : #muons
			#first bring the parameters back into range
			if eta<self.__muon_trk_eff_vs_eta_lowxs[0] : 
				eta = self.__muon_trk_eff_vs_eta_lowxs[0] 
				doubleetaerr=True
			elif eta>self.__muon_trk_eff_vs_eta_hixs[self.__muon_trk_eff_vs_eta_nbins-1] : 
				eta = self.__muon_trk_eff_vs_eta_hixs[self.__muon_trk_eff_vs_eta_nbins-1]- 0.0000001
				doubleetaerr=True
			pu = pileup
			if pu<self.__muon_trk_eff_vs_pu_lowxs[0] : 
				pu = self.__muon_trk_eff_vs_pu_lowxs[0] 
				doublepuerr=True
			elif pu>self.__muon_trk_eff_vs_pu_hixs[self.__muon_trk_eff_vs_pu_nbins-1] : 
				pu = self.__muon_trk_eff_vs_pu_hixs[self.__muon_trk_eff_vs_pu_nbins-1]- 0.0000001
				doublepuerr=True
			#get factors
			for i in range(self.__muon_trk_eff_vs_eta_nbins) :
				if eta>=self.__muon_trk_eff_vs_eta_lowxs[i] and eta<=self.__muon_trk_eff_vs_eta_hixs[i] :
					etafac 		= self.__muon_trk_eff_vs_eta_ys[i]
					etaerr_up 	= self.__muon_trk_eff_vs_eta_yerrups[i]
					etaerr_down = self.__muon_trk_eff_vs_eta_yerrdowns[i]
					if doubleetaerr :
						etaerr_up*=2; etaerr_down*=2
					break
			for i in range(self.__muon_trk_eff_vs_pu_nbins) :
				if pu>=self.__muon_trk_eff_vs_pu_lowxs[i] and pu<=self.__muon_trk_eff_vs_pu_hixs[i] :
					pufac 		= self.__muon_trk_eff_vs_pu_ys[i]
					puerr_up 	= self.__muon_trk_eff_vs_pu_yerrups[i]
					puerr_down = self.__muon_trk_eff_vs_pu_yerrdowns[i]
					if doublepuerr :
						puerr_up*=2; puerr_down*=2
					break
			#calculate total factors
			nomfac 	= etafac*pufac
			uperr 	= nomfac*sqrt((etaerr_up/etafac)**2+(puerr_up/pufac)**2)
			downerr = nomfac*sqrt((etaerr_down/etafac)**2+(puerr_down/pufac)**2)
		elif lepflav=='el' : #electrons
			#NOT IMPLEMENTED YET!!
			nomfac = 1.; uperr = 0.; downerr = 0.
		#return reweighting factors
		return nomfac,uperr,downerr

	def getJECforJet(self,jetvec,area,rho,npv,pp) :
		corrector=None
		if pp.find('AK4') != -1 :
			corrector = self.__ak4JetCorrector
		if pp.find('AK8') != -1 :
			corrector = self.__ak8JetCorrector
		corrector.setJetEta(jetvec.Eta())
		corrector.setJetPt(jetvec.Pt())
		corrector.setJetE(jetvec.E())
		corrector.setJetA(area)
		corrector.setRho(rho)
		corrector.setNPV(npv)
		newJEC = corrector.getCorrection()
		return newJEC

	def getJECuncForJet(self,jetvec,pp) :
		uncCorrector=None
		if pp.find('AK4') != -1 :
			uncCorrector = self.__ak4JecUncertainty
		if pp.find('AK8') != -1 :
			uncCorrector = self.__ak8JecUncertainty
		uncCorrector.setJetPhi(jetvec.Phi())
		uncCorrector.setJetEta(jetvec.Eta())
		uncCorrector.setJetPt(jetvec.Pt())
		uncDown = uncCorrector.getUncertainty(0)
		uncCorrector.setJetPhi(jetvec.Phi())
		uncCorrector.setJetEta(jetvec.Eta())
		uncCorrector.setJetPt(jetvec.Pt())
		uncUp   = uncCorrector.getUncertainty(1)
		return uncDown, uncUp

	def smearJet(self,jetvec,jer,genJetVec,ptres,dRCheck) :
		ptsmear=1.
		eta = jetvec.Eta()
		#get the smearing factors
		if jer=='nominal' :
			ptsmearfac = getJER(eta,0)
		elif jer=='down' :
			ptsmearfac = getJER(eta,-1)
		elif jer=='up' :
			ptsmearfac = getJER(eta,1)
		else :
			print 'WARNING: JER OPTION '+str(jer)+' NOT RECOGNIZED!!'
			return None
		#see which smearing method we should use based on MC matching
		#print '	ptsmearfac=%.3f, reco/gen dR = %.2f, abs(dpT)=%.3f, 3*ptres*pT=%.3f'%(ptsmearfac,jetvec.DeltaR(genJetVec),abs(jetvec.Pt()-genJetVec.Pt()),3.*ptres*jetvec.Pt()) #DEBUG
		if genJetVec!=None and jetvec.DeltaR(genJetVec)<dRCheck/2. and abs(jetvec.Pt()-genJetVec.Pt())<3.*ptres*jetvec.Pt() : #scaling
			#smear the pt
			recopt = jetvec.Pt()
			deltapt = (recopt-genJetVec.Pt())*(ptsmearfac-1.0)
			ptsmear = max(0.0, (recopt+deltapt)/recopt)
		else : #gaussian smearing
			sigma = sqrt(abs(ptsmearfac**2-1.))*ptres #Note the absolute value here, this is to handle cases where ptsmearfac<1., which are always out of selection range
			ptsmear = 1.+gauss(0.,sigma)
			#print '	sigma = %.4f, ptsmear = %.4f'%(sigma,ptsmear) #DEBUG
		ptsmearedjet = jetvec*ptsmear
		#return the smeared fourvector
		#print '	old jet = (%.1f,%.1f,%.1f,%.1f), newjet = (%.1f,%.1f,%.1f,%.1f)'%(jetvec.Pt(),jetvec.Eta(),jetvec.Phi(),jetvec.M(),ptsmearedjet.Pt(),ptsmearedjet.Eta(),ptsmearedjet.Phi(),ptsmearedjet.M()) #DEBUG
		return ptsmearedjet

	def getPileupReweight(self,pu_value) :
		mcbincontent = self.__MC_pu_histo.GetBinContent(self.__MC_pu_histo.FindFixBin(pu_value))
		sf = self.__data_pu_histo_nom.GetBinContent(self.__data_pu_histo_nom.FindFixBin(pu_value))/mcbincontent
		sf_up = self.__data_pu_histo_up.GetBinContent(self.__data_pu_histo_up.FindFixBin(pu_value))/mcbincontent
		sf_down = self.__data_pu_histo_down.GetBinContent(self.__data_pu_histo_down.FindFixBin(pu_value))/mcbincontent
		#print 'pileup = %d,	bin=%d,	mcbincontent=%6f,	sf=%.6f'%(pu_value,self.__MC_pu_histo.FindFixBin(pu_value),mcbincontent,sf) #DEBUG
		return sf, sf_up, sf_down

	def getTrigEff(self,topology,lepflav_or_lep,pt=-111.,eta=-111.) :
		nomfac_BtoF = 1.; upfac_BtoF = 1.; downfac_BtoF = 1.
		nomfac_GH = 1.;   upfac_GH = 1.;   downfac_GH = 1.
		doubleetaerr=False; doublepterr=False; doubleerr = False
		lepflav = lepflav_or_lep
		if pt==-111. :
			lepflav=lepflav_or_lep.getType()
			pt = lepflav_or_lep.getPt()
			if lepflav=='mu' :
				eta = lepflav_or_lep.getEta()
			elif lepflav=='el' :
				eta = lepflav_or_lep.getEtaSC()
		if lepflav=='mu' : #muons
			#first bring the parameters back into range and set the histograms to use
			if topology==1 or topology==2 :
				if eta<self.__muon_trig_b_eta_low : 
					eta = self.__muon_trig_b_eta_low
					doubleetaerr=True
				elif eta>self.__muon_trig_b_eta_hi : 
					eta = self.__muon_trig_b_eta_hi- 0.0000001
					doubleetaerr=True
				if pt<self.__muon_trig_b_pt_low : 
					pt = self.__muon_trig_b_pt_low 
					doublepterr=True
				elif pt>self.__muon_trig_b_pt_hi : 
					pt = self.__muon_trig_b_pt_hi - 0.0000001
					doublepterr=True
				btof_histo_vs_pt=self.__muon_trig_eff_boosted_BtoF_vs_pt
				btof_histo_vs_eta=self.__muon_trig_eff_boosted_BtoF_vs_eta
				gh_histo_vs_pt=self.__muon_trig_eff_boosted_GH_vs_pt
				gh_histo_vs_eta=self.__muon_trig_eff_boosted_GH_vs_eta
			elif topology==3 :
				if eta<self.__muon_trig_r_eta_low : 
					eta = self.__muon_trig_r_eta_low
					doubleetaerr=True
				elif eta>self.__muon_trig_r_eta_hi : 
					eta = self.__muon_trig_r_eta_hi- 0.0000001
					doubleetaerr=True
				if pt<self.__muon_trig_r_pt_low : 
					pt = self.__muon_trig_r_pt_low 
					doublepterr=True
				elif pt>self.__muon_trig_r_pt_hi : 
					pt = self.__muon_trig_r_pt_hi - 0.0000001
					doublepterr=True
				btof_histo_vs_pt=self.__muon_trig_eff_resolved_BtoF_vs_pt
				btof_histo_vs_eta=self.__muon_trig_eff_resolved_BtoF_vs_eta
				gh_histo_vs_pt=self.__muon_trig_eff_resolved_GH_vs_pt
				gh_histo_vs_eta=self.__muon_trig_eff_resolved_GH_vs_eta
			#get factors
			ptbin_BtoF 		 = btof_histo_vs_pt.FindFixBin(pt)
			ptfac_BtoF 		 = btof_histo_vs_pt.GetBinContent(ptbin_BtoF)
			pterr_up_BtoF 	 = btof_histo_vs_pt.GetBinErrorUp(ptbin_BtoF)
			pterr_down_BtoF  = btof_histo_vs_pt.GetBinErrorLow(ptbin_BtoF)
			etabin_BtoF 	 = btof_histo_vs_eta.FindFixBin(eta)
			etafac_BtoF 	 = btof_histo_vs_eta.GetBinContent(etabin_BtoF)
			etaerr_up_BtoF 	 = btof_histo_vs_eta.GetBinErrorUp(etabin_BtoF)
			etaerr_down_BtoF = btof_histo_vs_eta.GetBinErrorLow(etabin_BtoF)
			ptbin_GH 		 = gh_histo_vs_pt.FindFixBin(pt)
			ptfac_GH 		 = gh_histo_vs_pt.GetBinContent(ptbin_GH)
			pterr_up_GH 	 = gh_histo_vs_pt.GetBinErrorUp(ptbin_GH)
			pterr_down_GH 	 = gh_histo_vs_pt.GetBinErrorLow(ptbin_GH)
			etabin_GH 		 = gh_histo_vs_eta.FindFixBin(eta)
			etafac_GH 		 = gh_histo_vs_eta.GetBinContent(etabin_GH)
			etaerr_up_GH 	 = gh_histo_vs_eta.GetBinErrorUp(etabin_GH)
			etaerr_down_GH 	 = gh_histo_vs_eta.GetBinErrorLow(etabin_GH)
			if doublepterr :
				pterr_up_BtoF*=2; pterr_down_BtoF*=2; pterr_up_GH*=2; pterr_down_GH*=2
			if doubleetaerr :
				etaerr_up_BtoF*=2; etaerr_down_BtoF*=2; etaerr_up_GH*=2; etaerr_down_GH*=2
			#calculate total factors
			nomfac_BtoF  = ptfac_BtoF*etafac_BtoF
			upfac_BtoF 	 = nomfac_BtoF*(1.+sqrt((pterr_up_BtoF/ptfac_BtoF)**2+(etaerr_up_BtoF/etafac_BtoF)**2))
			downfac_BtoF = nomfac_BtoF*(1.-sqrt((pterr_down_BtoF/ptfac_BtoF)**2+(etaerr_down_BtoF/etafac_BtoF)**2))
			nomfac_GH 	 = ptfac_GH*etafac_GH
			upfac_GH 	 = nomfac_GH*(1.+sqrt((pterr_up_GH/ptfac_GH)**2+(etaerr_up_GH/etafac_GH)**2))
			downfac_GH 	 = nomfac_GH*(1.-sqrt((pterr_down_GH/ptfac_GH)**2+(etaerr_down_GH/etafac_GH)**2))
		elif lepflav=='el' : #electrons
			thiseta=eta; thispt=pt
			#first bring the parameters back into range and set the histograms to use
			if topology==1 or topology==2 :
				thiseta = abs(thiseta)
				if thiseta<self.__ele_trig_eff_r_eta_low : 
					thiseta = self.__ele_trig_eff_r_eta_low 
					doubleerr=True
				elif thiseta>self.__ele_trig_eff_r_eta_hi : 
					thiseta = self.__ele_trig_eff_r_eta_hi - 0.0000001
					doubleerr=True
				if thispt<self.__ele_trig_eff_r_pt_low : 
					thispt = self.__ele_trig_eff_r_pt_low 
					doubleerr=True
				elif thispt>self.__ele_trig_eff_r_pt_hi : 
					thispt = self.__ele_trig_eff_r_pt_hi - 0.0000001
					doubleerr=True
				btofhisto = self.__ele_trig_eff_boosted
				ghhisto = self.__ele_trig_eff_boosted
				#get factors
				binx_BtoF = btofhisto.GetXaxis().FindFixBin(thiseta)
				biny_BtoF = btofhisto.GetYaxis().FindFixBin(thispt)
				binx_GH = ghhisto.GetXaxis().FindFixBin(thiseta)
				biny_GH = ghhisto.GetYaxis().FindFixBin(thispt)
			elif topology==3 :
				if thiseta<self.__ele_trig_eff_r_eta_low : 
					thiseta = self.__ele_trig_eff_r_eta_low 
					doubleerr=True
				elif thiseta>self.__ele_trig_eff_r_eta_hi : 
					thiseta = self.__ele_trig_eff_r_eta_hi - 0.0000001
					doubleerr=True
				if thispt<self.__ele_trig_eff_r_pt_low : 
					thispt = self.__ele_trig_eff_r_pt_low 
					doubleerr=True
				elif thispt>self.__ele_trig_eff_r_pt_hi : 
					thispt = self.__ele_trig_eff_r_pt_hi - 0.0000001
					doubleerr=True
				btofhisto = self.__ele_trig_eff_resolved_BtoF
				ghhisto = self.__ele_trig_eff_resolved_GH
				#get factors
				binx_BtoF = btofhisto.GetXaxis().FindFixBin(thispt)
				biny_BtoF = btofhisto.GetYaxis().FindFixBin(thiseta)
				binx_GH = ghhisto.GetXaxis().FindFixBin(thispt)
				biny_GH = ghhisto.GetYaxis().FindFixBin(thiseta)
			nomfac_BtoF = btofhisto.GetBinContent(btofhisto.GetBin(binx_BtoF,biny_BtoF))
			err_up_BtoF = btofhisto.GetBinErrorUp(binx_BtoF,biny_BtoF)
			err_dn_BtoF = btofhisto.GetBinErrorLow(binx_BtoF,biny_BtoF)
			nomfac_GH = ghhisto.GetBinContent(ghhisto.GetBin(binx_GH,biny_GH))
			err_up_GH = ghhisto.GetBinErrorUp(binx_GH,biny_GH)
			err_dn_GH = ghhisto.GetBinErrorLow(binx_GH,biny_GH)
			if doubleerr :
				err_up_BtoF*=2; err_dn_BtoF*=2; err_up_GH*=2; err_dn_GH*=2
			#calculate total factors
			upfac_BtoF 	 = nomfac_BtoF+err_up_BtoF; upfac_GH   = nomfac_GH+err_up_GH
			downfac_BtoF = nomfac_BtoF-err_dn_BtoF; downfac_GH = nomfac_GH-err_dn_GH
			#print ( 'electron trig. SF_nom=%.3f, err=%.3f (%.4f%%), pt=%.1f(%.1f), abseta=%.2f(%.2f), binx=%d, biny=%d, doubleerr=%s'
			#		%(nomfac_BtoF, err_up, 100.*(err_up/nomfac_BtoF), thispt, pt, thiseta, eta, binx, biny, doubleerr) ) #DEBUG
		return nomfac_BtoF,upfac_BtoF,downfac_BtoF,nomfac_GH,upfac_GH,downfac_GH

	def getIDEff(self,pileup,lepflav_or_lep,pt=-111.,eta=-111.) :
		nomfac_BtoF = 1.; upfac_BtoF = 1.; downfac_BtoF = 1.
		nomfac_GH = 1.; upfac_GH = 1.; downfac_GH = 1.
		doublepuerr=False; doubleerr=False
		lepflav = lepflav_or_lep
		if pt==-111. :
			lepflav=lepflav_or_lep.getType()
			pt = lepflav_or_lep.getPt()
			eta = lepflav_or_lep.getEtaSC() if lepflav=='el' else lepflav_or_lep.getEta()
		if lepflav=='mu' : #muons
			eta=abs(eta)
			#first bring the parameters back into range
			if eta<self.__muon_id_eff_eta_low : 
				eta = self.__muon_id_eff_eta_low 
				doubleerr=True
			elif eta>self.__muon_id_eff_eta_hi : 
				eta = self.__muon_id_eff_eta_hi - 0.0000001
				doubleerr=True
			if pt<self.__muon_id_eff_pt_low : 
				pt = self.__muon_id_eff_pt_low 
				doubleerr=True
			elif pt>self.__muon_id_eff_pt_hi : 
				pt = self.__muon_id_eff_pt_hi - 0.0000001
				doubleerr=True
			pu = pileup
			if pu<self.__muon_id_eff_pu_low : 
				pu = self.__muon_id_eff_pu_low 
				doublepuerr=True
			elif pu>self.__muon_id_eff_pu_hi : 
				pu = self.__muon_id_eff_pu_hi - 0.0000001
				doublepuerr=True
			#get factors
			binx = self.__muon_id_eff_BtoF_abseta_vs_pt.GetXaxis().FindFixBin(pt)
			biny = self.__muon_id_eff_BtoF_abseta_vs_pt.GetYaxis().FindFixBin(eta)
			idfac_BtoF = self.__muon_id_eff_BtoF_abseta_vs_pt.GetBinContent(self.__muon_id_eff_BtoF_abseta_vs_pt.GetBin(binx,biny))
			iderr_up_BtoF = self.__muon_id_eff_BtoF_abseta_vs_pt.GetBinErrorUp(binx,biny)
			iderr_dn_BtoF = self.__muon_id_eff_BtoF_abseta_vs_pt.GetBinErrorLow(binx,biny)
			binx = self.__muon_id_eff_GH_abseta_vs_pt.GetXaxis().FindFixBin(pt)
			biny = self.__muon_id_eff_GH_abseta_vs_pt.GetYaxis().FindFixBin(eta)
			idfac_GH = self.__muon_id_eff_GH_abseta_vs_pt.GetBinContent(self.__muon_id_eff_GH_abseta_vs_pt.GetBin(binx,biny))
			iderr_up_GH = self.__muon_id_eff_GH_abseta_vs_pt.GetBinErrorUp(binx,biny)
			iderr_dn_GH = self.__muon_id_eff_GH_abseta_vs_pt.GetBinErrorLow(binx,biny)
			if idfac_BtoF==0.:
				idfac_BtoF=0.0000001
			if idfac_GH==0.:
				idfac_GH=0.0000001
			if doubleerr :
				iderr_up_BtoF*=2; iderr_dn_BtoF*=2; iderr_up_GH*=2; iderr_dn_GH*=2
			pubin_BtoF 		 = self.__muon_id_eff_BtoF_vs_pu.FindFixBin(pu)
			pufac_BtoF 		 = self.__muon_id_eff_BtoF_vs_pu.GetBinContent(pubin_BtoF)
			puerr_up_BtoF 	 = self.__muon_id_eff_BtoF_vs_pu.GetBinErrorUp(pubin_BtoF)
			puerr_down_BtoF  = self.__muon_id_eff_BtoF_vs_pu.GetBinErrorLow(pubin_BtoF)
			pubin_GH 		 = self.__muon_id_eff_GH_vs_pu.FindFixBin(pu)
			pufac_GH 		 = self.__muon_id_eff_GH_vs_pu.GetBinContent(pubin_GH)
			puerr_up_GH 	 = self.__muon_id_eff_GH_vs_pu.GetBinErrorUp(pubin_GH)
			puerr_down_GH 	 = self.__muon_id_eff_GH_vs_pu.GetBinErrorLow(pubin_GH)
			if doublepuerr :
				puerr_up_BtoF*=2; puerr_down_BtoF*=2; puerr_up_GH*=2; puerr_down_GH*=2
			#get the tracking efficiency and errors
			trkfac, trkerrup, trkerrdown = self.__getTrkEff__(pileup,lepflav,pt,eta)
			#calculate total factors
			nomfac_BtoF  = idfac_BtoF*pufac_BtoF*trkfac
			#print 'idfac_BtoF = %.4f, pufac_BtoF = %.4f, trkfac = %.4f'%(idfac_BtoF,pufac_BtoF,trkfac) #DEBUG
			upfac_BtoF 	 = nomfac_BtoF*(1.+sqrt((iderr_up_BtoF/idfac_BtoF)**2+(puerr_up_BtoF/pufac_BtoF)**2+(trkerrup/trkfac)**2))
			downfac_BtoF = nomfac_BtoF*(1.-sqrt((iderr_dn_BtoF/idfac_BtoF)**2+(puerr_down_BtoF/pufac_BtoF)**2+(trkerrdown/trkfac)**2))
			nomfac_GH  = idfac_GH*pufac_GH*trkfac
			upfac_GH   = nomfac_GH*(1.+sqrt((iderr_up_GH/idfac_GH)**2+(puerr_up_GH/pufac_GH)**2+(trkerrup/trkfac)**2))
			downfac_GH = nomfac_GH*(1.-sqrt((iderr_dn_GH/idfac_GH)**2+(puerr_down_GH/pufac_GH)**2+(trkerrdown/trkfac)**2))
		elif lepflav=='el' : #electrons
			thiseta=abs(eta); thispt=pt
			#just need to get two numbers with errors for id and reco efficiency
			#ID:
			#first bring the parameters back into range
			if thiseta<self.__ele_id_eff_eta_low : 
				thiseta = self.__ele_id_eff_eta_low 
				doubleerr=True
			elif thiseta>self.__ele_id_eff_eta_hi : 
				thiseta = self.__ele_id_eff_eta_hi - 0.0000001
				doubleerr=True
			if pt<self.__ele_id_eff_pt_low : 
				thispt = self.__ele_id_eff_pt_low 
				doubleerr=True
			elif pt>self.__ele_id_eff_pt_hi : 
				thispt = self.__ele_id_eff_pt_hi - 0.0000001
				doubleerr=True
			#get factors
			binx = self.__ele_id_eff_eta_vs_pt.GetXaxis().FindFixBin(thispt)
			biny = self.__ele_id_eff_eta_vs_pt.GetYaxis().FindFixBin(thiseta)
			idfac = self.__ele_id_eff_eta_vs_pt.GetBinContent(self.__ele_id_eff_eta_vs_pt.GetBin(binx,biny))
			iderr_up = self.__ele_id_eff_eta_vs_pt.GetBinErrorUp(binx,biny)
			iderr_dn = self.__ele_id_eff_eta_vs_pt.GetBinErrorLow(binx,biny)
			if doubleerr :
				iderr_up*=2; iderr_dn*=2
			doubleerr=False
			#Reco:
			#first bring the parameters back into range
			thiseta=eta; thispt=pt
			#print 'pt=%.1f, eta=%.2f, %.1f<pT<%.1f, %.2f<eta<%.2f, thispt=%.1f, thiseta=%.2f'%(pt,eta,self.__ele_reco_eff_pt_low,self.__ele_reco_eff_pt_hi,self.__ele_reco_eff_eta_low,self.__ele_reco_eff_eta_hi,thispt,thiseta) #DEBUG
			if eta<self.__ele_reco_eff_eta_low : 
				thiseta = self.__ele_reco_eff_eta_low 
				doubleerr=True
			elif eta>self.__ele_reco_eff_eta_hi : 
				thiseta = self.__ele_reco_eff_eta_hi - 0.0000001
				doubleerr=True
			if pt<self.__ele_reco_eff_pt_low : 
				thispt = self.__ele_reco_eff_pt_low 
				doubleerr=True
			elif pt>self.__ele_reco_eff_pt_hi : 
				thispt = self.__ele_reco_eff_pt_hi - 0.0000001
				doubleerr=True
			#get factors
			binx = self.__ele_reco_eff_pt_vs_eta.GetXaxis().FindFixBin(thiseta)
			biny = self.__ele_reco_eff_pt_vs_eta.GetYaxis().FindFixBin(thispt)
			recofac = self.__ele_reco_eff_pt_vs_eta.GetBinContent(self.__ele_reco_eff_pt_vs_eta.GetBin(binx,biny))
			recoerr_up = self.__ele_reco_eff_pt_vs_eta.GetBinErrorUp(binx,biny)
			recoerr_dn = self.__ele_reco_eff_pt_vs_eta.GetBinErrorLow(binx,biny)
			if recofac==0. :
				recofac=1.; recoerr_up=0.; recoerr_dn=0.
			if doubleerr :
				recoerr_up*=2; recoerr_dn*=2
			#print 'recofac=%.3f, errup=%.3f, errdown=%.3f, pt=%f, eta = %f (bin=%d, binx=%d, biny=%d)'%(recofac,recoerr_up,recoerr_dn,thispt,thiseta,self.__ele_reco_eff_pt_vs_eta.GetBin(binx,biny),binx,biny) #DEBUG
			if recofac==0. : recofac=1.
			#calculate total factors
			nomfac_BtoF  = idfac*recofac
			upfac_BtoF 	 = nomfac_BtoF*(1.+sqrt((iderr_up/idfac)**2+(recoerr_up/recofac)**2))
			downfac_BtoF = nomfac_BtoF*(1.-sqrt((iderr_dn/idfac)**2+(recoerr_dn/recofac)**2))
			nomfac_GH 	 = nomfac_BtoF
			upfac_GH 	 = upfac_BtoF
			downfac_GH 	 = downfac_BtoF
		#return reweighting factors
		return nomfac_BtoF,upfac_BtoF,downfac_BtoF,nomfac_GH,upfac_GH,downfac_GH

	def getIsoEff(self,pileup,topology,lepflav_or_lep,pt=-111.,eta=-111.) :
		nomfac_BtoF = 1.; upfac_BtoF = 1.; downfac_BtoF = 1.
		nomfac_GH = 1.; upfac_GH = 1.; downfac_GH = 1.
		doubleerr=False; doublepuerr=False
		lepflav = lepflav_or_lep
		if pt==-111. :
			lepflav=lepflav_or_lep.getType()
			pt = lepflav_or_lep.getPt()
			eta = lepflav_or_lep.getEta()
		if topology==3 and lepflav=='mu' : #muons in the resolved channel
			#first bring the parameters back into range
			thiseta=abs(eta)
			if thiseta<self.__muon_iso_eff_eta_low : 
				thiseta = self.__muon_iso_eff_eta_low 
				doubleerr=True
			elif thiseta>self.__muon_iso_eff_eta_hi : 
				thiseta = self.__muon_iso_eff_eta_hi - 0.0000001
				doubleerr=True
			if pt<self.__muon_iso_eff_pt_low : 
				pt = self.__muon_iso_eff_pt_low 
				doubleerr=True
			elif pt>self.__muon_iso_eff_pt_hi : 
				pt = self.__muon_iso_eff_pt_hi - 0.0000001
				doubleerr=True
			pu = pileup
			if pu<self.__muon_iso_eff_pu_low : 
				pu = self.__muon_iso_eff_pu_low 
				doublepuerr=True
			elif pu>self.__muon_iso_eff_pu_hi : 
				pu = self.__muon_iso_eff_pu_hi - 0.0000001
				doublepuerr=True
			#get factors
			binx = self.__muon_iso_eff_BtoF_abseta_vs_pt.GetXaxis().FindFixBin(pt)
			biny = self.__muon_iso_eff_BtoF_abseta_vs_pt.GetYaxis().FindFixBin(thiseta)
			isofac_BtoF = self.__muon_iso_eff_BtoF_abseta_vs_pt.GetBinContent(self.__muon_iso_eff_BtoF_abseta_vs_pt.GetBin(binx,biny))
			isoerr_up_BtoF = self.__muon_iso_eff_BtoF_abseta_vs_pt.GetBinErrorUp(binx,biny)
			isoerr_dn_BtoF = self.__muon_iso_eff_BtoF_abseta_vs_pt.GetBinErrorLow(binx,biny)
			binx = self.__muon_iso_eff_GH_abseta_vs_pt.GetXaxis().FindFixBin(pt)
			biny = self.__muon_iso_eff_GH_abseta_vs_pt.GetYaxis().FindFixBin(thiseta)
			isofac_GH = self.__muon_iso_eff_GH_abseta_vs_pt.GetBinContent(self.__muon_iso_eff_GH_abseta_vs_pt.GetBin(binx,biny))
			isoerr_up_GH = self.__muon_iso_eff_GH_abseta_vs_pt.GetBinErrorUp(binx,biny)
			isoerr_dn_GH = self.__muon_iso_eff_GH_abseta_vs_pt.GetBinErrorLow(binx,biny)
			if isofac_BtoF==0.:
				isofac_BtoF=0.0000001
			if isofac_GH==0.:
				isofac_GH=0.0000001
			if doubleerr :
				isoerr_up_BtoF*=2; isoerr_dn_BtoF*=2; isoerr_up_GH*=2; isoerr_dn_GH*=2
			pubin_BtoF 		 = self.__muon_iso_eff_BtoF_vs_pu.FindFixBin(pu)
			pufac_BtoF 		 = self.__muon_iso_eff_BtoF_vs_pu.GetBinContent(pubin_BtoF)
			puerr_up_BtoF 	 = self.__muon_iso_eff_BtoF_vs_pu.GetBinErrorUp(pubin_BtoF)
			puerr_down_BtoF  = self.__muon_iso_eff_BtoF_vs_pu.GetBinErrorLow(pubin_BtoF)
			pubin_GH 		 = self.__muon_iso_eff_GH_vs_pu.FindFixBin(pu)
			pufac_GH 		 = self.__muon_iso_eff_GH_vs_pu.GetBinContent(pubin_GH)
			puerr_up_GH 	 = self.__muon_iso_eff_GH_vs_pu.GetBinErrorUp(pubin_GH)
			puerr_down_GH 	 = self.__muon_iso_eff_GH_vs_pu.GetBinErrorLow(pubin_GH)
			if doublepuerr :
				puerr_up_BtoF*=2; puerr_down_BtoF*=2; puerr_up_GH*=2; puerr_down_GH*=2
			#get the tracking efficiency and errors
			trkfac, trkerrup, trkerrdown = self.__getTrkEff__(pileup,lepflav,pt,eta)
			#calculate total factors
			nomfac_BtoF = isofac_BtoF*pufac_BtoF*trkfac
			upfac_BtoF 	 = nomfac_BtoF*(1.+sqrt((isoerr_up_BtoF/isofac_BtoF)**2+(puerr_up_BtoF/pufac_BtoF)**2+(trkerrup/trkfac)**2))
			downfac_BtoF = nomfac_BtoF*(1.-sqrt((isoerr_dn_BtoF/isofac_BtoF)**2+(puerr_down_BtoF/pufac_BtoF)**2+(trkerrdown/trkfac)**2))
			nomfac_GH = isofac_GH*pufac_GH*trkfac
			upfac_GH 	 = nomfac_GH*(1.+sqrt((isoerr_up_GH/isofac_GH)**2+(puerr_up_GH/pufac_GH)**2+(trkerrup/trkfac)**2))
			downfac_GH = nomfac_GH*(1.-sqrt((isoerr_dn_GH/isofac_GH)**2+(puerr_down_GH/pufac_GH)**2+(trkerrdown/trkfac)**2))
		elif topology==3 and lepflav=='el' : #electrons
			#NOT IMPLEMENTED YET!!
			nomfac_BtoF = 1.; upfac_BtoF = 1.; downfac_BtoF = 1.
			nomfac_GH = 1.; upfac_GH = 1.; downfac_GH = 1.
		#return reweighting factors
		return nomfac_BtoF,upfac_BtoF,downfac_BtoF,nomfac_GH,upfac_GH,downfac_GH

	def getMiniIsoEff(self,pileup,lepflav_or_lep,pt=-111.,eta=-111.) :
		nomfac = 1.; upfac = 1.; downfac = 1.
		doubleerr=False
		lepflav = lepflav_or_lep
		if pt==-111. :
			lepflav=lepflav_or_lep.getType()
			pt = lepflav_or_lep.getPt()
			eta = lepflav_or_lep.getEtaSC() if lepflav=='el' else lepflav_or_lep.getEta()
		if lepflav=='mu' : #muons
			#NOT IMPLEMENTED YET!!
			nomfac = 1.; upfac = 1.; downfac = 1.
		elif lepflav=='el' : #electrons
			thiseta=abs(eta); thispt=pt
			#just need to get two numbers with errors for id and reco efficiency
			if thiseta<self.__ele_miniiso_eff_eta_low : 
				thiseta = self.__ele_miniiso_eff_eta_low 
				doubleerr=True
			elif thiseta>self.__ele_miniiso_eff_eta_hi : 
				thiseta = self.__ele_miniiso_eff_eta_hi - 0.0000001
				doubleerr=True
			if pt<self.__ele_miniiso_eff_pt_low : 
				thispt = self.__ele_miniiso_eff_pt_low 
				doubleerr=True
			elif pt>self.__ele_miniiso_eff_pt_hi : 
				thispt = self.__ele_miniiso_eff_pt_hi - 0.0000001
				doubleerr=True
			#get factors
			binx = self.__ele_miniiso_eff_abseta_vs_pt.GetXaxis().FindFixBin(thispt)
			biny = self.__ele_miniiso_eff_abseta_vs_pt.GetYaxis().FindFixBin(thiseta)
			fac = self.__ele_miniiso_eff_abseta_vs_pt.GetBinContent(self.__ele_miniiso_eff_abseta_vs_pt.GetBin(binx,biny))
			err_up = self.__ele_miniiso_eff_abseta_vs_pt.GetBinErrorUp(binx,biny)
			err_dn = self.__ele_miniiso_eff_abseta_vs_pt.GetBinErrorLow(binx,biny)
			if doubleerr :
				err_up*=2; err_dn*=2
			#calculate total factors
			nomfac  = fac
			upfac 	 = nomfac*(1.+err_up)
			downfac = nomfac*(1.-err_dn)
		#return reweighting factors
		return nomfac,upfac,downfac

	def getBTagEff(self,topology,jets) :
		nomfac = 1.; upfac = 0.; downfac = 0.
		#set the histograms and reader to use based on topology (loose or medium b-tags?)
		udsg_histo = self.__btag_MC_eff_udsg_histo_l
		cjet_histo = self.__btag_MC_eff_cjet_histo_l
		bjet_histo = self.__btag_MC_eff_bjet_histo_l
		reader = self.__reader_l
		if topology==3 :
			udsg_histo = self.__btag_MC_eff_udsg_histo_m
			cjet_histo = self.__btag_MC_eff_cjet_histo_m
			bjet_histo = self.__btag_MC_eff_bjet_histo_m
			reader = self.__reader_m
		#loop over all the jets
		for jet in jets :
			#get the MC btagging efficiency and scalefactors https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration
			flav = jet.getFlavor()
			pt = jet.getPt()
			apt = max(min(pt,self.__btag_MC_eff_pt_hi),self.__btag_MC_eff_pt_low)
			eta = jet.getEta()
			aeta = max(min(abs(eta),self.__btag_MC_eff_eta_hi),self.__btag_MC_eff_eta_low)
			jetistagged = jet.isMbTagged() if topology==3 else jet.isLbTagged()
			mc_eff = 1.; sf=1.; sfup=1.; sfdown=1.;
			if flav==0 :
				mc_eff = udsg_histo.GetBinContent(udsg_histo.FindFixBin(aeta,apt))
				if mc_eff==0. : mc_eff=udsg_histo.GetMaximum()
				sf = reader.eval_auto_bounds('central', # systematic (here also 'up'/'down' possible)
												2, 	   # jet flavor
												eta, 	   # eta
												pt) 	   # pt
				sfup = reader.eval_auto_bounds('up',2,eta,pt)
				sfdown = reader.eval_auto_bounds('down',2,eta,pt)
			elif flav==4 :
				mc_eff = cjet_histo.GetBinContent(cjet_histo.FindFixBin(aeta,apt))
				if mc_eff==0. : mc_eff=cjet_histo.GetMaximum()
				sf = reader.eval_auto_bounds('central', # systematic (here also 'up'/'down' possible)
												1, 	   # jet flavor
												eta, 	   # eta
												pt) 	   # pt
				sfup = reader.eval_auto_bounds('up',1,eta,pt)
				sfdown = reader.eval_auto_bounds('down',1,eta,pt)
			elif flav==5 :
				mc_eff = bjet_histo.GetBinContent(bjet_histo.FindFixBin(aeta,apt))
				if mc_eff==0. : mc_eff=bjet_histo.GetMaximum()
				sf = reader.eval_auto_bounds('central', # systematic (here also 'up'/'down' possible)
												0, 	   # jet flavor
												eta, 	   # eta
												pt) 	   # pt
				sfup = reader.eval_auto_bounds('up',0,eta,pt)
				sfdown = reader.eval_auto_bounds('down',0,eta,pt)
			#if the jet is b-tagged
			if jetistagged :
				#multiply weights just by the sf
				nomfac*=sf; 
				#up/down fractional errors add in quadrature
				uperr = sfup-sf; downerr=sf-sfdown
				upfac+=(uperr/sf)**2; downfac+=(downerr/sf)**2
			#and if not
			else :
				#multiply weights by ratio of (1-sf*eff)/(1-eff)
				newfac = (1.-sf*mc_eff)/(1.-mc_eff)
				nomfac*=newfac
				#up/down fractional errors add in quadrature
				uperr = ((sfup-sf)*mc_eff)/(1.-mc_eff); downerr = ((sf-sfdown)*mc_eff)/(1.-mc_eff)
				upfac+=(uperr/newfac)**2; downfac+=(downerr/newfac)**2
		#adjust the final weights
		upfac = nomfac*(1.+sqrt(upfac)); downfac=nomfac*(1.-sqrt(downfac))
		return nomfac,upfac,downfac

	def getGenReweights(self,branches) :
		returnlist = []
		#mu_R, from scale_weights
		returnlist.append(1.) #nominal
		returnlist.append(branches['scale_Weights'].getReadValue(2)/self.__renormdict['muRup']) #mu_R up
		returnlist.append(branches['scale_Weights'].getReadValue(4)/self.__renormdict['muRdown']) #mu_R down
		#mu_F
		returnlist.append(1.) #nominal
		returnlist.append(branches['scale_Weights'].getReadValue(0)/self.__renormdict['muFup']) #mu_F up
		returnlist.append(branches['scale_Weights'].getReadValue(1)/self.__renormdict['muFdown']) #mu_F down
		#combined mu_R/mu_F ("scale_comb")
		returnlist.append(1.) #nominal
		returnlist.append(branches['scale_Weights'].getReadValue(3)/self.__renormdict['scup']) #mu_R/F up
		returnlist.append(branches['scale_Weights'].getReadValue(5)/self.__renormdict['scdown']) #mu_R/F down
		#PDF and alpha_s ("pdf_alphas")
		#get the average and std dev of the pdf replica reweights
		pdfvaluearray = branches['pdf_Weights'].getReadArray()
		npdfweights = len(pdfvaluearray)
		pdfmean = sum(pdfvaluearray)/npdfweights
		pdfvalues2array = []
		for i in range(npdfweights) :
			pdfvalues2array.append(pdfvaluearray[i]**2)
		pdfunc = sqrt(abs((sum(pdfvalues2array)/npdfweights)-(pdfmean**2)))
		#get the alpha_s values
		alphas_up_unc 	= abs(branches['alphas_Weights'].getReadValue(1)-1.)*0.75 if len(branches['alphas_Weights'].getReadArray())>1 else 0.
		alphas_down_unc = abs(branches['alphas_Weights'].getReadValue(0)-1.)*0.75 if len(branches['alphas_Weights'].getReadArray())>1 else 0.
		#nominal value
		returnlist.append(pdfmean/self.__renormdict['pdfas'])
		#up/down with pdf and alpha_s added in quadrature
		returnlist.append((pdfmean+sqrt(pdfunc**2+alphas_up_unc**2))/self.__renormdict['pdfasup'])
		returnlist.append((pdfmean-sqrt(pdfunc**2+alphas_down_unc**2))/self.__renormdict['pdfasdown'])
		#print 'alphas_Weights[0]=%.4f, alphas_Weights[1]=%.4f'%(branches['alphas_Weights'].getReadValue(0),branches['alphas_Weights'].getReadValue(1))
		#print 'alphas_Weights_size=%d, pdfmean=%.2f, pdfas=%.2f, pdfunc=%.2f, alphas_up_unc=%.2f, alphas_down_unc=%.2f, nom=%.2f, up=%.2f, donw=%.2f'%(branches['alphas_size'].getReadValue(),pdfmean,self.__renormdict['pdfas'],pdfunc,alphas_up_unc,alphas_down_unc,pdfmean/self.__renormdict['pdfas'],(pdfmean+sqrt(pdfunc**2+alphas_up_unc**2))/self.__renormdict['pdfasup'],(pdfmean-sqrt(pdfunc**2+alphas_down_unc**2))/self.__renormdict['pdfasdown'])
		#return the list, made into a tuple so it's immutable
		return tuple(returnlist)

	def getAlphaEpsilon(self,beta) :
		alpha = 0.0; epsilon = 0.0
		for i in range(len(self.__alphalist)) :
			if (i==0 and beta<self.__alphalist[i][0]) or (i==len(self.__alphalist)-1 and beta>self.__alphalist[i][0]) or (beta>self.__alphalist[i][0] and beta<self.__alphalist[i+1][0]) :
				alpha = self.__alphalist[i][1]
				break
		for i in range(len(self.__epsilonlist)) :
			if (i==0 and beta<self.__epsilonlist[i][0]) or (i==len(self.__epsilonlist)-1 and beta>self.__epsilonlist[i][0]) or (beta>self.__epsilonlist[i][0] and beta<self.__epsilonlist[i+1][0]) :
				epsilon = self.__epsilonlist[i][1]
				break
		#print 'beta = %.4f, alpha = %.4f, epsilon = %.4f'%(beta,alpha,epsilon) #DEBUG
		return alpha, epsilon

	def __del__(self) :
		pass

	#-------------------------------- object functions to set up corrector histograms -------------------------------- 

	def setupPileupHistos(self,onGrid,pu_histo) :
		if pu_histo.Integral() > 0 : pu_histo.Scale(1./pu_histo.Integral())
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		self._data_nom_file  = TFile.Open(pp+DATA_PU_HISTO_NOMINAL_FILENAME)
		data_nom_histo  = self._data_nom_file.Get('pileup')
		data_nom_histo.Scale(1./data_nom_histo.Integral())
		self._data_up_file   = TFile.Open(pp+DATA_PU_HISTO_XS_UP_FILENAME)
		data_up_histo   = self._data_up_file.Get('pileup')
		data_up_histo.Scale(1./data_up_histo.Integral())
		self._data_down_file = TFile.Open(pp+DATA_PU_HISTO_XS_DOWN_FILENAME)
		data_down_histo = self._data_down_file.Get('pileup')
		data_down_histo.Scale(1./data_down_histo.Integral())
		#pu_histo.SetDirectory(0); data_nom_histo.SetDirectory(0); data_up_histo.SetDirectory(0); data_down_histo.SetDirectory(0)
		#data_nom_file.Close(); data_up_file.Close(); data_down_file.Close()
		return pu_histo, data_nom_histo, data_up_histo, data_down_histo

	def setupMuonTriggerHistos(self,onGrid) :
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		self._muon_trig_BtoF_file = TFile.Open(pp+MUON_TRIG_EFF_BTOF_ROOT_FILENAME)
		self._muon_trig_GH_file 	= TFile.Open(pp+MUON_TRIG_EFF_GH_ROOT_FILENAME)
		BtoF_pt_histo_b  = self._muon_trig_BtoF_file.Get(MUON_TRIG_EFF_BOOSTED_PT_HISTONAME)
		BtoF_eta_histo_b = self._muon_trig_BtoF_file.Get(MUON_TRIG_EFF_BOOSTED_ETA_HISTONAME)
		GH_pt_histo_b 	 = self._muon_trig_GH_file.Get(MUON_TRIG_EFF_BOOSTED_PT_HISTONAME)
		GH_eta_histo_b  = self._muon_trig_GH_file.Get(MUON_TRIG_EFF_BOOSTED_ETA_HISTONAME)
		BtoF_pt_histo_r  = self._muon_trig_BtoF_file.Get(MUON_TRIG_EFF_RESOLVED_PT_HISTONAME)
		BtoF_eta_histo_r = self._muon_trig_BtoF_file.Get(MUON_TRIG_EFF_RESOLVED_ETA_HISTONAME)
		GH_pt_histo_r 	 = self._muon_trig_GH_file.Get(MUON_TRIG_EFF_RESOLVED_PT_HISTONAME)
		GH_eta_histo_r  = self._muon_trig_GH_file.Get(MUON_TRIG_EFF_RESOLVED_ETA_HISTONAME)
		#BtoF_pt_histo_b.SetDirectory(0); BtoF_eta_histo_b.SetDirectory(0);
		#GH_pt_histo_b.SetDirectory(0); GH_eta_histo_b.SetDirectory(0);
		#BtoF_pt_histo_r.SetDirectory(0); BtoF_eta_histo_r.SetDirectory(0);
		#GH_pt_histo_r.SetDirectory(0); GH_eta_histo_r.SetDirectory(0);
		#muon_trig_BtoF_file.Close(); muon_trig_GH_file.Close()
		return BtoF_pt_histo_b, BtoF_eta_histo_b, GH_pt_histo_b, GH_eta_histo_b, BtoF_pt_histo_r, BtoF_eta_histo_r, GH_pt_histo_r, GH_eta_histo_r

	def setupEleTriggerHistos(self,onGrid) :
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		self._ele_trig_b_file = TFile.Open(pp+ELE_TRIG_EFF_BOOSTED_ROOT_FILENAME)
		self._ele_trig_r_BtoF_file = TFile.Open(pp+ELE_TRIG_EFF_RESOLVED_BTOF_ROOT_FILENAME)
		self._ele_trig_r_GH_file = TFile.Open(pp+ELE_TRIG_EFF_RESOLVED_GH_ROOT_FILENAME)
		histo_b = self._ele_trig_b_file.Get(ELE_TRIG_EFF_BOOSTED_HISTONAME)
		btof_histo_r = self._ele_trig_r_BtoF_file.Get(ELE_TRIG_EFF_RESOLVED_HISTONAME)
		gh_histo_r = self._ele_trig_r_GH_file.Get(ELE_TRIG_EFF_RESOLVED_HISTONAME)
		#histo_b.SetDirectory(0); btof_histo_r.SetDirectory(0); gh_histo_r.SetDirectory(0)
		#ele_trig_b_file.Close(); ele_trig_r_BtoF_file.Close(); ele_trig_r_GH_file.Close()
		return histo_b, btof_histo_r, gh_histo_r

	def setupMuonIDHistos(self,onGrid) :
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		self._muon_id_BtoF_file = TFile.Open(pp+MUON_ID_EFF_BTOF_ROOT_FILENAME)
		self._muon_id_GH_file   = TFile.Open(pp+MUON_ID_EFF_GH_ROOT_FILENAME)
		BtoF_abseta_vs_pt_histo  = self._muon_id_BtoF_file.Get(MUON_ID_EFF_BTOF_PT_ETA_HISTONAME)
		BtoF_vs_pu_histo  = self._muon_id_BtoF_file.Get(MUON_ID_EFF_BTOF_PU_HISTONAME)
		GH_abseta_vs_pt_histo  = self._muon_id_GH_file.Get(MUON_ID_EFF_GH_PT_ETA_HISTONAME)
		GH_vs_pu_histo  = self._muon_id_GH_file.Get(MUON_ID_EFF_GH_PU_HISTONAME)
		#BtoF_abseta_vs_pt_histo.SetDirectory(0); BtoF_vs_pu_histo.SetDirectory(0)
		#GH_abseta_vs_pt_histo.SetDirectory(0); GH_vs_pu_histo.SetDirectory(0)
		#muon_id_BtoF_file.Close()
		#muon_id_GH_file.Close()
		return BtoF_abseta_vs_pt_histo, BtoF_vs_pu_histo, GH_abseta_vs_pt_histo, GH_vs_pu_histo

	def setupEleIDHistos(self,onGrid) :
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		self._ele_id_file = TFile.Open(pp+ELE_ID_EFF_ROOT_FILENAME)
		self._ele_reco_file   = TFile.Open(pp+ELE_RECO_EFF_ROOT_FILENAME)
		ele_id_histo = self._ele_id_file.Get(ELE_ID_EFF_2D_HISTONAME)
		ele_reco_histo = self._ele_reco_file.Get(ELE_RECO_EFF_2D_HISTONAME)
		#ele_id_histo.SetDirectory(0); ele_reco_histo.SetDirectory(0)
		#ele_id_file.Close(); ele_reco_file.Close()
		return ele_id_histo, ele_reco_histo

	def setupMuonIsoHistos(self,onGrid) :
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		self._muon_iso_BtoF_file = TFile.Open(pp+MUON_ISO_EFF_BTOF_ROOT_FILENAME)
		self._muon_iso_GH_file   = TFile.Open(pp+MUON_ISO_EFF_GH_ROOT_FILENAME)
		BtoF_abseta_vs_pt_histo  = self._muon_iso_BtoF_file.Get(MUON_ISO_EFF_PT_ETA_HISTONAME)
		BtoF_vs_pu_histo  = self._muon_iso_BtoF_file.Get(MUON_ISO_EFF_PU_HISTONAME)
		GH_abseta_vs_pt_histo  = self._muon_iso_GH_file.Get(MUON_ISO_EFF_PT_ETA_HISTONAME)
		GH_vs_pu_histo  = self._muon_iso_GH_file.Get(MUON_ISO_EFF_PU_HISTONAME)
		#BtoF_abseta_vs_pt_histo.SetDirectory(0); BtoF_vs_pu_histo.SetDirectory(0)
		#GH_abseta_vs_pt_histo.SetDirectory(0); GH_vs_pu_histo.SetDirectory(0)
		#muon_iso_BtoF_file.Close()
		#muon_iso_GH_file.Close()
		return BtoF_abseta_vs_pt_histo, BtoF_vs_pu_histo, GH_abseta_vs_pt_histo, GH_vs_pu_histo

	def setupEleMiniIsoHistos(self,onGrid) :
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		self._ele_miniiso_file = TFile.Open(pp+ELE_MINIISO_EFF_ROOT_FILENAME)
		ele_miniiso_histo = self._ele_miniiso_file.Get(ELE_MINIISO_EFF_2D_HISTONAME)
		#ele_miniiso_histo.SetDirectory(0)
		#ele_miniiso_file.Close()
		return ele_miniiso_histo

	def setupMuonTrackingGraphs(self,onGrid) :
		pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
		self._muon_trk_file = TFile.Open(pp+MUON_TRK_EFF_ROOT_FILENAME)
		vs_eta_graph  = self._muon_trk_file.Get(MUON_TRK_EFF_ETA_GRAPHNAME)
		vs_pu_graph   = self._muon_trk_file.Get(MUON_TRK_EFF_PU_GRAPHNAME)
		#muon_trk_file.Close()
		return vs_eta_graph, vs_pu_graph


def setupJECCorrector(onGrid,isdata,jetType,runera) :
	#Define the JEC  parameters
	L1JetPar  = None
	L2JetPar  = None
	L3JetPar  = None
	pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
	#Get the files linked below from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
	if not isdata :
		L1JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016V4_MC_L1FastJet_'+jetType+'.txt')
		L2JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016V4_MC_L2Relative_'+jetType+'.txt')
		L3JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016V4_MC_L3Absolute_'+jetType+'.txt')
		jetUncertainty = JetCorrectionUncertainty(pp+'Summer16_23Sep2016V4_MC_Uncertainty_'+jetType+'.txt')
	else :
		if runera=='B' or runera=='C' or runera=='D' :
			L1JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016BCDV4_DATA_L1FastJet_'+jetType+'.txt')
			L2JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016BCDV4_DATA_L2Relative_'+jetType+'.txt')
			L3JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016BCDV4_DATA_L3Absolute_'+jetType+'.txt')
			ResJetPar = JetCorrectorParameters(pp+'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_'+jetType+'.txt')
			jetUncertainty = JetCorrectionUncertainty(pp+'Summer16_23Sep2016BCDV4_DATA_Uncertainty_'+jetType+'.txt')
		elif runera=='E' or runera=='F' :
			L1JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016EFV4_DATA_L1FastJet_'+jetType+'.txt')
			L2JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016EFV4_DATA_L2Relative_'+jetType+'.txt')
			L3JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016EFV4_DATA_L3Absolute_'+jetType+'.txt')
			ResJetPar = JetCorrectorParameters(pp+'Summer16_23Sep2016EFV4_DATA_L2L3Residual_'+jetType+'.txt')
			jetUncertainty = JetCorrectionUncertainty(pp+'Summer16_23Sep2016EFV4_DATA_Uncertainty_'+jetType+'.txt')
		elif runera=='G' :
			L1JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016GV4_DATA_L1FastJet_'+jetType+'.txt')
			L2JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016GV4_DATA_L2Relative_'+jetType+'.txt')
			L3JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016GV4_DATA_L3Absolute_'+jetType+'.txt')
			ResJetPar = JetCorrectorParameters(pp+'Summer16_23Sep2016GV4_DATA_L2L3Residual_'+jetType+'.txt')
			jetUncertainty = JetCorrectionUncertainty(pp+'Summer16_23Sep2016GV4_DATA_Uncertainty_'+jetType+'.txt')
		elif runera=='H' :
			L1JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016HV4_DATA_L1FastJet_'+jetType+'.txt')
			L2JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016HV4_DATA_L2Relative_'+jetType+'.txt')
			L3JetPar  = JetCorrectorParameters(pp+'Summer16_23Sep2016HV4_DATA_L3Absolute_'+jetType+'.txt')
			ResJetPar = JetCorrectorParameters(pp+'Summer16_23Sep2016HV4_DATA_L2L3Residual_'+jetType+'.txt')
			jetUncertainty = JetCorrectionUncertainty(pp+'Summer16_23Sep2016HV4_DATA_Uncertainty_'+jetType+'.txt')
		else :
			print 'WARNING: Can\'t recognize Data Run Era based on filename (Run Era variable is '+str(runera)+')! This will crash I think'
			return None, None
	#Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
	vParJec = vector('JetCorrectorParameters')()
	vParJec.push_back(L1JetPar)
	vParJec.push_back(L2JetPar)
	vParJec.push_back(L3JetPar)
	#add the resolution for data
	if isdata : 
		vParJec.push_back(ResJetPar)
	#define the corrector
	jetCorrector = FactorizedJetCorrector(vParJec)
	return jetCorrector, jetUncertainty

def getJER(jetEta, sysType) :
    jerSF = 1.0
    if ( (sysType==0 or sysType==-1 or sysType==1) == False):
        print "ERROR: Can't get JER! use type=0 (nom), -1 (down), +1 (up)"
        return float(jerSF)
    # Values from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    etamin = 	   [0.000,  0.522,  0.783,  1.131,  1.305,  1.740,  1.930,  2.043,  2.322,  2.500,  2.853,  2.964,  3.139]
    etamax = 	   [0.522,  0.783,  1.131,  1.305,  1.740,  1.930,  2.043,  2.322,  2.500,  2.853,  2.964,  3.139,  5.191]
    scale_nom =    [1.1595, 1.1948, 1.1464, 1.1609, 1.1278, 1.1000, 1.1426, 1.1512, 1.2963, 1.3418, 1.7788, 1.1869, 1.1922]
    scale_uncert = [0.0645, 0.0652, 0.0632, 0.1025, 0.0986, 0.1079, 0.1214, 0.1140, 0.2371, 0.2091, 0.2008, 0.1243, 0.1488]
    for iSF in range(0,len(scale_nom)) :
        if abs(jetEta) >= etamin[iSF] and abs(jetEta) < etamax[iSF] :
            if sysType < 0 :
                jerSF = scale_nom[iSF] - scale_uncert[iSF]
            elif sysType > 0 :
                jerSF = scale_nom[iSF] + scale_uncert[iSF]
            else :
                jerSF = scale_nom[iSF]
            break
    return float(jerSF)

def getList(fulldict,pp) :
	#find all the relevant keys
	keys = []
	for key in fulldict :
		if key.find(pp)!=-1 :
			keys.append(key)
	#if there was only one just return this single value
	if len(keys)==1 :
		return [(0.5,fulldict[keys[0]])]
	#otherwise split off the number and add a tuple to the list
	returnlist = []
	for key in keys :
		num = key.split('_')[1]
		returnlist.append((float(num),fulldict[key]))
	#sort the list based on the numbers
	returnlist.sort(key=lambda x: x[0])
	return returnlist