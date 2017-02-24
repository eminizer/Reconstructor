#Global variables
#Pileup
DATA_PU_HISTO_NOMINAL_FILENAME = 'DataPileupHistogram_69200.root'
DATA_PU_HISTO_XS_UP_FILENAME   = 'DataPileupHistogram_72660.root'
DATA_PU_HISTO_XS_DOWN_FILENAME = 'DataPileupHistogram_65740.root'
#Muon ID Efficiency
MUON_ID_EFF_ROOT_FILENAME = 'MuonID_Z_RunBCD_prompt80X_7p65.root'
MUON_ID_VS_ETA_HISTONAME = 'MC_NUM_LooseID_DEN_genTracks_PAR_eta/eta_ratio'
MUON_ID_VS_PT_HISTONAME  = 'MC_NUM_LooseID_DEN_genTracks_PAR_pt_alleta_bin1/pt_ratio'
MUON_ID_VS_PU_HISTONAME  = 'MC_NUM_LooseID_DEN_genTracks_PAR_pt_vtx/tag_nVertices_ratio'
#Electron Tracking efficiency
ELE_TRK_EFF_ROOT_FILENAME = 'egammaEffi_SF2D.root'
ELE_TRK_EFF_HISTO_NAME = 'EGamma_SF2D'

#Imports
from ROOT import vector, FactorizedJetCorrector, JetCorrectorParameters, JetCorrectionUncertainty, TFile
from random import gauss
from math import sqrt
import pickle

class Corrector(object) :

	##################################  #__init__ function  ##################################
	def __init__(self,isdata,eventType,onGrid,pu_histo) :
		self.__ak4JetCorrector, self.__ak4JecUncertainty = setupJECCorrector(onGrid,isdata,'AK4PFchs')
		self.__ak8JetCorrector, self.__ak8JecUncertainty = setupJECCorrector(onGrid,isdata,'AK8PFchs')
		self.__MC_pu_histo, self.__data_pu_histo_nom, self.__data_pu_histo_up, self.__data_pu_histo_down = setupPileupHistos(onGrid,pu_histo)
		self.__muon_id_eff_vs_eta, self.__muon_id_eff_vs_pt, self.__muon_id_eff_vs_pu = setupMuonIDHistos(onGrid)
		self.__muon_id_eff_vs_eta_low = self.__muon_id_eff_vs_eta.GetBinCenter(1); self.__muon_id_eff_vs_eta_hi = self.__muon_id_eff_vs_eta.GetBinCenter(self.__muon_id_eff_vs_eta.GetNbinsX())
		self.__muon_id_eff_vs_pt_low  = self.__muon_id_eff_vs_pt.GetBinCenter(1);  self.__muon_id_eff_vs_pt_hi  = self.__muon_id_eff_vs_pt.GetBinCenter(self.__muon_id_eff_vs_pt.GetNbinsX())
		self.__muon_id_eff_vs_pu_low  = self.__muon_id_eff_vs_pu.GetBinCenter(1);  self.__muon_id_eff_vs_pu_hi  = self.__muon_id_eff_vs_pu.GetBinCenter(self.__muon_id_eff_vs_pu.GetNbinsX())
		self.__ele_trk_eff_2D_histo = setupEleTrkHisto(onGrid)
		self.__ele_trk_eff_pt_low = self.__ele_trk_eff_2D_histo.GetYaxis().GetBinCenter(0)
		self.__ele_trk_eff_pt_hi  = self.__ele_trk_eff_2D_histo.GetYaxis().GetBinCenter(self.__ele_trk_eff_2D_histo.GetNbinsY())
		self.__ele_trk_eff_sceta_low = self.__ele_trk_eff_2D_histo.GetXaxis().GetBinCenter(0)
		self.__ele_trk_eff_sceta_hi  = self.__ele_trk_eff_2D_histo.GetXaxis().GetBinCenter(self.__ele_trk_eff_2D_histo.GetNbinsX())

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
		if jetvec.DeltaR(genJetVec)<dRCheck/2. and abs(jetvec.Pt()-genJetVec.Pt())<3.*ptres : #scaling
			#smear the pt
			recopt = jetvec.Pt()
			deltapt = (recopt-genJetVec.Pt())*(ptsmearfac-1.0)
			ptsmear = max(0.0, (recopt+deltapt)/recopt)
		else : #gaussian smearing
			sigma = sqrt(abs(ptsmearfac**2-1.))*ptres #Note the absolute value here, this is to handle cases where ptsmearfac<1., which are always out of selection range
			ptsmear = 1.+gauss(0.,sigma)
		ptsmearedjet = jetvec*ptsmear
		#return the smeared fourvector
		return ptsmearedjet

	def getPileupReweight(self,pu_value) :
		mcbincontent = self.__MC_pu_histo.GetBinContent(self.__MC_pu_histo.FindFixBin(pu_value))
		sf = self.__data_pu_histo_nom.GetBinContent(self.__data_pu_histo_nom.FindFixBin(pu_value))/mcbincontent
		sf_up = self.__data_pu_histo_up.GetBinContent(self.__data_pu_histo_up.FindFixBin(pu_value))/mcbincontent
		sf_down = self.__data_pu_histo_down.GetBinContent(self.__data_pu_histo_down.FindFixBin(pu_value))/mcbincontent
		return sf, sf_up, sf_down

	def getIDEff(self,pileup,lepton) :
		nomfac = upfac = downfac = 1.
		lepflavor=lepton.getType()
		if lepflavor=='mu' : #muons
			#first bring the parameters back into range
			eta = lepton.getEta()
			if eta<self.__muon_id_eff_vs_eta_low : eta = self.__muon_id_eff_vs_eta_low 
			elif eta>self.__muon_id_eff_vs_eta_hi : eta = self.__muon_id_eff_vs_eta_hi 
			pt = lepton.getPt()
			if pt<self.__muon_id_eff_vs_pt_low : pt = self.__muon_id_eff_vs_pt_low 
			elif pt>self.__muon_id_eff_vs_pt_hi : pt = self.__muon_id_eff_vs_pt_hi 
			pu = pileup
			if pu<self.__muon_id_eff_vs_pu_low : pu = self.__muon_id_eff_vs_pu_low 
			elif pu>self.__muon_id_eff_vs_pu_hi : pu = self.__muon_id_eff_vs_pu_hi 
			#get factors
			etabin = self.__muon_id_eff_vs_eta.FindFixBin(eta)
			etafac = self.__muon_id_eff_vs_eta.GetBinContent(etabin)
			etaerr = self.__muon_id_eff_vs_eta.GetBinError(etabin)
			ptbin = self.__muon_id_eff_vs_pt.FindFixBin(pt)
			ptfac = self.__muon_id_eff_vs_pt.GetBinContent(ptbin)
			pterr = self.__muon_id_eff_vs_pt.GetBinError(ptbin)
			pubin = self.__muon_id_eff_vs_pu.FindFixBin(pu)
			pufac = self.__muon_id_eff_vs_pu.GetBinContent(pubin)
			puerr = self.__muon_id_eff_vs_pu.GetBinError(pubin)
			#calculate total factors
			nomfac 	= etafac*ptfac*pufac
			upfac 	= (etafac+etaerr)*(ptfac+pterr)*(pufac+puerr)
			downfac = (etafac-etaerr)*(ptfac-pterr)*(pufac-puerr)
		elif lepflavor=='el' : #electrons
			#first bring the parameters back into range
			sceta = lepton.getEtaSC()
			if sceta<self.__ele_trk_eff_sceta_low : eta = self.__ele_trk_eff_sceta_low
			elif sceta>self.__ele_trk_eff_sceta_hi : eta = self.__ele_trk_eff_sceta_hi
			pt = lepton.getPt()
			if pt<self.__ele_trk_eff_pt_low : pt = self.__ele_trk_eff_pt_low 
			elif pt>self.__ele_trk_eff_pt_hi : pt = self.__ele_trk_eff_pt_hi
			#get tracking efficiency value
			trkeffbin = self.__ele_trk_eff_2D_histo.FindFixBin(sceta,pt)
			trkeffnom = self.__ele_trk_eff_2D_histo.GetBinContent(trkeffbin)
			trkefferr = self.__ele_trk_eff_2D_histo.GetBinError(trkeffbin)
			#calculate total factors
			nomfac 	= trkeffnom
			upfac 	= trkeffnom+trkefferr
			downfac = trkeffnom-trkefferr
		#return reweighting factors
		return nomfac, upfac, downfac

	def getGenReweights(self,branches) :
		returnlist = []
		#mu_R, from scale_weights
		returnlist.append(1.) #nominal
		returnlist.append(branches['scale_Weights'].getReadValue(2)) #mu_R up
		returnlist.append(branches['scale_Weights'].getReadValue(4)) #mu_R down
		#mu_F
		returnlist.append(1.) #nominal
		returnlist.append(branches['scale_Weights'].getReadValue(0)) #mu_F up
		returnlist.append(branches['scale_Weights'].getReadValue(1)) #mu_F down
		#combined mu_R/mu_F ("scale_comb")
		returnlist.append(1.) #nominal
		returnlist.append(branches['scale_Weights'].getReadValue(3)) #mu_R/F up
		returnlist.append(branches['scale_Weights'].getReadValue(5)) #mu_R/F down
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
		alphas_up_unc 	= abs(branches['alphas_Weights'].getReadValue(1)*0.75-1.)
		alphas_down_unc = abs(branches['alphas_Weights'].getReadValue(0)*0.75-1.)
		#nominal value
		returnlist.append(pdfmean)
		#up/down with pdf and alpha_s added in quadrature
		returnlist.append(pdfmean+sqrt(pdfunc**2+alphas_up_unc**2))
		returnlist.append(pdfmean-sqrt(pdfunc**2+alphas_down_unc**2))
		#return the list, made into a tuple so it's immutable
		return tuple(returnlist)

	def __del__(self) :
		pass

def setupJECCorrector(onGrid,isdata,jetType) :
	#Define the JEC  parameters
	L1JetPar  = None
	L2JetPar  = None
	L3JetPar  = None
	pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
	if not isdata :
		L1JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_MC_L1FastJet_'+jetType+'.txt')
		L2JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_MC_L2Relative_'+jetType+'.txt')
		L3JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_MC_L3Absolute_'+jetType+'.txt')
		jetUncertainty = JetCorrectionUncertainty(pp+'Spring16_25nsV6_MC_Uncertainty_'+jetType+'.txt')
	else :
		L1JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_DATA_L1FastJet_'+jetType+'.txt')
		L2JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_DATA_L2Relative_'+jetType+'.txt')
		L3JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_DATA_L3Absolute_'+jetType+'.txt')
		ResJetPar = JetCorrectorParameters(pp+'Spring16_25nsV6_DATA_L2L3Residual_'+jetType+'.txt')
		jetUncertainty = JetCorrectionUncertainty(pp+'Spring16_25nsV6_DATA_Uncertainty_'+jetType+'.txt')
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
    etamin = [0.0,0.5,0.8,1.1,1.3,1.7,1.9,2.1,2.3,2.5,2.8,3.0,3.2]
    etamax = [0.5,0.8,1.1,1.3,1.7,1.9,2.1,2.3,2.5,2.8,3.0,3.2,5.0]
    scale_nom =    [1.122,1.167,1.168,1.029,1.115,1.041,1.167,1.094,1.168,1.266,1.595,0.998,1.226]
    scale_uncert = [0.026,0.048,0.046,0.066,0.030,0.062,0.086,0.093,0.120,0.132,0.175,0.066,0.145]
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

def setupPileupHistos(onGrid,pu_histo) :
	if pu_histo.Integral() > 0 : pu_histo.Scale(1./pu_histo.Integral())
	pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
	data_nom_file  = TFile(pp+DATA_PU_HISTO_NOMINAL_FILENAME)
	data_nom_histo  = data_nom_file.Get('pileup')
	data_nom_histo.Scale(1./data_nom_histo.Integral())
	data_up_file   = TFile(pp+DATA_PU_HISTO_XS_UP_FILENAME)
	data_up_histo   = data_up_file.Get('pileup')
	data_up_histo.Scale(1./data_up_histo.Integral())
	data_down_file = TFile(pp+DATA_PU_HISTO_XS_DOWN_FILENAME)
	data_down_histo = data_down_file.Get('pileup')
	data_down_histo.Scale(1./data_down_histo.Integral())
	pu_histo.SetDirectory(0); data_nom_histo.SetDirectory(0); data_up_histo.SetDirectory(0); data_down_histo.SetDirectory(0)
	return pu_histo, data_nom_histo, data_up_histo, data_down_histo

def setupMuonIDHistos(onGrid) :
	pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
	muon_id_file = TFile(pp+MUON_ID_EFF_ROOT_FILENAME)
	vs_eta_histo = muon_id_file.Get(MUON_ID_VS_ETA_HISTONAME)
	vs_pt_histo  = muon_id_file.Get(MUON_ID_VS_PT_HISTONAME)
	vs_pu_histo  = muon_id_file.Get(MUON_ID_VS_PU_HISTONAME)
	vs_eta_histo.SetDirectory(0); vs_pt_histo.SetDirectory(0); vs_pu_histo.SetDirectory(0)
	return vs_eta_histo, vs_pt_histo, vs_pu_histo

def setupEleTrkHisto(onGrid) :
	pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
	ele_id_file = TFile(pp+ELE_TRK_EFF_ROOT_FILENAME)
	histo = ele_id_file.Get(ELE_TRK_EFF_HISTO_NAME)
	histo.SetDirectory(0)
	return histo