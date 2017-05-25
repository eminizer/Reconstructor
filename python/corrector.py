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
	def __init__(self,isdata,onGrid,pu_histo,runera,renormdict) :
		self.__ak4JetCorrector, self.__ak4JecUncertainty = setupJECCorrector(onGrid,isdata,'AK4PFchs',runera)
		self.__ak8JetCorrector, self.__ak8JecUncertainty = setupJECCorrector(onGrid,isdata,'AK8PFchs',runera)
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
		self.__renormdict = renormdict
		self.__alphalist = getList(renormdict,'alpha')
		self.__epsilonlist = getList(renormdict,'epsilon')

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
		alphas_up_unc 	= abs(branches['alphas_Weights'].getReadValue(1)*0.75-1.) if len(branches['alphas_Weights'].getReadArray())>1 else 0.
		alphas_down_unc = abs(branches['alphas_Weights'].getReadValue(0)*0.75-1.) if len(branches['alphas_Weights'].getReadArray())>1 else 0.
		#nominal value
		returnlist.append(pdfmean/self.__renormdict['pdfas'])
		#up/down with pdf and alpha_s added in quadrature
		returnlist.append((pdfmean+sqrt(pdfunc**2+alphas_up_unc**2))/self.__renormdict['pdfasup'])
		returnlist.append(pdfmean-sqrt(pdfunc**2+alphas_down_unc**2)/self.__renormdict['pdfasdown'])
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
    etamin = 	   [0.0,  0.5,  0.8,  1.1,  1.3,  1.7,  1.9,  2.1,  2.3,  2.5,  2.8,  3.0,  3.2]
    etamax = 	   [0.5,  0.8,  1.1,  1.3,  1.7,  1.9,  2.1,  2.3,  2.5,  2.8,  3.0,  3.2,  5.0]
    scale_nom =    [1.109,1.138,1.114,1.123,1.084,1.082,1.140,1.067,1.177,1.364,1.857,1.328,1.16]
    scale_uncert = [0.008,0.013,0.013,0.024,0.011,0.035,0.047,0.053,0.041,0.039,0.071,0.022,0.029]
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
	data_nom_file  = TFile.Open(pp+DATA_PU_HISTO_NOMINAL_FILENAME)
	data_nom_histo  = data_nom_file.Get('pileup')
	data_nom_histo.Scale(1./data_nom_histo.Integral())
	data_up_file   = TFile.Open(pp+DATA_PU_HISTO_XS_UP_FILENAME)
	data_up_histo   = data_up_file.Get('pileup')
	data_up_histo.Scale(1./data_up_histo.Integral())
	data_down_file = TFile.Open(pp+DATA_PU_HISTO_XS_DOWN_FILENAME)
	data_down_histo = data_down_file.Get('pileup')
	data_down_histo.Scale(1./data_down_histo.Integral())
	pu_histo.SetDirectory(0); data_nom_histo.SetDirectory(0); data_up_histo.SetDirectory(0); data_down_histo.SetDirectory(0)
	data_nom_file.Close(); data_up_file.Close(); data_down_file.Close()
	return pu_histo, data_nom_histo, data_up_histo, data_down_histo

def setupMuonIDHistos(onGrid) :
	pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
	muon_id_file = TFile.Open(pp+MUON_ID_EFF_ROOT_FILENAME)
	vs_eta_histo = muon_id_file.Get(MUON_ID_VS_ETA_HISTONAME)
	vs_pt_histo  = muon_id_file.Get(MUON_ID_VS_PT_HISTONAME)
	vs_pu_histo  = muon_id_file.Get(MUON_ID_VS_PU_HISTONAME)
	vs_eta_histo.SetDirectory(0); vs_pt_histo.SetDirectory(0); vs_pu_histo.SetDirectory(0)
	muon_id_file.Close()
	return vs_eta_histo, vs_pt_histo, vs_pu_histo

def setupEleTrkHisto(onGrid) :
	pp = './tardir/' if onGrid == 'yes' else '../other_input_files/'
	ele_id_file = TFile.Open(pp+ELE_TRK_EFF_ROOT_FILENAME)
	histo = ele_id_file.Get(ELE_TRK_EFF_HISTO_NAME)
	histo.SetDirectory(0)
	ele_id_file.Close()
	return histo

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