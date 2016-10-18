#Global variables

#Imports
from ROOT import vector, FactorizedJetCorrector, JetCorrectorParameters, JetCorrectionUncertainty
from random import gauss
from math import sqrt

class Corrector(object) :

	##################################  #__init__ function  ##################################
	def __init__(self,isdata,generator,eventType,onGrid) :
		self.__ak4JetCorrector, self.__ak4JecUncertainty = setupJECCorrector(onGrid,isdata,'AK4PFchs')
		self.__ak8JetCorrector, self.__ak8JecUncertainty = setupJECCorrector(onGrid,isdata,'AK8PFchs')

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

	def __del__(self) :
		pass

def setupJECCorrector(onGrid,isdata,jetType) :
	#Define the JEC  parameters
	L1JetPar  = None
	L2JetPar  = None
	L3JetPar  = None
	pp = ''
	if onGrid == 'yes' :
		pp+='./tardir/'
	else :
		pp+='../other_input_files/'
	if not isdata :
		L1JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_MC_L1FastJet_'+jetType+'.txt')
		L2JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_MC_L2Relative_'+jetType+'.txt')
		L3JetPar  = JetCorrectorParameters(pp+'Spring16_25nsV6_MC_L3Absolute_'+jetType+'.txt')
		jetUncertainty = JetCorrectionUncertainty(pp+'Spring16_25nsV6_MC_Uncertainty_'+jetType+'.txt')
	else :
		print 'I do not know where to get the corrections for data right now; this will crash.'
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
