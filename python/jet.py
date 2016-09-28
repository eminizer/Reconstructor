#Global variables

#imports
from ROOT import TLorentzVector
from math import *

class Jet(object) :

	def __init__(self,branches,index,jes,jer,lep,corrector,isdata,pp) :
		self.__fourvec = getfourvec(branches,index,jes,jer,lep,corrector,isdata,pp)
		self.__pt = self.__fourvec.Pt() if self.__fourvec!=None else None
		self.__eta = self.__fourvec.Eta() if self.__fourvec!=None else None
		self.__isIDed = self.__checkID__(branches,index,pp)
		self.__csvv2 = branches[pp+'_CSVv2'].getReadValue(index)

	def __checkID__(self,branches,index,pp) :
		if self.__eta==None :
			return None
		neutralHadronFraction = branches[pp+'_neutralHadronEnergyFrac'].getReadValue(index)
		neutralEMFraction = branches[pp+'_neutralEmEnergyFrac'].getReadValue(index)
		neutralMultiplicity = branches[pp+'_neutralMultiplicity'].getReadValue(index)
		chargedMultiplicity = branches[pp+'_chargedMultiplicity'].getReadValue(index)
		numConstituents = neutralMultiplicity+chargedMultiplicity
		chargedHadronFraction = branches[pp+'_chargedHadronEnergyFrac'].getReadValue(index)
		chargedEMFraction = branches[pp+'_chargedEmEnergyFrac'].getReadValue(index)
		if abs(self.__eta)<2.7 :
			if neutralHadronFraction<0.99 and neutralEMFraction<0.99 and numConstituents>1 :
				if abs(self.__eta)<2.4 :
					if chargedHadronFraction>0. and chargedMultiplicity>0 and chargedEMFraction<0.99 :
						return True
				else :
					return True
			else :
				return False
		elif abs(self.__eta)<3.0 :
			if neutralEMFraction<0.9 and neutralMultiplicity>2 :
				return True
			else :
				return False
		else :
			if neutralEMFraction<0.90 and neutralMultiplicity>10 :
				return True
			else :
				return False
		return False

	def getFourVector(self) :
		return self.__fourvec
	def getPt(self) :
		return self.__pt
	def getEta(self) :
		return self.__eta
	def isIDed(self) :
		return self.__isIDed
	def getCSVv2(self) :
		return self.__csvv2

class AK4Jet(Jet) :

	def __init__(self,branches,index,jes,jer,lep,corrector,isdata) :
		Jet.__init__(self,branches,index,jes,jer,lep,corrector,isdata,'jetAK4')

class AK8Jet(Jet) :

	def __init__(self,branches,index,jes,jer,lep,corrector,isdata) :
		Jet.__init__(self,branches,index,jes,jer,lep,corrector,isdata,'jetAK8')
		self.__tau3 = branches['jetAK8_tau3'].getReadValue(index)
		self.__tau2 = branches['jetAK8_tau2'].getReadValue(index)
		self.__tau1 = branches['jetAK8_tau1'].getReadValue(index)
		self.__sdm  = branches['jetAK8_softDropMass'].getReadValue(index)
		#print 'Adding AK8 jet with soft drop mass %.4f'%(self.__sdm) #DEBUG
		self.__subjets = self.__getsubjets__(branches,index,jes,jer,lep,corrector,isdata)
		self.__n_subjets = len(self.__subjets)
		self.__isttagged = self.__sdm>105. and self.__sdm<220. and self.__tau3!=0. and self.__tau2!=0. and (self.__tau3/self.__tau2)<0.81

	def __getsubjets__(self,branches,index,jes,jer,lep,corrector,isdata) :
		#start by making a list of all the subjets for this jet
		allsubjets = []
		#get the subjet indices
		subjet0index = branches['jetAK8_vSubjetIndex0'].getReadValue(index)
		subjet1index = branches['jetAK8_vSubjetIndex1'].getReadValue(index)
		#add the subjets
		if subjet0index>-1 :
			newSubjet = Jet(branches,int(subjet0index),jes,jer,lep,corrector,isdata,'subjetAK8')
			if newSubjet.getFourVector()!=None :
				allsubjets.append(newSubjet)
		if subjet1index>-1 :
			newSubjet = Jet(branches,int(subjet1index),jes,jer,lep,corrector,isdata,'subjetAK8')
			if newSubjet.getFourVector()!=None :
				allsubjets.append(newSubjet)
		#order the list of subjets by pt
		allsubjets.sort(key=lambda x: x.getPt(), reverse=True)
		#print '----------------------DONE ADDING SUBJETS----------------------' #DEBUG
		return allsubjets

	def getTau32(self) :
		if self.__tau2!=0 :
			return self.__tau3/self.__tau2
		else :
			return 9999.
	def getTau21(self) :
		if self.__tau1!=0 :
			return self.__tau2/self.__tau1
		else :
			return 9999.
	def getSDM(self) :
		return self.__sdm
	def getSubjet(self,index) :
		return self.__subjets[index]
	def getNSubjets(self) :
		return self.__n_subjets
	def isTopTagged(self) :
		return self.__isttagged

def getfourvec(branches,index,jes,jer,lep,corrector,isdata,pp) :
	#get all the jet fourvectors, etc.
	pt  = branches[pp+'_Pt'].getReadValue(index)
	eta = branches[pp+'_Eta'].getReadValue(index)
	phi = branches[pp+'_Phi'].getReadValue(index)
	E   = branches[pp+'_E'].getReadValue(index)
	jec0 = branches[pp+'_jecFactor0'].getReadValue(index)
	jecunc = 0. if pp.find('subjet')!=-1 else branches[pp+'_jecUncertainty'].getReadValue(index)
	#get the jet keys and subjet keys if necessary
	jetkeys = branches[pp+'_Keys'].getReadValue(index)
	subjet1keys = []; subjet2keys = []
	if pp.find('jetAK8')!=-1 and pp.find('subjet')==-1 :
		subjet0index = branches[pp+'_vSubjetIndex0'].getReadValue(index)
		subjet1index = branches[pp+'_vSubjetIndex1'].getReadValue(index)
		if subjet0index>-1 :
			subjet1keys = branches['sub'+pp+'_Keys'].getReadValue(int(subjet0index))
		if subjet1index>-1 :
			subjet2keys = branches['sub'+pp+'_Keys'].getReadValue(int(subjet1index))
	#roll back the jec to get the raw jet
	rawjet = TLorentzVector(); rawjet.SetPtEtaPhiE(pt,eta,phi,E); 
	if rawjet.M()==-900 :
		return None
	rawjet*=(1./jec0)
	cleanjet, subtracted = cleanJet(rawjet,jetkeys,subjet1keys,subjet2keys,lep.getFourVector(),lep.getKey())
	if cleanjet == None :
	#	print 'CLEANED JET WAS "NONE"' #DEBUG
		return None
	nominalJet = cleanjet
	#if this is a subjet, regardless of the systematics, return it without correcting
	if pp.find('subjet')!=-1 :
		return nominalJet
	newJEC = jec0
	#if there was a lepton subtracted from the jet
	if subtracted :
		#get and apply the new jec
		jetArea = branches[pp+'_jetArea'].getReadValue(index)
		rho = branches['rho'].getReadValue()
		npv = branches['npv'].getReadValue()
		newJEC = corrector.getJECforJet(cleanjet,jetArea,rho,npv,pp)
	nominalJet = cleanjet*newJEC
	#If this is data, don't apply any smearing or systematics, just return the corrected, cleaned jet
	if isdata :
		return nominalJet
	#The rest depends on whether we're doing JEC systematics
	#Also we need the generated pt, eta, phi
	genPt  = branches[pp+'_GenJetPt'].getReadValue(index)
	genEta = branches[pp+'_GenJetEta'].getReadValue(index)
	genPhi = branches[pp+'_GenJetPhi'].getReadValue(index)
	newJet = None
	if jes=='nominal' :
		newJet=corrector.smearJet(nominalJet,jer,genPt,genEta,genPhi)
	else :
		#otherwise get the new jec uncertainty
		newJECuncDown, newJECuncUp = jecunc, jecunc
		if subtracted :
			newJECuncDown, newJECuncUp = corrector.getJECuncForJet(cleanjet,pp)
		#And scale the fourvector up or down if we're looking for JES corrections
		if jes=='up' :
			jecUpJet = cleanjet*(newJEC+newJECuncUp)
			newJet=corrector.smearJet(jecUpJet,jer,genPt,genEta,genPhi)
		elif jes=='down' :
			jecDownJet = cleanjet*(newJEC+newJECuncDown)
			newJet=corrector.smearJet(jecDownJet,jer,genPt,genEta,genPhi)
	return newJet

def cleanJet(jetvec,jetkeys,sj1keys,sj2keys,lepvec,lepkey) :
	subtracted = False
	#if the lepton is within the jet, remove it
	#print 'deltaR = %.4f'%(jetvec.DeltaR(lepvec)) #DEBUG
	if lepkey in jetkeys or lepkey in sj1keys or lepkey in sj2keys :
	#	print 'found matching key %d'%(lepkey) #DEBUG
		if lepvec.E() > jetvec.E() :
	#		print 'THIS JET WAS BASICALLY JUST A LEPTON!!' #DEBUG
			return None, subtracted
	#	print 'REMOVING LEPTON FROM JET' #DEBUG
	#	print 'Jet Before = (pT, eta, phi, M) = (%.2f, %.2f, %.2f, %.2f)'%(jetvec.Pt(),jetvec.Eta(),jetvec.Phi(),jetvec.M()) #DEBUG
		jetvec-=lepvec
		subtracted = True
	#	print 'Jet After = (pT, eta, phi, M) = (%.2f, %.2f, %.2f, %.2f)'%(jetvec.Pt(),jetvec.Eta(),jetvec.Phi(),jetvec.M()) #DEBUG
	#else : #DEBUG
	#	print 'NO LEPTON CLEANING NEEDED' #DEBUG
	if jetvec.Pt()==0. :
	#	print 'RESULTING JET HAS NO PT' #DEBUG
		return None, subtracted
	return jetvec, subtracted
