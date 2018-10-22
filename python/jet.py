#Global variables

#imports
from ROOT import TLorentzVector
from math import *

class Jet(object) :

	def __init__(self,branches,index,jec,leps,corrector,isdata,pp) :
		#print '----------------------- New Jet --------------------------' #DEBUG
		self.__fourvec, self.__cleanedLeptons, self.__metCorrVec = getfourvec(branches,index,jec,leps,corrector,isdata,pp)
		self.__pt = self.__fourvec.Pt() if self.__fourvec!=None else -900
		self.__eta = self.__fourvec.Eta() if self.__fourvec!=None else -900
		self.__isIDed = self.__checkID__(branches,index,pp)
		self.__csvv2 = branches[pp+'_CSVv2'].getReadValue(index)
		self.__isbtaggedL = self.__csvv2>0.5426 #loose working point https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
		self.__isbtaggedM = self.__csvv2>0.8484 #medium working point https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco

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
		#These are the "Loose" ID criteria from here: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
		abseta = abs(self.__eta)
		if abseta<2.7 :
			if neutralHadronFraction<0.99 and neutralEMFraction<0.99 and numConstituents>1 :
				if abseta<2.4 :
					if chargedHadronFraction>0. and chargedMultiplicity>0 and chargedEMFraction<0.99 :
						return True
				else :
					return True
		elif abseta<3.0 :
			if neutralEMFraction>0.1 and neutralHadronFraction<0.98 and neutralMultiplicity>2 :
				return True
		else :
			if neutralEMFraction<0.90 and neutralMultiplicity>10 :
				return True
		return False

	def getMETCorrectionVec(self) :
		#print '	met correction vector = (%.1f,%.1f,%.1f,%.1f)'%(self.__metCorrVec.Pt(),self.__metCorrVec.Eta(),self.__metCorrVec.Phi(),self.__metCorrVec.M()) #DEBUG
		return self.__metCorrVec
	def getFourVector(self) :
		return self.__fourvec
	def getPt(self) :
		return self.__pt
	def getEta(self) :
		return self.__eta
	def isIDed(self) :
		return self.__isIDed
	def setIsValid(self,iv) :
		self.__isValid = iv
	def isValid(self) :
		return self.__isValid
	def getCSVv2(self) :
		return self.__csvv2
	def isLbTagged(self) :
		return self.__isbtaggedL
	def isMbTagged(self) :
		return self.__isbtaggedM
	def getListOfCleanedLeptons(self) :
		return self.__cleanedLeptons

class AK4Jet(Jet) :

	def __init__(self,branches,index,jes,jer,leps,corrector,isdata) :
		Jet.__init__(self,branches,index,jes,jer,leps,corrector,isdata,'jetAK4CHS')
		self.__flavor = branches['jetAK4CHS_HadronFlavour'].getReadValue(index)
		self.setIsValid(self.getPt()>30. and abs(self.getEta())<2.4 and self.isIDed()==1)
		self.__isValidForIsoCalc = self.getPt()>15. and abs(self.getEta())<3.0 and self.isIDed()==1
		#print '	isvalid = %s, for iso calc = %s'%(self.__isValid,self.__isValidForIsoCalc) #DEBUG

	def getFlavor(self) :
		return self.__flavor
	def isValidForIsoCalc(self) :
		return self.__isValidForIsoCalc

class AK8Jet(Jet) :

	def __init__(self,branches,index,jes,jer,leps,corrector,isdata) :
		Jet.__init__(self,branches,index,jes,jer,leps,corrector,isdata,'jetAK8CHS')
		#print 'Adding AK8 jet with soft drop mass %.4f'%(self.__sdm) #DEBUG
		self.__subjets = self.__getsubjets__(branches,index,jes,jer,leps,corrector,isdata)
		self.__n_subjets = len(self.__subjets)
		self.setIsValid(self.getPt()>200. and abs(self.getEta())<2.4 and self.isIDed()==1 and self.__n_subjets>1)
		self.__tau3 = branches['jetAK8CHS_tau3CHS'].getReadValue(index)
		self.__tau2 = branches['jetAK8CHS_tau2CHS'].getReadValue(index)
		self.__tau1 = branches['jetAK8CHS_tau1CHS'].getReadValue(index)
		self.__sdm  = branches['jetAK8CHS_softDropMassCHS'].getReadValue(index)
		self.__isttagged = self.getPt()>400. and self.__sdm>105. and self.__sdm<220. and self.__tau3!=0. and self.__tau2!=0. and (self.__tau3/self.__tau2)<0.80
		self.__isWtagged = self.getPt()>200. and self.__sdm>65. and self.__sdm<105. and self.__tau2!=0. and self.__tau1!=0. and (self.__tau2/self.__tau1)<0.55

	def __getsubjets__(self,branches,index,jes,jer,leps,corrector,isdata) :
		#start by making a list of all the subjets for this jet
		allsubjets = []
		#get the subjet indices
		subjet0index = branches['jetAK8CHS_vSubjetIndex0'].getReadValue(index)
		subjet1index = branches['jetAK8CHS_vSubjetIndex1'].getReadValue(index)
		#add the subjets
		if subjet0index>-1 :
			newSubjet = Jet(branches,int(subjet0index),jes,jer,leps,corrector,isdata,'subjetAK8CHS')
			if newSubjet.getFourVector()!=None :
				allsubjets.append(newSubjet)
		if subjet1index>-1 :
			newSubjet = Jet(branches,int(subjet1index),jes,jer,leps,corrector,isdata,'subjetAK8CHS')
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
	def isWTagged(self) :
		return self.__isWtagged

def getfourvec(branches,index,jec,leps,corrector,isdata,pp) :
	#get all the jet fourvectors, etc.
	pt  = branches[pp+'_Pt'].getReadValue(index)
	eta = branches[pp+'_Eta'].getReadValue(index)
	phi = branches[pp+'_Phi'].getReadValue(index)
	E   = branches[pp+'_E'].getReadValue(index)
	jec0 = branches[pp+'_jecFactor0'].getReadValue(index)
	#jecunc = 0. if pp.find('subjet')!=-1 else branches[pp+'_jecUncertainty'].getReadValue(index)
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
	metcorrvec = TLorentzVector(); metcorrvec.SetPtEtaPhiE(pt,eta,phi,E)
	#print '		metcorrvec initial = (%.1f,%.1f,%.1f,%.1f)'%(metcorrvec.Pt(),metcorrvec.Eta(),metcorrvec.Phi(),metcorrvec.M()) #DEBUG
	if rawjet.M()==-900 :
		#print 'RAW JET HAD NONSENSE MASS' #DEBUG
		return None, None, TLorentzVector()
	#print 'new jet with initial pT = %.2f'%(pt) #DEBUG
	rawjet*=jec0
	#print '	raw pT = %.2f'%(rawjet.Pt()) #DEBUG
	jetradius = 0.8 if pp.find('jetAK8')!=-1 else 0.4
	cleanjet, subtractedleps = cleanJet(rawjet,jetkeys,subjet1keys,subjet2keys,leps,jetradius)
	for lep in subtractedleps :
		metcorrvec=metcorrvec-lep.getFourVector()
	#print '		metcorrvec after lep sub = (%.1f,%.1f,%.1f,%.1f)'%(metcorrvec.Pt(),metcorrvec.Eta(),metcorrvec.Phi(),metcorrvec.M()) #DEBUG
	if cleanjet == None :
		#print 'CLEANED JET WAS "NONE"' #DEBUG
		return None, None, metcorrvec
	#print '	cleaned pT = %.2f'%(cleanjet.Pt()) #DEBUG
	#if this is a subjet, regardless of the systematics, return it without correcting
	if pp.find('subjet')!=-1 :
		return cleanjet, None, (metcorrvec-cleanjet)
	newJEC = 1./jec0
	#if there was a lepton subtracted from the jet get the on-the-fly JEC and recorrect it
	if len(subtractedleps)>0 :
		#get and apply the new jec
		jetArea = branches[pp+'_jetArea'].getReadValue(index)
		rho = branches['rho'].getReadValue()
		npv = branches['npv'].getReadValue()
		newJEC = corrector.getJECforJet(cleanjet,jetArea,rho,npv,pp)
	#adjust the new JEC if this is a systematic-shifted sample
	if jec!='nominal' and ((pp.find('jetAK4')!=-1 and jec.find('AK4JES')!=-1) or (pp.find('jetAK8')!=-1 and jec.find('AK8JES')!=-1)) :
		jecunc_up, jecunc_down = corrector.getJECuncForJet(cleanjet,pp,jec)
		if jec.endswith('_up') :
			newJEC = newJEC+jecunc_up
		elif jec.endswith('_dn') :
			newJEC = newJEC-jecunc_down
	adjJet = cleanjet*newJEC
	#print '	readjusted pT = %.2f'%(adjJet.Pt()) #DEBUG
	#If this is data, don't apply any smearing. just return the corrected, cleaned jet
	if isdata :
		return adjJet, subtractedleps, metcorrvec-adjJet
	#Otherwise smear the jet pT (need the generated pt, eta, phi, and pT resolution)
	genJetVec = TLorentzVector()
	genPt  = branches[pp+'_GenJetPt'].getReadValue(index)
	genEta = branches[pp+'_GenJetEta'].getReadValue(index)
	genPhi = branches[pp+'_GenJetPhi'].getReadValue(index)
	genE   = branches[pp+'_GenJetE'].getReadValue(index)
	if genPt==-900. or genEta==-900. or genPhi==-900. or genE==-900. :
		genJetVec = None
	else :
		genJetVec.SetPtEtaPhiE(genPt,genEta,genPhi,genE)
	dRCheck = 0.8 if pp.find('jetAK8')!=-1 else 0.4
	ptres = branches[pp+'_PtResolution'].getReadValue(index) 
	newJet = None
	#for JER-wiggled systematics
	if jec!='nominal' and ((pp.find('jetAK4')!=-1 and jec.find('AK4JER')!=-1) or (pp.find('jetAK8')!=-1 and jec.find('AK8JER')!=-1)) :
		newJet=corrector.smearJet(adjJet,jec,genJetVec,ptres,dRCheck)
	#for nominal smearing (not wiggling JER systematics)
	else :
		newJet=corrector.smearJet(nominalJet,'nominal',genJetVec,ptres,dRCheck)
	#print '	final pT = %.2f'%(newJet.Pt()) #DEBUG
	return newJet, subtractedleps, metcorrvec-newJet

def cleanJet(jetvec,jetkeys,sj1keys,sj2keys,leps,jetradius) :
	#keystring = '[' #DEBUG
	#for key in jetkeys : #DEBUG
	#	keystring+=str(key)+',' #DEBUG
	#keystring+= ']' #DEBUG
	#print 'cleaning jet with pT=%.1f and keys=%s'%(jetvec.Pt(),keystring) #DEBUG
	subtractedleps = []
	for lep in leps :
		lepvec=lep.getFourVector(); lepkey=lep.getKey(); leptype = lep.getType()
		#if the lepton is within the jet, remove it
	#	print '	jet/lep deltaR = %.4f'%(jetvec.DeltaR(lepvec)) #DEBUG
		if (leptype=='el' and lepvec.DeltaR(jetvec)<jetradius) or (leptype=='mu' and (lepkey in jetkeys or lepkey in sj1keys or lepkey in sj2keys)) :
	#		print '	found matching key %d'%(lepkey) #DEBUG
			if lepvec.E() > jetvec.E() :
	#			print '	THIS JET WAS BASICALLY JUST A LEPTON!!' #DEBUG
				return None, subtractedleps
	#		print '	REMOVING LEPTON FROM JET (%s with pT=%.1f)'%(lep.getType(),lep.getPt()) #DEBUG
	#		print '	Jet Before = (pT, eta, phi, M) = (%.2f, %.2f, %.2f, %.2f)'%(jetvec.Pt(),jetvec.Eta(),jetvec.Phi(),jetvec.M()) #DEBUG
			jetvec-=lepvec
			subtractedleps.append(lep)
	#		print '	Jet After = (pT, eta, phi, M) = (%.2f, %.2f, %.2f, %.2f)'%(jetvec.Pt(),jetvec.Eta(),jetvec.Phi(),jetvec.M()) #DEBUG
	#	else : #DEBUG
	#		print 'NO LEPTON CLEANING NEEDED' #DEBUG
		if jetvec.Pt()==0. :
	#		print '	RESULTING JET HAS NO PT' #DEBUG
			return None, subtractedleps
	return jetvec, subtractedleps
