#Global variables

#imports
from ROOT import TLorentzVector
from math import *

class Jet(object) :

	def __init__(self,branches,index,jes,jer,lep,corrector,isdata,pp) :
		self.__fourvec = getfourvec(branches,index,jes,jer,lep,corrector,isdata,pp)
		if self.__fourvec!=None :
			self.__pt = self.__fourvec.Pt()

	def getFourVector(self) :
		return self.__fourvec
	def getPt(self) :
		return self.__pt

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

	def getTau32(self) :
		return self.__tau3/self.__tau2
	def getTau21(self) :
		return self.__tau2/self.__tau1
	def getSDM(self) :
		return self.__sdm

def getfourvec(branches,index,jes,jer,lep,corrector,isdata,pp) :
	#get all the uncertainties
	pt  = branches[pp+'_Pt'].getReadValue(index)
	eta = branches[pp+'_Eta'].getReadValue(index)
	phi = branches[pp+'_Phi'].getReadValue(index)
	E   = branches[pp+'_E'].getReadValue(index)
	rawjet = TLorentzVector(); rawjet.SetPtEtaPhiE(pt,eta,phi,E)
	jec0 = branches[pp+'_jecFactor0'].getReadValue(index)
	#roll back the jec to get the raw jet
	rawjet*=(1./jec0)
	cleanjet = cleanJet(rawjet,lep.getFourVector())
	if cleanjet == None :
	#	print 'CLEANED JET WAS "NONE"' #DEBUG
		return None
	#get and apply the new jec
	jetArea = branches[pp+'_jetArea'].getReadValue(index)
	rho = branches['rho'].getReadValue()
	npv = branches['npv'].getReadValue()
	newJEC = corrector.getJECforJet(cleanjet,jetArea,rho,npv,pp)
	#If this is data, don't apply any smearing or systematics, just return the corrected, cleaned jet
	nominalJet = cleanjet*newJEC
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
		newJECuncDown, newJECuncUp = corrector.getJECuncForJet(cleanjet,pp)
		#And scale the fourvector up or down if we're looking for JES corrections
		if jes=='up' :
			jecUpJet = cleanjet*(newJEC+newJECuncUp)
			newJet=corrector.smearJet(jecUpJet,jer,genPt,genEta,genPhi)
		elif jes=='down' :
			jecDownJet = cleanjet*(newJEC+newJECuncDown)
			newJet=corrector.smearJet(jecDownJet,jer,genPt,genEta,genPhi)
	return newJet

def cleanJet(jetvec,lepvec) :
	#if the lepton is within 0.4 of the jet, remove it
	if jetvec.DeltaR(lepvec) < 0.4 :
		if lepvec.E() > jetvec.E() :
	#		print 'THIS JET WAS BASICALLY JUST A LEPTON!!' #DEBUG
			return None
	#	print 'REMOVING LEPTON FROM JET' #DEBUG
	#	print 'Jet Before = (pT, eta, phi, M) = (%.2f, %.2f, %.2f, %.2f)'%(jetvec.Pt(),jetvec.Eta(),jetvec.Phi(),jetvec.M()) #DEBUG
		jetvec-=lepvec
	#	print 'Jet After = (pT, eta, phi, M) = (%.2f, %.2f, %.2f, %.2f)'%(jetvec.Pt(),jetvec.Eta(),jetvec.Phi(),jetvec.M()) #DEBUG
	#print 'NO LEPTON CLEANING NEEDED' #DEBUG
	if jetvec.Pt()==0. :
	#	print 'RESULTING JET HAS NO PT' #DEBUG
		return None
	return jetvec
