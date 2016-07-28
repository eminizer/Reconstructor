#Global variables

#imports
from ROOT import TLorentzVector

class Jet(object) :

	def __init__(self,branches,index,jec,pp) :
		self.__fourvec = getfourvec(branches,index,jec,pp)
		self.__pt = self.__fourvec.Pt()

	def getFourVector(self) :
		return self.__fourvec
	def getPt(self) :
		return self.__pt

class AK4Jet(Jet) :

	def __init__(self,branches,index,jec) :
		Jet.__init__(self,branches,index,jec,'jetAK4')

class AK8Jet(Jet) :

	def __init__(self,branches,index,jec) :
		Jet.__init__(self,branches,index,jec,'jetAK8')
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

def getfourvec(branches,index,jec,pp) :
	pt  = branches[pp+'_Pt'].getReadValue(index)
	eta = branches[pp+'_Eta'].getReadValue(index)
	phi = branches[pp+'_Phi'].getReadValue(index)
	E   = branches[pp+'_E'].getReadValue(index)
	nomvec = TLorentzVector(); nomvec.SetPtEtaPhiE(pt,eta,phi,E)
	if jec == 'nominal' :
		return nomvec
	smpt  = branches[pp+'_SmearedPt'].getReadValue(index)
	smeta = branches[pp+'_SmearedPEta'].getReadValue(index)
	smphi = branches[pp+'_SmearedPhi'].getReadValue(index)
	smE   = branches[pp+'_SmearedE'].getReadValue(index)
	pterr  = (smpt-pt)/pt
	etaerr = (smeta-eta)/eta
	phierr = (smphi-phi)/phi
	Eerr   = (smE-E)/E
	jecerr = branches[pp+'_jecUncertainty'].getReadValue(index)
	totalerr = sqrt(jecerr**2+pterr**2+etaerr**2+phierr**2+Eerr**2)
	if jec == 'up' :
		return (1.+totalerr)*nomvec
	elif jec == 'down' :
		return (1.-totalerr)*nomvec
	else :
		print 'WARNING: JEC OPTION '+str(jec)+' NOT RECOGNIZED'
		return nomvec
