#Imports
from ROOT import TLorentzVector

class Lepton(object) :

	def __init__(self,branches,index,pp) :
		self.__pt = branches[pp+'_Pt'].getReadValue(index)
		eta = branches[pp+'_Eta'].getReadValue(index)
		phi = branches[pp+'_Phi'].getReadValue(index)
		E = branches[pp+'_E'].getReadValue(index)
		self.__fourvec = TLorentzVector(); self.__fourvec.SetPtEtaPhiE(self.__pt,eta,phi,E)
		self.__Q = branches[pp+'_Charge'].getReadValue(index)

	def calculateIsolation(self,jets) :
		nearestJet = findNearestJet(self.__fourvec,jets)
		nearestJetvec = nearestJet.getFourVector()
		self.__relPt = nearestJetvec.Pt(self.__fourvec.Vect())
		self.__dR = nearestJetvec.DeltaR(self.__fourvec)

	def getPt(self) :
		return self.__pt
	def getFourVector(self) :
		return self.__fourvec
	def getQ(self) :
		return self.__Q
	def getRelPt(self) :
		return self.__relPt
	def getDR(self) :
		return self.__dR

class Muon(Lepton) :

	def __init__(self,branches,index) :
		Lepton.__init__(self,branches,index,'mu')
		self.__ID = branches['mu_IsLooseMuon'].getReadValue(index)

	def getID(self) :
		return self.__ID

class Electron(Lepton) :

	def __init__(self,branches,index) :
		Lepton.__init__(self,branches,index,'el')
		self.__ID = branches['el_vidLoose'].getReadValue(index)

	def getID(self) :
		return self.__ID

def findNearestJet(lepvec,jets) :
	closestDR = jets[0].getFourVector().DeltaR(lepvec)
	closestJet = jets[0]
	for i in range(1,len(jets)) :
		jet = jets[i]
		if jet.getPt()<25 :
			break
		checkDR = jet.getFourVector().DeltaR(lepvec)
		if checkDR < closestDR :
			closestDR = checkDR
			closestJet = jet
	return closestJet
