#Imports
from ROOT import TLorentzVector

class Lepton(object) :

	def __init__(self,branches,index,pp) :
		self.__pt = branches[pp+'_Pt'].getReadValue(index)
		self.__eta = branches[pp+'_Eta'].getReadValue(index)
		phi = branches[pp+'_Phi'].getReadValue(index)
		E = branches[pp+'_E'].getReadValue(index)
		self.__fourvec = TLorentzVector(); self.__fourvec.SetPtEtaPhiE(self.__pt,self.__eta,phi,E)
		self.__Q = branches[pp+'_Charge'].getReadValue(index)
		self.__Key = branches[pp+'_Key'].getReadValue(index)

	def calculateIsolation(self,jets) :
		nearestJetVec = findNearestJet(self.__fourvec,jets).getFourVector()
		self.__dR = self.__fourvec.DeltaR(nearestJetVec)
		self.__relPt = self.__fourvec.Pt(nearestJetVec.Vect())

	def getPt(self) :
		return self.__pt
	def getEta(self) :
		return self.__eta
	def getFourVector(self) :
		return self.__fourvec
	def getQ(self) :
		return self.__Q
	def getKey(self) :
		return self.__Key
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
		self.__ID = branches['el_IDTight_NoIso'].getReadValue(index)
		self.__scEta = branches['el_SCEta'].getReadValue(index)

	def getID(self) :
		return self.__ID
	def getEtaSC(self) :
		return self.__scEta

def findNearestJet(lepvec,jets) :
	closestJet = jets[0]
	closestDR = closestJet.getFourVector().DeltaR(lepvec)
	for i in range(1,len(jets)) :
		jet = jets[i]
		if jet.getPt()<15. or abs(jet.getEta())>3. :
			break
		checkDR = jet.getFourVector().DeltaR(lepvec)
		if checkDR < closestDR :
			closestDR = checkDR
			closestJet = jet
	return closestJet
