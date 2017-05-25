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
		self.__type = pp

	def calculateIsolation(self,jets) :
		nearestJet = findNearestJet(self.__fourvec,jets)
		nearestJetVec = nearestJet.getFourVector()
		self.__dR = self.__fourvec.DeltaR(nearestJetVec)
		self.__relPt = self.__fourvec.Pt(nearestJetVec.Vect())
		self.__nearestJetPt = nearestJetVec.Pt()
		cleanedleplist = nearestJet.getListOfCleanedLeptons()
		self.__wasCleanedFromNearestJet = 1 if cleanedleplist!=None and self in cleanedleplist else 0

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
	def getNearestJetPt(self) :
		return self.__nearestJetPt
	def getType(self) :
		return self.__type
	def wasCleanedFromNearestJet(self) :
		return self.__wasCleanedFromNearestJet

class Muon(Lepton) :

	def __init__(self,branches,index,runera) :
		Lepton.__init__(self,branches,index,'mu')
		if runera=='B' or runera=='C' or runera=='D' or runera=='E' or runera=='F' :
			self.__ID = branches['mu_IsMediumMuon2016'].getReadValue(index)
		else :
			self.__ID = branches['mu_IsMediumMuon'].getReadValue(index)
		self.__iso = branches['mu_Iso04'].getReadValue(index)
		self.__miniIso = branches['mu_MiniIso'].getReadValue(index)
		self.__isValid = self.getPt()>55. and abs(self.getEta())<2.5 and self.__ID==1 #my selection
		#self.__isValid = self.getPt()>50. and abs(self.getEta())<2.1 and self.__ID==1 #Susan's selection

	def isValid(self) :
		return self.__isValid
	def getID(self) :
		return self.__ID
	def isIso(self) :
		return (self.__iso/Lepton.getPt(self))<0.1
	def isMiniIso(self) :
		return (self.__miniIso/Lepton.getPt(self))<0.2

class Electron(Lepton) :

	def __init__(self,branches,index) :
		Lepton.__init__(self,branches,index,'el')
		self.__ID = branches['el_IDMedium_NoIso'].getReadValue(index)
		self.__scEta = branches['el_SCEta'].getReadValue(index)
		self.__iso = branches['el_Iso03'].getReadValue(index)
		self.__miniIso = branches['el_MiniIso'].getReadValue(index)
		self.__isValid = self.getPt()>50. and abs(self.__scEta)<2.5 and self.__ID==1 #my selection
		#self.__isValid = self.getPt()>50. and abs(self.getEta())<2.5 and self.__ID==1 #Susan's selection

	def isValid(self) :
		return self.__isValid
	def getID(self) :
		return self.__ID
	def getEtaSC(self) :
		return self.__scEta
	def isIso(self) :
		return (self.__iso/Lepton.getPt(self))<0.12
	def isMiniIso(self) :
		return (self.__miniIso/Lepton.getPt(self))<0.2

def findNearestJet(lepvec,jets) :
	closestJet = jets[0]
	closestDR = lepvec.DeltaR(closestJet.getFourVector())
	for i in range(1,len(jets)) :
		jet = jets[i]
		if jet.getPt()<15. or abs(jet.getEta())>3. :
			break
		checkDR = lepvec.DeltaR(jet.getFourVector())
		if checkDR < closestDR :
			closestDR = checkDR
			closestJet = jet
	return closestJet
