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
		self.__miniIso = None
		self.__dR = None; self.__relPt = None; self.__nearestJetPt = None

	def calculateIsolation(self,jets) :
		nearestJet = findNearestJet(self.__fourvec,jets)
		nearestJetVec = nearestJet.getFourVector()
		self.__dR = self.__fourvec.DeltaR(nearestJetVec)
		self.__relPt = self.__fourvec.Pt(nearestJetVec.Vect())
		self.__nearestJetPt = nearestJetVec.Pt()
		cleanedleplist = nearestJet.getListOfCleanedLeptons()
		self.__wasCleanedFromNearestJet = 1 if cleanedleplist!=None and self in cleanedleplist else 0

	def isTightMiniIso(self) :
		return self.__miniIso<0.1
	def isMedMiniIso(self) :
		return self.__miniIso<0.2
	def isLooseMiniIso(self) :
		return self.__miniIso<0.4

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
	def setIsValid(self,iv) :
		self.__isValid = iv
	def isValid(self) :
		return self.__isValid
	def setID(self,thisid) :
		self.__ID = thisid
	def getID(self) :
		return self.__ID
	def setIso(self,thisiso) :
		self.__iso = thisiso
	def getIso(self) :
		return self.__iso
	def setMiniIso(self,miniiso) :
		self.__miniIso = miniiso
	def getMiniIso(self) :
		return self.__miniIso

class Muon(Lepton) :

	def __init__(self,branches,index,runera) :
		Lepton.__init__(self,branches,index,'mu')
		if runera=='B' or runera=='C' or runera=='D' or runera=='E' or runera=='F' :
			self.setID(branches['mu_IsMediumMuon2016'].getReadValue(index))
		else :
			self.setID(branches['mu_IsMediumMuon'].getReadValue(index))
		self.setIso(branches['mu_Iso04'].getReadValue(index))
		self.setMiniIso(branches['mu_MiniIso'].getReadValue(index))
		self.setIsValid(self.getPt()>55. and abs(self.getEta())<2.5 and self.getID()==1)

	def isLooseIso(self) :
		return self.getIso<0.25 #loose WP https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
	def isTightIso(self) :
		return self.getIso<0.15 #tight WP, already divided by pt https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
	def isIso(self) : #really just to match with electron functions
		return self.getIso<0.15 #tight WP, already divided by pt https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
	def is2DIso(self,eventtopology) :
		if eventtopology<3 :
			return self.getDR()>0.4 or self.getRelPt()>30.
		elif eventtopology==3 :
			return self.getDR()>0.4 or self.getRelPt()>30.

class Electron(Lepton) :

	def __init__(self,branches,index) :
		Lepton.__init__(self,branches,index,'el')
		self.setID(branches['el_IDMedium_NoIso'].getReadValue(index))
		self.__tightID = branches['el_IDTight_NoIso'].getReadValue(index)
		self.__scEta = branches['el_SCEta'].getReadValue(index)
		self.setIso(branches['el_Iso03'].getReadValue(index))
		self.setMiniIso(branches['el_MiniIso'].getReadValue(index))
		self.setIsValid(self.getPt()>55. and abs(self.__scEta)<2.5 and self.getID()==1)

	def getEtaSC(self) :
		return self.__scEta
	def isLooseIso(self) :
		return self.getIso<0.0695 #what was removed from the Medium ID https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	def isIso(self) :
		return self.getIso<0.0695 #what was removed from the Medium ID https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
	def is2DIso(self,eventtopology) :
		if eventtopology<3 :
			return self.getDR()>0.4 or self.getRelPt()>30.
		elif eventtopology==3 :
			return self.getDR()>0.4 or self.getRelPt()>20.
	def getTightID(self) :
		return self.__tightID

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
