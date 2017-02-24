#Global variables
#PDG ID Numbering scheme http://pdg.lbl.gov/2002/montecarlorpp.pdf
PROTON_ID = 2212
B_ID 	  = 5
TOP_ID    = 6
GLUON_ID  = 21
W_ID      = 24
ELECTRON_ID = 11
TAU_NEUTRINO_ID = 18
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0

#Imports
from ROOT import TLorentzVector
from math import *

#################################  REFERENCED FUNCTIONS  ##################################

#given the MC gen info branches and the event type, returns a tuple of (keep [bool], addtwice [0 or 1])
def keepEventType(branches,event_type) :
	#event types:
	#	0 = semileptonic qqbar
	#	1 = semileptonic gg, etc.
	#	2 = dileptonic
	#	3 = hadronic
	#	4 = everything else
	#no background is a wrong background : )
	if event_type==4 :
		return True,False #yes keep it, no don't add it twice
	else :
		#return typeCheckMCAtNLO(branches,event_type)
	#	print '-------------------------------------' #DEBUG
		addTwice = False
		#check the initial state partons if necessary
		if event_type<2 :
			ip_ids_t = []
			ip_ids_tbar = []
			#loop through and find the particles whose daughters include the ttbar pair
			for i in range(branches['gen_size'].getReadValue()) :
				status = branches['gen_Status'].getReadValue(i)
				thisID = branches['gen_ID'].getReadValue(i)
				Mom0ID = branches['gen_Mom0ID'].getReadValue(i)
				Mom1ID = branches['gen_Mom1ID'].getReadValue(i)
				Dau0ID = branches['gen_Dau0ID'].getReadValue(i)
				#find the tops that have two mothers and decayed only into themselves
				if abs(thisID)==TOP_ID and Mom0ID!=-900 and Mom1ID!=-900 and Dau0ID==thisID :
	#				print '[%d + %d] = %d -> (%d + %d)'%(Mom0ID,Mom1ID,thisID,Dau0ID,branches['gen_Dau1ID'].getReadValue(i)) #DEBUG
					if thisID>0 :
						ip_ids_t.append(Mom0ID); ip_ids_t.append(Mom1ID)
					else :
						ip_ids_tbar.append(Mom0ID); ip_ids_tbar.append(Mom1ID)
			#is it a qqbar event?
			is_qq = len(ip_ids_t)==2 and len(ip_ids_tbar)==2 and ip_ids_t[0]+ip_ids_t[1]==0 and ip_ids_tbar[0]+ip_ids_tbar[1]==0 
	#		print 'ISQQ = '+str(is_qq) #DEBUG
			#regardless, did it have an initially symmetric state
			addTwice = event_type==0 or (len(ip_ids_t)==2 and len(ip_ids_tbar)==2 and ip_ids_t[0]==GLUON_ID and ip_ids_t[1]==GLUON_ID and ip_ids_tbar[0]==GLUON_ID and ip_ids_tbar[1]==GLUON_ID)
	#		print 'ADDTWICE = '+str(addTwice) #DEBUG
			#return false if the event type is incorrect
			if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
				return False,addTwice
		return semilepCheck(branches,event_type),addTwice

#finds which "down-the-beampipe" fourvector corresponds to the initial quark
def findInitialPartons(branches) :
	#print '------------------------------------------------' #DEBUG
	#first algorithm for simple first-order processes (also goes through once so we use it to assign particles for the higher-order algorithm)
	inip1id = 0; inip2id = 0
	for i in range(branches['gen_size'].getReadValue()) :
		status = branches['gen_Status'].getReadValue(i)
		thisID = branches['gen_ID'].getReadValue(i)
		Mom0ID = branches['gen_Mom0ID'].getReadValue(i)
		Mom1ID = branches['gen_Mom1ID'].getReadValue(i)
		Dau0ID = branches['gen_Dau0ID'].getReadValue(i)
		Dau1ID = branches['gen_Dau1ID'].getReadValue(i)
	#	print '[%d + %d] = %d -> (%d + %d)'%(Mom0ID,Mom1ID,thisID,Dau0ID,branches['gen_Dau1ID'].getReadValue(i)) #DEBUG
		if abs(Mom0ID)==PROTON_ID and Mom1ID==-900. and abs(Dau0ID)==TOP_ID and abs(Dau1ID)==TOP_ID :
	#		print '	eta = %.4f'%(branches['gen_Eta'].getReadValue(i)) #DEBUG
			factor = 1.0 if branches['gen_Eta'].getReadValue(i) > 0 else -1.0
			factor*=abs(thisID)/thisID
			return TLorentzVector(1.0,0.0,factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY), TLorentzVector(1.0,0.0,-1.*factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
		elif abs(Mom0ID)!=TOP_ID and abs(Mom1ID)!=TOP_ID and abs(thisID)==TOP_ID and Dau0ID==thisID and Dau1ID==-900 and inip1id==0 and inip2id==0 : 
			inip1id = Mom0ID; inip2id=Mom1ID
	#second algorithm for higher order events
	if inip1id!=0 and inip2id!=0 :
	#	print 'USING HIGHER ORDER ALGORITHM' #DEBUG
		for i in range(branches['gen_size'].getReadValue()) :
			status = branches['gen_Status'].getReadValue(i)
			thisID = branches['gen_ID'].getReadValue(i)
			Mom0ID = branches['gen_Mom0ID'].getReadValue(i)
			Mom1ID = branches['gen_Mom1ID'].getReadValue(i)
			Dau0ID = branches['gen_Dau0ID'].getReadValue(i)
			Dau1ID = branches['gen_Dau1ID'].getReadValue(i)
			if  abs(Mom0ID)==PROTON_ID and Mom1ID==-900. and (thisID==inip1id or thisID==inip2id) and Dau0ID==thisID and Dau1ID==-900 :
				factor = 1.0 if branches['gen_Eta'].getReadValue(i) > 0 else -1.0
				factor*=abs(thisID)/thisID
				return TLorentzVector(1.0,0.0,factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY), TLorentzVector(1.0,0.0,-1.*factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
	#and just in case neither of those worked, here's this default!
	return TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY), TLorentzVector(1.0,0.0,-1.*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)

#returns a bunch of fourvectors from Monte Carlo for MC truth reconstruction and matching
def findMCParticles(branches) :
	tvec = TLorentzVector()
	tbarvec = TLorentzVector()
	lepvec = TLorentzVector()
	vvec = TLorentzVector()
	hadWvec = TLorentzVector()
	bs = []
	lepcharge = 0
	#print '---------------------------------------------------------------------' #DEBUG
	#loop over the particles and assign them
	for i in range(branches['gen_size'].getReadValue()) :
		status = branches['gen_Status'].getReadValue(i)
		thisID = branches['gen_ID'].getReadValue(i)
		Mom0ID = branches['gen_Mom0ID'].getReadValue(i)
		Mom1ID = branches['gen_Mom1ID'].getReadValue(i)
		Dau0ID = branches['gen_Dau0ID'].getReadValue(i)
		Dau1ID = branches['gen_Dau1ID'].getReadValue(i)
		dau_IDs = (abs(Dau0ID), abs(Dau1ID))
	#	s = '[%d + %d] = %d -> (%d + %d) '%(Mom0ID,Mom1ID,thisID,Dau0ID,branches['gen_Dau1ID'].getReadValue(i)) #DEBUG
		#top
		if thisID == TOP_ID and Mom0ID==thisID and B_ID in dau_IDs and W_ID in dau_IDs :
	#		s+='top' #DEBUG
			tvec.SetPtEtaPhiE(branches['gen_Pt'].getReadValue(i),branches['gen_Eta'].getReadValue(i),branches['gen_Phi'].getReadValue(i),branches['gen_E'].getReadValue(i))
		#antitop
		elif thisID == -1*TOP_ID and Mom0ID==thisID and B_ID in dau_IDs and W_ID in dau_IDs :
	#		s+='antitop' #DEBUG
			tbarvec.SetPtEtaPhiE(branches['gen_Pt'].getReadValue(i),branches['gen_Eta'].getReadValue(i),branches['gen_Phi'].getReadValue(i),branches['gen_E'].getReadValue(i))
		#lepton and neutrino
		elif abs(Mom0ID)==W_ID and abs(thisID) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) :
			if abs(thisID)%2==1 : #electron, muon, or tau
	#			s+='lepton' #DEBUG
				lepvec.SetPtEtaPhiE(branches['gen_Pt'].getReadValue(i),branches['gen_Eta'].getReadValue(i),branches['gen_Phi'].getReadValue(i),branches['gen_E'].getReadValue(i))
				lepcharge = 1 if thisID>0 else -1
			else : #neutrino
	#			s+='neutrino' #DEBUG
				vvec.SetPtEtaPhiE(branches['gen_Pt'].getReadValue(i),branches['gen_Eta'].getReadValue(i),branches['gen_Phi'].getReadValue(i),branches['gen_E'].getReadValue(i))
		#hadronic W
		elif abs(Mom0ID)==W_ID and thisID==Mom0ID and Dau0ID not in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) and Dau1ID not in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) :
	#		s+='hadronic W' #DEBUG
			hadWvec.SetPtEtaPhiE(branches['gen_Pt'].getReadValue(i),branches['gen_Eta'].getReadValue(i),branches['gen_Phi'].getReadValue(i),branches['gen_E'].getReadValue(i))
		#bs from tops
		elif abs(Mom0ID)==TOP_ID and abs(thisID)==B_ID :
	#		s+='decaying b' #DEBUG
			newb = TLorentzVector(); newb.SetPtEtaPhiE(branches['gen_Pt'].getReadValue(i),branches['gen_Eta'].getReadValue(i),branches['gen_Phi'].getReadValue(i),branches['gen_E'].getReadValue(i))
			bs.append((newb,thisID))
	#	print s #DEBUG
	#find out which was the hadronic/leptonic b
	hadb = [TLorentzVector(1.,0.,0.,1.),0]; lepb = [TLorentzVector(1.,0.,0.,1.),0]
	if len(bs)==2 and bs[0][1]==-1*bs[1][1] :
		hadb = bs[0] if bs[0][1]/lepcharge>0 else bs[1]
		lepb = bs[0] if bs[0][1]/lepcharge<0 else bs[1]
	#print 'hadronic b ID = %d, leptonic b ID = %d'%(hadb[1],lepb[1]) #DEBUG
	return (tvec,tbarvec,lepvec,lepcharge,vvec,lepb[0],hadWvec,hadb[0])

###################################  HELPER FUNCTIONS  ##################################

#semileptonic checker
def semilepCheck(branches,event_type) :
	#check the decay type
	leptons_from_Ws_IDs = []
	for i in range(branches['gen_size'].getReadValue()) :
		status = branches['gen_Status'].getReadValue(i)
		thisID = branches['gen_ID'].getReadValue(i)
		Mom0ID = branches['gen_Mom0ID'].getReadValue(i)
		Mom1ID = branches['gen_Mom1ID'].getReadValue(i)
		Dau0ID = branches['gen_Dau0ID'].getReadValue(i)
		Dau1ID = branches['gen_Dau1ID'].getReadValue(i)
		#look for the decaying Ws that have themselves (or their top) as their only mother and daughters that are not themselves
		if abs(thisID)==W_ID and (Mom0ID==thisID or abs(Mom0ID)==TOP_ID) and Mom1ID==-900 and Dau0ID!=thisID and Dau1ID!=-900 :
	#		print '[%d + %d] = %d -> (%d + %d)'%(Mom0ID,Mom1ID,thisID,Dau0ID,Dau1ID) #DEBUG
			if abs(Dau0ID) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) : leptons_from_Ws_IDs.append(Dau0ID)
			if abs(Dau1ID) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) : leptons_from_Ws_IDs.append(Dau1ID)
	if ( (len(leptons_from_Ws_IDs) == 2 and event_type < 2) or (len(leptons_from_Ws_IDs) == 4 and event_type == 2) or 
		(len(leptons_from_Ws_IDs) == 0 and event_type == 3) ) :
	#	print 'LEP RIGHT' #DEBUG
		return True
	#print 'LEP WRONG' #DEBUG
	return False
