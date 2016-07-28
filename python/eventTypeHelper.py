#Global variables
#PDG ID Numbering scheme http://pdg.lbl.gov/2002/montecarlorpp.pdf
PROTON_ID = 2212
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
def keepEventType(branches,event_type,generator) :
	#event types:
	#	0 = semileptonic qqbar
	#	1 = semileptonic gg, etc.
	#	2 = dileptonic
	#	3 = hadronic
	#	4 = everything else
	if generator == 'powheg' or generator == 'madgraph' or generator == 'mg5' or generator == 'mg' :
		return typeCheckPowheg(branches,event_type)
	elif generator == 'mcatnlo' :
		return typeCheckMCAtNLO(branches,event_type)
	elif generator == 'pythia8' :
		print 'PYTHIA NOT YET SUPPORTED'
		return None,None
	else :
		print 'ERROR: GENERATOR "'+generator+'"" NOT RECOGNIZED!!!'
		return None,None

#finds which "down-the-beampipe" fourvector corresponds to the initial quark
def findInitialPartons(branches,generator) :
	if generator == 'powheg' or generator == 'madgraph' or generator == 'mg5' or generator == 'mg' :
		return findInitialPartonsPowheg(branches)
	elif generator == 'mcatnlo' :
		return findInitialPartonsMCAtNLO(branches)
	elif generator == 'pythia8' :
		print 'PYTHIA NOT YET SUPPORTED'
		return None
	else :
		print 'ERROR: GENERATOR '+generator+' NOT RECOGNIZED!!!'
		return None

#returns a tuple of the (top, antitop) fourvectors from Monte Carlo
def findMCTops(branches,generator) :
	found_t    = False; found_tbar = False
	tvec, tbarvec = (TLorentzVector(1.0,0.0,0.0,1.0),TLorentzVector(1.0,0.0,0.0,1.0))
	for i in range(branches['gen_size'].getReadValue()) :
		if found_t and found_tbar :
			break
		if branches['gen_ID'].getReadValue(i) == TOP_ID and not found_t :
			tvec.SetPtEtaPhiE(branches['gen_Pt'].getReadValue(i),branches['gen_Eta'].getReadValue(i),branches['gen_Phi'].getReadValue(i),branches['gen_E'].getReadValue(i))
			found_t = True
		elif branches['gen_ID'].getReadValue(i) == -1*TOP_ID and not found_tbar :
			tbarvec.SetPtEtaPhiE(branches['gen_Pt'].getReadValue(i),branches['gen_Eta'].getReadValue(i),branches['gen_Phi'].getReadValue(i),branches['gen_E'].getReadValue(i))
			found_tbar = True
	return (tvec,tbarvec)

###################################  HELPER FUNCTIONS  ##################################

#Powheg typeCheck function
def typeCheckPowheg(branches,event_type) :
	addTwice = False
	#check the initial state partons if necessary
	if event_type<2 :
		initial_parton_ids = []
		#loop through and find the pdg IDs of the first daughters of the protons
		for i in range(branches['gen_size'].getReadValue()) :
			if branches['gen_Status'].getReadValue(i)!=3 :
				continue
			if branches['gen_ID'].getReadValue(i)==PROTON_ID :
				initial_parton_ids.append(branches['gen_Dau0ID'].getReadValue(i))
		#is it a qqbar event?
		is_qq = len(initial_parton_ids) == 2 and initial_parton_ids[0]+initial_parton_ids[1]==0
		#regardless, did it have an initially symmetric state
		addTwice = event_type==0 or (len(initial_parton_ids)==2 and initial_parton_ids[0] == initial_parton_ids[1])
		#return false if the event type is incorrect
		if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
			return False,addTwice
	return semilepCheck(branches,event_type),addTwice

#MC@NLO eventTypeCheck function
def typeCheckMCAtNLO(branches,event_type) :
	#print '-------------------------------------' #DEBUG
	addTwice = False
	#check the initial state partons if necessary
	if event_type<2 :
		ip_ids_t = []
		ip_ids_tbar = []
		#loop through and find the particles whose daughters include the ttbar pair
		for i in range(branches['gen_size'].getReadValue()) :
			thisID = branches['gen_ID'].getReadValue(i)
			Mom0ID = branches['gen_Mom0ID'].getReadValue(i)
			Mom1ID = branches['gen_Mom1ID'].getReadValue(i)
			Dau0ID = branches['gen_Dau0ID'].getReadValue(i)
			#find the tops that have two mothers and decayed only into themselves
			if abs(thisID)!=TOP_ID or Mom0ID==-900 or Mom1ID==-900 or Dau0ID!=thisID :
				continue
	#		print '[%d + %d] = %d -> (%d + %d)'%(Mom0ID,Mom1ID,thisID,Dau0ID,branches['gen_Dau1ID'].getReadValue(i)) #DEBUG
			if thisID>0 :
				ip_ids_t.append(Mom0ID); ip_ids_t.append(Mom1ID)
			else :
				ip_ids_tbar.append(Mom0ID); ip_ids_tbar.append(Mom1ID)
		#is it a qqbar event?
		is_qq = len(ip_ids_t)==2 and len(ip_ids_tbar)==2 and ip_ids_t[0]+ip_ids_t[1]==0 and ip_ids_tbar[0]+ip_ids_tbar[1]==0 and GLUON_ID not in ip_ids_t and GLUON_ID not in ip_ids_tbar
	#	print 'ISQQ = '+str(is_qq) #DEBUG
		#regardless, did it have an initially symmetric state
		addTwice = event_type==0 or (len(ip_ids_t)==2 and len(ip_ids_tbar)==2 and ip_ids_t[0]==GLUON_ID and ip_ids_t[1]==GLUON_ID and ip_ids_tbar[0]==GLUON_ID and ip_ids_tbar[1]==GLUON_ID)
	#	print 'ADDTWICE = '+str(addTwice) #DEBUG
		#return false if the event type is incorrect
		if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
			return False,addTwice
	return semilepCheck(branches,event_type),addTwice

#semileptonic checker
def semilepCheck(branches,event_type) :
	#check the decay type
	leptons_from_Ws_IDs = []
	for i in range(branches['gen_size'].getReadValue()) :
		thisID = branches['gen_ID'].getReadValue(i)
		Mom0ID = branches['gen_Mom0ID'].getReadValue(i)
		Mom1ID = branches['gen_Mom1ID'].getReadValue(i)
		Dau0ID = branches['gen_Dau0ID'].getReadValue(i)
		Dau1ID = branches['gen_Dau1ID'].getReadValue(i)
		#look for the decaying Ws that have themselves as their only mother and daughters that are not themselves
		if abs(thisID)!=W_ID or Mom0ID!=thisID or Mom1ID!=-900 or Dau0ID==thisID or Dau1ID==thisID or Dau0ID==-900 or Dau1ID==-900 :
			continue
	#	print '[%d + %d] = %d -> (%d + %d)'%(Mom0ID,Mom1ID,thisID,Dau0ID,Dau1ID) #DEBUG
		if abs(Dau0ID) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) :
			leptons_from_Ws_IDs.append(Dau0ID)
		if abs(Dau1ID) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) :
			leptons_from_Ws_IDs.append(Dau1ID)
	if ( (len(leptons_from_Ws_IDs) == 2 and event_type > 1) or (len(leptons_from_Ws_IDs) == 4 and event_type != 2) or 
		(len(leptons_from_Ws_IDs) == 0 and event_type != 3) ) :
	#	print 'LEP WRONG' #DEBUG
		return False
	#print 'LEP RIGHT' #DEBUG
	return True

#Powheg findInitialPartons function
def findInitialPartonsPowheg(branches) :
	#for i in range(branches['gen_size'].getReadValue()) :
	#	if branches['gen_Status'].getReadValue(i)!=3 :
	#		continue
	#	if branches['gen_ID'].getReadValue(i) == PROTON_ID and branches['gen_Dau0ID'].getReadValue(i) > 0:
	#		factor = 0.0
	#		if branches['gen_Dau0Eta'].getReadValue(i) > 0 : 
	#			factor = 1.0
	#		else :
	#			factor = -1.0
	#		return TLorentzVector(1.0,0.0,factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
	return TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY), TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN UNTIL NEW NTUPLES

#MC@NLO findInitialPartons function
def findInitialPartonsMCAtNLO(branches) :
	#print '------------------------------------------------' #DEBUG
	for i in range(branches['gen_size'].getReadValue()) :
	#	thisID = branches['gen_ID'].getReadValue(i) #DEBUG
	#	Mom0ID = branches['gen_Mom0ID'].getReadValue(i) #DEBUG
	#	Mom1ID = branches['gen_Mom1ID'].getReadValue(i) #DEBUG
	#	Dau0ID = branches['gen_Dau0ID'].getReadValue(i) #DEBUG
	#	print '[%d + %d] = %d -> (%d + %d)'%(Mom0ID,Mom1ID,thisID,Dau0ID,branches['gen_Dau1ID'].getReadValue(i)) #DEBUG
		if abs(branches['gen_Dau0ID'].getReadValue(i)) == TOP_ID and abs(branches['gen_Dau1ID'].getReadValue(i)) == TOP_ID and branches['gen_ID'].getReadValue(i) > 0 :
			factor = 0.0
			if branches['gen_Eta'] > 0 :
				factor = 1.0
			else :
				factor = -1.0
			return TLorentzVector(1.0,0.0,factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
	return TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY), TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN UNTIL NEW NTUPLES
