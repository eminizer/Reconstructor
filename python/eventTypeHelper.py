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
from memsniffer import checkmem

#################################  REFERENCED FUNCTIONS  ##################################

#given the MC gen info branches and the event type, returns a tuple of (keep [bool], addtwice [0 or 1])
def getEventType(branches) :
	#event types:
	#	0 = semileptonic qqbar
	#	1 = semileptonic gg, etc.
	#	2 = dileptonic
	#	3 = hadronic
	#	4 = everything else
	etype = -1
	p1id = branches['MC_part1_ID'].getReadValue(); p2id = branches['MC_part2_ID'].getReadValue()
	#first check the initial parton IDs to assign the production mechanism
	if p1id==p2id==-9999. : #neither initial parton ID was filled because nothing made a ttbar pair; this is "other" background so just return
		return 4
	elif p1id+p2id==0 : #parton IDs are opposites of each other, qqbar
		etype=0
	else :
		etype=1 #otherwise we'll say it's qg/gg for now
	#next check whether it's background
	isbkg = branches['MC_lepb_pt'].getReadValue()==-9999. and branches['MC_hadb_pt'].getReadValue()==-9999.
	if not isbkg : #if it's signal return its current event type
		return etype
	haslepside = branches['MC_lep_pt'].getReadValue()!=-9999.
	hashadside = branches['MC_hadW_pt'].getReadValue()!=-9999.
	if haslepside : #if it's ttbar background with a leptonic side it's dileptonic
		return 2
	elif hashadside : #if it's ttbar background with a hadronic side it's fully hadronic
		return 3
	else : #if we made it here then something is screwy
		print 'WARNING: COULD NOT IDENTIFY EVENT TYPE!! (p1ID = %d, p2ID = %d, isbkg = %s, haslepside = %s, hashadside = %s)'%(p1id,p2id,isbkg,haslepside,hashadside)
		return 5

#finds which "down-the-beampipe" fourvector corresponds to the initial quark
def findInitialPartons(branches) :
	checkmem('&&&','step_into_findInitialPartons')
	factor = branches['MC_part1_factor'].getReadValue()
	ID = branches['MC_part1_ID'].getReadValue()
	qvec = TLorentzVector(1.0,0.0,(abs(ID)/ID)*factor*sqrt(BEAM_ENERGY*BEAM_ENERGY - 1.*1.),BEAM_ENERGY)
	qbar_vec = TLorentzVector(1.0,0.0,-1.*(abs(ID)/ID)*factor*sqrt(BEAM_ENERGY*BEAM_ENERGY - 1.*1.),BEAM_ENERGY)
	checkmem('&&&','ret_from_findInitialPartons')
	return qvec,qbar_vec

#returns a bunch of fourvectors from Monte Carlo for MC truth reconstruction and matching
def findMCParticles(branches) :
	checkmem('&&&','step_into_findMCParticles')
	returnlist = []
	names = ['t','tbar','lep','nu','lepb','hadW','hadb'] 
	for name in names :
		checkmem('&&&&','read_branchvals_for_'+name)
		thispt = branches['MC_'+name+'_pt'].getReadValue()
		thiseta = branches['MC_'+name+'_eta'].getReadValue()
		thisphi = branches['MC_'+name+'_phi'].getReadValue()
		thisE = branches['MC_'+name+'_E'].getReadValue()
		thisvec = None
		if thispt!=-9999. :
			checkmem('&&&&','make_MCvec_for_'+name)
			thisvec = TLorentzVector()
			thisvec.SetPtEtaPhiE(thispt,thiseta,thisphi,thisE)
		checkmem('&&&&','append_MCvec_for_'+name)
		returnlist.append(thisvec)
		checkmem('&&&&','done_with_'+name)
	lepID = branches['MC_lep_ID'].getReadValue()
	returnlist.append(-1.*abs(lepID)/lepID)
	checkmem('&&&','ret_from_findMCParticles')
	return returnlist
