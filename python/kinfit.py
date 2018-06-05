#Imports
from ROOT import TMinuit, Long, Double, TLorentzVector
from math import *
from array import array
import multiprocessing

#Global variables
SIGMAJ  = 0.10 #jet momentum resolution
SIGMAL  = 0.02 #lepton momentum resolution
MW = 80.4 #W mass
MT = 172.5 #top mass
MT_l1 = 172.6 #leptonic top mass for type 1 tops
MT_h1 = 182.6 #hadronic top mass for type 1 tops
MT_l2 = 172.9 #leptonic top mass for type 2 tops
MT_h2 = 173.7 #hadronic top mass for type 2 tops
MT_l3 = 172.2 #leptonic top mass for type 3 tops
MT_h3 = 168.1 #hadronic top mass for type 3 tops
WW = 2.0 #W width
WT = 1.4 #top width
WT_l1 = 15.9 #leptonic top width for type 1 tops
WT_h1 = 17.3 #hadronic top width for type 1 tops
WT_l2 = 15.7 #leptonic top width for type 2 tops
WT_h2 = 17.7 #hadronic top width for type 2 tops
WT_l3 = 16.0 #leptonic top width for type 3 tops
WT_h3 = 19.3 #hadronic top width for type 3 tops
lep_gvs=[]
met_gvs=[]
blep_gvs=[]
hadt_gvs=[]
had1_gvs=[]
had2_gvs=[]
had3_gvs=[]
allveclists = [lep_gvs,met_gvs,blep_gvs,hadt_gvs,had1_gvs,had2_gvs,had3_gvs]

class KinFit(object) :

	def __init__(self,topology,fitindex,hypothesis) :
		self.topology = topology
		self.fitindex = fitindex
		#list of final best fit parameters
		self.bestParValues = []
		#fit error flag
		self.errflag = -1
		#fit final chi2
		self.chi2value = 1000000.
		#load this fit's fourvectors to the shared fourvector lists
		KinFit.__addHypothesisVectors__(topology,hypothesis,self.fitindex==0)
		#print 'just initialized fit with index %d, length of list of global lepton vectors = %d'%(self.fitindex,len(KinFit.lep_gvs))#DEBUG

	def dofit(self) :
		#print 'starting fit for hypothesis index %d'%(self.fitindex)#DEBUG
		#set up the minuit object, etc.
		parNames = ['fitindex','pZv','scalelep','scaleblep','scalehad1','scalehad2','scalehad3']
		parinivals = [self.fitindex,0.,1.,1.,1.,1.,1.]
		parerrs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
		bestParValues = []
		nPars = 6 if self.topology==1 else 7
		for i in range(1,nPars) :
			bestParValues.append(parinivals[i])
		ierflag = Long(1)
		arglist = array( 'd', [-1.0] )
		minuit = TMinuit(nPars)
		minuit.mnexcm('SET PRINT', arglist, 1,ierflag)
		minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
		arglist[0] = 100000.
		#add parameters to the fitter
		for i in range(nPars) :
			minuit.mnparm(i,parNames[i],parinivals[i],1.0,0.,0.,ierflag)
		#fix the first 'fit index' parameter
		minuit.mnfixp(0,ierflag)
		#set the minimization function to use
		if self.topology==1 : minuit.SetFCN(fcnt1)
		elif self.topology==2 : minuit.SetFCN(fcnt2)
		elif self.topology==3 : minuit.SetFCN(fcnt3)
		#minimize and get the error flag
		minuit.mnexcm('MIGRAD', arglist, 1, ierflag)
		#print 'index = %d, errflag = %s'%(self.fitindex,ierflag) #DEBUG
		errflag = int(ierflag)
		#Check fit Chi2 
		tmp1 = Double(1.0); tmp2 = Double(1.0); tmp3 = Double(1.0)
		minuit.mnstat(tmp1,tmp2,tmp3,Long(1),Long(1),Long(1))
		#print 'chi2 = %.4f'%(tmp1) #DEBUG
		chi2value = float(tmp1)
		#Get the bestfit parameters back from minuit
		for j in range(1,nPars) :
			tmp = Double(1.0)
			minuit.GetParameter(j,tmp,Double(parerrs[j]))
			bestParValues[j-1] = float(tmp)
		return errflag, chi2value, bestParValues
		#print '	done fit for hypothesis index %d, pZv = %.3f, had1scale = %.3f, chi2 = %.2f'%(self.fitindex,self.bestParValues[0],self.bestParValues[3],self.chi2value)#DEBUG

	def getErrFlag(self) :
		return self.errflag
	def setErrFlag(self,nef) :
		self.errflag=nef
	def getFitChi2(self) :
		return self.chi2value
	def setFitChi2(self,nc2v) :
		self.chi2value=nc2v
	def getIndex(self) :
		return self.fitindex
	def getBestParValues(self) :
		return self.bestParValues
	def setBestParValues(self,nbpv) :
		self.bestParValues=nbpv

	@staticmethod
	def __addHypothesisVectors__(t,h,clearvecs) :
		if clearvecs :
			for veclist in allveclists :
				while len(veclist)>0 :
					veclist.pop()
		lep_gvs.append(h[0].getFourVector())
		met_gvs.append(h[1])
		blep_gvs.append(h[2].getFourVector())
		if t==1 :
			hadt_gvs.append(h[3][0].getFourVector())
			had1_gvs.append(h[3][0].getSubjet(0).getFourVector())
			had2_gvs.append(h[3][0].getSubjet(1).getFourVector())
		else :
			had1_gvs.append(h[3][0].getFourVector())
			had2_gvs.append(h[3][1].getFourVector())
			had3_gvs.append(h[3][2].getFourVector())

	##################################  Static minuit functions  ##################################

#type 1 (fully merged) top minimization function
def fcnt1(npar, deriv, f, par, flag) :
	ind = int(par[0])
	#Build rescaled versions of the vectors involved
	l = rescale(lep_gvs[ind],par[2])
	bl = rescale(blep_gvs[ind],par[3])
	hs1 = rescale(had1_gvs[ind],par[4])
	hs2 = rescale(had2_gvs[ind],par[5])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_gvs[ind].Px()+(1.0-par[2])*lep_gvs[ind].Px()+(1.0-par[3])*blep_gvs[ind].Px()
					+(1.0-par[4])*had1_gvs[ind].Px()+(1.0-par[5])*had2_gvs[ind].Px() )
	newmety = ( met_gvs[ind].Py()+(1.0-par[2])*lep_gvs[ind].Py()+(1.0-par[3])*blep_gvs[ind].Py()
					+(1.0-par[4])*had1_gvs[ind].Py()+(1.0-par[5])*had2_gvs[ind].Py() )
	v = rescale(met_gvs[ind],1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[1])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl; th = hadt_gvs[ind]-had1_gvs[ind]-had2_gvs[ind]+hs1+hs2
	pdf = getPDF(tl.M(),MT_l1,WT_l1,th.M(),MT_h1,WT_h1,MT,wl.M2(),MW,WW)
	f[0] = ( pdf+(par[2]-1.)*(par[2]-1.)/(SIGMAL*SIGMAL)+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)
				+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ)+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )

#type 2 (boosted untagged) top minimization function
def fcnt2(npar, deriv, f, par, flag) :
	ind = int(par[0])
	#Build rescaled versions of the vectors involved
	l = rescale(lep_gvs[ind],par[2])
	bl = rescale(blep_gvs[ind],par[3])
	h1 = rescale(had1_gvs[ind],par[4])
	h2 = rescale(had2_gvs[ind],par[5])
	h3 = rescale(had3_gvs[ind],par[6])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_gvs[ind].Px()+(1.0-par[2])*lep_gvs[ind].Px()+(1.0-par[3])*blep_gvs[ind].Px()
					+(1.0-par[4])*had1_gvs[ind].Px()+(1.0-par[5])*had2_gvs[ind].Px()+(1.0-par[6])*had3_gvs[ind].Px() )
	newmety = ( met_gvs[ind].Py()+(1.0-par[2])*lep_gvs[ind].Py()+(1.0-par[3])*blep_gvs[ind].Py()
					+(1.0-par[4])*had1_gvs[ind].Py()+(1.0-par[5])*had2_gvs[ind].Py()+(1.0-par[6])*had3_gvs[ind].Py() )
	v = rescale(met_gvs[ind],1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[1])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl; 
	th = h1+h2+h3
	pdf = getPDF(tl.M(),MT_l2,WT_l2,th.M(),MT_h2,WT_h2,MT,wl.M2(),MW,WW)
	f[0] = ( pdf+(par[2]-1.)*(par[2]-1.)/(SIGMAL*SIGMAL)+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)
				+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ)+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) 
				+(par[6]-1.)*(par[6]-1.)/(SIGMAJ*SIGMAJ) )

#type 3 (resolved) kinematic fitting function
def fcnt3(npar, deriv, f, par, flag) :
	ind = int(par[0])
	#Build rescaled versions of the vectors involved
	l = rescale(lep_gvs[ind],par[2])
	bl = rescale(blep_gvs[ind],par[3])
	h1 = rescale(had1_gvs[ind],par[4])
	h2 = rescale(had2_gvs[ind],par[5])
	h3 = rescale(had3_gvs[ind],par[6])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_gvs[ind].Px()+(1.0-par[2])*lep_gvs[ind].Px()+(1.0-par[3])*blep_gvs[ind].Px()
					+(1.0-par[4])*had1_gvs[ind].Px()+(1.0-par[5])*had2_gvs[ind].Px()+(1.0-par[6])*had3_gvs[ind].Px() )
	newmety = ( met_gvs[ind].Py()+(1.0-par[2])*lep_gvs[ind].Py()+(1.0-par[3])*blep_gvs[ind].Py()
					+(1.0-par[4])*had1_gvs[ind].Py()+(1.0-par[5])*had2_gvs[ind].Py()+(1.0-par[6])*had3_gvs[ind].Py() )
	v = rescale(met_gvs[ind],1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[1])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl; 
	th = h1+h2+h3
	pdf = getPDF(tl.M(),MT_l3,WT_l3,th.M(),MT_h3,WT_h3,MT,wl.M2(),MW,WW)
	f[0] = ( pdf+(par[2]-1.)*(par[2]-1.)/(SIGMAL*SIGMAL)+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)
				+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ)+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) 
			+(par[6]-1.)*(par[6]-1.)/(SIGMAJ*SIGMAJ) )

#PDF calculator function
def getPDF(mtl,MT_l,WT_l,mth,MT_h,WT_h,MT,mwl2,MW,WW) :
	pdftl = (mtl-MT_l)**2/WT_l**2
	pdfth = (mth-MT_h)**2/WT_h**2
	lwl = ((MT**2-mwl2)*(MT**2-mwl2)*(2*MT**2+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW)
	pdf = pdftl+pdfth
	if lwl > 0.0 :
		pdf += -2.0*log(lwl)    #need positive f
	else :
		#print 'WARNING -- pdf is negative!!!'
		pdf += -2.0*log(1.e-50)
	#print 'mtl = %.3f, mth = %.3f, pdftl = %.3f, pdfth = %.3f, lwl = %.3f, pdf = %.3f'%(mtl,mth,pdftl,pdfth,lwl,pdf) #DEBUG
	return pdf

#fourvector rescaling function
def rescale(vec,fac) :
	p2 = fac*fac*vec.Vect().Mag2()
	m2 = vec.M()*vec.M()
	newE = sqrt(p2+m2)
	#note that the function returns a NEW TLorentzVector, meaning the original is unaltered
	return TLorentzVector(fac*vec.Px(),fac*vec.Py(),fac*vec.Pz(),newE)

def reconstructParallel(kfo,c_fpt) :
	fitinfotuple = kfo.dofit()
	c_fpt.send(fitinfotuple)
	c_fpt.close()


#reconstruct takes in the list of hypotheses and returns a tuple of best fit information invluded the corrected fourvectors and chi2 value
def reconstruct(hypotheses,topology,ongrid) :
	#print '-----------------------------'#DEBUG
	#First declare and set up the kinematic fit objects (they'll be done in parallel in batches)
	allkinfitobjlists = []; allkinfitobjs = []
	j=-1
	batchsize = 1 if ongrid=='yes' else 5
	for i in range(len(hypotheses)) :
		if i%batchsize==0 :
			j+=1
			allkinfitobjlists.append([])
		newfit = KinFit(topology,i,hypotheses[i])
		allkinfitobjlists[j].append(newfit)
		allkinfitobjs.append(newfit)
	#Now do all the fits in parallel
	for kinfitobjlist in allkinfitobjlists :
		procs = []
		for kfo in kinfitobjlist :
			p_fpt, c_fpt = multiprocessing.Pipe()
			#kfo.dofit()
			p = multiprocessing.Process(target=reconstructParallel,args=(kfo,c_fpt))
			p.start()
			procs.append(p)
			newfitinfotuple = p_fpt.recv()
			kfo.setErrFlag(newfitinfotuple[0])
			kfo.setFitChi2(newfitinfotuple[1])
			kfo.setBestParValues(newfitinfotuple[2])
		for p in procs :
			p.join()
	#find the index of the best fit hypothesis
	bestfitchi2 = 1000000.; bestfitindex = -1; bestParValues = None
	for kfo in allkinfitobjs :
		if kfo.getErrFlag()==0 and kfo.getFitChi2()<bestfitchi2 :
			bestfitchi2 = kfo.getFitChi2(); bestfitindex = kfo.getIndex(); bestParValues = kfo.getBestParValues()
	#print 'bestfitindex = %d, bestfitchi2 = %.4f'%(bestfitindex,bestfitchi2) #DEBUG
	#if no fits converged, return garbage
	if bestfitindex==-1 :
		return -1, None, None, None, None, None, None, None, 1000000000., []
	#get the best fit hypothesis
	besthypothesis = hypotheses[bestfitindex]
	#start with the MET as given from the hypothesis
	final_met = TLorentzVector(); final_met.SetPtEtaPhiE(besthypothesis[1].Pt(),besthypothesis[1].Eta(),besthypothesis[1].Phi(),besthypothesis[1].E())
	#rescale the lepton and jet four vectors based on the final parameters
	hadt_return = None; had1_return = None; had2_return = None; had3_return = None
	orig_lep = besthypothesis[0].getFourVector()
	lep_return = rescale(orig_lep,bestParValues[1])
	orig_lepb = besthypothesis[2].getFourVector()
	lepb_return = rescale(orig_lepb,bestParValues[2])
	if topology==1 :
		orig_hadt = besthypothesis[3][0].getFourVector()
		orig_had1 = besthypothesis[3][0].getSubjet(0).getFourVector()
		had1_return = rescale(orig_had1,bestParValues[3])
		orig_had2 = besthypothesis[3][0].getSubjet(1).getFourVector()
		had2_return = rescale(orig_had2,bestParValues[4])
		hadt_return = orig_hadt-orig_had1+had1_return-orig_had2+had2_return
	else :
		orig_had1 = besthypothesis[3][0].getFourVector()
		had1_return = rescale(orig_had1,bestParValues[3])
		orig_had2 = besthypothesis[3][1].getFourVector()
		had2_return = rescale(orig_had2,bestParValues[4])
		orig_had3 = besthypothesis[3][2].getFourVector()
		had3_return = rescale(orig_had3,bestParValues[5])
		hadt_return=had1_return+had2_return+had3_return
	#rebuild the neutrino post-rescaling
	newmetx = final_met.Px()+ (1.0-bestParValues[1])*orig_lep.Px()
	newmety = final_met.Py()+ (1.0-bestParValues[1])*orig_lep.Py()
	newmetx += (1.0-bestParValues[2])*orig_lepb.Px()
	newmety += (1.0-bestParValues[2])*orig_lepb.Py()
	newmetx += (1.0-bestParValues[3])*orig_had1.Px()
	newmety += (1.0-bestParValues[3])*orig_had1.Py()
	newmetx += (1.0-bestParValues[4])*orig_had2.Px()
	newmety += (1.0-bestParValues[4])*orig_had2.Py()
	if had3_return!=None :
		newmetx += (1.0-bestParValues[5])*orig_had3.Px()
		newmety += (1.0-bestParValues[5])*orig_had3.Py()
	final_met.SetPx(newmetx); final_met.SetPy(newmety); final_met.SetPz(bestParValues[0])
	final_met.SetE(final_met.Vect().Mag())
	#return everything
	return (bestfitindex,lep_return,final_met,lepb_return,hadt_return,had1_return,had2_return,had3_return,bestfitchi2,bestParValues)