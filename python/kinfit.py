#Imports
from ROOT import TMinuit, Long, Double, TLorentzVector
from math import *
from array import array
import multiprocessing

#Global variables
SIGMAJ  = 0.10 #jet momentum resolution
SIGMAL  = 0.02 #lepton momentum resolution

class KinFit(object) :

	lep_global_vecs = []
	met_global_vecs = []
	blep_global_vecs = []
	hadt_global_vecs = []
	had1_global_vecs = []
	had2_global_vecs = []
	had3_global_vecs = []

	def __init__(self,topology,fitindex,hypothesis) :
		self.fitindex = fitindex
		#Set up Minuit
		self.parNames = ['fitindex','pZv','scalelep','scaleblep','scalehad1','scalehad2','scalehad3']
		self.parinivals = [fitindex,0.,1.,1.,1.,1.,1.]
		self.parerrs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
		#fit setup
		self.nPars = 6 if topology==1 else 7
		self.minuit = TMinuit(self.nPars)
		self.ierflag = Long(1)
		self.arglist = array( 'd', [-1.0] )
		self.minuit.mnexcm('SET PRINT', self.arglist, 1,self.ierflag)
		self.minuit.mnexcm('SET NOWARNINGS',self.arglist,1,self.ierflag)
		self.arglist[0] = 100000.
		#add parameters to the fitter
		for i in range(self.nPars) :
			self.minuit.mnparm(i,self.parNames[i],self.parinivals[i],1.0,0.,0.,self.ierflag)
		#fix the first 'fit index' parameter
		self.minuit.mnfixp(0,self.ierflag)
		#set the minimization function to use
		if topology==1 : self.minuit.SetFCN(KinFit.fcnt1)
		elif topology==2 : self.minuit.SetFCN(KinFit.fcnt2)
		elif topology==3 : self.minuit.SetFCN(KinFit.fcnt3)
		#list of final best fit parameters
		self.bestParValues = []
		for i in range(1,self.nPars) :
			self.bestParValues.append(self.parinivals[i])
		#fit error flag
		self.errflag = -1
		#fit final chi2
		self.chi2value = 1000000.
		#load this fit's fourvectors to the shared fourvector lists
		KinFit.__addHypothesisVectors__(topology,hypothesis)

	def dofit(self) :
		#minimize and get the error flag
		self.minuit.mnexcm('MIGRAD', self.arglist, 1, self.ierflag)
		#print 'errflag = '+str(self.ierflag) #DEBUG
		self.errflag = self.ierflag
		#Check fit Chi2 
		tmp1 = Double(1.0); tmp2 = Double(1.0); tmp3 = Double(1.0)
		self.minuit.mnstat(tmp1,tmp2,tmp3,Long(1),Long(1),Long(1))
		#print 'chi2 = %.4f'%(tmp1) #DEBUG
		self.chi2value = tmp1
		#Get the bestfit parameters back from minuit
		for j in range(1,self.nPars) :
			tmp = Double(1.0)
			self.minuit.GetParameter(j,tmp,Double(self.parerrs[j]))
			self.bestParValues[j] = tmp

	def getErrFlag(self) :
		return self.errflag
	def getFitChi2(self) :
		return self.chi2value
	def getIndex(self) :
		return self.fitindex
	def getBestParValues(self) :
		return self.bestParValues

	@staticmethod
	def __addHypothesisVectors__(t,h) :
		KinFit.lep_global_vecs.append(h[0].getFourVector())
		KinFit.met_global_vecs.append(h[1])
		KinFit.blep_global_vecs.append(h[2].getFourVector())
		if t==1 :
			KinFit.hadt_global_vecs.append(h[3][0].getFourVector())
			KinFit.had1_global_vecs.append(h[3][0].getSubjet(0).getFourVector())
			KinFit.had2_global_vecs.append(h[3][0].getSubjet(1).getFourVector())
		else :
			KinFit.had1_global_vecs.append(h[3][0].getFourVector())
			KinFit.had2_global_vecs.append(h[3][1].getFourVector())
			KinFit.had3_global_vecs.append(h[3][2].getFourVector())

	##################################  Static minuit functions  ##################################

	#type 1 (fully merged) top minimization function
	@staticmethod
	def fcnt1(npar, deriv, f, par, flag) :
		ind = int(par[0])
		#Build rescaled versions of the vectors involved
		l = rescale(KinFit.lep_global_vecs[ind],par[2])
		bl = rescale(KinFit.blep_global_vecs[ind],par[3])
		hs1 = rescale(KinFit.had1_global_vecs[ind],par[4])
		hs2 = rescale(KinFit.had2_global_vecs[ind],par[5])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( KinFit.met_global_vecs[ind].Px()+(1.0-par[2])*KinFit.lep_global_vecs[ind].Px()+(1.0-par[3])*KinFit.blep_global_vecs[ind].Px()
						+(1.0-par[4])*KinFit.had1_global_vecs[ind].Px()+(1.0-par[5])*KinFit.had2_global_vecs[ind].Px() )
		newmety = ( KinFit.met_global_vecs[ind].Py()+(1.0-par[2])*KinFit.lep_global_vecs[ind].Py()+(1.0-par[3])*KinFit.blep_global_vecs[ind].Py()
						+(1.0-par[4])*KinFit.had1_global_vecs[ind].Py()+(1.0-par[5])*KinFit.had2_global_vecs[ind].Py() )
		v = rescale(KinFit.met_global_vecs[ind],1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[1])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; th = KinFit.hadt_global_vecs[ind]-KinFit.had1_global_vecs[ind]-KinFit.had2_global_vecs[ind]+hs1+hs2
		pdf = getPDF(tl.M(),KinFit.MT_l1,KinFit.WT_l1,th.M(),KinFit.MT_h1,
					 KinFit.WT_h1,KinFit.MT,wl.M2(),KinFit.MW,KinFit.WW)
		f[0] = ( pdf+(par[2]-1.)*(par[2]-1.)/(SIGMAL*SIGMAL)+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)
					+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ)+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )

	#type 2 (boosted untagged) top minimization function
	@staticmethod
	def fcnt2(npar, deriv, f, par, flag) :
		ind = int(par[0])
		#Build rescaled versions of the vectors involved
		l = rescale(KinFit.lep_global_vecs[ind],par[2])
		bl = rescale(KinFit.blep_global_vecs[ind],par[3])
		h1 = rescale(KinFit.had1_global_vecs[ind],par[4])
		h2 = rescale(KinFit.had2_global_vecs[ind],par[5])
		h3 = rescale(KinFit.had3_global_vecs[ind],par[6])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( KinFit.met_global_vecs[ind].Px()+(1.0-par[2])*KinFit.lep_global_vecs[ind].Px()+(1.0-par[3])*KinFit.blep_global_vecs[ind].Px()
						+(1.0-par[4])*KinFit.had1_global_vecs[ind].Px()+(1.0-par[5])*KinFit.had2_global_vecs[ind].Px()+(1.0-par[6])*KinFit.had3_global_vecs[ind].Px() )
		newmety = ( KinFit.met_global_vecs[ind].Py()+(1.0-par[2])*KinFit.lep_global_vecs[ind].Py()+(1.0-par[3])*KinFit.blep_global_vecs[ind].Py()
						+(1.0-par[4])*KinFit.had1_global_vecs[ind].Py()+(1.0-par[5])*KinFit.had2_global_vecs[ind].Py()+(1.0-par[6])*KinFit.had3_global_vecs[ind].Py() )
		v = rescale(KinFit.met_global_vecs[ind],1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[1])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; 
		th = h1+h2+h3
		pdf = getPDF(tl.M(),KinFit.MT_l2,KinFit.WT_l2,th.M(),KinFit.MT_h2,
					 KinFit.WT_h2,KinFit.MT,wl.M2(),KinFit.MW,KinFit.WW)
		f[0] = ( pdf+(par[2]-1.)*(par[2]-1.)/(SIGMAL*SIGMAL)+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)
					+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ)+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) 
					+(par[6]-1.)*(par[6]-1.)/(SIGMAJ*SIGMAJ) )

	#type 3 (resolved) kinematic fitting function
	@staticmethod
	def fcnt3(npar, deriv, f, par, flag) :
		ind = int(par[0])
		#Build rescaled versions of the vectors involved
		l = rescale(KinFit.lep_global_vecs[ind],par[2])
		bl = rescale(KinFit.blep_global_vecs[ind],par[3])
		h1 = rescale(KinFit.had1_global_vecs[ind],par[4])
		h2 = rescale(KinFit.had2_global_vecs[ind],par[5])
		h3 = rescale(KinFit.had3_global_vecs[ind],par[6])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( KinFit.met_global_vecs[ind].Px()+(1.0-par[2])*KinFit.lep_global_vecs[ind].Px()+(1.0-par[3])*KinFit.blep_global_vecs[ind].Px()
						+(1.0-par[4])*KinFit.had1_global_vecs[ind].Px()+(1.0-par[5])*KinFit.had2_global_vecs[ind].Px()+(1.0-par[6])*KinFit.had3_global_vecs[ind].Px() )
		newmety = ( KinFit.met_global_vecs[ind].Py()+(1.0-par[2])*KinFit.lep_global_vecs[ind].Py()+(1.0-par[3])*KinFit.blep_global_vecs[ind].Py()
						+(1.0-par[4])*KinFit.had1_global_vecs[ind].Py()+(1.0-par[5])*KinFit.had2_global_vecs[ind].Py()+(1.0-par[6])*KinFit.had3_global_vecs[ind].Py() )
		v = rescale(KinFit.met_global_vecs[ind],1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[1])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; 
		th = h1+h2+h3
		pdf = getPDF(tl.M(),KinFit.MT_l3,KinFit.WT_l3,th.M(),KinFit.MT_h3,
					 KinFit.WT_h3,KinFit.MT,wl.M2(),KinFit.MW,KinFit.WW)
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

def reconstructParallel(kfo) :
	kfo.dofit()

#reconstruct takes in the list of hypotheses and returns a tuple of best fit information invluded the corrected fourvectors and chi2 value
def reconstruct(hypotheses,topology) :
	#First declare and set up the kinematic fit objects
	allkinfitobjs = []
	for i in range(len(hypotheses)) :
		allkinfitobjs.append(KinFit(topology,i,hypotheses[i]))
	#Now do all the fits in parallel
	procs = []
	for kfo in allkinfitobjs :
		p = multiprocessing.Process(target=reconstructParallel,args=(kfo,))
		p.start()
		procs.append(p)
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