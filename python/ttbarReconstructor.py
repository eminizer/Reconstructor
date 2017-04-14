#Imports
from ROOT import TMinuit, Long, Double, TLorentzVector
from math import *
from array import array

#Global variables
SIGMAJ  = 0.10 #jet momentum resolution
SIGMAL  = 0.02 #lepton momentum resolution

class TTBarReconstructor(object) :
	#Class variables
	MW = 80.4 #W mass
	MT = 172.5 #top mass
	MT_l1 = 178.9 #leptonic top mass for type 1 tops
	MT_h1 = 180.0 #hadronic top mass for type 1 tops
	MT_l2 = 172.5 #leptonic top mass for type 2 tops
	MT_h2 = 142.0 #hadronic top mass for type 2 tops
	MT_l3 = 173.5 #leptonic top mass for type 3 tops
	MT_h3 = 153.8 #hadronic top mass for type 3 tops
	MT_l4 = 172.8 #leptonic top mass for type 4 tops
	MT_h4 = 171.3 #hadronic top mass for type 4 tops
	WW = 2.0 #W width
	WT = 1.4 #top width
	WT_l1 = 24.7 #leptonic top width for type 1 tops
	WT_h1 = 24.0 #hadronic top width for type 1 tops
	WT_l2 = 20.9 #leptonic top width for type 2 tops
	WT_h2 = 40.9 #hadronic top width for type 2 tops
	WT_l3 = 21.3 #leptonic top width for type 3 tops
	WT_h3 = 32.5 #hadronic top width for type 3 tops
	WT_l4 = 19.9 #leptonic top width for type 4 tops
	WT_h4 = 28.6 #hadronic top width for type 4 tops
	#global fourvectors for the fit
	lep_global_vec = TLorentzVector()
	met_global_vec = TLorentzVector()
	blep_global_vec = TLorentzVector()
	bigjet_global_vec = TLorentzVector()
	had1_global_vec = TLorentzVector()
	had2_global_vec = TLorentzVector()
	had3_global_vec = TLorentzVector()

	def __init__(self) :
		#Set up Minuit
		self.parNames = ['pZv','scalelep','scaleblep','scalehad1','scalehad2','scalehad3']
		self.parinivals = [0.,1.,1.,1.,1.,1.]
		self.parerrs = [0.0,0.0,0.0,0.0,0.0,0.0]
		#fit setup stuff common to every iteration
		self.minuit = TMinuit(6)
		self.ierflag = Long(1)
		self.arglist = array( 'd', [-1.0] )
		self.minuit.mnexcm('SET PRINT', self.arglist, 1,self.ierflag)
		self.minuit.mnexcm('SET NOWARNINGS',self.arglist,1,self.ierflag)
		self.arglist[0] = 100000.
		#add parameters to the fitter
		for i in range(len(self.parNames)) :
			self.minuit.mnparm(i,self.parNames[i],self.parinivals[i],1.0,0.,0.,self.ierflag)

	#reconstruct takes in the list of hypotheses and returns a tuple of fit information invluded the corrected fourvectors and chi2 value
	def reconstruct(self,hypotheses,topology) :
		#lists of final parameters and chi2 values
		bestParValues = []; finalChi2s = []; errflags = []
		for i in range(len(hypotheses)): 
			bestParValues.append([hypotheses[i][1].Pz(),1.0,1.0,1.0,1.0,1.0])
			finalChi2s.append(1000000000.)
			errflags.append(-1)
		#perform the fit once for each jet/MET hypothesis
		#print '------------------------------------------------------------------------' #DEBUG
		for i in range(len(hypotheses)) :
			freeparameters = False
			#print 'Running with hypothesis number = %d'%(i) #DEBUG
			hypothesis = hypotheses[i]
			nhsjets = len(hypothesis[3])
			#pop unneeded parameter values
			if topology==3 and nhsjets==1 :
				bestParValues[i].pop(5); bestParValues[i].pop(4)
				self.minuit.mnfixp(4,self.ierflag); self.minuit.mnfixp(4,self.ierflag)
				freeparameters=True
			elif topology==1 or (topology==3 and nhsjets==2) :
				bestParValues[i].pop(5)
				self.minuit.mnfixp(5,self.ierflag)
				freeparameters=True
			#set parameter values
			for j in range(len(bestParValues[i])) :
				self.minuit.mnparm(j,self.parNames[j],bestParValues[i][j],1.0,0,0,self.ierflag)
			#set constants
			lep_vec = hypothesis[0].getFourVector()
			met = hypothesis[1]
			lepb_vec = hypothesis[2].getFourVector()
			#set the global fourvector variables
			self.lep_global_vec.SetPtEtaPhiE(lep_vec.Pt(),lep_vec.Eta(),lep_vec.Phi(),lep_vec.E())
			self.met_global_vec.SetPtEtaPhiE(met.Pt(),met.Eta(),met.Phi(),met.E())
			self.blep_global_vec.SetPtEtaPhiE(lepb_vec.Pt(),lepb_vec.Eta(),lepb_vec.Phi(),lepb_vec.E())
			#the rest are topology-dependent
			if topology<3 :
				hypak8 = hypothesis[3][0] 
				hypak8v = hypak8.getFourVector()
				self.bigjet_global_vec.SetPtEtaPhiE(hypak8v.Pt(),hypak8v.Eta(),hypak8v.Phi(),hypak8v.E())
				h1v = hypak8.getSubjet(0).getFourVector() #had1 is the first subjet of the merged W
				self.had1_global_vec.SetPtEtaPhiE(h1v.Pt(),h1v.Eta(),h1v.Phi(),h1v.E())
				h2v = hypak8.getSubjet(1).getFourVector() #had2 is the second subjet of the merged W
				self.had2_global_vec.SetPtEtaPhiE(h2v.Pt(),h2v.Eta(),h2v.Phi(),h2v.E())
				if topology==1 : 
					self.minuit.SetFCN(TTBarReconstructor.fcnt1)
				elif topology==2 :
					h3v = hypothesis[3][1].getFourVector()
					self.had3_global_vec.SetPtEtaPhiE(h3v.Pt(),h3v.Eta(),h3v.Phi(),h3v.E())
					self.minuit.SetFCN(TTBarReconstructor.fcnt2)
			else :
				h1v = hypothesis[3][0].getFourVector()
				self.had1_global_vec.SetPtEtaPhiE(h1v.Pt(),h1v.Eta(),h1v.Phi(),h1v.E())
				if nhsjets==1 :
					self.minuit.SetFCN(TTBarReconstructor.fcnt31)
				else :
					h2v = hypothesis[3][1].getFourVector()
					self.had2_global_vec.SetPtEtaPhiE(h2v.Pt(),h2v.Eta(),h2v.Phi(),h2v.E())
					if nhsjets==2 :
						self.minuit.SetFCN(TTBarReconstructor.fcnt32)
					elif nhsjets==3 :
						h3v = hypothesis[3][1].getFourVector()
						self.had3_global_vec.SetPtEtaPhiE(h3v.Pt(),h3v.Eta(),h3v.Phi(),h3v.E())
						if topology==3 : self.minuit.SetFCN(TTBarReconstructor.fcnt33)
						elif topology==4 : self.minuit.SetFCN(TTBarReconstructor.fcnt4)
			#print 'self.lep_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.lep_global_vec.Pt(),self.lep_global_vec.Eta(),self.lep_global_vec.Phi(),self.lep_global_vec.M()) #DEBUG
			#print 'self.met_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.met_global_vec.Pt(),self.met_global_vec.Eta(),self.met_global_vec.Phi(),self.met_global_vec.M()) #DEBUG
			#print 'self.blep_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.blep_global_vec.Pt(),self.blep_global_vec.Eta(),self.blep_global_vec.Phi(),self.blep_global_vec.M()) #DEBUG
			#print 'self.had1_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.had1_global_vec.Pt(),self.had1_global_vec.Eta(),self.had1_global_vec.Phi(),self.had1_global_vec.M()) #DEBUG
			#print 'self.had2_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.had2_global_vec.Pt(),self.had2_global_vec.Eta(),self.had2_global_vec.Phi(),self.had2_global_vec.M()) #DEBUG
			#print 'self.had3_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.had3_global_vec.Pt(),self.had3_global_vec.Eta(),self.had3_global_vec.Phi(),self.had3_global_vec.M()) #DEBUG
			#print 'self.bigjet_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.bigjet_global_vec.Pt(),self.bigjet_global_vec.Eta(),self.bigjet_global_vec.Phi(),self.bigjet_global_vec.M()) #DEBUG
			#minimize
			self.minuit.mnexcm('MIGRAD', self.arglist, 1, self.ierflag)
			#print 'errflag = '+str(ierflag) #DEBUG
			errflags[i] = self.ierflag
			#Get the bestfit parameters back from minuit
			for j in range(len(bestParValues[i])) :
				tmp = Double(1.0)
				self.minuit.GetParameter(j,tmp,Double(self.parerrs[j]))
				bestParValues[i][j] = tmp
			#check to make sure it didn't just get stuck in the "negative LL" trap
			gotstuck = True 
			for j in range(1,len(bestParValues[i])) : 
				if bestParValues[i][j]!=1.0 :  
					gotstuck = False; break 
			if gotstuck :
				print 'This fit never had positive likelihood' #DEBUG
				errflags[i] = -1
			#Set fit Chi2 for this pZ solution
			tmp1 = Double(1.0); tmp2 = Double(1.0); tmp3 = Double(1.0)
			self.minuit.mnstat(tmp1,tmp2,tmp3,Long(1),Long(1),Long(1))
			#print 'chi2 = %.4f'%(tmp1) #DEBUG
			finalChi2s[i] = tmp1
			if freeparameters :
				self.minuit.mnfree(0)
				freeparameters=False
		#if no fits converged, return garbage
		noneconverged = True
		for errflag in errflags :
			if errflag==0 :
				noneconverged = False; break
		if noneconverged :
			return -1, None, None, None, None, None, None, None, 1000000000., []
		#find which fit gave the best results and record best parameter values
		final_pars = []
		bestfitchi2 = min(finalChi2s)
		bestfitindex = finalChi2s.index(bestfitchi2)
		for i in range(len(bestParValues[bestfitindex])) :
			final_pars.append(bestParValues[bestfitindex][i])
		besthypothesis = hypotheses[bestfitindex]
		final_met = TLorentzVector(); final_met.SetPtEtaPhiE(besthypothesis[1].Pt(),besthypothesis[1].Eta(),besthypothesis[1].Phi(),besthypothesis[1].E())
		#print 'bestfitindex = %d, bestfitchi2 = %.4f'%(bestfitindex,bestfitchi2) #DEBUG
		#rescale the lepton and jet four vectors based on the final parameters
		bigjet_return = None; had1_return = None; had2_return = None; had3_return = None
		orig_lep = besthypothesis[0].getFourVector()
		lep_return = rescale(orig_lep,final_pars[1])
		orig_lepb = besthypothesis[2].getFourVector()
		lepb_return = rescale(orig_lepb,final_pars[2])
		if topology<3 :
			orig_bigjet = besthypothesis[3][0].getFourVector()
			orig_had1 = besthypothesis[3][0].getSubjet(0).getFourVector()
			had1_return = rescale(orig_had1,final_pars[3])
			orig_had2 = besthypothesis[3][0].getSubjet(1).getFourVector()
			had2_return = rescale(orig_had2,final_pars[4])
			bigjet_return = orig_bigjet-orig_had1+had1_return-orig_had2+had2_return
			if topology==2 :
				orig_had3 = besthypothesis[3][1].getFourVector()
				had3_return = rescale(orig_had3,final_pars[5])
		else :
			orig_had1 = besthypothesis[3][0].getFourVector()
			had1_return = rescale(orig_had1,final_pars[3])
			if len(besthypothesis[3])>1 :
				orig_had2 = besthypothesis[3][1].getFourVector()
				had2_return = rescale(orig_had2,final_pars[4])
				if len(besthypothesis[3])>2 :
					orig_had3 = besthypothesis[3][2].getFourVector()
					had3_return = rescale(orig_had3,final_pars[5])
		#rebuild the neutrino post-rescaling
		newmetx = final_met.Px()+ (1.0-final_pars[1])*orig_lep.Px()
		newmety = final_met.Py()+ (1.0-final_pars[1])*orig_lep.Py()
		newmetx += (1.0-final_pars[2])*orig_lepb.Px()
		newmety += (1.0-final_pars[2])*orig_lepb.Py()
		newmetx += (1.0-final_pars[3])*orig_had1.Px()
		newmety += (1.0-final_pars[3])*orig_had1.Py()
		if had2_return!=None :
			newmetx += (1.0-final_pars[4])*orig_had2.Px()
			newmety += (1.0-final_pars[4])*orig_had2.Py()
		if had3_return!=None :
			newmetx += (1.0-final_pars[5])*orig_had3.Px()
			newmety += (1.0-final_pars[5])*orig_had3.Py()
		final_met.SetPx(newmetx); final_met.SetPy(newmety); final_met.SetPz(final_pars[0])
		final_met.SetE(final_met.Vect().Mag())
		#return everything
		return (bestfitindex,lep_return,final_met,lepb_return,bigjet_return,had1_return,had2_return,had3_return,bestfitchi2,final_pars)

	#type 1 (fully merged) top minimization function
	@staticmethod
	def fcnt1(npar, deriv, f, par, flag) :
		#Build rescaled versions of the vectors involved
		l = rescale(TTBarReconstructor.lep_global_vec,par[1])
		bl = rescale(TTBarReconstructor.blep_global_vec,par[2])
		hs1 = rescale(TTBarReconstructor.had1_global_vec,par[3])
		hs2 = rescale(TTBarReconstructor.had2_global_vec,par[4])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( TTBarReconstructor.met_global_vec.Px()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Px()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Px()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Px()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Px() )
		newmety = ( TTBarReconstructor.met_global_vec.Py()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Py()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Py()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Py()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Py() )
		v = rescale(TTBarReconstructor.met_global_vec,1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; th = TTBarReconstructor.bigjet_global_vec-TTBarReconstructor.had1_global_vec-TTBarReconstructor.had2_global_vec+hs1+hs2
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l1,TTBarReconstructor.WT_l1,th.M(),TTBarReconstructor.MT_h1,TTBarReconstructor.WT_h1,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
		f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
					+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) )

	#type 2 (partially merged) top minimization function
	@staticmethod
	def fcnt2(npar, deriv, f, par, flag) :
		#Build rescaled versions of the vectors involved
		l = rescale(TTBarReconstructor.lep_global_vec,par[1])
		bl = rescale(TTBarReconstructor.blep_global_vec,par[2])
		hs1 = rescale(TTBarReconstructor.had1_global_vec,par[3])
		hs2 = rescale(TTBarReconstructor.had2_global_vec,par[4])
		bh = rescale(TTBarReconstructor.had3_global_vec,par[5])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( TTBarReconstructor.met_global_vec.Px()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Px()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Px()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Px()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Px()+(1.0-par[5])*TTBarReconstructor.had3_global_vec.Px() )
		newmety = ( TTBarReconstructor.met_global_vec.Py()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Py()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Py()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Py()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Py()+(1.0-par[5])*TTBarReconstructor.had3_global_vec.Py() )
		v = rescale(TTBarReconstructor.met_global_vec,1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; 
		wh = TTBarReconstructor.bigjet_global_vec-TTBarReconstructor.had1_global_vec-TTBarReconstructor.had2_global_vec+hs1+hs2
		th = wh+bh
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l2,TTBarReconstructor.WT_l2,th.M(),TTBarReconstructor.MT_h2,TTBarReconstructor.WT_h2,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
		f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
					+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) 
					+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )

	#type 3 (boosted untagged) top minimization functions for various numbers of hadronic-side jets
	@staticmethod
	def fcnt33(npar, deriv, f, par, flag) :
		#Build rescaled versions of the vectors involved
		l = rescale(TTBarReconstructor.lep_global_vec,par[1])
		bl = rescale(TTBarReconstructor.blep_global_vec,par[2])
		h1 = rescale(TTBarReconstructor.had1_global_vec,par[3])
		h2 = rescale(TTBarReconstructor.had2_global_vec,par[4])
		h3 = rescale(TTBarReconstructor.had3_global_vec,par[5])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( TTBarReconstructor.met_global_vec.Px()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Px()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Px()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Px()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Px()+(1.0-par[5])*TTBarReconstructor.had3_global_vec.Px() )
		newmety = ( TTBarReconstructor.met_global_vec.Py()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Py()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Py()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Py()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Py()+(1.0-par[5])*TTBarReconstructor.had3_global_vec.Py() )
		v = rescale(TTBarReconstructor.met_global_vec,1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; 
		th = h1+h2+h3
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l3,TTBarReconstructor.WT_l3,th.M(),TTBarReconstructor.MT_h3,TTBarReconstructor.WT_h3,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
		f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
					+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) 
					+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )
	@staticmethod
	def fcnt32(npar, deriv, f, par, flag) :
		#Build rescaled versions of the vectors involved
		l = rescale(TTBarReconstructor.lep_global_vec,par[1])
		bl = rescale(TTBarReconstructor.blep_global_vec,par[2])
		h1 = rescale(TTBarReconstructor.had1_global_vec,par[3])
		h2 = rescale(TTBarReconstructor.had2_global_vec,par[4])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( TTBarReconstructor.met_global_vec.Px()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Px()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Px()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Px()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Px() )
		newmety = ( TTBarReconstructor.met_global_vec.Py()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Py()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Py()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Py()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Py() )
		v = rescale(TTBarReconstructor.met_global_vec,1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; 
		th = h1+h2
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l3,TTBarReconstructor.WT_l3,th.M(),TTBarReconstructor.MT_h3,TTBarReconstructor.WT_h3,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
		f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
					+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) 
					+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )
	@staticmethod
	def fcnt31(npar, deriv, f, par, flag) :
		#Build rescaled versions of the vectors involved
		l = rescale(TTBarReconstructor.lep_global_vec,par[1])
		bl = rescale(TTBarReconstructor.blep_global_vec,par[2])
		h1 = rescale(TTBarReconstructor.had1_global_vec,par[3])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( TTBarReconstructor.met_global_vec.Px()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Px()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Px()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Px() )
		newmety = ( TTBarReconstructor.met_global_vec.Py()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Py()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Py()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Py() )
		v = rescale(TTBarReconstructor.met_global_vec,1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; 
		th = h1
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l3,TTBarReconstructor.WT_l3,th.M(),TTBarReconstructor.MT_h3,TTBarReconstructor.WT_h3,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
		f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
					+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) 
					+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )

	#type 4 (resolved) kinematic fitting function
	@staticmethod
	def fcnt4(npar, deriv, f, par, flag) :
		#Build rescaled versions of the vectors involved
		l = rescale(TTBarReconstructor.lep_global_vec,par[1])
		bl = rescale(TTBarReconstructor.blep_global_vec,par[2])
		h1 = rescale(TTBarReconstructor.had1_global_vec,par[3])
		h2 = rescale(TTBarReconstructor.had2_global_vec,par[4])
		h3 = rescale(TTBarReconstructor.had3_global_vec,par[5])
		#rebuild the neutrino from the met post-rescaling
		newmetx = ( TTBarReconstructor.met_global_vec.Px()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Px()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Px()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Px()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Px()+(1.0-par[5])*TTBarReconstructor.had3_global_vec.Px() )
		newmety = ( TTBarReconstructor.met_global_vec.Py()+(1.0-par[1])*TTBarReconstructor.lep_global_vec.Py()+(1.0-par[2])*TTBarReconstructor.blep_global_vec.Py()
						+(1.0-par[3])*TTBarReconstructor.had1_global_vec.Py()+(1.0-par[4])*TTBarReconstructor.had2_global_vec.Py()+(1.0-par[5])*TTBarReconstructor.had3_global_vec.Py() )
		v = rescale(TTBarReconstructor.met_global_vec,1.0)
		v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
		v.SetE(v.Vect().Mag())
		wl = v + l; tl = wl + bl; 
		th = h1+h2+h3
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l4,TTBarReconstructor.WT_l4,th.M(),TTBarReconstructor.MT_h4,TTBarReconstructor.WT_h4,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
		f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
					+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) 
					+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )

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