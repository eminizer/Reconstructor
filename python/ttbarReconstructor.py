#reconstruct takes in the list of hypotheses and returns a tuple of best fit information invluded the corrected fourvectors and chi2 value
def reconstruct(self,hypotheses,topology) :
	#set the number of fit parameters based on topology
	fitpars = 6
	if topology==1 :
		self.minuit.mnfixp(5,self.ierflag)
		fitpars = 5
	#list of best parameter values
	bestParValues = []
	for i in range(fitpars) :
		bestParValues.append(self.parinivals[i])
	#list of flags for each fit
	errflags = []
	for i in range(len(hypotheses)): 
		errflags.append(-1)
	#perform the fit once for each jet/MET hypothesis
	#print '------------------------------------------------------------------------' #DEBUG
	bestfitchi2 = 10000000.; bestfitindex = -1; 
	for i in range(len(hypotheses)) :
		#print 'Running with hypothesis number = %d'%(i) #DEBUG
		hypothesis = hypotheses[i]
		#set parameter values
		for j in range(fitpars) :
			self.minuit.mnparm(j,self.parNames[j],self.parinivals[j],1.0,0,0,self.ierflag)
		#set constants
		lep_vec = hypothesis[0].getFourVector()
		met = hypothesis[1]
		lepb_vec = hypothesis[2].getFourVector()
		#set the global fourvector variables
		self.lep_global_vec.SetPtEtaPhiE(lep_vec.Pt(),lep_vec.Eta(),lep_vec.Phi(),lep_vec.E())
		self.met_global_vec.SetPtEtaPhiE(met.Pt(),met.Eta(),met.Phi(),met.E())
		self.blep_global_vec.SetPtEtaPhiE(lepb_vec.Pt(),lepb_vec.Eta(),lepb_vec.Phi(),lepb_vec.E())
		#the rest are topology-dependent
		if topology==1 :
			hypak8 = hypothesis[3][0] 
			hypak8v = hypak8.getFourVector()
			self.hadt_global_vec.SetPtEtaPhiE(hypak8v.Pt(),hypak8v.Eta(),hypak8v.Phi(),hypak8v.E())
			h1v = hypak8.getSubjet(0).getFourVector() #had1 is the first subjet of the merged W
			self.had1_global_vec.SetPtEtaPhiE(h1v.Pt(),h1v.Eta(),h1v.Phi(),h1v.E())
			h2v = hypak8.getSubjet(1).getFourVector() #had2 is the second subjet of the merged W
			self.had2_global_vec.SetPtEtaPhiE(h2v.Pt(),h2v.Eta(),h2v.Phi(),h2v.E())
			self.minuit.SetFCN(TTBarReconstructor.fcnt1)
		else :
			h1v = hypothesis[3][0].getFourVector()
			self.had1_global_vec.SetPtEtaPhiE(h1v.Pt(),h1v.Eta(),h1v.Phi(),h1v.E())
			h2v = hypothesis[3][1].getFourVector()
			self.had2_global_vec.SetPtEtaPhiE(h2v.Pt(),h2v.Eta(),h2v.Phi(),h2v.E())
			h3v = hypothesis[3][1].getFourVector()
			self.had3_global_vec.SetPtEtaPhiE(h3v.Pt(),h3v.Eta(),h3v.Phi(),h3v.E())
			if topology==2 : self.minuit.SetFCN(TTBarReconstructor.fcnt2)
			elif topology==3 : self.minuit.SetFCN(TTBarReconstructor.fcnt3)
		#print 'self.lep_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.lep_global_vec.Pt(),self.lep_global_vec.Eta(),self.lep_global_vec.Phi(),self.lep_global_vec.M()) #DEBUG
		#print 'self.met_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.met_global_vec.Pt(),self.met_global_vec.Eta(),self.met_global_vec.Phi(),self.met_global_vec.M()) #DEBUG
		#print 'self.blep_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.blep_global_vec.Pt(),self.blep_global_vec.Eta(),self.blep_global_vec.Phi(),self.blep_global_vec.M()) #DEBUG
		#print 'self.had1_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.had1_global_vec.Pt(),self.had1_global_vec.Eta(),self.had1_global_vec.Phi(),self.had1_global_vec.M()) #DEBUG
		#print 'self.had2_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.had2_global_vec.Pt(),self.had2_global_vec.Eta(),self.had2_global_vec.Phi(),self.had2_global_vec.M()) #DEBUG
		#print 'self.had3_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.had3_global_vec.Pt(),self.had3_global_vec.Eta(),self.had3_global_vec.Phi(),self.had3_global_vec.M()) #DEBUG
		#print 'self.hadt_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(self.hadt_global_vec.Pt(),self.hadt_global_vec.Eta(),self.hadt_global_vec.Phi(),self.hadt_global_vec.M()) #DEBUG
		#minimize
		self.minuit.mnexcm('MIGRAD', self.arglist, 1, self.ierflag)
		#print 'errflag = '+str(ierflag) #DEBUG
		errflags[i] = self.ierflag
		#Check fit Chi2 for this hypothesis
		tmp1 = Double(1.0); tmp2 = Double(1.0); tmp3 = Double(1.0)
		self.minuit.mnstat(tmp1,tmp2,tmp3,Long(1),Long(1),Long(1))
		#print 'chi2 = %.4f'%(tmp1) #DEBUG
		if tmp1<bestfitchi2 :
			bestfitchi2=tmp1
			bestfitindex=i
			#Get the bestfit parameters back from minuit
			for j in range(fitpars) :
				tmp = Double(1.0)
				self.minuit.GetParameter(j,tmp,Double(self.parerrs[j]))
				bestParValues[j] = tmp
	#undo fixing parameters if necessary
	if fitpars==5 :
		self.minuit.mnfree(0)
	#if no fits converged, return garbage
	noneconverged = True
	for errflag in errflags :
		if errflag==0 :
			noneconverged = False; break
	if noneconverged :
		return -1, None, None, None, None, None, None, None, 1000000000., []
	#find which fit gave the best results and record best parameter values
	besthypothesis = hypotheses[bestfitindex]
	final_met = TLorentzVector(); final_met.SetPtEtaPhiE(besthypothesis[1].Pt(),besthypothesis[1].Eta(),besthypothesis[1].Phi(),besthypothesis[1].E())
	#print 'bestfitindex = %d, bestfitchi2 = %.4f'%(bestfitindex,bestfitchi2) #DEBUG
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
