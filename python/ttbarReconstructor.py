#Imports
from ROOT import TMinuit, Long, Double, TLorentzVector
from math import *
from array import array

#Global variables
MW = 80.4 #W mass
MW_ht2 = 80.4 #hadronic W mass for type 2 tops
MT = 172.5 #top mass
MT_lt1 = 174.9 #leptonic top mass for type 1 tops
MT_lt2 = 173.1 #leptonic top mass for type 2 tops
MT_ht1 = 187.4 #hadronic top mass for type 1 tops
MT_ht2 = 173.8 #hadronic top mass for type 2 tops
WW = 2.0 #W width
WW_ht2 = 2.0 #hadronic W width for type 2 tops
WT = 1.4 #top width
WT_lt1 = 33.17 #leptonic top width for type 1 tops
WT_lt2 = 26.64 #leptonic top width for type 2 tops
WT_ht1 = 19.54 #hadronic top width for type 1 tops
WT_ht2 = 16.39 #hadronic top width for type 2 tops
SIGMAJ  = 0.10 #jet momentum resolution
SIGMAL  = 0.03 #lepton momentum resolution
#global fourvectors for the fit
lep_global_vec = TLorentzVector()
met_global_vec = TLorentzVector()
blep_global_vec = TLorentzVector()
had1_global_vec = TLorentzVector()
had2_global_vec = TLorentzVector()
had3_global_vec = TLorentzVector()
hadt_global_vec = TLorentzVector()
hadW_global_vec = TLorentzVector()

#reconstruct takes in the list of hypotheses and returns a tuple of fit information invluded the corrected fourvectors and chi2 value
def reconstruct(hypotheses) :
	#lists of final parameters and chi2 values
	bestParValues = []; finalChi2s = []; errflags = []
	for i in range(len(hypotheses)): 
		if len(hypotheses[i])==4 :
			bestParValues.append([hypotheses[i][1].Pz(),1.0,1.0,1.0,1.0])
		elif len(hypotheses[i])==5 :
			bestParValues.append([hypotheses[i][1].Pz(),1.0,1.0,1.0,1.0,1.0])
		finalChi2s.append(1000000000.)
		errflags.append(-1)
	parNames = ['pZv','scalelep','scaleblep','scalehad1','scalehad2']
	parerrs = [0.0,0.0,0.0,0.0,0.0]
	if len(hypotheses[0])==5 :
		parNames.append('scalehad3'); parerrs.append(0.0)
	#fit setup stuff common to every iteration
	if len(hypotheses[0])==4 :
		minuit = TMinuit(5)
		minuit.SetFCN(fcnt1)
	elif len(hypotheses[0])==5 :
		minuit = TMinuit(6)
		minuit.SetFCN(fcnt2)
	ierflag = Long(1)
	arglist = array( 'd', [-1.0] )
	minuit.mnexcm('SET PRINT', arglist, 1,ierflag)
	minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
	arglist[0] = 100000.
	#perform the fit once for each jet/MET hypothesis
	#print '------------------------------------------------------------------------' #DEBUG
	for i in range(len(hypotheses)) :
	#	print 'Running with hypothesis number = %d'%(i) #DEBUG
		hypothesis = hypotheses[i]
		#set the global fourvector variables
		lep_vec = hypothesis[0].getFourVector()
		lep_global_vec.SetPtEtaPhiE(lep_vec.Pt(),lep_vec.Eta(),lep_vec.Phi(),lep_vec.E())
		met = hypothesis[1]
		met_global_vec.SetPtEtaPhiE(met.Pt(),met.Eta(),met.Phi(),met.E())
		lepb_vec = hypothesis[2].getFourVector()
		blep_global_vec.SetPtEtaPhiE(lepb_vec.Pt(),lepb_vec.Eta(),lepb_vec.Phi(),lepb_vec.E())
		h1v = hypothesis[3].getSubjet(0).getFourVector()
		h2v = hypothesis[3].getSubjet(1).getFourVector()
		had1_global_vec.SetPtEtaPhiE(h1v.Pt(),h1v.Eta(),h1v.Phi(),h1v.E())
		had2_global_vec.SetPtEtaPhiE(h2v.Pt(),h2v.Eta(),h2v.Phi(),h2v.E())
		if len(hypothesis)==4 :
			htv = hypothesis[3].getFourVector()
			hadt_global_vec.SetPtEtaPhiE(htv.Pt(),htv.Eta(),htv.Phi(),htv.E())
		elif len(hypothesis)==5 :
			hWv = hypothesis[3].getFourVector()
			hadW_global_vec.SetPtEtaPhiE(hWv.Pt(),hWv.Eta(),hWv.Phi(),hWv.E())
			h3v = hypothesis[4].getFourVector()
			had3_global_vec.SetPtEtaPhiE(h3v.Pt(),h3v.Eta(),h3v.Phi(),h3v.E())
		#add parameters to the fitter
		for j in range(len(bestParValues[0])) :
			minuit.mnparm(j,parNames[j],bestParValues[i][j],1.0,0,0,ierflag)
	#	print 'lep_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(lep_global_vec.Pt(),lep_global_vec.Eta(),lep_global_vec.Phi(),lep_global_vec.M()) #DEBUG
	#	print 'met_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(met_global_vec.Pt(),met_global_vec.Eta(),met_global_vec.Phi(),met_global_vec.M()) #DEBUG
	#	print 'blep_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(blep_global_vec.Pt(),blep_global_vec.Eta(),blep_global_vec.Phi(),blep_global_vec.M()) #DEBUG
	#	print 'had1_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(had1_global_vec.Pt(),had1_global_vec.Eta(),had1_global_vec.Phi(),had1_global_vec.M()) #DEBUG
	#	print 'had2_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(had2_global_vec.Pt(),had2_global_vec.Eta(),had2_global_vec.Phi(),had2_global_vec.M()) #DEBUG
	#	print 'had3_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(had3_global_vec.Pt(),had3_global_vec.Eta(),had3_global_vec.Phi(),had3_global_vec.M()) #DEBUG
		#minimize
		minuit.mnexcm('MIGRAD', arglist, 1,ierflag)
		errflags[i] = ierflag
		#Get the bestfit parameters back from minuit
		for j in range(len(bestParValues[0])) :
			tmp = Double(1.0)
			minuit.GetParameter(j,tmp,Double(parerrs[j]))
			bestParValues[i][j] = tmp
		#check to make sure it didn't just get stuck in the "negative LL" trap
		gotstuck = True
		for j in range(1,len(bestParValues[i])) :
			if bestParValues[i][j]!=1.0 : 
				gotstuck = False; break
		if gotstuck :
	#		print 'This fit never had positive likelihood' #DEBUG
			errflags[i] = -1
		#Set fit Chi2 for this pZ solution
		tmp1 = Double(1.0); tmp2 = Double(1.0); tmp3 = Double(1.0)
		minuit.mnstat(tmp1,tmp2,tmp3,Long(1),Long(1),Long(1))
		finalChi2s[i] = tmp1
	#if no fits converged, return garbage
	noneconverged = True
	for errflag in errflags :
		if errflag==0 :
			noneconverged = False; break
	if noneconverged :
		return -1, None, None, None, None, None, None, 1000000000., []
	#find which fit gave the best results and record best parameter values
	final_pars = []
	bestfitchi2 = min(finalChi2s)
	bestfitindex = finalChi2s.index(bestfitchi2)
	for i in range(len(bestParValues[0])) :
		final_pars.append(bestParValues[bestfitindex][i])
	besthypothesis = hypotheses[bestfitindex]
	final_met = TLorentzVector(); final_met.SetPtEtaPhiE(besthypothesis[1].Pt(),besthypothesis[1].Eta(),besthypothesis[1].Phi(),besthypothesis[1].E())
	#print 'bestfitindex = %d, bestfitchi2 = %.4f'%(bestfitindex,bestfitchi2) #DEBUG
	#rescale the lepton and jet four vectors based on the final parameters
	orig_lep = besthypothesis[0].getFourVector()
	lep_return = rescale(orig_lep,final_pars[1])
	orig_lepb = besthypothesis[2].getFourVector()
	lepb_return = rescale(orig_lepb,final_pars[2])
	orig_subj1 = besthypothesis[3].getSubjet(0).getFourVector()
	orig_subj2 = besthypothesis[3].getSubjet(1).getFourVector()
	newsubj1 = rescale(orig_subj1,final_pars[3])
	newsubj2 = rescale(orig_subj2,final_pars[4])
	if len(besthypothesis)==4 :
		hadw_return = None
		hadb_return = None
		hadt_return = besthypothesis[3].getFourVector()-orig_subj1-orig_subj2+newsubj1+newsubj2
	elif len(besthypothesis)==5 :
		hadw_return = besthypothesis[3].getFourVector()-orig_subj1-orig_subj2+newsubj1+newsubj2
		hadb_return = rescale(besthypothesis[4].getFourVector(),final_pars[5])
		hadt_return = hadw_return+hadb_return
	#rebuild the neutrino post-rescaling
	newmetx = final_met.Px()+ (1.0-final_pars[1])*orig_lep.Px()
	newmety = final_met.Py()+ (1.0-final_pars[1])*orig_lep.Py()
	newmetx += (1.0-final_pars[2])*orig_lepb.Px()
	newmety += (1.0-final_pars[2])*orig_lepb.Py()
	newmetx += (1.0-final_pars[3])*orig_subj1.Px()
	newmety += (1.0-final_pars[3])*orig_subj1.Py()
	newmetx += (1.0-final_pars[4])*orig_subj2.Px()
	newmety += (1.0-final_pars[4])*orig_subj2.Py()
	if len(besthypothesis)==5 :
		newmetx += (1.0-final_pars[5])*besthypothesis[4].getFourVector().Px()
		newmety += (1.0-final_pars[5])*besthypothesis[4].getFourVector().Py()
	final_met.SetPx(newmetx); final_met.SetPy(newmety); final_met.SetPz(final_pars[0])
	final_met.SetE(final_met.Vect().Mag())
	#return everything
	return (bestfitindex,lep_return,final_met,lepb_return,hadw_return,hadb_return,hadt_return,bestfitchi2,final_pars)

#type 1 top fitting function
def fcnt1(npar, deriv, f, par, flag) :
	#Build rescaled versions of the vectors involved
	l = rescale(lep_global_vec,par[1])
	bl = rescale(blep_global_vec,par[2])
	hs1 = rescale(had1_global_vec,par[3])
	hs2 = rescale(had2_global_vec,par[4])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_global_vec.Px()+(1.0-par[1])*lep_global_vec.Px()+(1.0-par[2])*blep_global_vec.Px()
					+(1.0-par[3])*had1_global_vec.Px()+(1.0-par[4])*had2_global_vec.Px() )
	newmety = ( met_global_vec.Py()+(1.0-par[1])*lep_global_vec.Py()+(1.0-par[2])*blep_global_vec.Py()
					+(1.0-par[3])*had1_global_vec.Py()+(1.0-par[4])*had2_global_vec.Py() )
	v = rescale(met_global_vec,1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl; th = hadt_global_vec-had1_global_vec-had2_global_vec+hs1+hs2
	mwl2 = wl.M2(); mtl2 = tl.M2(); mth2 = th.M2()
	pdftl = 1./((mtl2-MT_lt1*MT_lt1)**2+MT_lt1*MT_lt1*WT_lt1*WT_lt1)
	pdfth = 1./((mth2-MT_ht1*MT_ht1)**2+MT_ht1*MT_ht1*WT_ht1*WT_ht1)
	pdfwl  = ((MT_lt1*MT_lt1-mwl2)*(MT_lt1*MT_lt1-mwl2)*(2*MT_lt1*MT_lt1+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW)
	pdf = pdfwl*pdftl*pdfth #vanilla
	lnL = 0.
	if pdf > 0.0 :
		lnL += log(pdf)    #need positive f
	else :
		#print 'WARNING -- pdf is negative!!!'
		pdf = 1.e-50
		lnL += log(pdf)
	f[0] = ( -2.0*lnL+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
				+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) )
	#and we don't need to return anything because minuit

#type 2 top fitting function
def fcnt2(npar, deriv, f, par, flag) :
	#Build rescaled versions of the vectors involved
	l = rescale(lep_global_vec,par[1])
	bl = rescale(blep_global_vec,par[2])
	hs1 = rescale(had1_global_vec,par[3])
	hs2 = rescale(had2_global_vec,par[4])
	hs3 = rescale(had3_global_vec,par[5])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_global_vec.Px()+(1.0-par[1])*lep_global_vec.Px()+(1.0-par[2])*blep_global_vec.Px()
					+(1.0-par[3])*had1_global_vec.Px()+(1.0-par[4])*had2_global_vec.Px()+(1.0-par[5])*had3_global_vec.Px() )
	newmety = ( met_global_vec.Py()+(1.0-par[1])*lep_global_vec.Py()+(1.0-par[2])*blep_global_vec.Py()
					+(1.0-par[3])*had1_global_vec.Py()+(1.0-par[4])*had2_global_vec.Py()+(1.0-par[5])*had3_global_vec.Py() )
	v = rescale(met_global_vec,1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl; wh = hadW_global_vec-had1_global_vec-had2_global_vec+hs1+hs2; th = wh+hs3
	mwl2 = wl.M2(); mtl2 = tl.M2(); mwh2 = wh.M2(); mth2 = th.M2()
	pdftl = 1./((mtl2-MT_lt2*MT_lt2)**2+MT_lt2*MT_lt2*WT_lt2*WT_lt2)
	pdfth = 1./((mth2-MT_ht2*MT_ht2)**2+MT_ht2*MT_ht2*WT_ht2*WT_ht2)
	pdfwl  = ((MT_lt2*MT_lt2-mwl2)*(MT_lt2*MT_lt2-mwl2)*(2*MT_lt2*MT_lt2+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW)
	pdfwh  = ((MT_ht2*MT_ht2-mwh2)*(MT_ht2*MT_ht2-mwh2)*(2*MT_ht2*MT_ht2+mwh2))/((mwh2-MW_ht2*MW_ht2)*(mwh2-MW_ht2*MW_ht2)+MW_ht2*MW_ht2*WW_ht2*WW_ht2)
	pdf = pdfwl*pdftl*pdfwh*pdfth #vanilla
	lnL = 0.
	if pdf > 0.0 :
		lnL += log(pdf)    #need positive f
	else :
		#print 'WARNING -- pdf is negative!!!'
		pdf = 1.e-50
		lnL += log(pdf)
	f[0] = ( -2.0*lnL+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
				+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ)+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )
	#and we don't need to return anything because minuit

#fourvector rescaling function
def rescale(vec,fac) :
	p2 = fac*fac*vec.Vect().Mag2()
	m2 = vec.M()*vec.M()
	newE = sqrt(p2+m2)
	#note that the function returns a NEW TLorentzVector, meaning the original is unaltered
	return TLorentzVector(fac*vec.Px(),fac*vec.Py(),fac*vec.Pz(),newE)