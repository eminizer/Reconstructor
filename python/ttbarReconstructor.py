#Imports
from ROOT import TMinuit, Long, Double, TLorentzVector
from math import *
from array import array

#Global variables
MW = 80.4 #W mass
MT = 172.5 #top mass
MT_l1 = 171.2 #leptonic top mass for type 1 tops
MT_h1 = 190.9 #hadronic top mass for type 1 tops
MT_l2 = 171.2 #leptonic top mass for type 2 tops
MW_h2 = 80.4  #hadronic W mass for type 2 tops
MT_h2 = 190.9 #hadronic top mass for type 2 tops
MT_l3 = 171.2 #leptonic top mass for type 3 tops
MW_h3 = 80.4  #hadronic W mass for type 3 tops
MT_h3 = 190.9 #hadronic top mass for type 3 tops
WW = 2.0 #W width
WT = 1.4 #top width
WT_l1 = 29.3 #leptonic top width for type 1 tops
WT_h1 = 19.4 #hadronic top width for type 1 tops
WT_l2 = 29.3 #leptonic top width for type 2 tops
WW_h2 = 2.0  #hadronic W width for type 2 tops
WT_h2 = 19.4 #hadronic top width for type 2 tops
WT_l3 = 29.3 #leptonic top width for type 3 tops
WW_h3 = 2.0  #hadronic W width for type 3 tops
WT_h3 = 19.4 #hadronic top width for type 3 tops
SIGMAJ  = 0.10 #jet momentum resolution
SIGMAL  = 0.02 #lepton momentum resolution
#global fourvectors for the fit
lep_global_vec = TLorentzVector()
met_global_vec = TLorentzVector()
blep_global_vec = TLorentzVector()
hads1_global_vec = TLorentzVector()
hads2_global_vec = TLorentzVector()
hadW_global_vec = TLorentzVector()
hadb_global_vec = TLorentzVector()
hadt_global_vec = TLorentzVector()

#reconstruct takes in the list of hypotheses and returns a tuple of fit information invluded the corrected fourvectors and chi2 value
def reconstruct(hypotheses) :
	#first find the topology based on the number of parameters in the hypothesis list
	topology=len(hypotheses[0])-3
	#lists of final parameters and chi2 values
	bestParValues = []; finalChi2s = []; errflags = []
	for i in range(len(hypotheses)): 
		bestParValues.append([hypotheses[i][1].Pz(),1.0,1.0,1.0,1.0])
		if topology==2 or topology==3 : bestParValues[i].append(1.0)
		finalChi2s.append(1000000000.)
		errflags.append(-1)
	parNames = ['pZv','scalelep','scaleblep','scalehad1','scalehad2']
	parerrs = [0.0,0.0,0.0,0.0,0.0]
	if topology==2 or topology==3 : 
		parNames.append('scalehadb'); parerrs.append(0.0)
	#fit setup stuff common to every iteration
	minuit=None
	if topology==1 : minuit = TMinuit(5)
	elif topology==2 or topology==3 : minuit = TMinuit(6)
	if topology==1 : minuit.SetFCN(fcn1)
	elif topology==2 : minuit.SetFCN(fcn2)
	elif topology==3 : minuit.SetFCN(fcn3)
	#if minuit == None : print 'topology=%d'%(topology) #DEBUG
	ierflag = Long(1)
	arglist = array( 'd', [-1.0] )
	minuit.mnexcm('SET PRINT', arglist, 1,ierflag)
	minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
	arglist[0] = 100000.
	#perform the fit once for each jet/MET hypothesis
	#print '------------------------------------------------------------------------' #DEBUG
	for i in range(len(hypotheses)) :
		#print 'Running with hypothesis number = %d'%(i) #DEBUG
		hypothesis = hypotheses[i]
		#set constants
		lep_vec = hypothesis[0].getFourVector()
		met = hypothesis[1]
		lepb_vec = hypothesis[2].getFourVector()
		#set the global fourvector variables
		lep_global_vec.SetPtEtaPhiE(lep_vec.Pt(),lep_vec.Eta(),lep_vec.Phi(),lep_vec.E())
		met_global_vec.SetPtEtaPhiE(met.Pt(),met.Eta(),met.Phi(),met.E())
		blep_global_vec.SetPtEtaPhiE(lepb_vec.Pt(),lepb_vec.Eta(),lepb_vec.Phi(),lepb_vec.E())
		#the rest are topology-dependent
		if topology==1 :
			ht = hypothesis[3] #fully-merged, hadt is the single merged top
			htv = ht.getFourVector()
			hadt_global_vec.SetPtEtaPhiE(htv.Pt(),htv.Eta(),htv.Phi(),htv.E())
			h1v = ht.getSubjet(0).getFourVector() #had1 is the first subjet of the merged top
			hads1_global_vec.SetPtEtaPhiE(h1v.Pt(),h1v.Eta(),h1v.Phi(),h1v.E())
			h2v = ht.getSubjet(1).getFourVector() #had2 is the second subjet of the merged top
			hads2_global_vec.SetPtEtaPhiE(h2v.Pt(),h2v.Eta(),h2v.Phi(),h2v.E())
		elif topology==2 :
			hb = hypothesis[4] #partially-merged, hadb is the hadronic b-vector
			hbv = hb.getFourVector()
			hadb_global_vec.SetPtEtaPhiE(hbv.Pt(),hbv.Eta(),hbv.Phi(),hbv.E())
			hW = hypothesis[3] 
			hWv = hW.getFourVector()
			hadW_global_vec.SetPtEtaPhiE(hWv.Pt(),hWv.Eta(),hWv.Phi(),hWv.E())
			h1v = hW.getSubjet(0).getFourVector() #had1 is the first subjet of the merged W
			hads1_global_vec.SetPtEtaPhiE(h1v.Pt(),h1v.Eta(),h1v.Phi(),h1v.E())
			h2v = hW.getSubjet(1).getFourVector() #had2 is the second subjet of the merged W
			hads2_global_vec.SetPtEtaPhiE(h2v.Pt(),h2v.Eta(),h2v.Phi(),h2v.E())
		elif topology==3 :
			hb = hypothesis[5] #resolved, hadb is the hadronic b-vector
			hbv = hb.getFourVector()
			hadb_global_vec.SetPtEtaPhiE(hbv.Pt(),hbv.Eta(),hbv.Phi(),hbv.E())
			h1v = hypothesis[3].getFourVector() #had1 is the first hadronic W-jet
			hads1_global_vec.SetPtEtaPhiE(h1v.Pt(),h1v.Eta(),h1v.Phi(),h1v.E())
			h2v = hypothesis[4].getFourVector() #had2 is the second hadronic W-jet
			hads2_global_vec.SetPtEtaPhiE(h2v.Pt(),h2v.Eta(),h2v.Phi(),h2v.E())
		#add parameters to the fitter
		for j in range(len(bestParValues[i])) :
			minuit.mnparm(j,parNames[j],bestParValues[i][j],1.0,0,0,ierflag)
	#	print 'lep_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(lep_global_vec.Pt(),lep_global_vec.Eta(),lep_global_vec.Phi(),lep_global_vec.M()) #DEBUG
	#	print 'met_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(met_global_vec.Pt(),met_global_vec.Eta(),met_global_vec.Phi(),met_global_vec.M()) #DEBUG
	#	print 'blep_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(blep_global_vec.Pt(),blep_global_vec.Eta(),blep_global_vec.Phi(),blep_global_vec.M()) #DEBUG
	#	print 'hads1_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(hads1_global_vec.Pt(),hads1_global_vec.Eta(),hads1_global_vec.Phi(),hads1_global_vec.M()) #DEBUG
	#	print 'hads2_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(hads2_global_vec.Pt(),hads2_global_vec.Eta(),hads2_global_vec.Phi(),hads2_global_vec.M()) #DEBUG
	#	print 'hadW_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(hadW_global_vec.Pt(),hadW_global_vec.Eta(),hadW_global_vec.Phi(),hadW_global_vec.M()) #DEBUG
	#	print 'hadb_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(hadb_global_vec.Pt(),hadb_global_vec.Eta(),hadb_global_vec.Phi(),hadb_global_vec.M()) #DEBUG
	#	print 'hadt_global_vec = (pt,eta,phi,M) = (%.2f,%.2f,%.2f,%.2f)'%(hadt_global_vec.Pt(),hadt_global_vec.Eta(),hadt_global_vec.Phi(),hadt_global_vec.M()) #DEBUG
		#minimize
		minuit.mnexcm('MIGRAD', arglist, 1,ierflag)
	#	print 'errflag = '+str(ierflag) #DEBUG
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
			print 'This fit never had positive likelihood' #DEBUG
			errflags[i] = -1
		#Set fit Chi2 for this pZ solution
		tmp1 = Double(1.0); tmp2 = Double(1.0); tmp3 = Double(1.0)
		minuit.mnstat(tmp1,tmp2,tmp3,Long(1),Long(1),Long(1))
		#print 'chi2 = %.4f'%(tmp1) #DEBUG
		finalChi2s[i] = tmp1
	#if no fits converged, return garbage
	noneconverged = True
	for errflag in errflags :
		if errflag==0 :
			noneconverged = False; break
	if noneconverged :
		return -1, None, None, None, None, None, None, None, None, 1000000000., []
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
	hadWs1_return=None; hadWs2_return=None; hadW_return=None; hadb_return=None; hadt_return=None
	orig_lep = besthypothesis[0].getFourVector()
	lep_return = rescale(orig_lep,final_pars[1])
	orig_lepb = besthypothesis[2].getFourVector()
	lepb_return = rescale(orig_lepb,final_pars[2])
	if topology==1 :
		orig_hadt = besthypothesis[3].getFourVector()
		orig_subj1 = besthypothesis[3].getSubjet(0).getFourVector()
		newsubj1 = rescale(orig_subj1,final_pars[3])
		orig_subj2 = besthypothesis[3].getSubjet(1).getFourVector()
		newsubj2 = rescale(orig_subj2,final_pars[4])
		hadt_return = orig_hadt-orig_subj1+newsubj1-orig_subj2+newsubj2
	elif topology==2 :
		orig_hadW = besthypothesis[3].getFourVector()
		orig_subj1 = besthypothesis[3].getSubjet(0).getFourVector()
		newsubj1 = rescale(orig_subj1,final_pars[3])
		orig_subj2 = besthypothesis[3].getSubjet(1).getFourVector()
		newsubj2 = rescale(orig_subj2,final_pars[4])
		hadW_return = orig_hadW-orig_subj1+newsubj1-orig_subj2+newsubj2
		orig_hadb = besthypothesis[4].getFourVector()
		hadb_return = rescale(orig_hadb,final_pars[5])
		hadt_return = hadW_return+hadb_return
	elif topology==3 :
		orig_hadWs1 = besthypothesis[3].getFourVector()
		hadWs1_return = rescale(orig_hadWs1,final_pars[3])
		orig_hadWs2 = besthypothesis[4].getFourVector()
		hadWs2_return = rescale(orig_hadWs2,final_pars[4])
		orig_hadb = besthypothesis[5].getFourVector()
		hadb_return = rescale(orig_hadb,final_pars[5])
		hadW_return = hadWs1_return+hadWs2_return
		hadt_return = hadW_return+hadb_return
	#rebuild the neutrino post-rescaling
	newmetx = final_met.Px()+ (1.0-final_pars[1])*orig_lep.Px()
	newmety = final_met.Py()+ (1.0-final_pars[1])*orig_lep.Py()
	newmetx += (1.0-final_pars[2])*orig_lepb.Px()
	newmety += (1.0-final_pars[2])*orig_lepb.Py()
	if topology==1 or topology==2 :
		newmetx += (1.0-final_pars[3])*orig_subj1.Px()
		newmety += (1.0-final_pars[3])*orig_subj1.Py()
		newmetx += (1.0-final_pars[4])*orig_subj2.Px()
		newmety += (1.0-final_pars[4])*orig_subj2.Py()
	elif topology==3 :
		newmetx += (1.0-final_pars[3])*orig_hadWs1.Px()
		newmety += (1.0-final_pars[3])*orig_hadWs1.Py()
		newmetx += (1.0-final_pars[4])*orig_hadWs2.Px()
		newmety += (1.0-final_pars[4])*orig_hadWs2.Py()
	if topology==2 or topology==3 :
		newmetx += (1.0-final_pars[5])*orig_hadb.Px()
		newmety += (1.0-final_pars[5])*orig_hadb.Py()
	final_met.SetPx(newmetx); final_met.SetPy(newmety); final_met.SetPz(final_pars[0])
	final_met.SetE(final_met.Vect().Mag())
	#return everything
	return (bestfitindex,lep_return,final_met,lepb_return,hadWs1_return,hadWs2_return,hadW_return,hadb_return,hadt_return,bestfitchi2,final_pars)

#type 1 (fully merged) top minimization function
def fcn1(npar, deriv, f, par, flag) :
	#Build rescaled versions of the vectors involved
	l = rescale(lep_global_vec,par[1])
	bl = rescale(blep_global_vec,par[2])
	hs1 = rescale(hads1_global_vec,par[3])
	hs2 = rescale(hads2_global_vec,par[4])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_global_vec.Px()+(1.0-par[1])*lep_global_vec.Px()+(1.0-par[2])*blep_global_vec.Px()
					+(1.0-par[3])*hads1_global_vec.Px()+(1.0-par[4])*hads2_global_vec.Px() )
	newmety = ( met_global_vec.Py()+(1.0-par[1])*lep_global_vec.Py()+(1.0-par[2])*blep_global_vec.Py()
					+(1.0-par[3])*hads1_global_vec.Py()+(1.0-par[4])*hads2_global_vec.Py() )
	v = rescale(met_global_vec,1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl; th = hadt_global_vec-hads1_global_vec-hads2_global_vec+hs1+hs2
	mwl2 = wl.M2(); mtl = tl.M(); mth = th.M()
	pdftl = (mtl-MT_l1)**2/WT_l1**2
	pdfth = (mth-MT_h1)**2/WT_h1**2
	lwl = ((MT**2-mwl2)*(MT**2-mwl2)*(2*MT**2+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW)
	pdf = pdftl+pdfth
	if lwl > 0.0 :
		pdf += -2.0*log(lwl)    #need positive f
	else :
		#print 'WARNING -- pdf is negative!!!'
		pdf += -2.0*log(1.e-50)
	f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
				+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) )
	#print 'mtl = %.3f, mth = %.3f, vpz = %.3f, pdftl = %.3f, pdfth = %.3f, lwl = %.3f, pdf = %.3f, f[0] = %.3f'%(mtl,mth,par[0],pdftl,pdfth,lwl,pdf,f[0]) #DEBUG

#type 2 (partially merged) top minimization function
def fcn2(npar, deriv, f, par, flag) :
	#Build rescaled versions of the vectors involved
	l = rescale(lep_global_vec,par[1])
	bl = rescale(blep_global_vec,par[2])
	hs1 = rescale(hads1_global_vec,par[3])
	hs2 = rescale(hads2_global_vec,par[4])
	bh = rescale(hadb_global_vec,par[5])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_global_vec.Px()+(1.0-par[1])*lep_global_vec.Px()+(1.0-par[2])*blep_global_vec.Px()
					+(1.0-par[3])*hads1_global_vec.Px()+(1.0-par[4])*hads2_global_vec.Px()+(1.0-par[5])*hadb_global_vec.Px() )
	newmety = ( met_global_vec.Py()+(1.0-par[1])*lep_global_vec.Py()+(1.0-par[2])*blep_global_vec.Py()
					+(1.0-par[3])*hads1_global_vec.Py()+(1.0-par[4])*hads2_global_vec.Py()+(1.0-par[5])*hadb_global_vec.Py() )
	v = rescale(met_global_vec,1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl; 
	wh = hadW_global_vec-hads1_global_vec-hads2_global_vec+hs1+hs2
	th = wh+bh
	mwl2 = wl.M2(); mtl = tl.M(); mwh = wh.M(); mth = th.M()
	pdftl = (mtl-MT_l2)**2/WT_l2**2
	pdfwh = (mwh-MW_h2)**2/WW_h2**2
	pdfth = (mth-MT_h2)**2/WT_h2**2
	lwl = ((MT**2-mwl2)*(MT**2-mwl2)*(2*MT**2+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW)
	pdf = pdftl+pdfwh+pdfth
	if lwl > 0.0 :
		pdf += -2.0*log(lwl)    #need positive f
	else :
		#print 'WARNING -- pdf is negative!!!'
		pdf += -2.0*log(1.e-50)
	f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
				+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) 
				+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )
	#print '	mtl = %.3f, mwh = %.3f, mth = %.3f, vpz = %.3f, pdftl = %.3f, pdfwh = %.3f, pdfth = %.3f, lwl = %.3f, pdf = %.3f, f[0] = %.3f'%(mtl,mwh,mth,par[0],pdftl,pdfwh,pdfth,lwl,pdf,f[0]) #DEBUG

#type 3 (resolved) top minimization function
def fcn3(npar, deriv, f, par, flag) :
	#Build rescaled versions of the vectors involved
	l = rescale(lep_global_vec,par[1])
	bl = rescale(blep_global_vec,par[2])
	hWs1 = rescale(hads1_global_vec,par[3])
	hWs2 = rescale(hads2_global_vec,par[4])
	bh = rescale(hadb_global_vec,par[5])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_global_vec.Px()+(1.0-par[1])*lep_global_vec.Px()+(1.0-par[2])*blep_global_vec.Px()
					+(1.0-par[3])*hads1_global_vec.Px()+(1.0-par[4])*hads2_global_vec.Px()+(1.0-par[5])*hadb_global_vec.Px() )
	newmety = ( met_global_vec.Py()+(1.0-par[1])*lep_global_vec.Py()+(1.0-par[2])*blep_global_vec.Py()
					+(1.0-par[3])*hads1_global_vec.Py()+(1.0-par[4])*hads2_global_vec.Py()+(1.0-par[5])*hadb_global_vec.Py() )
	v = rescale(met_global_vec,1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl; 
	wh = hWs1+hWs2
	th = wh+bh
	mwl2 = wl.M2(); mtl = tl.M(); mwh = wh.M(); mth = th.M()
	pdftl = (mtl-MT_l3)**2/WT_l3**2
	pdfwh = (mwh-MW_h3)**2/WW_h3**2
	pdfth = (mth-MT_h3)**2/WT_h3**2
	lwl = ((MT**2-mwl2)*(MT**2-mwl2)*(2*MT**2+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW)
	pdf = pdftl+pdfwh+pdfth
	if lwl > 0.0 :
		pdf += -2.0*log(lwl)    #need positive f
	else :
		#print 'WARNING -- pdf is negative!!!'
		pdf += -2.0*log(1.e-50)
	f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
				+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) 
				+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )
	#print 'mtl = %.3f, mth = %.3f, vpz = %.3f, pdftl = %.3f, pdfth = %.3f, lwl = %.3f, pdf = %.3f, f[0] = %.3f'%(mtl,mth,par[0],pdftl,pdfth,lwl,pdf,f[0]) #DEBUG

#fourvector rescaling function
def rescale(vec,fac) :
	p2 = fac*fac*vec.Vect().Mag2()
	m2 = vec.M()*vec.M()
	newE = sqrt(p2+m2)
	#note that the function returns a NEW TLorentzVector, meaning the original is unaltered
	return TLorentzVector(fac*vec.Px(),fac*vec.Py(),fac*vec.Pz(),newE)