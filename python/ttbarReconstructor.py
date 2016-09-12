#Imports
from ROOT import TMinuit, Long, Double, TLorentzVector
from math import *
from array import array

#Global variables
MW = 80.4 #W mass
MT = 172.5 #top mass
WW = 2.0 #W width
WT = 1.4 #top width
SIGMAJ  = 0.10 #jet momentum resolution
SIGMAL  = 0.03 #lepton momentum resolution
#global fourvectors for the fit
lep_global_vec = TLorentzVector()
met_global_vec = TLorentzVector()
blep_global_vec = TLorentzVector()
thad_global_vec = TLorentzVector()

#reconstruct takes in the lepton, 2 met, lepb, and hadt fourvectors
#and returns a 5-tuple of:
#	1) the corrected lepton fourvector
#	2) the corrected (and selected) met fourvector
#	3) the corrected leptonic b jet fourvector
#	4) the corrected hadronic top jet fourvector
#	5) the final Chi2 value from the kinematic fit
def reconstruct(lepton,met1,met2,lepb,hadt) :
	mets = [met1,met2]
	#lists of final parameters and chi2 values
	bestParValues = [[met1.Pz(),1.0,1.0,1.0],[met2.Pz(),1.0,1.0,1.0]]
	parNames = ['pZv','scalelep','scaleblep','scalethad']
	parerrs = [0.0,0.0,0.0,0.0]
	finalChi2s = [1000000000.,1000000000.]
	errflags = [-1,-1]
	#see whether we have to fit twice based on how many solutions for the neutrino Pz we have
	nFits = 2
	if met1.Pz() == met2.Pz() :
		nFits = 1
	#fit setup stuff common to both iterations
	minuit = TMinuit(4)
	minuit.SetFCN(fcn)
	ierflag = Long(1)
	arglist = array( 'd', [-1.0] )
	minuit.mnexcm('SET PRINT', arglist, 1,ierflag)
	minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
	arglist[0] = 100000.
	#perform the fit once for each unique neutrino solution
	for iFit in range(nFits) :
		#Set which met solution we're looking at
		met = mets[iFit]
		#set the parameters in minuit
		for i in range(len(bestParValues[0])) :
			minuit.mnparm(i,parNames[i],bestParValues[iFit][i],1.0,0,0,ierflag)
		#set the global fourvector variables
		lep_global_vec.SetPtEtaPhiM(lepton.Pt(),lepton.Eta(),lepton.Phi(),lepton.M())
		met_global_vec.SetPtEtaPhiM(met.Pt(),met.Eta(),met.Phi(),met.M())
		blep_global_vec.SetPtEtaPhiM(lepb.Pt(),lepb.Eta(),lepb.Phi(),lepb.M())
		thad_global_vec.SetPtEtaPhiM(hadt.Pt(),hadt.Eta(),hadt.Phi(),hadt.M())
		#minimize
		minuit.mnexcm('MIGRAD', arglist, 1,ierflag)
		errflags[iFit] = ierflag
		#Get the best parameters back from minuit
		for i in range(len(bestParValues[0])) :
			tmp = Double(1.0)
			minuit.GetParameter(i,tmp,Double(parerrs[i]))
			bestParValues[iFit][i] = tmp
		#Set fit Chi2 for this pZ solution
		tmp1 = Double(1.0); tmp2 = Double(1.0); tmp3 = Double(1.0)
		minuit.mnstat(tmp1,tmp2,tmp3,Long(1),Long(1),Long(1))
		finalChi2s[iFit] = tmp1
	#if neither fit converged, return garbage
	if errflags[0]!=0 and errflags[1]!=0 :
		return None, None, None, None, 1000000000., -900., -900., -900., -900.
	#find which pZ solution gave better results and record best parameter values
#	print 'finalChi2s = '+str(finalChi2s)+'' #DEBUGGING
	final_par_vals = []
	final_met = TLorentzVector()
	if finalChi2s[0] < finalChi2s[1] :
		for i in range(len(bestParValues[0])) :
			final_par_vals.append(bestParValues[0][i])
		final_met.SetPtEtaPhiM(met1.Pt(),met1.Eta(),met1.Phi(),met1.M())
	else :
		for i in range(len(bestParValues[1])) :
			final_par_vals.append(bestParValues[1][i])
		final_met.SetPtEtaPhiM(met2.Pt(),met2.Eta(),met2.Phi(),met2.M())
		finalChi2s[0] = finalChi2s[1]
	#rescale the lepton and jet four vectors based on the final parameters
	lep_return = rescale(lepton,final_par_vals[1])
	lepb_return = rescale(lepb,final_par_vals[2])
	hadt_return = rescale(hadt,final_par_vals[3])
	#rebuild the neutrino post-rescaling
	newmetx = final_met.Px()+ (1.0-final_par_vals[1])*lep_return.Px()
	newmety = final_met.Py()+ (1.0-final_par_vals[1])*lep_return.Py()
	newmetx += (1.0-final_par_vals[2])*lepb_return.Px()
	newmety += (1.0-final_par_vals[2])*lepb_return.Py()
	newmetx += (1.0-final_par_vals[3])*hadt_return.Px()
	newmety += (1.0-final_par_vals[3])*hadt_return.Py()
	final_met.SetPx(newmetx); final_met.SetPy(newmety); final_met.SetPz(final_par_vals[0])
	final_met.SetE(final_met.Vect().Mag())
	#return everything
	return (lep_return,final_met,lepb_return,hadt_return,finalChi2s[0],final_par_vals[0],final_par_vals[1],final_par_vals[2],final_par_vals[3])

#top fitting function
def fcn(npar, deriv, f, par, flag) :
	#Build rescaled versions of the vectors involved
	l = rescale(lep_global_vec,par[1])
	bl = rescale(blep_global_vec,par[2])
	th = rescale(thad_global_vec,par[3])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_global_vec.Px()+(1.0-par[1])*lep_global_vec.Px()+(1.0-par[2])*blep_global_vec.Px()
					+(1.0-par[3])*thad_global_vec.Px() )
	newmety = ( met_global_vec.Py()+(1.0-par[1])*lep_global_vec.Py()+(1.0-par[2])*blep_global_vec.Py()
					+(1.0-par[3])*thad_global_vec.Py() )
	v = rescale(met_global_vec,1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl
	mwl2 = wl.M2(); mtl2 = tl.M2(); mth2 = th.M2();
	pdftl = 1./((mtl2-MT*MT)**2+MT*MT*WT*WT)
	pdfth = 1./((mth2-MT*MT)**2+MT*MT*WT*WT)
	pdfw  = ((MT*MT-mwl2)*(MT*MT-mwl2)*(2*MT*MT+mwl2))/((mwl2-MW*MW)*(mwl2-MW*MW)+MW*MW*WW*WW)
	pdf = pdftl*pdfth*pdfw
	lnL = 0.
	if pdf > 0.0 :
		lnL += log(pdf)    #need positive f
	else :
		print 'WARNING -- pdf is negative!!!'
		pdf = 1.e-50
		lnL += log(pdf)
	f[0] = ( -2.0*lnL+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
				+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ) )
	#and we don't need to return anything because minuit

#fourvector rescaling function
def rescale(vec,fac) :
	p2 = fac*fac*vec.Vect().Mag2()
	m2 = vec.M()*vec.M()
	newE = sqrt(p2+m2)
	#note that the function returns a NEW TLorentzVector, meaning the original is unaltered
	return TLorentzVector(fac*vec.Px(),fac*vec.Py(),fac*vec.Pz(),newE)