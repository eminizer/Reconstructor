#Imports
from ROOT import TMinuit, Long, Double, TLorentzVector
from math import *
from array import array

#Global variables
SIGMAJ  = 0.10 #jet momentum resolution
SIGMAL  = 0.02 #lepton momentum resolution

class KinFit(object) :

	def __init__(self,topology) :
		#Set up Minuit
		self.parNames = ['pZv','scalelep','scaleblep','scalehad1','scalehad2','scalehad3']
		self.parinivals = [0.,1.,1.,1.,1.,1.]
		self.parerrs = [0.0,0.0,0.0,0.0,0.0,0.0]
		#fit setup
		nPars = 5 if topology==1 else 6
		self.minuit = TMinuit(nPars)
		self.ierflag = Long(1)
		self.arglist = array( 'd', [-1.0] )
		self.minuit.mnexcm('SET PRINT', self.arglist, 1,self.ierflag)
		self.minuit.mnexcm('SET NOWARNINGS',self.arglist,1,self.ierflag)
		self.arglist[0] = 100000.
		#add parameters to the fitter
		for i in range(nPars) :
			self.minuit.mnparm(i,self.parNames[i],self.parinivals[i],1.0,0.,0.,self.ierflag)

	##################################  Static minuit functions  ##################################

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
		wl = v + l; tl = wl + bl; th = TTBarReconstructor.hadt_global_vec-TTBarReconstructor.had1_global_vec-TTBarReconstructor.had2_global_vec+hs1+hs2
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l1,TTBarReconstructor.WT_l1,th.M(),TTBarReconstructor.MT_h1,
					 TTBarReconstructor.WT_h1,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
		f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
					+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) )

	#type 2 (boosted untagged) top minimization function
	@staticmethod
	def fcnt2(npar, deriv, f, par, flag) :
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
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l2,TTBarReconstructor.WT_l2,th.M(),TTBarReconstructor.MT_h2,
					 TTBarReconstructor.WT_h2,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
		f[0] = ( pdf+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
					+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ)+(par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) 
					+(par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ) )

	#type 3 (resolved) kinematic fitting function
	@staticmethod
	def fcnt3(npar, deriv, f, par, flag) :
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
		pdf = getPDF(tl.M(),TTBarReconstructor.MT_l3,TTBarReconstructor.WT_l3,th.M(),TTBarReconstructor.MT_h3,
					 TTBarReconstructor.WT_h3,TTBarReconstructor.MT,wl.M2(),TTBarReconstructor.MW,TTBarReconstructor.WW)
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