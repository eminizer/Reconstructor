#Global variables
#constants
MW = 80.4

#Imports
import copy
from math import *

#finds the two possible neutrino fourvectors based on the lep and met without pZ
def setupMET(lep_vec,metvec) :
	met1 = copy.deepcopy(metvec); met2 = copy.deepcopy(metvec)
	pTv    = metvec.Pt()
	phivec = [cos(metvec.Phi()),sin(metvec.Phi())]
	Elep   = lep_vec.E()
	plep   = lep_vec.Vect().Mag()
	pZlep  = lep_vec.Pz()
	pPhi   = lep_vec.Px()*phivec[0]+lep_vec.Py()*phivec[1]
	arg0   = MW*MW+plep*plep-Elep*Elep+2.*pTv*pPhi
	arg    = Elep*Elep*(4.*pTv*pTv*(pZlep*pZlep-Elep*Elep)+arg0*arg0) #discriminant in the quadratic equation solution
#	print ' arg = %.4f = (%.4f)^2*(4*(%.4f)^2*((%.4f)^2-(%.4f)^2)+(%.4f)^2'%(arg,Elep,pTv,pZlep,Elep,arg0) #DEBUGGING
	if not arg > 0 : #If discriminant is imaginary
		pzv1 = pZlep*arg0/(2.*(Elep*Elep-pZlep*pZlep))
		met1.SetPz(pzv1)
		met1.SetE(sqrt(met1.Px()*met1.Px()+met1.Py()*met1.Py()+met1.Pz()*met1.Pz()))
		met2.SetPtEtaPhiM(met1.Pt(),met1.Eta(),met1.Phi(),met1.M())
	else : #have two choices for the neutrino Pz from the quadratic equation
		pzv1 = (pZlep*arg0+sqrt(arg))/(2.*(Elep*Elep-pZlep*pZlep))
		met1.SetPz(pzv1)
		met1.SetE(sqrt(met1.Px()*met1.Px()+met1.Py()*met1.Py()+met1.Pz()*met1.Pz()))
		pzv2 = (pZlep*arg0-sqrt(arg))/(2.*(Elep*Elep-pZlep*pZlep))
		met2.SetPz(pzv2)
		met2.SetE(sqrt(met2.Px()*met2.Px()+met2.Py()*met2.Py()+met2.Pz()*met2.Pz()))
	return (met1,met2)
