from ROOT import *
from math import *
import os, glob
from optparse import OptionParser

parser = OptionParser()
#Run options
parser.add_option('--generator', 	  type='string', action='store', default='powheg', dest='generator')
(options, args) = parser.parse_args()

generator=options.generator.lower()

#directory with input files
inputdir = '../'+generator+'_TT'

#output file
outfilename = 'ttbar_sample_fraction_table_'+generator+'.txt'

#Chain up the files
chain = TChain('tree')
for f in glob.glob(inputdir+'/aggregated_*.root') :
	chain.Add(f)

#Cut details
cutnames = []; cutstrings = []; eventnumbers = []

cutnames.append('type-1 qqbar semilep e+jets'); cutstrings.append('eventType==0 && lepflavor==2 && eventTopology==1')
cutnames.append('type-1 gg semilep e+jets'); cutstrings.append('eventType==1 && lepflavor==2 && eventTopology==1')
cutnames.append('type-1 qqbar semilep mu+jets'); cutstrings.append('eventType==0 && lepflavor==1 && eventTopology==1')
cutnames.append('type-1 gg semilep mu+jets'); cutstrings.append('eventType==1 && lepflavor==1 && eventTopology==1')

cutnames.append('type-2 qqbar semilep e+jets'); cutstrings.append('eventType==0 && lepflavor==2 && eventTopology==2')
cutnames.append('type-2 gg semilep e+jets'); cutstrings.append('eventType==1 && lepflavor==2 && eventTopology==2')
cutnames.append('type-2 qqbar semilep mu+jets'); cutstrings.append('eventType==0 && lepflavor==1 && eventTopology==2')
cutnames.append('type-2 gg semilep mu+jets'); cutstrings.append('eventType==1 && lepflavor==1 && eventTopology==2')

cutnames.append('type-3 qqbar semilep e+jets'); cutstrings.append('eventType==0 && lepflavor==2 && eventTopology==3')
cutnames.append('type-3 gg semilep e+jets'); cutstrings.append('eventType==1 && lepflavor==2 && eventTopology==3')
cutnames.append('type-3 qqbar semilep mu+jets'); cutstrings.append('eventType==0 && lepflavor==1 && eventTopology==3')
cutnames.append('type-3 gg semilep mu+jets'); cutstrings.append('eventType==1 && lepflavor==1 && eventTopology==3')

#weight string
weightstring = '(((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_btag_eff*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'

os.system('echo "Event Numbers:">'+outfilename)

#plot into plots and get event numbers
for i in range(len(cutnames)) :
	print 'drawing %d of %d....'%(i+1,len(cutnames))
	chain.Draw('lep_Q>>h'+str(i)+'(4,-2.,2.)','('+weightstring+')*(weight!=0 && '+cutstrings[i]+')')
	thishist = gROOT.FindObject('h'+str(i))
	eventnumbers.append(thishist.Integral())
	os.system('echo "'+cutnames[i]+': '+str(eventnumbers[i])+'">>'+outfilename)

os.system('echo "\nFractions:">>'+outfilename)

total = sum(eventnumbers)

for i in range(len(cutnames)) :
	this = eventnumbers[i]
	frac = this/total
	frac_err = frac*sqrt(1./this+1./total)
	os.system('echo "'+cutnames[i]+': '+str(frac)+' +/- '+str(frac_err)+'">>'+outfilename)