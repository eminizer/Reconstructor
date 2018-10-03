from ROOT import *
from array import array
from math import *
import os, glob
from optparse import OptionParser

def float_to_str(f):
    float_string = repr(f)
    if 'e' in float_string:  # detect scientific notation
        digits, exp = float_string.split('e')
        digits = digits.replace('.', '').replace('-', '')
        exp = int(exp)
        zero_padding = '0' * (abs(int(exp)) - 1)  # minus 1 for decimal point in the sci notation
        sign = '-' if f < 0 else ''
        if exp > 0:
            float_string = '{}{}{}.0'.format(sign, digits, zero_padding)
        else:
            float_string = '{}0.{}{}'.format(sign, zero_padding, digits)
    return float_string

parser = OptionParser()
#Run options
parser.add_option('--leptype', 	  type='string', action='store', default='all_leptons', dest='leptype',	   	  
	help='Use SingleMu or SingleEl data? ("muons" or "electrons")')
parser.add_option('--topologies', 	  type='string', action='store', default='t1_t2_t3', dest='topologies',	   	  
	help='Which topologies should we make tables for? (default: t1_t2_t3)')
(options, args) = parser.parse_args()

leptype=options.leptype.lower()
topologies = options.topologies.split('_')

cutflow_filename = 'cutflow_count'
cutflow_filename+='_'+leptype+'_'+options.topologies+'.csv'

tex_filename = 'cutflow_count'
tex_filename+='_'+leptype+'_'+options.topologies+'_texlines.txt'

filenames = []
shortnames = []
weights = []
#POWHEG TT
filenames.append('powheg_TT');					shortnames.append('powheg semilep')
filenames.append('powheg_TT');					shortnames.append('powheg dilep')
filenames.append('powheg_TT');					shortnames.append('powheg hadronic')
#MCATNLO TT
filenames.append('mcatnlo_TT');					shortnames.append('mcatnlo semilep')
filenames.append('mcatnlo_TT');					shortnames.append('mcatnlo dilep')
filenames.append('mcatnlo_TT');					shortnames.append('mcatnlo hadronic')
#Single top
filenames.append('ST_s-c');						shortnames.append('Single T')
filenames.append('ST_t-c_top');					shortnames.append('Single T')
filenames.append('ST_t-c_antitop');				shortnames.append('Single T')
filenames.append('ST_tW-c_top');				shortnames.append('Single T')
filenames.append('ST_tW-c_antitop');			shortnames.append('Single T')
#WJets
filenames.append('WJets_HT-200to400');			shortnames.append('WJets')
filenames.append('WJets_HT-400to600');			shortnames.append('WJets')
filenames.append('WJets_HT-600to800');			shortnames.append('WJets')
filenames.append('WJets_HT-800to1200');			shortnames.append('WJets')
filenames.append('WJets_HT-1200to2500');		shortnames.append('WJets')
filenames.append('WJets_HT-2500toInf');			shortnames.append('WJets')
#DYJets
filenames.append('DYJets_M-50_HT-70to100');		shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-100to200');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-200to400');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-400to600');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-600to800');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-800to1200');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-1200to2500');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-2500toInf');	shortnames.append('DYJets')
###QCD
#filenames.append('QCD_HT-100to200');			shortnames.append('QCD')
#filenames.append('QCD_HT-200to300');			shortnames.append('QCD')
#filenames.append('QCD_HT-300to500');			shortnames.append('QCD')
#filenames.append('QCD_HT-500to700');			shortnames.append('QCD')
#filenames.append('QCD_HT-700to1000');			shortnames.append('QCD')
#filenames.append('QCD_HT-1000to1500');			shortnames.append('QCD')
#filenames.append('QCD_HT-1500to2000');			shortnames.append('QCD')
#filenames.append('QCD_HT-2000toInf');			shortnames.append('QCD')
##Multiboson
#filenames.append('WW_to_L_Nu_2Q');				shortnames.append('Multiboson')
#filenames.append('WW_to_2L_2Nu');				shortnames.append('Multiboson')
#filenames.append('WZ_to_L_Nu_2Q');				shortnames.append('Multiboson')
#filenames.append('WZ_to_L_3Nu');				shortnames.append('Multiboson')
#filenames.append('WZ_to_2L_2Q');				shortnames.append('Multiboson')
#filenames.append('WZ_to_3L_Nu');				shortnames.append('Multiboson')
#filenames.append('ZZ_to_2L_2Nu');				shortnames.append('Multiboson')
#filenames.append('ZZ_to_2L_2Q');				shortnames.append('Multiboson')
#filenames.append('ZZ_to_4L');					shortnames.append('Multiboson')
#data
muon_data_filenames = []
muon_data_filenames.append('SingleMu_Run2016Bv2')
muon_data_filenames.append('SingleMu_Run2016C')
muon_data_filenames.append('SingleMu_Run2016D')
muon_data_filenames.append('SingleMu_Run2016E')
muon_data_filenames.append('SingleMu_Run2016F')
muon_data_filenames.append('SingleMu_Run2016G')
muon_data_filenames.append('SingleMu_Run2016Hv2')
muon_data_filenames.append('SingleMu_Run2016Hv3')
ele_data_filenames = []
ele_data_filenames.append('SingleEl_Run2016Bv2')
ele_data_filenames.append('SingleEl_Run2016C')
ele_data_filenames.append('SingleEl_Run2016D')
ele_data_filenames.append('SingleEl_Run2016E')
ele_data_filenames.append('SingleEl_Run2016F')
ele_data_filenames.append('SingleEl_Run2016G')
ele_data_filenames.append('SingleEl_Run2016Hv2')
ele_data_filenames.append('SingleEl_Run2016Hv3')

#Chain up the files
MC_chains = []
muon_data_chain = TChain('tree')
ele_data_chain = TChain('tree')
shortnames_done = []
for i in range(len(filenames)) :
	index = 0
	if shortnames[i] not in shortnames_done :
		MC_chains.append(TChain('tree'))
		index = len(MC_chains)-1
		shortnames_done.append(shortnames[i])
	else :
		index = shortnames_done.index(shortnames[i])
	filenamelist = glob.glob('../'+filenames[i]+'/aggregated_'+filenames[i]+'_*.root')
	for filename in filenamelist :
		if filename.find('JES')==-1 and filename.find('JER')==-1 :
			MC_chains[index].Add(filename)
for muon_data_filename in muon_data_filenames :
	filenamelist = glob.glob('../'+muon_data_filename+'/aggregated_'+muon_data_filename+'_*.root')
	for filename in filenamelist :
		muon_data_chain.Add(filename)
for ele_data_filename in ele_data_filenames :
	filenamelist = glob.glob('../'+ele_data_filename+'/aggregated_'+ele_data_filename+'_*.root')
	for filename in filenamelist :
		ele_data_chain.Add(filename)

#Cut details
cutnames = []; cutstrings = []; prior_cutstrings = []
goodpv = 'goodpv'
metfilters = 'metfilters==1'
METcuts = 'METcuts==1'
trigger = 'trigger==1'
isolepton = 'isolepton==1'
onelepton = 'onelepton==1'
ak4jetmult = 'ak4jetmult==1'
btags = 'btags==1'
jetcuts = 'jetcuts==1'
lepcuts = 'lepcuts==1'
validminimization = 'validminimization==1'
kinfitchi2 = 'kinfitchi2==1'
recoleptM = 'recoleptM==1'
fullselection = 'fullselection==1'
wjetscrselection = 'wjets_cr_selection==1'
qcdASRselection = 'qcd_A_SR_selection==1'
qcdBSRselection = 'qcd_B_SR_selection==1'
qcdCSRselection = 'qcd_C_SR_selection==1'
qcdACRselection = 'qcd_A_CR_selection==1'
qcdBCRselection = 'qcd_B_CR_selection==1'
qcdCCRselection = 'qcd_C_CR_selection==1'

if 't1' in topologies :
	cutnames.append('t1 skim'); 				  cutstrings.append('weight!=0.'); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 good PV');				  cutstrings.append(goodpv); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 MET filters'); 			  cutstrings.append(metfilters); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 trigger'); 				  cutstrings.append(trigger); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 MET cuts'); 			  cutstrings.append(METcuts); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 lep cuts'); 			  cutstrings.append(lepcuts); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 lepton 2D isolation'); 	  cutstrings.append(isolepton); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 additional lepton veto'); cutstrings.append(onelepton); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 nbtags'); 				  cutstrings.append(btags); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 jet cuts'); 			  cutstrings.append(jetcuts); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 valid minimization'); 	  cutstrings.append(validminimization); prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 reco lept mass'); 	  	  cutstrings.append(recoleptM); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 kinfit chi2'); 	  	  	  cutstrings.append(kinfitchi2); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 full selection'); 		  cutstrings.append(fullselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 WJets CR selection'); 	  cutstrings.append(wjetscrselection); 	prior_cutstrings.append('weight!=0.')
	#if leptype=='electrons' :
	cutnames.append('t1 QCD A SR selection'); 	  cutstrings.append(qcdASRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 QCD B SR selection'); 	  cutstrings.append(qcdBSRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 QCD C SR selection'); 	  cutstrings.append(qcdCSRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 QCD A CR selection'); 	  cutstrings.append(qcdACRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 QCD B CR selection'); 	  cutstrings.append(qcdBCRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t1 QCD C CR selection'); 	  cutstrings.append(qcdCCRselection); 	prior_cutstrings.append('weight!=0.')


if 't2' in topologies :
	cutnames.append('t2 skim'); 				  cutstrings.append('weight!=0.'); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 good PV');				  cutstrings.append(goodpv); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 MET filters'); 			  cutstrings.append(metfilters); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 trigger'); 				  cutstrings.append(trigger); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 MET cuts'); 			  cutstrings.append(METcuts); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 lep cuts'); 			  cutstrings.append(lepcuts); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 lepton 2D isolation'); 	  cutstrings.append(isolepton); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 additional lepton veto'); cutstrings.append(onelepton); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 nbtags'); 				  cutstrings.append(btags); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 jet cuts'); 			  cutstrings.append(jetcuts); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 valid minimization'); 	  cutstrings.append(validminimization); prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 reco lept mass'); 	  	  cutstrings.append(recoleptM); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 kinfit chi2'); 	  	  	  cutstrings.append(kinfitchi2); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 full selection'); 		  cutstrings.append(fullselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 WJets CR selection'); 	  cutstrings.append(wjetscrselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 QCD A SR selection'); 	  cutstrings.append(qcdASRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 QCD B SR selection'); 	  cutstrings.append(qcdBSRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 QCD C SR selection'); 	  cutstrings.append(qcdCSRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 QCD A CR selection'); 	  cutstrings.append(qcdACRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 QCD B CR selection'); 	  cutstrings.append(qcdBCRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t2 QCD C CR selection'); 	  cutstrings.append(qcdCCRselection); 	prior_cutstrings.append('weight!=0.')

if 't3' in topologies :
	cutnames.append('t3 skim'); 				  cutstrings.append('weight!=0.'); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 good PV');				  cutstrings.append(goodpv); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 MET filters'); 			  cutstrings.append(metfilters); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 trigger'); 				  cutstrings.append(trigger); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 MET cuts'); 			  cutstrings.append(METcuts); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 lepton isolation'); 	  cutstrings.append(isolepton); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 additional lepton veto'); cutstrings.append(onelepton); 		prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 nbtags'); 				  cutstrings.append(btags); 			prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 valid minimization'); 	  cutstrings.append(validminimization); prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 full selection'); 		  cutstrings.append(fullselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 QCD A SR selection'); 	  cutstrings.append(qcdASRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 QCD B SR selection'); 	  cutstrings.append(qcdBSRselection); 	prior_cutstrings.append('weight!=0.')
	cutnames.append('t3 QCD C SR selection'); 	  cutstrings.append(qcdCSRselection); 	prior_cutstrings.append('weight!=0.')

#weight string
weightstring = '(((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_btag_eff*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'

#print out the number of events in data and the efficiencies for the MC samples
#first line is just table headings for each cutflow and each sample type
first_line = 'Cut,Data Events,Data Eff, Data Eff Err,'
first_tex_line = 'Cut & Data events '
for shortname in shortnames_done :
	first_line += shortname+' events,'+shortname+' eff,'+shortname+' eff err,'
	first_tex_line+='& '+shortname+' eff '
print first_line
first_tex_line+=' \\\\'
os.system('echo "'+first_line+'" > '+cutflow_filename+'')
os.system('echo "'+first_tex_line+'" > '+tex_filename+'')

#main loop
data_events_at_cut = []; data_events_at_prior_cut = []
events_at_cut = []; events_at_prior_cut = []
dist = TH1F('dist','distribution',20,-1.0,1.0)
for i in range(len(cutnames)) :
	if cutnames[i].startswith('t1') :
		cutstrings[i]+=' && eventTopology==1'
		prior_cutstrings[i]+=' && eventTopology==1'
	elif cutnames[i].startswith('t2') :
		cutstrings[i]+=' && eventTopology==2'
		prior_cutstrings[i]+=' && eventTopology==2'
	elif cutnames[i].startswith('t3') :
		cutstrings[i]+=' && eventTopology==3'
		prior_cutstrings[i]+=' && eventTopology==3'	
	if (leptype=='muons' and len(muon_data_filenames)>0) or (leptype=='electrons' and len(ele_data_filenames)>0) :
		print 'Getting numbers of data events for cut '+cutnames[i]+' ('+str(i+1)+' out of '+str(len(cutnames))+')'
		if leptype=='muons' :
			tmp = dist.Clone('tmp')
			muon_data_chain.Draw('cstar>>tmp','(weight!=0.)*('+cutstrings[i]+' && lepflavor==1)')
			data_events_at_cut.append(tmp.Integral())
			tmp = dist.Clone('tmp')
			muon_data_chain.Draw('cstar>>tmp','(weight!=0.)*('+prior_cutstrings[i]+' && lepflavor==1)')
			data_events_at_prior_cut.append(tmp.Integral())
		elif leptype=='electrons' :
			tmp = dist.Clone('tmp')
			ele_data_chain.Draw('cstar>>tmp','(weight!=0.)*('+cutstrings[i]+' && lepflavor==2)')
			data_events_at_cut.append(tmp.Integral())
			tmp = dist.Clone('tmp')
			ele_data_chain.Draw('cstar>>tmp','(weight!=0.)*('+prior_cutstrings[i]+' && lepflavor==2)')
			data_events_at_prior_cut.append(tmp.Integral())
	events_at_cut.append([]); events_at_prior_cut.append([])
	print 'Getting numbers of MC events for cut '+cutnames[i]+' ('+str(i+1)+' out of '+str(len(cutnames))+')'
	for j in range(len(MC_chains)) :
		thiscut = cutstrings[i]
		priorcut = prior_cutstrings[i]
		if shortnames_done[j].find('semilep')!=-1 :
			thiscut+=' && eventType<2'
			priorcut+=' && eventType<2'
		elif shortnames_done[j].find('dilep')!=-1 :
			thiscut+=' && eventType==2'
			priorcut+=' && eventType==2'
		elif shortnames_done[j].find('hadronic')!=-1 :
			thiscut+=' && eventType==3'
			priorcut+=' && eventType==3'
		print '	Doing '+shortnames_done[j]+' ('+str(j+1)+' out of '+str(len(MC_chains))+')'
		tmp = dist.Clone('tmp')
		if leptype=='muons' :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+thiscut+' && lepflavor==1)')
		elif leptype=='electrons' :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+thiscut+' && lepflavor==2)')
		else :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+thiscut+')')
		events_at_cut[i].append(tmp.Integral())
		tmp = dist.Clone('tmp')
		if leptype=='muons' :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+priorcut+' && lepflavor==1)')
		elif leptype=='electrons' :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+priorcut+' && lepflavor==2)')
		else :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+priorcut+')')
		events_at_prior_cut[i].append(tmp.Integral())
	#Each line is the cutflow number, then the number of events in data, then the eff. with uncertainty for each sample type
	#Begin with the cutflow number and the number of events in data
	data_eff=900.; data_eff_err=900.
	if data_events_at_prior_cut[i]!=0. :
		data_eff = data_events_at_cut[i]/data_events_at_prior_cut[i]
		if data_events_at_cut[i]!=0. :
			data_eff_err = data_eff*sqrt(1./data_events_at_cut[i]+1./data_events_at_prior_cut[i])
	next_line = '%s,%d,%.4f,%.4f'%(cutnames[i],data_events_at_cut[i],data_eff,data_eff_err)
	next_tex_line = '%s & %d'%(cutnames[i],data_events_at_cut[i])
	for j in range(len(MC_chains)) :
		#calculate the efficiency for this sample type
		eff=900.; eff_err=900.
		if events_at_prior_cut[i][j]!=0. :
			eff = abs(events_at_cut[i][j]/events_at_prior_cut[i][j])
			if events_at_cut[i][j]!=0. :
				eff_err = abs(eff*sqrt(abs(1./events_at_cut[i][j]+1./events_at_prior_cut[i][j])))
		#adjust the values for sig figs
		if len(str(eff_err).replace('.','').strip('0'))>2 :
			eff_err_str='('
			eff_err_digit1 = float_to_str(eff_err).replace('.','').strip('0')[0]
			eff_err_digit2 = float_to_str(eff_err).replace('.','').strip('0')[1]
			eff_err_digit3 = float_to_str(eff_err).replace('.','').strip('0')[2]
			#print 'eff = %s, eff_err = %s, dig1 = %s, dig2 = %s, dig3 = %s'%(str(eff),str(eff_err),str(eff_err_digit1),str(eff_err_digit2),str(eff_err_digit3)) #DEBUG
			if int(eff_err_digit1)<3 :
				eff_err_str+=eff_err_digit1
				if int(eff_err_digit3)<5 or int(eff_err_digit2)==9 :
					eff_err_str+=eff_err_digit2
				else :
					eff_err_str+=str(int(eff_err_digit2)+1)
			else :
				if int(eff_err_digit2)<5 or int(eff_err_digit1)==9 :
					eff_err_str+=eff_err_digit1
				else :
					eff_err_str+=str(int(eff_err_digit1)+1)
			eff_err_str+=')'
			#print 'str(eff)=%s, str(eff_err).find(eff_err_digit1)=%d, len(eff_err_str)=%d, str(eff_err)=%s, eff_err_digit1=%s, eff_err_digit2=%s'%(str(eff),str(eff_err).find(eff_err_digit1),len(eff_err_str),str(eff_err),eff_err_digit1,eff_err_digit2) #DEBUG
			eff_str_plus_one = float_to_str(eff)[:float_to_str(eff_err).find(eff_err_digit1)+len(eff_err_str)-1]
			eff_str = ''
			#print 'eff_err_str=%s, eff_str_plus_one=%s, eff=%.4f, eff_err=%.4f'%(eff_err_str,eff_str_plus_one,eff,eff_err) #DEBUG
			if eff_str_plus_one.endswith('.') :
				if int(float_to_str(eff)[:float_to_str(eff_err).find(eff_err_digit1)+len(eff_err_str)][-1])<5 or int(eff_str_plus_one[-2])==9 :
					eff_str+=eff_str_plus_one[:-1]+eff_err_str
				else :
					eff_str+=eff_str_plus_one[:-2]+str(int(eff_str_plus_one[-2])+1)+eff_err_str
			else :
				if int(eff_str_plus_one[-1])<5 or int(eff_str_plus_one[-2])==9 :
					eff_str+=eff_str_plus_one[:len(eff_str_plus_one)-1]+eff_err_str
				else :
					eff_str+=eff_str_plus_one[:len(eff_str_plus_one)-2]+str(int(eff_str_plus_one[-2])+1)+eff_err_str
		else :
			eff_str = float_to_str(eff)+'('+float_to_str(eff_err)+')'
		#print 'eff_err = %s, d1 = %s, d2 = %s, d3 = %s, eff = %s, eff_str = %s'%(str(eff_err),eff_err_digit1,eff_err_digit2,eff_err_digit3,str(eff),eff_str) #DEBUG
		#add to the line
		next_line+=',%.1f,%.4f,%.4f'%(events_at_cut[i][j],eff,eff_err)
		next_tex_line+=' & %s'%eff_str
	next_tex_line+=' \\\\'
	print next_line
	print next_tex_line
	os.system('echo "'+next_line+'" >> '+cutflow_filename+'')
	os.system('echo "'+next_tex_line+'" >> '+tex_filename+'')