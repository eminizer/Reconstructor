from ROOT import *
from array import array
from math import *
import os, glob
from optparse import OptionParser

parser = OptionParser()
#Run options
parser.add_option('--leptype', 	  type='string', action='store', default='all_leptons', dest='leptype',	   	  
	help='Use SingleMu or SingleEl data? ("muons" of "electrons")')
(options, args) = parser.parse_args()

leptype=options.leptype.lower()

cutflow_filename = 'cutflow_count'
cutflow_filename+='_'+leptype+'.csv'

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
#DYJets
filenames.append('DYJets_M-50_HT-70to100');		shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-100to200');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-200to400');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-400to600');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-600to800');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-800to1200');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-1200to2500');	shortnames.append('DYJets')
filenames.append('DYJets_M-50_HT-2500toInf');	shortnames.append('DYJets')
#WJets
filenames.append('WJets_HT-200to400');			shortnames.append('WJets')
filenames.append('WJets_HT-400to600');			shortnames.append('WJets')
filenames.append('WJets_HT-600to800');			shortnames.append('WJets')
filenames.append('WJets_HT-800to1200');			shortnames.append('WJets')
filenames.append('WJets_HT-1200to2500');		shortnames.append('WJets')
filenames.append('WJets_HT-2500toInf');			shortnames.append('WJets')
##QCD
filenames.append('QCD_HT-100to200');			shortnames.append('QCD')
filenames.append('QCD_HT-200to300');			shortnames.append('QCD')
filenames.append('QCD_HT-300to500');			shortnames.append('QCD')
filenames.append('QCD_HT-500to700');			shortnames.append('QCD')
filenames.append('QCD_HT-700to1000');			shortnames.append('QCD')
filenames.append('QCD_HT-1000to1500');			shortnames.append('QCD')
filenames.append('QCD_HT-1500to2000');			shortnames.append('QCD')
filenames.append('QCD_HT-2000toInf');			shortnames.append('QCD')
#Multiboson
filenames.append('WW_to_L_Nu_2Q');				shortnames.append('Multiboson')
filenames.append('WW_to_2L_2Nu');				shortnames.append('Multiboson')
filenames.append('WZ_to_L_Nu_2Q');				shortnames.append('Multiboson')
filenames.append('WZ_to_L_3Nu');				shortnames.append('Multiboson')
filenames.append('WZ_to_2L_2Q');				shortnames.append('Multiboson')
filenames.append('WZ_to_3L_Nu');				shortnames.append('Multiboson')
filenames.append('ZZ_to_2L_2Nu');				shortnames.append('Multiboson')
filenames.append('ZZ_to_2L_2Q');				shortnames.append('Multiboson')
filenames.append('ZZ_to_4L');					shortnames.append('Multiboson')
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
metfilters = 'metfilters==1'
trigger = 'trigger==1'
onelepton = 'onelepton==1'
isolepton = 'isolepton==1'
jetcuts = 'jetcuts==1'
fullselection = 'fullselection==1'
preselection = metfilters+' && '+trigger
lepside = onelepton+' && '+isolepton

cutnames.append('t1 skim'); 				  cutstrings.append('weight!=0.'); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t1 MET filters'); 			  cutstrings.append(metfilters); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t1 trigger'); 				  cutstrings.append(trigger); 			prior_cutstrings.append('weight!=0.')
cutnames.append('t1 additional lepton veto'); cutstrings.append(onelepton); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t1 lepton 2D isolation'); 	  cutstrings.append(isolepton); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t1 jet cuts'); 			  cutstrings.append(jetcuts); 			prior_cutstrings.append('weight!=0.')
cutnames.append('t1 preselection'); 		  cutstrings.append(preselection); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t1 lepcuts given precuts');  cutstrings.append(preselection+' && '+lepside); 			prior_cutstrings.append(preselection)
cutnames.append('t1 hadcuts given precuts');  cutstrings.append(preselection+' && '+jetcuts); 			prior_cutstrings.append(preselection)
cutnames.append('t1 hadcuts given lepcuts');  cutstrings.append(preselection+' && '+lepside+' && '+jetcuts); 			prior_cutstrings.append(preselection+' && '+lepside)
cutnames.append('t1 full selection'); 		  cutstrings.append(fullselection); 	prior_cutstrings.append('weight!=0.')

cutnames.append('t2 skim'); 				  cutstrings.append('weight!=0.'); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t2 MET filters'); 			  cutstrings.append(metfilters); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t2 trigger'); 				  cutstrings.append(trigger); 			prior_cutstrings.append('weight!=0.')
cutnames.append('t2 additional lepton veto'); cutstrings.append(onelepton); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t2 lepton 2D isolation'); 	  cutstrings.append(isolepton); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t2 jet cuts'); 			  cutstrings.append(jetcuts); 			prior_cutstrings.append('weight!=0.')
cutnames.append('t2 preselection'); 		  cutstrings.append(preselection); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t2 lepcuts given precuts');  cutstrings.append(preselection+' && '+lepside); 			prior_cutstrings.append(preselection)
cutnames.append('t2 hadcuts given precuts');  cutstrings.append(preselection+' && '+jetcuts); 			prior_cutstrings.append(preselection)
cutnames.append('t2 hadcuts given lepcuts');  cutstrings.append(preselection+' && '+lepside+' && '+jetcuts); 			prior_cutstrings.append(preselection+' && '+lepside)
cutnames.append('t2 full selection'); 		  cutstrings.append(fullselection); 	prior_cutstrings.append('weight!=0.')

cutnames.append('t3 skim'); 				  cutstrings.append('weight!=0.'); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t3 MET filters'); 			  cutstrings.append(metfilters); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t3 trigger'); 				  cutstrings.append(trigger); 			prior_cutstrings.append('weight!=0.')
cutnames.append('t3 additional lepton veto'); cutstrings.append(onelepton); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t3 lepton 2D isolation'); 	  cutstrings.append(isolepton); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t3 jet cuts'); 			  cutstrings.append(jetcuts); 			prior_cutstrings.append('weight!=0.')
cutnames.append('t3 preselection'); 		  cutstrings.append(preselection); 		prior_cutstrings.append('weight!=0.')
cutnames.append('t3 lepcuts given precuts');  cutstrings.append(preselection+' && '+lepside); 			prior_cutstrings.append(preselection)
cutnames.append('t3 hadcuts given precuts');  cutstrings.append(preselection+' && '+jetcuts); 			prior_cutstrings.append(preselection)
cutnames.append('t3 hadcuts given lepcuts');  cutstrings.append(preselection+' && '+lepside+' && '+jetcuts); 			prior_cutstrings.append(preselection+' && '+lepside)
cutnames.append('t3 full selection'); 		  cutstrings.append(fullselection); 	prior_cutstrings.append('weight!=0.')

#weight string
weightstring = '35867.*weight*sf_pileup*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'

#print out the number of events in data and the efficiencies for the MC samples
#first line is just table headings for each cutflow and each sample type
first_line = 'Cut,Data Events,Data Eff, Data Eff Err,'
for shortname in shortnames_done :
	first_line += shortname+' events,'+shortname+' eff,'+shortname+' eff err,'
print first_line
os.system('echo "'+first_line+'" > '+cutflow_filename+'')

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
			muon_data_chain.Draw('cstar>>tmp','(weight!=0.)*('+cutstrings[i]+')')
			data_events_at_cut.append(tmp.Integral())
			tmp = dist.Clone('tmp')
			muon_data_chain.Draw('cstar>>tmp','(weight!=0.)*('+prior_cutstrings[i]+')')
			data_events_at_prior_cut.append(tmp.Integral())
		elif leptype=='electrons' :
			tmp = dist.Clone('tmp')
			ele_data_chain.Draw('cstar>>tmp','(weight!=0.)*('+cutstrings[i]+')')
			data_events_at_cut.append(tmp.Integral())
			tmp = dist.Clone('tmp')
			ele_data_chain.Draw('cstar>>tmp','(weight!=0.)*('+prior_cutstrings[i]+')')
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
	#Each line after that is the cutflow number, then the number of events in data, then the eff. with uncertainty for each sample type
	#Begin with the cutflow number and the number of events in data
	data_eff=900.; data_eff_err=900.
	if data_events_at_prior_cut[i]!=0. :
		data_eff = data_events_at_cut[i]/data_events_at_prior_cut[i]
		if data_events_at_cut[i]!=0. :
			data_eff_err = data_eff*sqrt(1./data_events_at_cut[i]+1./data_events_at_prior_cut[i])
	next_line = '%s,%d,%.4f,%.4f'%(cutnames[i],data_events_at_cut[i],data_eff,data_eff_err)
	for j in range(len(MC_chains)) :
		#calculate the efficiency for this sample type
		eff=900.; eff_err=900.
		if events_at_prior_cut[i][j]!=0. :
			eff = events_at_cut[i][j]/events_at_prior_cut[i][j]
			if events_at_cut[i][j]!=0. :
				eff_err = eff*sqrt(1./events_at_cut[i][j]+1./events_at_prior_cut[i][j])
		#add to the line
		next_line+=',%.1f,%.4f,%.4f'%(events_at_cut[i][j],eff,eff_err)
	print next_line
	os.system('echo "'+next_line+'" >> '+cutflow_filename+'')