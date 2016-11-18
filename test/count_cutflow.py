from ROOT import *
from array import array
from math import *
import os, glob

leptype='muons'
#leptype='electrons'
#leptype='all_leptons'

cutflow_filename = 'cutflow_count'
cutflow_filename+='_'+leptype+'.txt'

filenames = []
shortnames = []
weights = []
#MCATNLO TTJETS
#semileptonic qq
filenames.append('mcatnlo_qq_semilep_TT');					shortnames.append('Semileptonic TTBar')
#semileptonic gg
filenames.append('mcatnlo_gg_semilep_TT');					shortnames.append('Semileptonic TTBar')
#dileptonic 
filenames.append('mcatnlo_dilep_TT');						shortnames.append('Dileptonic TTBar')
#hadronic
filenames.append('mcatnlo_had_TT');							shortnames.append('Hadronic TTBar')
#Single top
filenames.append('ST_s-c');							shortnames.append('Single Top')
filenames.append('ST_t-c_top');							shortnames.append('Single Top')
filenames.append('ST_t-c_antitop');							shortnames.append('Single Top')
filenames.append('ST_tW-c_top');							shortnames.append('Single Top')
filenames.append('ST_tW-c_antitop');							shortnames.append('Single Top')

muon_data_filenames = []; ele_data_filenames = []
#data
muon_data_filenames.append('SingleMu_Run2016Bv2')
muon_data_filenames.append('SingleMu_Run2016C')
muon_data_filenames.append('SingleMu_Run2016D')
ele_data_filenames.append('SingleEl_Run2016Bv2')
ele_data_filenames.append('SingleEl_Run2016C')
ele_data_filenames.append('SingleEl_Run2016D')

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
	filenamelist = glob.glob('../total_ttree_files/'+filenames[i]+'_skim_all.root')
	#filenamelist = glob.glob('../total_ttree_files/'+filenames[i]+'_all.root')
	for filename in filenamelist :
		if filename.find('JES')==-1 and filename.find('JER')==-1 :
			MC_chains[index].Add(filename)
for muon_data_filename in muon_data_filenames :
	filenamelist = glob.glob('../total_ttree_files/'+muon_data_filename+'_skim_all.root')
	#filenamelist = glob.glob('../total_ttree_files/'+muon_data_filename+'_all.root')
	for filename in filenamelist :
		muon_data_chain.Add(filename)
for ele_data_filename in ele_data_filenames :
	filenamelist = glob.glob('../total_ttree_files/'+ele_data_filename+'_skim_all.root')
	#filenamelist = glob.glob('../total_ttree_files/'+ele_data_filename+'_all.root')
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
had_pretag = fullselection+' && hadt_SDM>25. && hadt_SDM<750. && M>750.'

cutnames.append('preselection'); 	  cutstrings.append('weight!=0.'); 		prior_cutstrings.append('weight!=0.')
cutnames.append('MET filters'); 		  cutstrings.append(metfilters); 		prior_cutstrings.append('weight!=0.')
cutnames.append('trigger'); 		  cutstrings.append(trigger); 			prior_cutstrings.append('weight!=0.')
cutnames.append('additional lepton veto'); 		  cutstrings.append(onelepton); 		prior_cutstrings.append('weight!=0.')
cutnames.append('lepton 2D isolation'); 		  cutstrings.append(isolepton); 		prior_cutstrings.append('weight!=0.')
cutnames.append('jet cuts'); 		  cutstrings.append(jetcuts); 			prior_cutstrings.append('weight!=0.')
cutnames.append('hadronic pretag'); 	  	  cutstrings.append(had_pretag); 		prior_cutstrings.append('weight!=0.')
cutnames.append('signal mass'); 	  	  cutstrings.append(had_pretag+' && hadt_SDM>105. && hadt_SDM<220.'); 		prior_cutstrings.append(had_pretag)
cutnames.append('signal $\\tau_{32}$'); 	  	  cutstrings.append(had_pretag+' && hadt_tau32<0.69'); 		prior_cutstrings.append(had_pretag)
cutnames.append('antitag $\\tau_{32}$'); 	  	  cutstrings.append(had_pretag+' && hadt_tau32>0.69'); 		prior_cutstrings.append(had_pretag)
cutnames.append('full signal selection'); 	  	  cutstrings.append(had_pretag+' && hadt_SDM>105. && hadt_SDM<220. && hadt_tau32<0.69'); 		prior_cutstrings.append('weight!=0.')
cutnames.append('full antitag selection'); 	  cutstrings.append(had_pretag+' && hadt_SDM>105. && hadt_SDM<220. && hadt_tau32>0.69'); 		prior_cutstrings.append('weight!=0.')
cutnames.append('full sideband selection'); 	  	  cutstrings.append(had_pretag+' && (hadt_SDM<105. || hadt_SDM>220.)'); 		prior_cutstrings.append('weight!=0.')

#weight string
weightstring = '12917.*weight*sf_pileup*sf_lep_ID*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'

data_events_at_cut = []; data_events_at_prior_cut = []
events_at_cut = []; events_at_prior_cut = []
dist = TH1F('dist','distribution',20,-1.0,1.0)
for i in range(len(cutnames)) :
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
		print '	Doing '+shortnames_done[j]+' ('+str(j+1)+' out of '+str(len(MC_chains))+')'
		tmp = dist.Clone('tmp')
		if leptype=='muons' :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+cutstrings[i]+' && lepflavor==1)')
		elif leptype=='electrons' :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+cutstrings[i]+' && lepflavor==2)')
		else :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+cutstrings[i]+')')
		events_at_cut[i].append(tmp.Integral())
		tmp = dist.Clone('tmp')
		if leptype=='muons' :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+prior_cutstrings[i]+' && lepflavor==1)')
		elif leptype=='electrons' :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+prior_cutstrings[i]+' && lepflavor==2)')
		else :
			MC_chains[j].Draw('cstar>>tmp','('+weightstring+')*('+prior_cutstrings[i]+')')
		events_at_prior_cut[i].append(tmp.Integral())

#print out the number of events in data and the efficiencies for the MC samples
#first line is just table headings for each cutflow and each sample type
first_line = 'Cut 		 	Data Events 		'
for shortname in shortnames_done :
	first_line += shortname+' eff		'
print first_line
os.system('echo "'+first_line+'" > '+cutflow_filename+'')
#second line is just the total number of events in data
second_line = 'Total events in data     & '+str(data_events_at_cut[0])+' 			'
for shortname in shortnames_done :
	second_line += ' &            '
second_line+='\\\\'
print second_line
os.system('echo "'+second_line+'" >> '+cutflow_filename+'')
#Each line after that is the cutflow number, then the number of events in data, then the eff. with uncertainty for each sample type
for i in range(1,len(cutnames)) :
	#Begin with the cutflow number and the number of events in data
	next_line = ''
	next_line +='%-28s'%(cutnames[i])
	next_line+='& %-8d'%(data_events_at_cut[i])
	for j in range(len(MC_chains)) :
		#calculate the efficiency for this sample type
		eff=900.; eff_err=900.
		if events_at_prior_cut[i][j]!=0. :
			eff = events_at_cut[i][j]/events_at_prior_cut[i][j]
			if events_at_cut[i][j]!=0. :
				eff_err = eff*sqrt(1./events_at_cut[i][j]+1./events_at_prior_cut[i][j])
		#add to the line
		next_line+=' & %-6.4f(%-6.4f)'%(eff,eff_err)
	next_line+='\\\\'
	print next_line
	os.system('echo "'+next_line+'" >> '+cutflow_filename+'')