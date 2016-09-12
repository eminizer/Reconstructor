from ROOT import *
from math import *

#input filename
#single_filename = 'qq_semilep_TT_TEST.root'
single_filename = '../total_ttree_files/mcatnlo_semilep_TT_skim_all.root'

#output filename
outfilename = 'kinfit_debug_plots.root'

#plot class
class Plot(object) :
	def __init__(self,name_,object_,canvas_) :
		self._name = name_
		self._object = object_
		self._canvas = canvas_
	def getName(self) :
		return self._name
	def getObject(self) :
		return self._object
	def getCanvas(self) :
		return self._canvas
#SingleHist class (1- or 2-D histogram of a single file and a single set of variables given a cut for each)
class SingleHist(Plot) :
	def __init__(self,filename_=single_filename,treename_='tree',name_='new_1D_plot',title_='',bintuple_=None,plotstrings_='',cutstrings_='',weightstring_='',optionstring_='') :
		if weightstring_=='' :
			weightstring_=common_weights
		if optionstring_=='' :
			optionstring_ = 'COLZ' if len(bintuple_)==6 else 'HIST'
		canvas = TCanvas(name_+'_c',name_+' canvas',1100,900)
		f = TFile(filename_)
		tree = f.Get(treename_)
		if len(bintuple_)==3 :
			histo = TH1F(name_,title_,bintuple_[0],bintuple_[1],bintuple_[2])
		if len(bintuple_)==6 :
			histo = TH2F(name_,title_,bintuple_[0],bintuple_[1],bintuple_[2],bintuple_[3],bintuple_[4],bintuple_[5])
		histo.SetDirectory(0)
		for i in range(len(plotstrings_)) :
			if len(bintuple_)==3 :
				fulldrawstring='%s>>%s(%s,%s,%s)'%(plotstrings_[i],name_+'_'+str(i),bintuple_[0],bintuple_[1],bintuple_[2])
			elif len(bintuple_)==6 :
				fulldrawstring='%s>>%s(%s,%s,%s,%s,%s,%s)'%(plotstrings_[i],name_+'_'+str(i),bintuple_[0],bintuple_[1],bintuple_[2],bintuple_[3],bintuple_[4],bintuple_[5])
			fullweightcutstring = ''
			fullweightcutstring+='('+weightstring_+')' if weightstring_!='' else '(1.)'
			fullweightcutstring+='*('+cutstrings_[i]+')' if cutstrings_[i]!='' else '*(weight!=0.)'
			tree.Draw(fulldrawstring,fullweightcutstring,optionstring_)
			histo.Add(gROOT.FindObject(name_+'_'+str(i)))
		histo.SetTitle(title_)
		if optionstring_=='HIST' :
			histo.SetLineWidth(3)
		histo.Draw(optionstring_)
		Plot.__init__(self,name_,histo,canvas)

#start up output file
outfile = TFile(outfilename,'recreate')

#lists of all plots to save
all_plots = []

#weight strings
common_weights = '12900.*weight'
#cutstrings
preselection = '(hadt_pt>300. && hadt_SDM>50.)'
mu_selection = '(muTrig==1 && muon1_ID==1 && muon1_pt>40. && muon1_pt>ele1_pt && abs(muon1_eta)<2.4)'
el_selection = '(elTrig==1 && ele1_ID==1 && ele1_pt>40. && ele1_pt>muon1_pt && abs(ele1_eta)<2.4)'
lep_side_selection = '(scaled_lept_M>140. && scaled_lept_M<250.)'
had_side_pretag = '(hadt_SDM>50. && M_scaled>750. && hadt_tau21>0.1)'
antitag = '(hadt_SDM>140. && hadt_SDM<250. && hadt_tau32>0.55)'
full  = '(hadt_SDM>140. && hadt_SDM<250. && hadt_tau32<0.55)'
lowmass_sideband = '(hadt_SDM>50. && hadt_SDM<140.)'
highmass_sideband = '(hadt_SDM>250.)'
had_pretag_selection = '((%s) && ((%s) || (%s)) && (%s) && (%s))'%(preselection,mu_selection,el_selection,lep_side_selection,had_side_pretag)
antitag_selection 	 = '((%s) && ((%s) || (%s)) && (%s) && (%s) && (%s))'%(preselection,mu_selection,el_selection,lep_side_selection,had_side_pretag,antitag)
full_selection 	 = '((%s) && ((%s) || (%s)) && (%s) && (%s) && (%s))'%(preselection,mu_selection,el_selection,lep_side_selection,had_side_pretag,full)
lowmass_sideband_selection = '((%s) && ((%s) || (%s)) && (%s) && (%s) && (%s))'%(preselection,mu_selection,el_selection,lep_side_selection,had_side_pretag,lowmass_sideband)
highmass_sideband_selection = '((%s) && ((%s) || (%s)) && (%s) && (%s) && (%s))'%(preselection,mu_selection,el_selection,lep_side_selection,had_side_pretag,highmass_sideband)

#cutoption dictionary
#cutoptions = {'raw':'(weight!=0.)', 'pretag':had_pretag_selection, 'antitag':antitag_selection, 'full':full_selection, 'lowmass':lowmass_sideband_selection, 'highmass':highmass_sideband_selection}
cutoptions = {'raw':'(weight!=0.)', 'pretag':had_pretag_selection, 'full':full_selection}

#declare and draw the varietal plots
for cut in cutoptions :
	#kinematic fit chi2
	all_plots.append(SingleHist(name_=cut+'_chi2',
								title_='Kinematic fit #Chi^{2} ('+cut+' events); #Chi^{2}; events',
								bintuple_=(20,0.,100.),
								plotstrings_=['chi2'],
								cutstrings_=[cutoptions[cut]]))
	#Mass resolution
	all_plots.append(SingleHist(name_=cut+'_M_res',
								title_='Mass resolution ('+cut+' events); M_{r}-M_{MC} (GeV); events',
								bintuple_=(25,-500,500.),
								plotstrings_=['M_scaled-M_MC'],
								cutstrings_=[cutoptions[cut]]))
	#chi2 vs Mass resolution
	all_plots.append(SingleHist(name_=cut+'_chi2_vs_M_res',
								title_='Kinematic Fit #Chi^{2} vs. Mass Resolution ('+cut+' events); M_{r}-M_{MC} (GeV); #Chi^{2}',
								bintuple_=(25,-500,500.,20,0.,100.),
								plotstrings_=['chi2:(M_scaled-M_MC)'],
								cutstrings_=[cutoptions[cut]]))
	#par0 vs Chi2
	all_plots.append(SingleHist(name_=cut+'_par0_vs_chi2',
								title_='Neutrino p_{z} vs. kinematic fit #Chi^{2} ('+cut+' events); #Chi^{2}; p^{#nu}_{z} (GeV)',
								bintuple_=(20,0.,100.,25,-500,500.),
								plotstrings_=['par_0:chi2'],
								cutstrings_=[cutoptions[cut]]))
	#par1 vs Chi2
	all_plots.append(SingleHist(name_=cut+'_par1_vs_chi2',
								title_='Lepton scalefactor vs. kinematic fit #Chi^{2} ('+cut+' events); #Chi^{2}; #lambda_{l} (GeV)',
								bintuple_=(20,0.,100.,25,0.95,1.05),
								plotstrings_=['par_1:chi2'],
								cutstrings_=[cutoptions[cut]]))
	#par2 vs Chi2
	all_plots.append(SingleHist(name_=cut+'_par2_vs_chi2',
								title_='Leptonic b jet scalefactor vs. kinematic fit #Chi^{2} ('+cut+' events); #Chi^{2}; #lambda_{bl} (GeV)',
								bintuple_=(20,0.,100.,25,0.75,1.25),
								plotstrings_=['par_2:chi2'],
								cutstrings_=[cutoptions[cut]]))
	#par3 vs Chi2
	all_plots.append(SingleHist(name_=cut+'_par3_vs_chi2',
								title_='Hadronic t jet scalefactor vs. kinematic fit #Chi^{2} ('+cut+' events); #Chi^{2}; #lambda_{th} (GeV)',
								bintuple_=(20,0.,100.,25,0.75,1.25),
								plotstrings_=['par_3:chi2'],
								cutstrings_=[cutoptions[cut]]))
	#par0 vs Mass resolution
	all_plots.append(SingleHist(name_=cut+'_par0_vs_M_res',
								title_='Neutrino p_{z} vs. Mass Resolution ('+cut+' events); M_{r}-M_{MC} (GeV); p^{#nu}_{z} (GeV)',
								bintuple_=(25,-500,500.,25,-500,500.),
								plotstrings_=['par_0:(M_scaled-M_MC)'],
								cutstrings_=[cutoptions[cut]]))
	#par1 vs Mass resolution
	all_plots.append(SingleHist(name_=cut+'_par1_vs_M_res',
								title_='Lepton scalefactor vs. Mass Resolution ('+cut+' events); M_{r}-M_{MC} (GeV); #lambda_{l} (GeV)',
								bintuple_=(25,-500,500.,25,0.95,1.05),
								plotstrings_=['par_1:(M_scaled-M_MC)'],
								cutstrings_=[cutoptions[cut]]))
	#par2 vs Mass resolution
	all_plots.append(SingleHist(name_=cut+'_par2_vs_M_res',
								title_='Leptonic b jet scalefactor vs. Mass Resolution ('+cut+' events); M_{r}-M_{MC} (GeV); #lambda_{bl} (GeV)',
								bintuple_=(25,-500,500.,25,0.75,1.25),
								plotstrings_=['par_2:(M_scaled-M_MC)'],
								cutstrings_=[cutoptions[cut]]))
	#par3 vs Mass resolution
	all_plots.append(SingleHist(name_=cut+'_par3_vs_M_res',
								title_='Hadronic t jet scalefactor vs. Mass Resolution ('+cut+' events); M_{r}-M_{MC} (GeV); #lambda_{th} (GeV)',
								bintuple_=(25,-500,500.,25,0.75,1.25),
								plotstrings_=['par_3:(M_scaled-M_MC)'],
								cutstrings_=[cutoptions[cut]]))

	######################################################################################################################################
	###################################### 				HEMISPHERE CUT DEBUGGING PLOTS 				######################################
	######################################################################################################################################
	##NON-candidate AK4 jet delta Phi
	#all_plots.append(SingleHist(name_=cut+'_nonlepb_dPhi',
	#							title_='#Delta #phi between hadronic top and hard non-b AK4 jets; #Delta #phi; events',
	#							bintuple_=(20,0.,2.),
	#							plotstrings_=['ak4_1_dPhi',
	#										  'ak4_2_dPhi'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index>0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>1)']))
	##NON-candidate AK4 jet delta R
	#all_plots.append(SingleHist(name_=cut+'_nonlepb_dR',
	#							title_='#Delta R between hadronic top and hard non-b AK4 jets; #Delta R; events',
	#							bintuple_=(20,0.,6.),
	#							plotstrings_=['ak4_1_dR',
	#										  'ak4_2_dR'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index>0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>1)']))
	##NON-candidate AK4 jet delta R (version 2)
	#all_plots.append(SingleHist(name_=cut+'_nonlepb_dR_v2',
	#							title_='#Delta R (version two) between hadronic top and hard non-b AK4 jets; #Delta R; events',
	#							bintuple_=(20,0.,6.),
	#							plotstrings_=['ak4_1_dR_v2',
	#										  'ak4_2_dR_v2'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index>0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>1)']))
	##Leptonic b candidate delta R
	#all_plots.append(SingleHist(name_=cut+'_lepb_dR',
	#							title_='#Delta R between hadronic top and leptonic b candidate; #Delta R; events',
	#							bintuple_=(20,0.,6.),
	#							plotstrings_=['ak4_1_dR',
	#										  'ak4_2_dR',
	#										  'ak4_3_dR'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate delta Phi
	#all_plots.append(SingleHist(name_=cut+'_lepb_dPhi',
	#							title_='#Delta #phi between hadronic top and leptonic b candidate; #Delta #phi; events',
	#							bintuple_=(20,1.,4.),
	#							plotstrings_=['ak4_1_dPhi',
	#										  'ak4_2_dPhi',
	#										  'ak4_3_dPhi'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate delta R (Version 2)
	#all_plots.append(SingleHist(name_=cut+'_lepb_dR_v2',
	#							title_='#Delta R between hadronic top and leptonic b candidate; #Delta R; events',
	#							bintuple_=(20,0.,6.),
	#							plotstrings_=['ak4_1_dR_v2',
	#										  'ak4_2_dR_v2',
	#										  'ak4_3_dR_v2'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate delta Phi (Version 2)
	#all_plots.append(SingleHist(name_=cut+'_lepb_dPhi_v2',
	#							title_='#Delta #phi between hadronic top and leptonic b candidate; #Delta #phi; events',
	#							bintuple_=(20,1.,4.),
	#							plotstrings_=['ak4_1_dPhi_v2',
	#										  'ak4_2_dPhi_v2',
	#										  'ak4_3_dPhi_v2'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##AK4 jet delta R (Version 1 vs. Version 2)
	#all_plots.append(SingleHist(name_=cut+'_ak4_dR_v2_vs_v1',
	#							title_='#Delta R between hadronic top and AK4 jets, v2 vs. v1; #Delta R (v1); #Delta R (v2)',
	#							bintuple_=(20,0.,6.,20,0.,6.),
	#							plotstrings_=['ak4_1_dR_v2:ak4_1_dR',
	#										  'ak4_2_dR_v2:ak4_2_dR',
	#										  'ak4_3_dR_v2:ak4_3_dR'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index>=0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>=1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>=2)']))
	##AK4 jet delta phi (Version 1 vs. Version 2)
	#all_plots.append(SingleHist(name_=cut+'_ak4_dPhi_v2_vs_v1',
	#							title_='#Delta #phi between hadronic top and AK4 jets, v2 vs. v1; #Delta #phi (v1); #Delta #phi (v2)',
	#							bintuple_=(20,pi/2,pi,20,pi/2,pi),
	#							plotstrings_=['ak4_1_dPhi_v2:ak4_1_dPhi',
	#										  'ak4_2_dPhi_v2:ak4_2_dPhi',
	#										  'ak4_3_dPhi_v2:ak4_3_dPhi'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##AK4 jet delta R (Version 2 - Version 1)
	#all_plots.append(SingleHist(name_=cut+'_ak4_dR_v2_minus_v1',
	#							title_='#Delta R between hadronic top and AK4 jets, v2 - v1; #Delta R (v2) - #Delta R (v1); events',
	#							bintuple_=(20,-5.,5.),
	#							plotstrings_=['ak4_1_dR_v2-ak4_1_dR',
	#										  'ak4_2_dR_v2-ak4_2_dR',
	#										  'ak4_3_dR_v2-ak4_3_dR'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index>=0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>=1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>=2)']))
	##AK4 jet delta Phi (Version 2 - Version 1)
	#all_plots.append(SingleHist(name_=cut+'_ak4_dPhi_v2_minus_v1',
	#							title_='#Delta #phi between hadronic top and AK4 jets, v2 - v1; #Delta #phi (v2) - #Delta #phi (v1); events',
	#							bintuple_=(20,-5.,5.),
	#							plotstrings_=['ak4_1_dPhi_v2-ak4_1_dPhi',
	#										  'ak4_2_dPhi_v2-ak4_2_dPhi',
	#										  'ak4_3_dPhi_v2-ak4_3_dPhi'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index>=0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>=1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index>=2)']))
	##Leptonic b candidate Pt
	#all_plots.append(SingleHist(name_=cut+'_lepb_pt',
	#							title_='leptonic b candidate p_{T}; p_{T} (GeV); events',
	#							bintuple_=(30,0.,600.),
	#							plotstrings_=['ak4_1_pt',
	#										  'ak4_2_pt',
	#										  'ak4_3_pt'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate Delta R vs Pt
	#all_plots.append(SingleHist(name_=cut+'_lepb_dR_vs_pt',
	#							title_='leptonic b candidate #Delta R vs. p_{T}; p_{T} (GeV); #Delta R',
	#							bintuple_=(30,0.,600.,20,0.,10.),
	#							plotstrings_=['ak4_1_dR:ak4_1_pt',
	#										  'ak4_2_dR:ak4_2_pt',
	#										  'ak4_3_dR:ak4_3_pt'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate Delta Phi vs Pt
	#all_plots.append(SingleHist(name_=cut+'_lepb_dPhi_vs_pt',
	#							title_='leptonic b candidate #Delta #phi vs. p_{T}; p_{T} (GeV); #Delta #phi',
	#							bintuple_=(30,0.,600.,20,1.,4.),
	#							plotstrings_=['ak4_1_dPhi:ak4_1_pt',
	#										  'ak4_2_dPhi:ak4_2_pt',
	#										  'ak4_3_dPhi:ak4_3_pt'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate Delta R vs. cstar
	#all_plots.append(SingleHist(name_=cut+'_lepb_dR_vs_cstar',
	#							title_='leptonic b candidate #Delta R vs. c*; c*; #Delta R',
	#							bintuple_=(20,-1.,1.,20,0.,6.),
	#							plotstrings_=['ak4_1_dR:cstar_scaled',
	#										  'ak4_2_dR:cstar_scaled',
	#										  'ak4_3_dR:cstar_scaled'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate Delta R v2 vs. cstar
	#all_plots.append(SingleHist(name_=cut+'_lepb_dR_v2_vs_cstar',
	#							title_='leptonic b candidate #Delta R (version 2) vs. c*; c*; #Delta R',
	#							bintuple_=(20,-1.,1.,20,0.,6.),
	#							plotstrings_=['ak4_1_dR_v2:cstar_scaled',
	#										  'ak4_2_dR_v2:cstar_scaled',
	#										  'ak4_3_dR_v2:cstar_scaled'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate Delta R v2 vs. cstar
	#all_plots.append(SingleHist(name_=cut+'_lepb_dR_v2_vs_x_F',
	#							title_='leptonic b candidate #Delta R (version 2) vs. |x_{F}|; |x_{F}|; #Delta R',
	#							bintuple_=(20,0.,0.7,20,0.,6.),
	#							plotstrings_=['ak4_1_dR_v2:x_F_scaled',
	#										  'ak4_2_dR_v2:x_F_scaled',
	#										  'ak4_3_dR_v2:x_F_scaled'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))
	##Leptonic b candidate Delta R v2 vs. M
	#all_plots.append(SingleHist(name_=cut+'_lepb_dR_v2_vs_M',
	#							title_='leptonic b candidate #Delta R (version 2) vs. M; M (GeV); #Delta R',
	#							bintuple_=(20,750.,2750.,20,0.,6.),
	#							plotstrings_=['ak4_1_dR_v2:M_scaled',
	#										  'ak4_2_dR_v2:M_scaled',
	#										  'ak4_3_dR_v2:M_scaled'],
	#							cutstrings_=['('+cutoptions[cut]+' && lepb_cand_index==0)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==1)',
	#										 '('+cutoptions[cut]+' && lepb_cand_index==2)']))

#write plots and canvases to output file
outfile.cd()
for plot in all_plots :
	plot.getObject().Write()
for plot in all_plots :
	plot.getCanvas().Write()
outfile.Write()
outfile.Close()