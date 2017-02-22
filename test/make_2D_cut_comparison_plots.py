from ROOT import *

#declare output file
outputfile = TFile('lepton_2D_cut_plots.root','recreate')

#declare plots
sig_h = TH2D('sig_h','lepton 2D cut for semileptonic t#bar{t} signal; #Delta R; p_{T}^{rel} (GeV)',30,0.,1.5,30,0.,150.)
bkg_h = TH2D('bkg_h','lepton 2D cut for various backgrounds; #Delta R; p_{T}^{rel} (GeV)',30,0.,1.5,30,0.,150.)

#open files
sig_filenames = ['mcatnlo_semilep_TT']
bkg_filenames = ['mcatnlo_dilep_TT',
				 'mcatnlo_had_TT',
				 'DYJets_M-50_HT-100to200',
				 'DYJets_M-50_HT-200to400',
				 'DYJets_M-50_HT-400to600',
				 'ST_s-c',
				 'ST_t-c_antitop',
				 'ST_t-c_top',
				 'ST_tW-c_antitop',
				 'ST_tW-c_top',
				 'WJets_HT-200to400',
				 'WJets_HT-400to600',
				 'WW_to_2L_2Nu',
				 'WW_to_L_Nu_2Q',
				 'WZ_to_2L_2Q',
				 'WZ_to_3L_Nu',
				 'WZ_to_L_3Nu',
				 'WZ_to_L_Nu_2Q',
				 'ZZ_to_2L_2Nu',
				 'ZZ_to_2L_2Q',
				 'ZZ_to_4L']
for i in range(len(sig_filenames)) :
	sig_filenames[i]='../total_ttree_files/'+sig_filenames[i]+'_all.root'
for i in range(len(bkg_filenames)) :
	bkg_filenames[i]='../total_ttree_files/'+bkg_filenames[i]+'_skim_all.root'

#plot into plots
#weights and cuts
weights = '12917.*weight*sf_pileup*sf_lep_ID*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'
cuts = 'metfilters==1 && trigger==1 && onelepton==1 && jetcuts==1'
#signal plots
sig_chain = TChain('tree')
for n in sig_filenames :
	sig_chain.Add(n)
#sig_chain.Draw('muon1_relPt:muon1_dR>>sig_h_1(30,0.,1.5,30,0.,150.)','('+weights+')*('+cuts+' && lepflavor==1)','COLZ')
sig_chain.Draw('muon1_relPt:muon1_dR>>sig_h_1(30,0.,1.5,30,0.,150.)','('+weights+')*('+cuts+' && lepflavor==1 && (muon1_relPt>20. || muon1_dR>0.4))','COLZ')
#sig_chain.Draw('ele1_relPt:ele1_dR>>sig_h_2(30,0.,1.5,30,0.,150.)','('+weights+')*('+cuts+' && lepflavor==2)','COLZ')
sig_chain.Draw('ele1_relPt:ele1_dR>>sig_h_2(30,0.,1.5,30,0.,150.)','('+weights+')*('+cuts+' && lepflavor==2 && (ele1_relPt>20. || ele1_dR>0.4) )','COLZ')
sig_h.Add(gROOT.FindObject('sig_h_1')); sig_h.Add(gROOT.FindObject('sig_h_2'))
#background plots
bkg_chain = TChain('tree')
for n in bkg_filenames :
	bkg_chain.Add(n)
#bkg_chain.Draw('muon1_relPt:muon1_dR>>bkg_h_1(30,0.,1.5,30,0.,150.)','('+weights+')*('+cuts+' && lepflavor==1)','COLZ')
bkg_chain.Draw('muon1_relPt:muon1_dR>>bkg_h_1(30,0.,1.5,30,0.,150.)','('+weights+')*('+cuts+' && lepflavor==1 && (muon1_relPt>20. || muon1_dR>0.4))','COLZ')
#bkg_chain.Draw('ele1_relPt:ele1_dR>>bkg_h_2(30,0.,1.5,30,0.,150.)','('+weights+')*('+cuts+' && lepflavor==2)','COLZ')
bkg_chain.Draw('ele1_relPt:ele1_dR>>bkg_h_2(30,0.,1.5,30,0.,150.)','('+weights+')*('+cuts+' && lepflavor==2 && (ele1_relPt>20. || ele1_dR>0.4) )','COLZ')
bkg_h.Add(gROOT.FindObject('bkg_h_1')); bkg_h.Add(gROOT.FindObject('bkg_h_2'))

#save plots to file
outputfile.cd()
sig_h.Write()
bkg_h.Write()

#close file
outputfile.Write()
outputfile.Close()