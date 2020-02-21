from ROOT import *
from math import sqrt
gROOT.SetBatch()

#data filenames
mu_data_fns = ['SingleMu_Run2016Bv2_skim_all.root',
				'SingleMu_Run2016C_skim_all.root',
				'SingleMu_Run2016D_skim_all.root',
				'SingleMu_Run2016E_skim_all.root',
				'SingleMu_Run2016F_skim_all.root',
				'SingleMu_Run2016G_skim_all.root',
				'SingleMu_Run2016Hv2_skim_all.root',
				'SingleMu_Run2016Hv3_skim_all.root']
el_data_fns = ['SingleEl_Run2016Bv2_skim_all.root',
				'SingleEl_Run2016C_skim_all.root',
				'SingleEl_Run2016D_skim_all.root',
				'SingleEl_Run2016E_skim_all.root',
				'SingleEl_Run2016F_skim_all.root',
				'SingleEl_Run2016G_skim_all.root',
				'SingleEl_Run2016Hv2_skim_all.root',
				'SingleEl_Run2016Hv3_skim_all.root']
#charge asymmetric background filenames
bkg_fns=['ST_s-c_skim_all.root',
		'ST_t-c_top_skim_all.root',
		'ST_t-c_antitop_skim_all.root',
		'ST_tW-c_top_skim_all.root',
		'ST_tW-c_antitop_skim_all.root',
		'WJets_HT-70to100_skim_all.root',
		'WJets_HT-100to200_skim_all.root',
		'WJets_HT-200to400_skim_all.root',
		'WJets_HT-400to600_skim_all.root',
		'WJets_HT-600to800_skim_all.root',
		'WJets_HT-800to1200_skim_all.root',
		'WJets_HT-1200to2500_skim_all.root',
		'WJets_HT-2500toInf_skim_all.root']
#selection & reweighting strings
skim_cut = 'fullselection==1'
muflav_cut = 'lepflavor==1'
elflav_cut = 'lepflavor==2'
lpos_cut = 'lep_Q>0'
lneg_cut = 'lep_Q<0'
t1top_cut = 'eventTopology==1'
t2top_cut = 'eventTopology==2'
t3top_cut = 'eventTopology==3'
def_weights  = '(((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_ttag_eff_merged*sf_ttag_eff_semimerged*sf_ttag_eff_notmerged*sf_btag_eff_heavy*sf_btag_eff_light*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'

#chain up files
mu_data_chain = TChain('tree')
for fn in mu_data_fns :
	mu_data_chain.Add(fn)
el_data_chain = TChain('tree')
for fn in el_data_fns :
	el_data_chain.Add(fn)
bkg_chain = TChain('tree')
for fn in bkg_fns :
	bkg_chain.Add(fn)

#skim chains
print 'skimming chains....'
print '	t1muplus_data: (%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut)
t1muplus_data_skimmed = mu_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut))
print '	t1muminus_data: (%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut)
t1muminus_data_skimmed = mu_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut))
print '	t1elplus_data: (%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut)
t1elplus_data_skimmed = el_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut))
print '	t1elminus_data: (%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut)
t1elminus_data_skimmed = el_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut))
print '	t2muplus_data: (%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut)
t2muplus_data_skimmed = mu_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut))
print '	t2muminus_data: (%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut)
t2muminus_data_skimmed = mu_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut))
print '	t2elplus_data: (%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut)
t2elplus_data_skimmed = el_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut))
print '	t2elminus_data: (%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut)
t2elminus_data_skimmed = el_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut))
print '	t3muplus_data: (%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut)
t3muplus_data_skimmed = mu_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut))
print '	t3muminus_data: (%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut)
t3muminus_data_skimmed = mu_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut))
print '	t3elplus_data: (%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut)
t3elplus_data_skimmed = el_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut))
print '	t3elminus_data: (%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut)
t3elminus_data_skimmed = el_data_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut))
print '	t1muplus_bkg: (%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut)
t1muplus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut))
print '	t1muminus_bkg: (%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut)
t1muminus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut))
print '	t1elplus_bkg: (%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut)
t1elplus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut))
print '	t1elminus_bkg: (%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut)
t1elminus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut))
print '	t2muplus_bkg: (%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut)
t2muplus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut))
print '	t2muminus_bkg: (%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut)
t2muminus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut))
print '	t2elplus_bkg: (%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut)
t2elplus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut))
print '	t2elminus_bkg: (%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut)
t2elminus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut))
print '	t3muplus_bkg: (%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut)
t3muplus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut))
print '	t3muminus_bkg: (%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut)
t3muminus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut))
print '	t3elplus_bkg: (%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut)
t3elplus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut))
print '	t3elminus_bkg: (%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut)
t3elminus_bkg_skimmed = bkg_chain.CopyTree('(%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut))
print 'Done.'

#get event numbers
print 'Drawing histograms....'
print '	t1muplus_data'
t1muplus_data_histo = TH1D('t1muplus_data_histo','',20,-1.,1.)
t1muplus_data_histo_f = TH1D('t1muplus_data_histo','',20,-1.,1.)
t1muplus_data_histo_b = TH1D('t1muplus_data_histo','',20,-1.,1.)
t1muplus_data_skimmed.Draw('cstar>>t1muplus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut))
t1muplus_data_skimmed.Draw('cstar>>t1muplus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut))
t1muplus_data_skimmed.Draw('cstar>>t1muplus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut))
t1muplus_data_histo.Add(gROOT.FindObject('t1muplus_data'))
t1muplus_data_histo_f.Add(gROOT.FindObject('t1muplus_data_f'))
t1muplus_data_histo_b.Add(gROOT.FindObject('t1muplus_data_b'))
t1muplus_data=t1muplus_data_histo.Integral()
t1muplus_data_f=t1muplus_data_histo_f.Integral()
t1muplus_data_b=t1muplus_data_histo_b.Integral()
print '	t1muplus_bkg'
t1muplus_bkg_histo = TH1D('t1muplus_bkg_histo','',20,-1.,1.)
t1muplus_bkg_histo_f = TH1D('t1muplus_bkg_histo_f','',20,-1.,1.)
t1muplus_bkg_histo_b = TH1D('t1muplus_bkg_histo_b','',20,-1.,1.)
t1muplus_bkg_skimmed.Draw('cstar>>t1muplus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut,def_weights))
t1muplus_bkg_skimmed.Draw('cstar>>t1muplus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut,def_weights))
t1muplus_bkg_skimmed.Draw('cstar>>t1muplus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t1top_cut,muflav_cut,lpos_cut,def_weights))
t1muplus_bkg_histo.Add(gROOT.FindObject('t1muplus_bkg'))
t1muplus_bkg_histo_f.Add(gROOT.FindObject('t1muplus_bkg_f'))
t1muplus_bkg_histo_b.Add(gROOT.FindObject('t1muplus_bkg_b'))
t1muplus_bkg=t1muplus_bkg_histo.Integral()
t1muplus_bkg_f=t1muplus_bkg_histo_f.Integral()
t1muplus_bkg_b=t1muplus_bkg_histo_b.Integral()
print '	t1muminus_data'
t1muminus_data_histo = TH1D('t1muminus_data_histo','',20,-1.,1.)
t1muminus_data_histo_f = TH1D('t1muminus_data_histo','',20,-1.,1.)
t1muminus_data_histo_b = TH1D('t1muminus_data_histo','',20,-1.,1.)
t1muminus_data_skimmed.Draw('cstar>>t1muminus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut))
t1muminus_data_skimmed.Draw('cstar>>t1muminus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut))
t1muminus_data_skimmed.Draw('cstar>>t1muminus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut))
t1muminus_data_histo.Add(gROOT.FindObject('t1muminus_data'))
t1muminus_data_histo_f.Add(gROOT.FindObject('t1muminus_data_f'))
t1muminus_data_histo_b.Add(gROOT.FindObject('t1muminus_data_b'))
t1muminus_data=t1muminus_data_histo.Integral()
t1muminus_data_f=t1muminus_data_histo_f.Integral()
t1muminus_data_b=t1muminus_data_histo_b.Integral()
print '	t1muminus_bkg'
t1muminus_bkg_histo = TH1D('t1muminus_bkg_histo','',20,-1.,1.)
t1muminus_bkg_histo_f = TH1D('t1muminus_bkg_histo_f','',20,-1.,1.)
t1muminus_bkg_histo_b = TH1D('t1muminus_bkg_histo_b','',20,-1.,1.)
t1muminus_bkg_skimmed.Draw('cstar>>t1muminus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut,def_weights))
t1muminus_bkg_skimmed.Draw('cstar>>t1muminus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut,def_weights))
t1muminus_bkg_skimmed.Draw('cstar>>t1muminus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t1top_cut,muflav_cut,lneg_cut,def_weights))
t1muminus_bkg_histo.Add(gROOT.FindObject('t1muminus_bkg'))
t1muminus_bkg_histo_f.Add(gROOT.FindObject('t1muminus_bkg_f'))
t1muminus_bkg_histo_b.Add(gROOT.FindObject('t1muminus_bkg_b'))
t1muminus_bkg=t1muminus_bkg_histo.Integral()
t1muminus_bkg_f=t1muminus_bkg_histo_f.Integral()
t1muminus_bkg_b=t1muminus_bkg_histo_b.Integral()
print '	t1elplus_data'
t1elplus_data_histo = TH1D('t1elplus_data_histo','',20,-1.,1.)
t1elplus_data_histo_f = TH1D('t1elplus_data_histo','',20,-1.,1.)
t1elplus_data_histo_b = TH1D('t1elplus_data_histo','',20,-1.,1.)
t1elplus_data_skimmed.Draw('cstar>>t1elplus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut))
t1elplus_data_skimmed.Draw('cstar>>t1elplus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut))
t1elplus_data_skimmed.Draw('cstar>>t1elplus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut))
t1elplus_data_histo.Add(gROOT.FindObject('t1elplus_data'))
t1elplus_data_histo_f.Add(gROOT.FindObject('t1elplus_data_f'))
t1elplus_data_histo_b.Add(gROOT.FindObject('t1elplus_data_b'))
t1elplus_data=t1elplus_data_histo.Integral()
t1elplus_data_f=t1elplus_data_histo_f.Integral()
t1elplus_data_b=t1elplus_data_histo_b.Integral()
print '	t1elplus_bkg'
t1elplus_bkg_histo = TH1D('t1elplus_bkg_histo','',20,-1.,1.)
t1elplus_bkg_histo_f = TH1D('t1elplus_bkg_histo_f','',20,-1.,1.)
t1elplus_bkg_histo_b = TH1D('t1elplus_bkg_histo_b','',20,-1.,1.)
t1elplus_bkg_skimmed.Draw('cstar>>t1elplus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut,def_weights))
t1elplus_bkg_skimmed.Draw('cstar>>t1elplus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut,def_weights))
t1elplus_bkg_skimmed.Draw('cstar>>t1elplus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t1top_cut,elflav_cut,lpos_cut,def_weights))
t1elplus_bkg_histo.Add(gROOT.FindObject('t1elplus_bkg'))
t1elplus_bkg_histo_f.Add(gROOT.FindObject('t1elplus_bkg_f'))
t1elplus_bkg_histo_b.Add(gROOT.FindObject('t1elplus_bkg_b'))
t1elplus_bkg=t1elplus_bkg_histo.Integral()
t1elplus_bkg_f=t1elplus_bkg_histo_f.Integral()
t1elplus_bkg_b=t1elplus_bkg_histo_b.Integral()
print '	t1elminus_data'
t1elminus_data_histo = TH1D('t1elminus_data_histo','',20,-1.,1.)
t1elminus_data_histo_f = TH1D('t1elminus_data_histo','',20,-1.,1.)
t1elminus_data_histo_b = TH1D('t1elminus_data_histo','',20,-1.,1.)
t1elminus_data_skimmed.Draw('cstar>>t1elminus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut))
t1elminus_data_skimmed.Draw('cstar>>t1elminus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut))
t1elminus_data_skimmed.Draw('cstar>>t1elminus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut))
t1elminus_data_histo.Add(gROOT.FindObject('t1elminus_data'))
t1elminus_data_histo_f.Add(gROOT.FindObject('t1elminus_data_f'))
t1elminus_data_histo_b.Add(gROOT.FindObject('t1elminus_data_b'))
t1elminus_data=t1elminus_data_histo.Integral()
t1elminus_data_f=t1elminus_data_histo_f.Integral()
t1elminus_data_b=t1elminus_data_histo_b.Integral()
print '	t1elminus_bkg'
t1elminus_bkg_histo = TH1D('t1elminus_bkg_histo','',20,-1.,1.)
t1elminus_bkg_histo_f = TH1D('t1elminus_bkg_histo_f','',20,-1.,1.)
t1elminus_bkg_histo_b = TH1D('t1elminus_bkg_histo_b','',20,-1.,1.)
t1elminus_bkg_skimmed.Draw('cstar>>t1elminus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut,def_weights))
t1elminus_bkg_skimmed.Draw('cstar>>t1elminus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut,def_weights))
t1elminus_bkg_skimmed.Draw('cstar>>t1elminus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t1top_cut,elflav_cut,lneg_cut,def_weights))
t1elminus_bkg_histo.Add(gROOT.FindObject('t1elminus_bkg'))
t1elminus_bkg_histo_f.Add(gROOT.FindObject('t1elminus_bkg_f'))
t1elminus_bkg_histo_b.Add(gROOT.FindObject('t1elminus_bkg_b'))
t1elminus_bkg=t1elminus_bkg_histo.Integral()
t1elminus_bkg_f=t1elminus_bkg_histo_f.Integral()
t1elminus_bkg_b=t1elminus_bkg_histo_b.Integral()
print '	t2muplus_data'
t2muplus_data_histo = TH1D('t2muplus_data_histo','',20,-1.,1.)
t2muplus_data_histo_f = TH1D('t2muplus_data_histo','',20,-1.,1.)
t2muplus_data_histo_b = TH1D('t2muplus_data_histo','',20,-1.,1.)
t2muplus_data_skimmed.Draw('cstar>>t2muplus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut))
t2muplus_data_skimmed.Draw('cstar>>t2muplus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut))
t2muplus_data_skimmed.Draw('cstar>>t2muplus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut))
t2muplus_data_histo.Add(gROOT.FindObject('t2muplus_data'))
t2muplus_data_histo_f.Add(gROOT.FindObject('t2muplus_data_f'))
t2muplus_data_histo_b.Add(gROOT.FindObject('t2muplus_data_b'))
t2muplus_data=t2muplus_data_histo.Integral()
t2muplus_data_f=t2muplus_data_histo_f.Integral()
t2muplus_data_b=t2muplus_data_histo_b.Integral()
print '	t2muplus_bkg'
t2muplus_bkg_histo = TH1D('t2muplus_bkg_histo','',20,-1.,1.)
t2muplus_bkg_histo_f = TH1D('t2muplus_bkg_histo_f','',20,-1.,1.)
t2muplus_bkg_histo_b = TH1D('t2muplus_bkg_histo_b','',20,-1.,1.)
t2muplus_bkg_skimmed.Draw('cstar>>t2muplus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut,def_weights))
t2muplus_bkg_skimmed.Draw('cstar>>t2muplus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut,def_weights))
t2muplus_bkg_skimmed.Draw('cstar>>t2muplus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t2top_cut,muflav_cut,lpos_cut,def_weights))
t2muplus_bkg_histo.Add(gROOT.FindObject('t2muplus_bkg'))
t2muplus_bkg_histo_f.Add(gROOT.FindObject('t2muplus_bkg_f'))
t2muplus_bkg_histo_b.Add(gROOT.FindObject('t2muplus_bkg_b'))
t2muplus_bkg=t2muplus_bkg_histo.Integral()
t2muplus_bkg_f=t2muplus_bkg_histo_f.Integral()
t2muplus_bkg_b=t2muplus_bkg_histo_b.Integral()
print '	t2muminus_data'
t2muminus_data_histo = TH1D('t2muminus_data_histo','',20,-1.,1.)
t2muminus_data_histo_f = TH1D('t2muminus_data_histo','',20,-1.,1.)
t2muminus_data_histo_b = TH1D('t2muminus_data_histo','',20,-1.,1.)
t2muminus_data_skimmed.Draw('cstar>>t2muminus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut))
t2muminus_data_skimmed.Draw('cstar>>t2muminus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut))
t2muminus_data_skimmed.Draw('cstar>>t2muminus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut))
t2muminus_data_histo.Add(gROOT.FindObject('t2muminus_data'))
t2muminus_data_histo_f.Add(gROOT.FindObject('t2muminus_data_f'))
t2muminus_data_histo_b.Add(gROOT.FindObject('t2muminus_data_b'))
t2muminus_data=t2muminus_data_histo.Integral()
t2muminus_data_f=t2muminus_data_histo_f.Integral()
t2muminus_data_b=t2muminus_data_histo_b.Integral()
print '	t2muminus_bkg'
t2muminus_bkg_histo = TH1D('t2muminus_bkg_histo','',20,-1.,1.)
t2muminus_bkg_histo_f = TH1D('t2muminus_bkg_histo_f','',20,-1.,1.)
t2muminus_bkg_histo_b = TH1D('t2muminus_bkg_histo_b','',20,-1.,1.)
t2muminus_bkg_skimmed.Draw('cstar>>t2muminus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut,def_weights))
t2muminus_bkg_skimmed.Draw('cstar>>t2muminus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut,def_weights))
t2muminus_bkg_skimmed.Draw('cstar>>t2muminus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t2top_cut,muflav_cut,lneg_cut,def_weights))
t2muminus_bkg_histo.Add(gROOT.FindObject('t2muminus_bkg'))
t2muminus_bkg_histo_f.Add(gROOT.FindObject('t2muminus_bkg_f'))
t2muminus_bkg_histo_b.Add(gROOT.FindObject('t2muminus_bkg_b'))
t2muminus_bkg=t2muminus_bkg_histo.Integral()
t2muminus_bkg_f=t2muminus_bkg_histo_f.Integral()
t2muminus_bkg_b=t2muminus_bkg_histo_b.Integral()
print '	t2elplus_data'
t2elplus_data_histo = TH1D('t2elplus_data_histo','',20,-1.,1.)
t2elplus_data_histo_f = TH1D('t2elplus_data_histo','',20,-1.,1.)
t2elplus_data_histo_b = TH1D('t2elplus_data_histo','',20,-1.,1.)
t2elplus_data_skimmed.Draw('cstar>>t2elplus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut))
t2elplus_data_skimmed.Draw('cstar>>t2elplus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut))
t2elplus_data_skimmed.Draw('cstar>>t2elplus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut))
t2elplus_data_histo.Add(gROOT.FindObject('t2elplus_data'))
t2elplus_data_histo_f.Add(gROOT.FindObject('t2elplus_data_f'))
t2elplus_data_histo_b.Add(gROOT.FindObject('t2elplus_data_b'))
t2elplus_data=t2elplus_data_histo.Integral()
t2elplus_data_f=t2elplus_data_histo_f.Integral()
t2elplus_data_b=t2elplus_data_histo_b.Integral()
print '	t2elplus_bkg'
t2elplus_bkg_histo = TH1D('t2elplus_bkg_histo','',20,-1.,1.)
t2elplus_bkg_histo_f = TH1D('t2elplus_bkg_histo_f','',20,-1.,1.)
t2elplus_bkg_histo_b = TH1D('t2elplus_bkg_histo_b','',20,-1.,1.)
t2elplus_bkg_skimmed.Draw('cstar>>t2elplus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut,def_weights))
t2elplus_bkg_skimmed.Draw('cstar>>t2elplus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut,def_weights))
t2elplus_bkg_skimmed.Draw('cstar>>t2elplus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t2top_cut,elflav_cut,lpos_cut,def_weights))
t2elplus_bkg_histo.Add(gROOT.FindObject('t2elplus_bkg'))
t2elplus_bkg_histo_f.Add(gROOT.FindObject('t2elplus_bkg_f'))
t2elplus_bkg_histo_b.Add(gROOT.FindObject('t2elplus_bkg_b'))
t2elplus_bkg=t2elplus_bkg_histo.Integral()
t2elplus_bkg_f=t2elplus_bkg_histo_f.Integral()
t2elplus_bkg_b=t2elplus_bkg_histo_b.Integral()
print '	t2elminus_data'
t2elminus_data_histo = TH1D('t2elminus_data_histo','',20,-1.,1.)
t2elminus_data_histo_f = TH1D('t2elminus_data_histo','',20,-1.,1.)
t2elminus_data_histo_b = TH1D('t2elminus_data_histo','',20,-1.,1.)
t2elminus_data_skimmed.Draw('cstar>>t2elminus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut))
t2elminus_data_skimmed.Draw('cstar>>t2elminus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut))
t2elminus_data_skimmed.Draw('cstar>>t2elminus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut))
t2elminus_data_histo.Add(gROOT.FindObject('t2elminus_data'))
t2elminus_data_histo_f.Add(gROOT.FindObject('t2elminus_data_f'))
t2elminus_data_histo_b.Add(gROOT.FindObject('t2elminus_data_b'))
t2elminus_data=t2elminus_data_histo.Integral()
t2elminus_data_f=t2elminus_data_histo_f.Integral()
t2elminus_data_b=t2elminus_data_histo_b.Integral()
print '	t2elminus_bkg'
t2elminus_bkg_histo = TH1D('t2elminus_bkg_histo','',20,-1.,1.)
t2elminus_bkg_histo_f = TH1D('t2elminus_bkg_histo_f','',20,-1.,1.)
t2elminus_bkg_histo_b = TH1D('t2elminus_bkg_histo_b','',20,-1.,1.)
t2elminus_bkg_skimmed.Draw('cstar>>t2elminus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut,def_weights))
t2elminus_bkg_skimmed.Draw('cstar>>t2elminus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut,def_weights))
t2elminus_bkg_skimmed.Draw('cstar>>t2elminus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t2top_cut,elflav_cut,lneg_cut,def_weights))
t2elminus_bkg_histo.Add(gROOT.FindObject('t2elminus_bkg'))
t2elminus_bkg_histo_f.Add(gROOT.FindObject('t2elminus_bkg_f'))
t2elminus_bkg_histo_b.Add(gROOT.FindObject('t2elminus_bkg_b'))
t2elminus_bkg=t2elminus_bkg_histo.Integral()
t2elminus_bkg_f=t2elminus_bkg_histo_f.Integral()
t2elminus_bkg_b=t2elminus_bkg_histo_b.Integral()
print '	t3muplus_data'
t3muplus_data_histo = TH1D('t3muplus_data_histo','',20,-1.,1.)
t3muplus_data_histo_f = TH1D('t3muplus_data_histo','',20,-1.,1.)
t3muplus_data_histo_b = TH1D('t3muplus_data_histo','',20,-1.,1.)
t3muplus_data_skimmed.Draw('cstar>>t3muplus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut))
t3muplus_data_skimmed.Draw('cstar>>t3muplus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut))
t3muplus_data_skimmed.Draw('cstar>>t3muplus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut))
t3muplus_data_histo.Add(gROOT.FindObject('t3muplus_data'))
t3muplus_data_histo_f.Add(gROOT.FindObject('t3muplus_data_f'))
t3muplus_data_histo_b.Add(gROOT.FindObject('t3muplus_data_b'))
t3muplus_data=t3muplus_data_histo.Integral()
t3muplus_data_f=t3muplus_data_histo_f.Integral()
t3muplus_data_b=t3muplus_data_histo_b.Integral()
print '	t3muplus_bkg'
t3muplus_bkg_histo = TH1D('t3muplus_bkg_histo','',20,-1.,1.)
t3muplus_bkg_histo_f = TH1D('t3muplus_bkg_histo_f','',20,-1.,1.)
t3muplus_bkg_histo_b = TH1D('t3muplus_bkg_histo_b','',20,-1.,1.)
t3muplus_bkg_skimmed.Draw('cstar>>t3muplus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut,def_weights))
t3muplus_bkg_skimmed.Draw('cstar>>t3muplus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut,def_weights))
t3muplus_bkg_skimmed.Draw('cstar>>t3muplus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t3top_cut,muflav_cut,lpos_cut,def_weights))
t3muplus_bkg_histo.Add(gROOT.FindObject('t3muplus_bkg'))
t3muplus_bkg_histo_f.Add(gROOT.FindObject('t3muplus_bkg_f'))
t3muplus_bkg_histo_b.Add(gROOT.FindObject('t3muplus_bkg_b'))
t3muplus_bkg=t3muplus_bkg_histo.Integral()
t3muplus_bkg_f=t3muplus_bkg_histo_f.Integral()
t3muplus_bkg_b=t3muplus_bkg_histo_b.Integral()
print '	t3muminus_data'
t3muminus_data_histo = TH1D('t3muminus_data_histo','',20,-1.,1.)
t3muminus_data_histo_f = TH1D('t3muminus_data_histo','',20,-1.,1.)
t3muminus_data_histo_b = TH1D('t3muminus_data_histo','',20,-1.,1.)
t3muminus_data_skimmed.Draw('cstar>>t3muminus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut))
t3muminus_data_skimmed.Draw('cstar>>t3muminus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut))
t3muminus_data_skimmed.Draw('cstar>>t3muminus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut))
t3muminus_data_histo.Add(gROOT.FindObject('t3muminus_data'))
t3muminus_data_histo_f.Add(gROOT.FindObject('t3muminus_data_f'))
t3muminus_data_histo_b.Add(gROOT.FindObject('t3muminus_data_b'))
t3muminus_data=t3muminus_data_histo.Integral()
t3muminus_data_f=t3muminus_data_histo_f.Integral()
t3muminus_data_b=t3muminus_data_histo_b.Integral()
print '	t3muminus_bkg'
t3muminus_bkg_histo = TH1D('t3muminus_bkg_histo','',20,-1.,1.)
t3muminus_bkg_histo_f = TH1D('t3muminus_bkg_histo_f','',20,-1.,1.)
t3muminus_bkg_histo_b = TH1D('t3muminus_bkg_histo_b','',20,-1.,1.)
t3muminus_bkg_skimmed.Draw('cstar>>t3muminus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut,def_weights))
t3muminus_bkg_skimmed.Draw('cstar>>t3muminus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut,def_weights))
t3muminus_bkg_skimmed.Draw('cstar>>t3muminus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t3top_cut,muflav_cut,lneg_cut,def_weights))
t3muminus_bkg_histo.Add(gROOT.FindObject('t3muminus_bkg'))
t3muminus_bkg_histo_f.Add(gROOT.FindObject('t3muminus_bkg_f'))
t3muminus_bkg_histo_b.Add(gROOT.FindObject('t3muminus_bkg_b'))
t3muminus_bkg=t3muminus_bkg_histo.Integral()
t3muminus_bkg_f=t3muminus_bkg_histo_f.Integral()
t3muminus_bkg_b=t3muminus_bkg_histo_b.Integral()
print '	t3elplus_data'
t3elplus_data_histo = TH1D('t3elplus_data_histo','',20,-1.,1.)
t3elplus_data_histo_f = TH1D('t3elplus_data_histo','',20,-1.,1.)
t3elplus_data_histo_b = TH1D('t3elplus_data_histo','',20,-1.,1.)
t3elplus_data_skimmed.Draw('cstar>>t3elplus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut))
t3elplus_data_skimmed.Draw('cstar>>t3elplus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut))
t3elplus_data_skimmed.Draw('cstar>>t3elplus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut))
t3elplus_data_histo.Add(gROOT.FindObject('t3elplus_data'))
t3elplus_data_histo_f.Add(gROOT.FindObject('t3elplus_data_f'))
t3elplus_data_histo_b.Add(gROOT.FindObject('t3elplus_data_b'))
t3elplus_data=t3elplus_data_histo.Integral()
t3elplus_data_f=t3elplus_data_histo_f.Integral()
t3elplus_data_b=t3elplus_data_histo_b.Integral()
print '	t3elplus_bkg'
t3elplus_bkg_histo = TH1D('t3elplus_bkg_histo','',20,-1.,1.)
t3elplus_bkg_histo_f = TH1D('t3elplus_bkg_histo_f','',20,-1.,1.)
t3elplus_bkg_histo_b = TH1D('t3elplus_bkg_histo_b','',20,-1.,1.)
t3elplus_bkg_skimmed.Draw('cstar>>t3elplus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut,def_weights))
t3elplus_bkg_skimmed.Draw('cstar>>t3elplus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut,def_weights))
t3elplus_bkg_skimmed.Draw('cstar>>t3elplus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t3top_cut,elflav_cut,lpos_cut,def_weights))
t3elplus_bkg_histo.Add(gROOT.FindObject('t3elplus_bkg'))
t3elplus_bkg_histo_f.Add(gROOT.FindObject('t3elplus_bkg_f'))
t3elplus_bkg_histo_b.Add(gROOT.FindObject('t3elplus_bkg_b'))
t3elplus_bkg=t3elplus_bkg_histo.Integral()
t3elplus_bkg_f=t3elplus_bkg_histo_f.Integral()
t3elplus_bkg_b=t3elplus_bkg_histo_b.Integral()
print '	t3elminus_data'
t3elminus_data_histo = TH1D('t3elminus_data_histo','',20,-1.,1.)
t3elminus_data_histo_f = TH1D('t3elminus_data_histo','',20,-1.,1.)
t3elminus_data_histo_b = TH1D('t3elminus_data_histo','',20,-1.,1.)
t3elminus_data_skimmed.Draw('cstar>>t3elminus_data(20,-1.,1.)','(%s && %s && %s && %s)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut))
t3elminus_data_skimmed.Draw('cstar>>t3elminus_data_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut))
t3elminus_data_skimmed.Draw('cstar>>t3elminus_data_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut))
t3elminus_data_histo.Add(gROOT.FindObject('t3elminus_data'))
t3elminus_data_histo_f.Add(gROOT.FindObject('t3elminus_data_f'))
t3elminus_data_histo_b.Add(gROOT.FindObject('t3elminus_data_b'))
t3elminus_data=t3elminus_data_histo.Integral()
t3elminus_data_f=t3elminus_data_histo_f.Integral()
t3elminus_data_b=t3elminus_data_histo_b.Integral()
print '	t3elminus_bkg'
t3elminus_bkg_histo = TH1D('t3elminus_bkg_histo','',20,-1.,1.)
t3elminus_bkg_histo_f = TH1D('t3elminus_bkg_histo_f','',20,-1.,1.)
t3elminus_bkg_histo_b = TH1D('t3elminus_bkg_histo_b','',20,-1.,1.)
t3elminus_bkg_skimmed.Draw('cstar>>t3elminus_bkg(20,-1.,1.)','(%s && %s && %s && %s)*(%s)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut,def_weights))
t3elminus_bkg_skimmed.Draw('cstar>>t3elminus_bkg_f(20,-1.,1.)','(%s && %s && %s && %s && cstar>0)*(%s)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut,def_weights))
t3elminus_bkg_skimmed.Draw('cstar>>t3elminus_bkg_b(20,-1.,1.)','(%s && %s && %s && %s && cstar<0)*(%s)'%(skim_cut,t3top_cut,elflav_cut,lneg_cut,def_weights))
t3elminus_bkg_histo.Add(gROOT.FindObject('t3elminus_bkg'))
t3elminus_bkg_histo_f.Add(gROOT.FindObject('t3elminus_bkg_f'))
t3elminus_bkg_histo_b.Add(gROOT.FindObject('t3elminus_bkg_b'))
t3elminus_bkg=t3elminus_bkg_histo.Integral()
t3elminus_bkg_f=t3elminus_bkg_histo_f.Integral()
t3elminus_bkg_b=t3elminus_bkg_histo_b.Integral()

#calculate and print lepton asymmetries
print '############### charge asymmetries ###############'
t1mu_assym = ((t1muplus_data-t1muplus_bkg)-(t1muminus_data-t1muminus_bkg))/((t1muplus_data-t1muplus_bkg)+(t1muminus_data-t1muminus_bkg))
t1mu_assymerr = (2./((t1muplus_data-t1muplus_bkg)+(t1muminus_data-t1muminus_bkg))**2)*sqrt((t1muplus_data+t1muplus_bkg)*(t1muminus_data-t1muminus_bkg)**2+(t1muminus_data+t1muminus_bkg)*(t1muplus_data-t1muplus_bkg)**2)
print 'type-1 muon charge asymmetry = %.6f +/- %.6f'%(t1mu_assym,t1mu_assymerr)
t1el_assym = ((t1elplus_data-t1elplus_bkg)-(t1elminus_data-t1elminus_bkg))/((t1elplus_data-t1elplus_bkg)+(t1elminus_data-t1elminus_bkg))
t1el_assymerr = (2./((t1elplus_data-t1elplus_bkg)+(t1elminus_data-t1elminus_bkg))**2)*sqrt((t1elplus_data+t1elplus_bkg)*(t1elminus_data-t1elminus_bkg)**2+(t1elminus_data+t1elminus_bkg)*(t1elplus_data-t1elplus_bkg)**2)
print 'type-1 electron charge asymmetry = %.6f +/- %.6f'%(t1el_assym,t1el_assymerr)
t2mu_assym = ((t2muplus_data-t2muplus_bkg)-(t2muminus_data-t2muminus_bkg))/((t2muplus_data-t2muplus_bkg)+(t2muminus_data-t2muminus_bkg))
t2mu_assymerr = (2./((t2muplus_data-t2muplus_bkg)+(t2muminus_data-t2muminus_bkg))**2)*sqrt((t2muplus_data+t2muplus_bkg)*(t2muminus_data-t2muminus_bkg)**2+(t2muminus_data+t2muminus_bkg)*(t2muplus_data-t2muplus_bkg)**2)
print 'type-2 muon charge asymmetry = %.6f +/- %.6f'%(t2mu_assym,t2mu_assymerr)
t2el_assym = ((t2elplus_data-t2elplus_bkg)-(t2elminus_data-t2elminus_bkg))/((t2elplus_data-t2elplus_bkg)+(t2elminus_data-t2elminus_bkg))
t2el_assymerr = (2./((t2elplus_data-t2elplus_bkg)+(t2elminus_data-t2elminus_bkg))**2)*sqrt((t2elplus_data+t2elplus_bkg)*(t2elminus_data-t2elminus_bkg)**2+(t2elminus_data+t2elminus_bkg)*(t2elplus_data-t2elplus_bkg)**2)
print 'type-2 electron charge asymmetry = %.6f +/- %.6f'%(t2el_assym,t2el_assymerr)
t3mu_assym = ((t3muplus_data-t3muplus_bkg)-(t3muminus_data-t3muminus_bkg))/((t3muplus_data-t3muplus_bkg)+(t3muminus_data-t3muminus_bkg))
t3mu_assymerr = (2./((t3muplus_data-t3muplus_bkg)+(t3muminus_data-t3muminus_bkg))**2)*sqrt((t3muplus_data+t3muplus_bkg)*(t3muminus_data-t3muminus_bkg)**2+(t3muminus_data+t3muminus_bkg)*(t3muplus_data-t3muplus_bkg)**2)
print 'type-3 muon charge asymmetry = %.6f +/- %.6f'%(t3mu_assym,t3mu_assymerr)
t3el_assym = ((t3elplus_data-t3elplus_bkg)-(t3elminus_data-t3elminus_bkg))/((t3elplus_data-t3elplus_bkg)+(t3elminus_data-t3elminus_bkg))
t3el_assymerr = (2./((t3elplus_data-t3elplus_bkg)+(t3elminus_data-t3elminus_bkg))**2)*sqrt((t3elplus_data+t3elplus_bkg)*(t3elminus_data-t3elminus_bkg)**2+(t3elminus_data+t3elminus_bkg)*(t3elplus_data-t3elplus_bkg)**2)
print 'type-3 electron charge asymmetry = %.6f +/- %.6f'%(t3el_assym,t3el_assymerr)

#calculate and print fb asymmetries
print '############### fb asymmetries ###############'
t1muf_data = t1muplus_data_f+t1muminus_data_f; t1mub_data = t1muplus_data_b+t1muminus_data_b;
t1muf_bkg = t1muplus_bkg_f+t1muminus_bkg_f; t1mub_bkg = t1muplus_bkg_b+t1muminus_bkg_b;
t1mu_fbassym = ((t1muf_data-t1muf_bkg)-(t1mub_data-t1mub_bkg))/((t1muf_data-t1muf_bkg)+(t1mub_data-t1mub_bkg))
t1mu_fbassymerr = (2./((t1muf_data-t1muf_bkg)+(t1mub_data-t1mub_bkg))**2)*sqrt((t1muf_data+t1muf_bkg)*(t1mub_data-t1mub_bkg)**2+(t1mub_data+t1mub_bkg)*(t1muf_data-t1muf_bkg)**2)
print 'type-1 muon fb asymmetry = %.6f +/- %.6f'%(t1mu_fbassym,t1mu_fbassymerr)
t1elf_data = t1elplus_data_f+t1elminus_data_f; t1elb_data = t1elplus_data_b+t1elminus_data_b;
t1elf_bkg = t1elplus_bkg_f+t1elminus_bkg_f; t1elb_bkg = t1elplus_bkg_b+t1elminus_bkg_b;
t1el_fbassym = ((t1elf_data-t1elf_bkg)-(t1elb_data-t1elb_bkg))/((t1elf_data-t1elf_bkg)+(t1elb_data-t1elb_bkg))
t1el_fbassymerr = (2./((t1elf_data-t1elf_bkg)+(t1elb_data-t1elb_bkg))**2)*sqrt((t1elf_data+t1elf_bkg)*(t1elb_data-t1elb_bkg)**2+(t1elb_data+t1elb_bkg)*(t1elf_data-t1elf_bkg)**2)
print 'type-1 electron fb asymmetry = %.6f +/- %.6f'%(t1el_fbassym,t1el_fbassymerr)
t2muf_data = t2muplus_data_f+t2muminus_data_f; t2mub_data = t2muplus_data_b+t2muminus_data_b;
t2muf_bkg = t2muplus_bkg_f+t2muminus_bkg_f; t2mub_bkg = t2muplus_bkg_b+t2muminus_bkg_b;
t2mu_fbassym = ((t2muf_data-t2muf_bkg)-(t2mub_data-t2mub_bkg))/((t2muf_data-t2muf_bkg)+(t2mub_data-t2mub_bkg))
t2mu_fbassymerr = (2./((t2muf_data-t2muf_bkg)+(t2mub_data-t2mub_bkg))**2)*sqrt((t2muf_data+t2muf_bkg)*(t2mub_data-t2mub_bkg)**2+(t2mub_data+t2mub_bkg)*(t2muf_data-t2muf_bkg)**2)
print 'type-2 muon fb asymmetry = %.6f +/- %.6f'%(t2mu_fbassym,t2mu_fbassymerr)
t2elf_data = t2elplus_data_f+t2elminus_data_f; t2elb_data = t2elplus_data_b+t2elminus_data_b;
t2elf_bkg = t2elplus_bkg_f+t2elminus_bkg_f; t2elb_bkg = t2elplus_bkg_b+t2elminus_bkg_b;
t2el_fbassym = ((t2elf_data-t2elf_bkg)-(t2elb_data-t2elb_bkg))/((t2elf_data-t2elf_bkg)+(t2elb_data-t2elb_bkg))
t2el_fbassymerr = (2./((t2elf_data-t2elf_bkg)+(t2elb_data-t2elb_bkg))**2)*sqrt((t2elf_data+t2elf_bkg)*(t2elb_data-t2elb_bkg)**2+(t2elb_data+t2elb_bkg)*(t2elf_data-t2elf_bkg)**2)
print 'type-2 electron fb asymmetry = %.6f +/- %.6f'%(t2el_fbassym,t2el_fbassymerr)
t3muf_data = t3muplus_data_f+t3muminus_data_f; t3mub_data = t3muplus_data_b+t3muminus_data_b;
t3muf_bkg = t3muplus_bkg_f+t3muminus_bkg_f; t3mub_bkg = t3muplus_bkg_b+t3muminus_bkg_b;
t3mu_fbassym = ((t3muf_data-t3muf_bkg)-(t3mub_data-t3mub_bkg))/((t3muf_data-t3muf_bkg)+(t3mub_data-t3mub_bkg))
t3mu_fbassymerr = (2./((t3muf_data-t3muf_bkg)+(t3mub_data-t3mub_bkg))**2)*sqrt((t3muf_data+t3muf_bkg)*(t3mub_data-t3mub_bkg)**2+(t3mub_data+t3mub_bkg)*(t3muf_data-t3muf_bkg)**2)
print 'type-3 muon fb asymmetry = %.6f +/- %.6f'%(t3mu_fbassym,t3mu_fbassymerr)
t3elf_data = t3elplus_data_f+t3elminus_data_f; t3elb_data = t3elplus_data_b+t3elminus_data_b;
t3elf_bkg = t3elplus_bkg_f+t3elminus_bkg_f; t3elb_bkg = t3elplus_bkg_b+t3elminus_bkg_b;
t3el_fbassym = ((t3elf_data-t3elf_bkg)-(t3elb_data-t3elb_bkg))/((t3elf_data-t3elf_bkg)+(t3elb_data-t3elb_bkg))
t3el_fbassymerr = (2./((t3elf_data-t3elf_bkg)+(t3elb_data-t3elb_bkg))**2)*sqrt((t3elf_data+t3elf_bkg)*(t3elb_data-t3elb_bkg)**2+(t3elb_data+t3elb_bkg)*(t3elf_data-t3elf_bkg)**2)
print 'type-3 electron fb asymmetry = %.6f +/- %.6f'%(t3el_fbassym,t3el_fbassymerr)

print '############### total A_FB effects ###############'
t1mu_AFB_effect = t1mu_assym*t1mu_fbassym
t1mu_AFB_effecterr = t1mu_AFB_effect*sqrt((t1mu_assymerr/t1mu_assym)**2+(t1mu_fbassymerr/t1mu_fbassym)**2)
print 'type-1 muon fb asymmetry impact = %.6f +/- %.6f'%(t1mu_AFB_effect,t1mu_AFB_effecterr)
t1el_AFB_effect = t1el_assym*t1el_fbassym
t1el_AFB_effecterr = t1el_AFB_effect*sqrt((t1el_assymerr/t1el_assym)**2+(t1el_fbassymerr/t1el_fbassym)**2)
print 'type-1 electron fb asymmetry impact = %.6f +/- %.6f'%(t1el_AFB_effect,t1el_AFB_effecterr)
t2mu_AFB_effect = t2mu_assym*t2mu_fbassym
t2mu_AFB_effecterr = t2mu_AFB_effect*sqrt((t2mu_assymerr/t2mu_assym)**2+(t2mu_fbassymerr/t2mu_fbassym)**2)
print 'type-2 muon fb asymmetry impact = %.6f +/- %.6f'%(t2mu_AFB_effect,t2mu_AFB_effecterr)
t2el_AFB_effect = t2el_assym*t2el_fbassym
t2el_AFB_effecterr = t2el_AFB_effect*sqrt((t2el_assymerr/t2el_assym)**2+(t2el_fbassymerr/t2el_fbassym)**2)
print 'type-2 electron fb asymmetry impact = %.6f +/- %.6f'%(t2el_AFB_effect,t2el_AFB_effecterr)
t3mu_AFB_effect = t3mu_assym*t3mu_fbassym
t3mu_AFB_effecterr = t3mu_AFB_effect*sqrt((t3mu_assymerr/t3mu_assym)**2+(t3mu_fbassymerr/t3mu_fbassym)**2)
print 'type-3 muon fb asymmetry impact = %.6f +/- %.6f'%(t3mu_AFB_effect,t3mu_AFB_effecterr)
t3el_AFB_effect = t3el_assym*t3el_fbassym
t3el_AFB_effecterr = t3el_AFB_effect*sqrt((t3el_assymerr/t3el_assym)**2+(t3el_fbassymerr/t3el_fbassym)**2)
print 'type-3 electron fb asymmetry impact = %.6f +/- %.6f'%(t3el_AFB_effect,t3el_AFB_effecterr)


