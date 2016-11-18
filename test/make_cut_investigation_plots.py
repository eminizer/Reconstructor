from ROOT import *
import glob

#output file name
OUTPUT_FILE_NAME = 'cut_investigation_plots.root'

#sample class
class Sample(object) :
	def __init__(self,uniqueName_=None,sType_='Signal',color_=kBlack,filePrep_='../',treeName_='tree') :
		self._uniqueName = uniqueName_
		self._sType = sType_
		self._color = color_
		filelist = glob.glob(filePrep_+uniqueName_+'/'+uniqueName_+'_*_tree.root')
		toremove = []
		for filename in filelist :
			if filename.find('JER_up')!=-1 or filename.find('JER_down')!=-1 or filename.find('JES_up')!=-1 or filename.find('JES_down')!=-1 or filename.find('_skim_')!=-1 :
				toremove.append(filename)
		for r in toremove :
			filelist.pop(filelist.index(r))
		self._chain = TChain(treeName_)
		for filename in filelist :
			self._chain.Add(filename)

	def getUniqueName(self) :
		return self._uniqueName
	def getSType(self) :
		return self._sType
	def getColor(self) :
		return self._color
	def getChain(self) :
		return self._chain

#plot class
class Plot(object) :
	def __init__(self,name_,title_,bintuple_,plotStrings_,cutStrings_=[],weightStrings_=[],optionStrings_=[]) :
		print 'Making plots with name '+name_
		#add defaults
		if optionStrings_==[] : 
			for p in plotStrings_ : 
				optionStrings_.append('HIST')
		if weightStrings_==[] : 
			for p in plotStrings_ : 
				weightStrings_.append('12917.*weight')
		if cutStrings_==[] : 
			for p in plotStrings_ : 
				cutStrings_.append('weight!=0.')
		#make sure the input doesn't suck
		if not (len(plotStrings_) == len(cutStrings_) == len(weightStrings_) == len(optionStrings_)) :
			print 'lists for plot named %s have different lengths, sorry, doublecheck!'%(name_)
			return
		self._allCanvases = []; self._allHistos = []
		#separated backgrounds
		#plot the event yields
		print '	Doing separated backgrounds, event yields plot'
		self._stype_histos_sep_ev = []
		self._allCanvases.append(TCanvas(name_+'_sep_ev_c',name_+' separated backgrounds, event yields'))
		self._sep_bck_leg = TLegend(0.5,0.1,0.9,0.5)
		for sType in allSTypes :
			print '		Doing sample type %s'%(sType)
			stypehistoname = ''
			for word in sType.split(' ') : 
				stypehistoname+=word.lower()+'_'
			stypehistoname+='h'
			stype_event_yield_histo = TH1D(stypehistoname,title_+'; Events',bintuple_[0],bintuple_[1],bintuple_[2])
			stype_event_yield_histo.SetDirectory(0)
			self._stype_histos_sep_ev.append(stype_event_yield_histo); self._allHistos.append(stype_event_yield_histo)
			for sample in allSamples :
				if not sample.getSType()==sType : continue
				color = sample.getColor()
				chain = sample.getChain()
				for i in range(len(plotStrings_)) :
					chain.Draw('('+plotStrings_[i]+')>>'+sample.getUniqueName()+'_'+str(i)+'('+str(bintuple_[0])+','+str(bintuple_[1])+','+str(bintuple_[2])+')','('+weightStrings_[i]+')*('+cutStrings_[i]+')',optionStrings_[i])
					stype_event_yield_histo.Add(gROOT.FindObject(sample.getUniqueName()+'_'+str(i)))
			stype_event_yield_histo.SetTitle(title_+'; Events')
			stype_event_yield_histo.SetLineWidth(3)
			stype_event_yield_histo.SetLineColor(color)
			self._sep_bck_leg.AddEntry(stype_event_yield_histo,sType,'L')
		self._allCanvases[len(self._allCanvases)-1].cd()
		self._stype_histos_sep_ev[0].Draw(optionStrings_[0])
		for i in range(1,len(self._stype_histos_sep_ev)) :
			self._stype_histos_sep_ev[i].Draw(optionStrings_[0]+' SAME')
		self._sep_bck_leg.Draw()
		#plot the fractional plots
		print '	Doing separated backgrounds, shapes plot'
		self._stype_histos_sep_fr = []
		self._allCanvases.append(TCanvas(name_+'_sep_fr_c',name_+' separated backgrounds, shapes'))
		for sType in allSTypes :
			print '		Doing sample type %s'%(sType)
			stypehistoname = ''
			for word in sType.split(' ') : 
				stypehistoname+=word.lower()+'_'
			stypehistoname+='h'
			stype_fractional_histo = TH1D(stypehistoname,title_+'; Fraction',bintuple_[0],bintuple_[1],bintuple_[2])
			stype_fractional_histo.SetDirectory(0)
			self._stype_histos_sep_fr.append(stype_fractional_histo); self._allHistos.append(stype_fractional_histo)
			for sample in allSamples :
				if not sample.getSType()==sType : continue
				color = sample.getColor()
				chain = sample.getChain()
				for i in range(len(plotStrings_)) :
					chain.Draw('('+plotStrings_[i]+')>>'+sample.getUniqueName()+'_'+str(i)+'('+str(bintuple_[0])+','+str(bintuple_[1])+','+str(bintuple_[2])+')','('+weightStrings_[i]+')*('+cutStrings_[i]+')',optionStrings_[i])
					stype_fractional_histo.Add(gROOT.FindObject(sample.getUniqueName()+'_'+str(i)))
			stype_fractional_histo.SetTitle(title_+'; Fraction')
			stype_fractional_histo.SetLineWidth(3)
			stype_fractional_histo.SetLineColor(color)
			stype_fractional_histo.Scale(1./stype_fractional_histo.Integral())
		self._allCanvases[len(self._allCanvases)-1].cd()
		self._stype_histos_sep_fr[0].Draw(optionStrings_[0])
		for i in range(1,len(self._stype_histos_sep_fr)) :
			self._stype_histos_sep_fr[i].Draw('SAME')
		self._sep_bck_leg.Draw()
		#combined backgrounds
			#plot the event yields
			#plot the fractional plots

	def getCanvases(self) :
		return self._allCanvases
	def getHistoList(self) :
		return self._allHistos

#countplot class
#plot class
class CountPlot(object) :
	def __init__(self,name_=None,title_=None,bintuple_=None,selecStrings_=[],weightStrings_=[],optionStrings_=[]) :
		#add defaults
		if optionStrings_==[] : 
			for s in selecStrings_ : 
				optionStrings_.append('HIST')
		if weightStrings_==[] : 
			for s in selecStrings_ : 
				weightStrings_.append('12917.*weight')
		#make sure the input doesn't suck
		if not (len(selecStrings_) == len(weightStrings_) == len(optionStrings_)) :
			print 'lists for counting plot named %s have different lengths, sorry, doublecheck!'%(name_)
			return
		#separated backgrounds
			#plot the event yields
			#plot the fractional plots
		#combined backgrounds
			#plot the event yields
			#plot the fractional plots

#2D plot class
class Plot2D(object) :
	def __init__(self,name_=None,title_=None,bintuple_=None,plotStrings_=None,cutStrings_=[],weightStrings_=[],optionStrings_=[]) :
		#add defaults
		if optionStrings_==[] : 
			for p in plotStrings_ : 
				optionStrings_.append('HIST')
		if weightStrings_==[] : 
			for p in plotStrings_ : 
				weightStrings_.append('12917.*weight')
		if cutStrings_==[] : 
			for p in plotStrings_ : 
				cutStrings_.append('weight!=0.')
		#make sure the input doesn't suck
		if not (len(plotStrings_) == len(cutStrings_) == len(weightStrings_) == len(optionStrings_)) :
			print 'lists for plot named %s have different lengths, sorry, doublecheck!'%(name_)
			return
		#separated backgrounds
			#plot the event yields
			#plot the fractional plots
		#combined backgrounds
			#plot the event yields
			#plot the fractional plots

#Declare Samples
allSamples = []
#signal
allSamples.append(Sample(uniqueName_='mcatnlo_qq_semilep_TT'))
allSamples.append(Sample(uniqueName_='mcatnlo_gg_semilep_TT'))
#dileptonic background
allSamples.append(Sample(uniqueName_='mcatnlo_dilep_TT',sType_='Dileptonic t#bar{t}',color_=kRed-7))
#hadronic background
allSamples.append(Sample(uniqueName_='mcatnlo_had_TT',sType_='Hadronic t#bar{t}',color_=kRed-5))
#single top background
allSamples.append(Sample(uniqueName_='ST_s-c',sType_='Single Top',color_=kMagenta))
allSamples.append(Sample(uniqueName_='ST_t-c_top',sType_='Single Top',color_=kMagenta))
allSamples.append(Sample(uniqueName_='ST_tW-c_top',sType_='Single Top',color_=kMagenta))
allSamples.append(Sample(uniqueName_='ST_t-c_antitop',sType_='Single Top',color_=kMagenta))
allSamples.append(Sample(uniqueName_='ST_tW-c_antitop',sType_='Single Top',color_=kMagenta))
#WJets background
allSamples.append(Sample(uniqueName_='WJets_HT-200to400',sType_='W + Jets',color_=kGreen-3))
allSamples.append(Sample(uniqueName_='WJets_HT-400to600',sType_='W + Jets',color_=kGreen-3))
allSamples.append(Sample(uniqueName_='WJets_HT-600to800',sType_='W + Jets',color_=kGreen-3))
allSamples.append(Sample(uniqueName_='WJets_HT-800to1200',sType_='W + Jets',color_=kGreen-3))
#allSamples.append(Sample(uniqueName_='WJets_HT-2500toInf',sType_='W + Jets',color_=kGreen-3))
#DYJets background
allSamples.append(Sample(uniqueName_='DYJets_M-50_HT-100to200',sType_='Z/#gamma + Jets',color_=kAzure-2))
allSamples.append(Sample(uniqueName_='DYJets_M-50_HT-200to400',sType_='Z/#gamma + Jets',color_=kAzure-2))
allSamples.append(Sample(uniqueName_='DYJets_M-50_HT-400to600',sType_='Z/#gamma + Jets',color_=kAzure-2))
allSamples.append(Sample(uniqueName_='DYJets_M-50_HT-600toInf',sType_='Z/#gamma + Jets',color_=kAzure-2))
#Diboson background
allSamples.append(Sample(uniqueName_='WW_to_2L_2Nu',sType_='Diboson',color_=kOrange-3))
allSamples.append(Sample(uniqueName_='WW_to_L_Nu_2Q',sType_='Diboson',color_=kOrange-3))
allSamples.append(Sample(uniqueName_='WZ_to_L_Nu_2Q',sType_='Diboson',color_=kOrange-3))
allSamples.append(Sample(uniqueName_='WZ_to_L_3Nu',sType_='Diboson',color_=kOrange-3))
allSamples.append(Sample(uniqueName_='WZ_to_2L_2Q',sType_='Diboson',color_=kOrange-3))
allSamples.append(Sample(uniqueName_='WZ_to_3L_Nu',sType_='Diboson',color_=kOrange-3))
allSamples.append(Sample(uniqueName_='ZZ_to_2L_2Nu',sType_='Diboson',color_=kOrange-3))
allSamples.append(Sample(uniqueName_='ZZ_to_2L_2Q',sType_='Diboson',color_=kOrange-3))
allSamples.append(Sample(uniqueName_='ZZ_to_4L',sType_='Diboson',color_=kOrange-3))
#aggregate the unique background types
allSTypes = []
for sample in allSamples :
	thistype = sample.getSType()
	if thistype!=None and thistype not in allSTypes :
		allSTypes.append(thistype)

#Set up output file
outfile = TFile(OUTPUT_FILE_NAME,'recreate')

#cutstrings to use
metfilters = 'metfilters==1'
trigger = 'trigger==1'
onelep = 'onelepton==1'
isolep = 'isolepton==1'
jetcuts = 'jetcuts==1'
lepHT = '(scaled_lep_pt + met_pt)>150.'
elemet = 'met_pt>120.'
mumet  = 'met_pt>50.'
othercuts = '(lepflavor==1 && (muon1_pt + met_pt)>150. && met_pt>50.) || (lepflavor==2 && met_pt>120.)'
fullvanilla = 'fullselection==1'

#Declare Plots
allPlots = []
#hardest AK4 jet pt
allPlots.append(Plot(name_='ak41_pt_all_other',
					 title_='hardest AK4 jet p_{T}, all other selection; p_{T} (GeV)',
					 bintuple_=(60,0.,1500.),
					 plotStrings_=['ak41_pt'],
					 cutStrings_=['('+metfilters+') && ('+trigger+') && ('+onelep+') && ('+isolep+') && ('+othercuts+') && ak42_pt>50.']))
#second-hardest AK4 jet pt
allPlots.append(Plot(name_='ak42_pt_all_other',
					 title_='second hardest AK4 jet p_{T}, all other selection; p_{T} (GeV)',
					 bintuple_=(60,0.,1500.),
					 plotStrings_=['ak42_pt'],
					 cutStrings_=['('+metfilters+') && ('+trigger+') && ('+onelep+') && ('+isolep+') && ('+othercuts+') && ak41_pt>150.']))
##leptonic HT
#allPlots.append(Plot(name_='lepside_HT_all_other',
#					 title_='leptonic side HT, all other selection; HT (GeV)',
#					 bintuple_=(60,0.,1500.),
#					 plotStrings_=['muon1_pt + met_pt', 
#					 			   'ele1_pt  + met_pt'],
#					 cutStrings_=['('+metfilters+') && ('+trigger+') && ('+onelep+') && ('+isolep+') && ('+jetcuts+') && lepflavor==1 && met_pt>50.',
#					 			  '('+metfilters+') && ('+trigger+') && ('+onelep+') && ('+isolep+') && ('+jetcuts+') && lepflavor==2 && met_pt>120.']))
##MET
#allPlots.append(Plot(name_='MET_all_other',
#					 title_='missing E_{T}, all other selection; E_{T}^{miss} (GeV)',
#					 bintuple_=(60,0.,1500.),
#					 plotStrings_=['muon1_pt + met_pt', 'ele1_pt + met_pt'],
#					 cutStrings_=['('+metfilters+') && ('+trigger+') && ('+onelep+') && ('+isolep+') && ('+jetcuts+') && lepflavor==1 && (muon1_pt+met_pt)>150.',
#					 			  '('+metfilters+') && ('+trigger+') && ('+onelep+') && ('+isolep+') && ('+jetcuts+') && lepflavor==2']))
##lepton multiplicity
#lepmult_stdcuts = '('+metfilters+') && ('+trigger+') && ('+isolep+') && ('+jetcuts+') && ('+othercuts+')'
#allPlots.append(CountPlot(name_='lep_mult_all_other',
#						  title_='lepton multiplicity, all other selection; lepton multiplicity',
#						  selecStrings_=['((muon1_pt>0. && muon2_pt<0. && ele1_pt<0. && ele2_pt<0.) || (ele1_pt>0. && ele2_pt<0. && muon1_pt<0. && muon2_pt<0.)) && '+lepmult_stdcuts],
#						  				['((muon1_pt>0. && muon2_pt>0. && ele1_pt<0. && ele2_pt<0.) || (muon1_pt>0. && muon2_pt<0. && ele1_pt>0. && ele2_pt<0.) || (muon1_pt<0. && muon2_pt<0. && ele1_pt>0. && ele2_pt>0.)) && '+lepmult_stdcuts],
#						  				['((muon1_pt>0. && muon2_pt>0. && ele1_pt>0. && ele2_pt<0.) || (muon1_pt>0. && muon2_pt<0. && ele1_pt>0. && ele2_pt>0.)) && '+lepmult_stdcuts],
#						  				['(muon1_pt>0. && muon2_pt>0. && ele1_pt>0. && ele2_pt>0.) && '+lepmult_stdcuts]))
##lepton isolation
#allPlots.append(Plot2D(name_='lep_iso_all_other',
#					   title_='lepton isolation cut, all other selection; dR(lep,jet); p_{T}^{rel}(lep,jet) (GeV)',
#					   bintuple_=(25,0.,2.5,25,0.,100.),
#					   plotStrings_=['muon1_relPt:muon1_dR','ele1_relPt:ele1_dR'],
#					   cutStrings_=['('+metfilters+') && ('+trigger+') && ('+onelep+') && ('+jetcuts+') && ('+othercuts+') && lepflavor==1',
#					 			  '('+metfilters+') && ('+trigger+') && ('+onelep+') && ('+jetlep+') && ('+othercuts+') && lepflavor==2']))

#Save plots
outfile.cd()
for plot in allPlots :
	for canv in plot.getCanvases() :
		canv.Write()
for plot in allPlots :
	for histo in plot.getHistoList() :
		histo.Write()
outfile.Close()