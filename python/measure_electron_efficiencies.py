#Trigger and ID/isolation efficiency calculator designed to measure efficiency 
#of the HLT_Ele45_WPLoose_Gsf trigger path and electron medium ID (no iso)
#using SingleMu or SingleEl data and ttbar MC events

##########								   Imports  								##########

from ROOT import *
from corrector import Corrector
from optparse import OptionParser
from array import array
from math import *
from glob import glob
import sys

#Global variables
##Histogram bins
ele_pt_bins 	 = array('d',[0.,40.,50.,95.,130.,165.,210.,300.,500.])
eta_bins 		 = array('d',[-2.5,-1.,-0.5,0.,0.5,0.75,1.,2.5])
ele_pt_bins_2D 	 = array('d',[50.,95.,130.,165.,210.,300.])
eta_bins_2D 	 = array('d',[0.,0.5,0.75,1.,2.5])
pu_bins 		 = array('d',[0.,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5,32.5,34.5,36.5,38.5,40.5,42.5,44.5])
n_ele_pt_bins 	 = len(ele_pt_bins)-1
n_eta_bins 		 = len(eta_bins)-1
n_ele_pt_bins_2D = len(ele_pt_bins_2D)-1
n_eta_bins_2D 	 = len(eta_bins_2D)-1
n_pu_bins 		 = len(pu_bins)-1
#constants
LUMINOSITYBTOFELE = 19171.010
LUMINOSITYGHELE   = 16214.862
LUMINOSITYBTOFMU  = 19690.184
LUMINOSITYGHMU 	  = 16226.452

##########								Parser Options								##########

parser = OptionParser()
#Run options
parser.add_option('--mode', 	  type='string', action='store', dest='mode',	   	  
	help='What would you like to measure with this run? ("trigger/trig", "ID", or "iso[lation]")')
parser.add_option('--ttbar_generator', 	  type='string', action='store', default='powheg', dest='ttbar_generator',	   	  
	help='Use "powheg" or "mcatnlo" ttbar MC?')
parser.add_option('--topology', 	  type='string', action='store', default='boosted', dest='topology',	   	  
	help='Measure for "boosted" or "resolved" events?')
parser.add_option('--use_preskimmed_files', type='string', action='store', default='no', dest='usepreskim')
parser.add_option('--print_every',type='int',    action='store', default=1000000,	  dest='print_every', 
	help='Print progress after how many events?')
parser.add_option('--max_events', type='int',    action='store', default=-1,	  dest='max_events',  
	help='Maximum number of events to process (default is -1 for "all")')
(options, args) = parser.parse_args()

#process some input options
generator = options.ttbar_generator.lower()
mode = options.mode.lower()
if mode in ['trig','trigger'] :
	mode = 't'
elif mode in ['id'] :
	mode = 'id'
elif mode in ['iso','isolation'] :
	mode = 'iso'
print 'Running electron %s efficiencies'%('trigger' if mode=='t' else mode)

##########							Sample List 								##########

samplenames = []

##Multiboson
#samplenames.append('WW_to_L_Nu_2Q')
#samplenames.append('WW_to_2L_2Nu')
#samplenames.append('WZ_to_L_Nu_2Q')
#samplenames.append('WZ_to_L_3Nu')
#samplenames.append('WZ_to_2L_2Q')
#samplenames.append('WZ_to_3L_Nu')
#samplenames.append('ZZ_to_2L_2Nu')
#samplenames.append('ZZ_to_2L_2Q')
#samplenames.append('ZZ_to_4L')
##QCD
#samplenames.append('QCD_HT-100to200')
#samplenames.append('QCD_HT-200to300')
#samplenames.append('QCD_HT-300to500')
#samplenames.append('QCD_HT-500to700')
#samplenames.append('QCD_HT-700to1000')
#samplenames.append('QCD_HT-1000to1500')
#samplenames.append('QCD_HT-1500to2000')
#samplenames.append('QCD_HT-2000toInf')
##WJets
#samplenames.append('WJets_HT-200to400')
#samplenames.append('WJets_HT-400to600')
#samplenames.append('WJets_HT-600to800')
#samplenames.append('WJets_HT-800to1200')
#samplenames.append('WJets_HT-1200to2500')
#samplenames.append('WJets_HT-2500toInf')
##DYJets
#samplenames.append('DYJets_M-50_HT-70to100')
#samplenames.append('DYJets_M-50_HT-100to200')
#samplenames.append('DYJets_M-50_HT-200to400')
#samplenames.append('DYJets_M-50_HT-400to600')
#samplenames.append('DYJets_M-50_HT-600to800')
#samplenames.append('DYJets_M-50_HT-800to1200')
#samplenames.append('DYJets_M-50_HT-1200to2500')
#samplenames.append('DYJets_M-50_HT-2500toInf')
##Single top
#samplenames.append('ST_s-c')
#samplenames.append('ST_t-c_top')
#samplenames.append('ST_t-c_antitop')
#samplenames.append('ST_tW-c_top')
#samplenames.append('ST_tW-c_antitop')
#POWHEG TT
if generator == 'powheg' :
	samplenames.append('powheg_TT')
#MCATNLO TT
elif generator == 'mcatnlo' :
	samplenames.append('mcatnlo_TT')
else :
	print 'generator type '+generator+' not recognized! Exiting.'
	exit()
#data
data_samplenames = []
#if mode=='t' :
#	samplenames.append('SingleMu_Run2016Bv2')
#	samplenames.append('SingleMu_Run2016C')
#	samplenames.append('SingleMu_Run2016D')
#	samplenames.append('SingleMu_Run2016E')
#	samplenames.append('SingleMu_Run2016F')
#	samplenames.append('SingleMu_Run2016G')
#	samplenames.append('SingleMu_Run2016Hv2')
#	samplenames.append('SingleMu_Run2016Hv3')
if mode=='t' or mode=='id' or mode=='iso' :
	samplenames.append('SingleEl_Run2016Bv2')
	samplenames.append('SingleEl_Run2016C')
	samplenames.append('SingleEl_Run2016D')
	samplenames.append('SingleEl_Run2016E')
	samplenames.append('SingleEl_Run2016F')
	samplenames.append('SingleEl_Run2016G')
	samplenames.append('SingleEl_Run2016Hv2')
	samplenames.append('SingleEl_Run2016Hv3')
else :
	print 'mode '+mode+' not recognized! Exiting.'
	exit()

##########							Filegroup Class 							##########

class filegroup : 
	#docstring
	"""filegroup class"""
	
	#__init__function
	def __init__(self,name,chain,outputfile,mode) :
		self.name = name
		isdata = name.find('SingleEl')!=-1 or name.find('SingleMu')!=-1
		self.chain = chain
		self.mode=mode
		self.corrector = Corrector(isdata,'no',TH1F('pileup','pileup',100,0.,100.),'B',{'alpha':1.0,'epsilon':1.0})
		#Get relevant branches
		#Physics objects
		self.phyobjbs = {}
		themuon_pt 	  = array('f',[-1.0]);  self.chain.SetBranchAddress('eltrig_themuon_pt',themuon_pt); 	  self.phyobjbs['themuon_pt'] = (themuon_pt,-1.0)
		themuon_eta   = array('f',[100.0]); self.chain.SetBranchAddress('eltrig_themuon_eta',themuon_eta); 	  self.phyobjbs['themuon_eta'] = (themuon_eta,100.0)
		theele_pt 	  = array('f',[-1.0]);  self.chain.SetBranchAddress('eltrig_theelectron_pt',theele_pt);   self.phyobjbs['theele_pt'] = (theele_pt,-1.0)
		theele_eta 	  = array('f',[100.0]); self.chain.SetBranchAddress('eltrig_theelectron_eta',theele_eta); self.phyobjbs['theele_eta'] = (theele_eta,100.0)
		tag_pt 		  = array('f',[-1.0]);  self.chain.SetBranchAddress('elID_tag_pt',tag_pt); 				  self.phyobjbs['tag_pt'] = (tag_pt,-1.0)
		tag_eta 	  = array('f',[100.0]); self.chain.SetBranchAddress('elID_tag_eta',tag_eta); 			  self.phyobjbs['tag_eta'] = (tag_eta,100.0)
		probe_pt 	  = array('f',[-1.0]);  self.chain.SetBranchAddress('elID_probe_pt',probe_pt); 			  self.phyobjbs['probe_pt'] = (probe_pt,-1.0)
		probe_eta 	  = array('f',[100.0]); self.chain.SetBranchAddress('elID_probe_eta',probe_eta); 		  self.phyobjbs['probe_eta'] = (probe_eta,100.0)
		evt_pu 		  = array('i',[-1]); 	self.chain.SetBranchAddress('ngoodvtx',evt_pu); 				  self.phyobjbs['evt_pu'] = (evt_pu,-1)
		evt_top 	  = array('i',[0]); 	self.chain.SetBranchAddress('eventTopology',evt_top); 			  self.phyobjbs['evt_top'] = (evt_top,0)
		nbTags 		  = array('i',[0]); 	self.chain.SetBranchAddress('nbTags',nbTags); 					  self.phyobjbs['nbTags'] = (nbTags,0)
		lepflavor 	  = array('I',[0]); 	self.chain.SetBranchAddress('lepflavor',lepflavor); 			  self.phyobjbs['lepflavor'] = (lepflavor,0)
		#selection and trigger variables
		self.selecpassbs = {}
		selectiontrig = array('I',[2]); self.chain.SetBranchAddress('eltrig_fullselection',selectiontrig); self.selecpassbs['selectiontrig'] = (selectiontrig,2)
		passtrig 	  = array('I',[2]); self.chain.SetBranchAddress('eltrig_eltrigger',passtrig); 		   self.selecpassbs['passtrig'] = (passtrig,2)
		selectionID   = array('I',[2]); self.chain.SetBranchAddress('elID_fullselection',selectionID); 	   self.selecpassbs['selectionID'] = (selectionID,2)
		passID 		  = array('I',[2]); self.chain.SetBranchAddress('elID_probeID',passID); 			   self.selecpassbs['passID'] = (passID,2)
		selectioniso  = array('I',[2]); self.chain.SetBranchAddress('elID_fullselection',selectioniso);    self.selecpassbs['selectioniso'] = (selectioniso,2)
		passiso 	  = array('I',[2]); self.chain.SetBranchAddress('elID_isoprobe',passiso); 			   self.selecpassbs['passiso'] = (passiso,2)
		jets_cut 	  = array('I',[2]); self.chain.SetBranchAddress('ak4jetcuts',jets_cut); 			   self.selecpassbs['jets_cut'] = (jets_cut,2)
		#reweighting factors we can read from the files without recalculating
		self.rwbs = {}
		weight 		  = array('f',[1.0]);  self.chain.SetBranchAddress('weight',weight); 			   self.rwbs['weight'] = (weight,1.0)
		sf_pileup 	  = array('f',[1.0]);  self.chain.SetBranchAddress('sf_pileup',sf_pileup); 		   self.rwbs['sf_pileup'] = (sf_pileup,1.0)
		sf_btag_eff   = array('f',[1.0]);  self.chain.SetBranchAddress('sf_btag_eff',sf_btag_eff); 	   self.rwbs['sf_btag_eff'] = (sf_btag_eff,1.0)
		sf_mu_R 	  = array('f',[1.0]);  self.chain.SetBranchAddress('sf_mu_R',sf_mu_R);			   self.rwbs['sf_mu_R'] = (sf_mu_R,1.0)
		sf_mu_F 	  = array('f',[1.0]);  self.chain.SetBranchAddress('sf_mu_F',sf_mu_F); 			   self.rwbs['sf_mu_F'] = (sf_mu_F,1.0)
		sf_scale_comb = array('f',[1.0]);  self.chain.SetBranchAddress('sf_scale_comb',sf_scale_comb); self.rwbs['sf_scale_comb'] = (sf_scale_comb,1.0)
		sf_pdf_alphas = array('f',[1.0]);  self.chain.SetBranchAddress('sf_pdf_alphas',sf_pdf_alphas); self.rwbs['sf_pdf_alphas'] = (sf_pdf_alphas,1.0)
		#Set up output TTree
		self.tree  = TTree('tree','tree')
		self.otbs = {}
		el_pt 	   = array('f',[-1.0]); self.tree.Branch('el_pt',el_pt,'el_pt/F'); self.otbs['el_pt'] = (el_pt,-1.0)
		el_eta 	   = array('f',[100.]); self.tree.Branch('el_eta',el_eta,'el_eta/F'); self.otbs['el_eta'] = (el_eta,100.)
		pu 		   = array('i',[-1]);   self.tree.Branch('pu',pu,'pu/I'); self.otbs['pu'] = (pu,-1)
		topology   = array('i',[0]); 	self.tree.Branch('topology',topology,'topology/I'); self.otbs['topology'] = (topology,0)
		pass_trig  = array('I',[2]);    self.tree.Branch('pass_trig',pass_trig,'pass_trig/i'); self.otbs['pass_trig'] = (pass_trig,2)
		pass_ID    = array('I',[2]);    self.tree.Branch('pass_ID',pass_ID,'pass_ID/i'); self.otbs['pass_ID'] = (pass_ID,2)
		pass_iso   = array('I',[2]);    self.tree.Branch('pass_iso',pass_iso,'pass_iso/i'); self.otbs['pass_iso'] = (pass_iso,2)
		evt_weight = array('f',[1.0]); 	self.tree.Branch('evt_weight',evt_weight,'evt_weight/F'); self.otbs['evt_weight'] = (evt_weight,1.0)
		self.allbranchdicts = [self.phyobjbs,self.selecpassbs,self.rwbs,self.otbs]
		#Set up histograms and efficiency graphs
		self.histos_and_graphs = []
		self.ele_pt_all   = TH1D(self.name+'_ele_pt_all',self.name+' electron p_{T} for all events; p_{T} (GeV)',n_ele_pt_bins,ele_pt_bins) 
		self.histos_and_graphs.append(self.ele_pt_all)
		self.ele_eta_all  = TH1D(self.name+'_ele_eta_all',self.name+' electron #eta for all events; #eta',n_eta_bins,eta_bins) 
		self.histos_and_graphs.append(self.ele_eta_all)
		self.evt_pu_all   = TH1D(self.name+'_evt_pu_all',self.name+' pileup for all events; # vertices',n_pu_bins,pu_bins) 
		self.histos_and_graphs.append(self.evt_pu_all)
		self.ele_pt_pass  = TH1D(self.name+'_ele_pt_pass',self.name+' electron p_{T} for passing events; p_{T} (GeV)',n_ele_pt_bins,ele_pt_bins) 
		self.histos_and_graphs.append(self.ele_pt_pass)
		self.ele_eta_pass = TH1D(self.name+'_ele_eta_pass',self.name+' electron #eta for passing events; #eta',n_eta_bins,eta_bins) 
		self.histos_and_graphs.append(self.ele_eta_pass)
		self.evt_pu_pass  = TH1D(self.name+'_evt_pu_pass',self.name+' pileup for passing events; # vertices',n_pu_bins,pu_bins) 
		self.histos_and_graphs.append(self.evt_pu_pass)
		self.histo_2d_all  = TH2D(self.name+'_histo_2D_all','',n_ele_pt_bins_2D,ele_pt_bins_2D,n_eta_bins_2D,eta_bins_2D); self.histos_and_graphs.append(self.histo_2d_all)
		self.histo_2d_pass = TH2D(self.name+'_histo_2D_pass','',n_ele_pt_bins_2D,ele_pt_bins_2D,n_eta_bins_2D,eta_bins_2D); self.histos_and_graphs.append(self.histo_2d_pass)
		ele_pt_x   = array('d',n_ele_pt_bins*[0.])
		ele_pt_xe  = array('d',n_ele_pt_bins*[0.])
		ele_pt_y   = array('d',n_ele_pt_bins*[0.])
		ele_pt_ye  = array('d',n_ele_pt_bins*[0.])
		ele_eta_x  = array('d',n_eta_bins*[0.])
		ele_eta_xe = array('d',n_eta_bins*[0.])
		ele_eta_y  = array('d',n_eta_bins*[0.])
		ele_eta_ye = array('d',n_eta_bins*[0.])
		pu_x 	   = array('d',n_pu_bins*[0.])
		pu_xe 	   = array('d',n_pu_bins*[0.])
		pu_y 	   = array('d',n_pu_bins*[0.])
		pu_ye 	   = array('d',n_pu_bins*[0.])
		self.ele_pt_gr    = TGraphErrors(n_ele_pt_bins,ele_pt_x,ele_pt_y,ele_pt_xe,ele_pt_ye); self.histos_and_graphs.append(self.ele_pt_gr)
		self.ele_pt_gr.SetName(self.name+'ele_pt_gr'); self.ele_pt_gr.SetTitle(self.name+' probe efficiency vs. electron p_{T}'); 
		self.ele_pt_gr.GetXaxis().SetName('electron p_{T} (GeV)'); self.ele_pt_gr.GetYaxis().SetName('Probe efficiency')
		self.ele_eta_gr   = TGraphErrors(n_eta_bins,ele_eta_x,ele_eta_y,ele_eta_xe,ele_eta_ye); self.histos_and_graphs.append(self.ele_eta_gr)
		self.ele_eta_gr.SetName(self.name+'ele_eta_gr'); self.ele_eta_gr.SetTitle(self.name+' probe efficiency vs. electron #eta'); 
		self.ele_eta_gr.GetXaxis().SetName('#eta'); self.ele_eta_gr.GetYaxis().SetName('Probe efficiency')
		self.pu_gr   = TGraphErrors(n_pu_bins,pu_x,pu_y,pu_xe,pu_ye); self.histos_and_graphs.append(self.pu_gr)
		self.pu_gr.SetName(self.name+'pu_gr'); self.pu_gr.SetTitle(self.name+' probe efficiency vs.pileup'); 
		self.pu_gr.GetXaxis().SetName('pileup'); self.pu_gr.GetYaxis().SetName('Probe efficiency')
		#Counter
		count = 0
		##########								Main Event Loop								##########
		print 'Filling trees for '+self.name+'. . .'
		nEntries = chain.GetEntries()
		for entry in range(nEntries) :
			#check the max events
			count+=1
			if count == options.max_events+1 :
				print 'Processed event number '+str(count-1)+', exiting'
				break
			#print progress
			if count % options.print_every == 0 or count == 1:
				print 'Count at '+str(count)+' out of '+str(nEntries)+', (%.4f%% complete)'%(float(count) / float(nEntries) * 100.0)
			
			chain.GetEntry(entry)
			
			cuts = []
			cuts.append(self.name.lower().find('data')!=-1 or weight[0]!=1.0)
			#cuts only dependent on topology
			if options.topology=='boosted' :
				cuts.append(evt_top[0]<3)
				cuts.append(nbTags[0]>0)
			elif options.topology=='resolved' :
				cuts.append(evt_top[0]==3)
				cuts.append(nbTags[0]>1)
			#cuts for trigger analysis
			if mode=='t' :
				cuts.append(selectiontrig[0]==1)
				cuts.append(lepflavor[0]==2)
			#cuts for ID analysis
			elif mode=='id' :
				#print 'selection = %d'%(selectionID[0]) #DEBUG
				cuts.append(selectionID[0]==1)
				cuts.append(lepflavor[0]==2)
			#cuts for isolation analysis
			elif mode=='iso' :
				cuts.append(selectioniso[0]==1)
			#miscellaneous
			cuts.append(jets_cut[0]==1)
			#check all cuts
			if cuts.count(False) > 0 :
				continue
			
			#fill the tree
			el_pt[0] = theele_pt[0] 
			el_eta[0] = theele_eta[0]
			pu[0] = evt_pu[0]
			topology[0] = evt_top[0]
			pass_trig[0] = passtrig[0]
			pass_ID[0] 	 = passID[0]
			pass_iso[0]  = passiso[0]
			evt_weight[0] = self.__getEvtWeight__()
			self.tree.Fill()
			
			#fill the histograms (bring parameters in range)
			self.ele_pt_all.Fill(max(min(el_pt[0],ele_pt_bins[n_ele_pt_bins-1]-0.0001),ele_pt_bins[0]+0.0001),evt_weight[0])
			self.ele_eta_all.Fill(max(min(el_eta[0],eta_bins[n_eta_bins-1]-0.0001),eta_bins[0]+0.0001),evt_weight[0])
			self.evt_pu_all.Fill(max(min(pu[0],pu_bins[n_pu_bins-1]-0.0001),pu_bins[0]+0.0001),evt_weight[0])
			self.histo_2d_all.Fill(max(min(el_pt[0],ele_pt_bins_2D[n_ele_pt_bins_2D-1]-0.0001),ele_pt_bins_2D[0]+0.0001),max(min(abs(el_eta[0]),eta_bins_2D[n_eta_bins_2D-1]-0.0001),eta_bins_2D[0]+0.0001),evt_weight[0])
			if mode=='id': #DEBUG
				print 'pass_ID[0]=%d'%(pass_ID[0]) #DEBUG
			if (mode=='t' and pass_trig[0]==1) or (mode=='id' and pass_ID[0]==1) or (mode=='iso' and pass_iso[0]==1) :
				self.ele_pt_pass.Fill(max(min(el_pt[0],ele_pt_bins[n_ele_pt_bins-1]-0.0001),ele_pt_bins[0]+0.0001),evt_weight[0])
				self.ele_eta_pass.Fill(max(min(el_eta[0],eta_bins[n_eta_bins-1]-0.0001),eta_bins[0]+0.0001),evt_weight[0])
				self.evt_pu_pass.Fill(max(min(pu[0],pu_bins[n_pu_bins-1]-0.0001),pu_bins[0]+0.0001),evt_weight[0])
				self.histo_2d_pass.Fill(max(min(el_pt[0],ele_pt_bins_2D[n_ele_pt_bins_2D-1]-0.0001),ele_pt_bins_2D[0]+0.0001),max(min(abs(el_eta[0]),eta_bins_2D[n_eta_bins_2D-1]-0.0001),eta_bins_2D[0]+0.0001),evt_weight[0])
		
			#reset all of the arrays holding tree branches
			for branchdict in self.allbranchdicts :
				for branchtuple in branchdict.values() :
					branchtuple[0][0] = branchtuple[1]

		print 'Done'
		
		#Make the graph y-values and errors
		for i in range(n_ele_pt_bins) :
			x_value = (ele_pt_bins[i+1]+ele_pt_bins[i])/2
			x_err   = (ele_pt_bins[i+1]-ele_pt_bins[i])/2
			passing_events = self.ele_pt_pass.GetBinContent(self.ele_pt_pass.FindBin(x_value))
			all_events = self.ele_pt_all.GetBinContent(self.ele_pt_all.FindBin(x_value))
			y_value = 0.; y_err = 0.
			if all_events > 0. :
				y_value = passing_events/all_events
				if passing_events>0. :
					y_err = y_value*sqrt(1./passing_events+1./all_events)
			self.ele_pt_gr.SetPoint(i,x_value,y_value)
			self.ele_pt_gr.SetPointError(i,x_err,y_err)
		for i in range(n_eta_bins) :
			x_value = (eta_bins[i+1]+eta_bins[i])/2
			x_err   = (eta_bins[i+1]-eta_bins[i])/2
			passing_events = self.ele_eta_pass.GetBinContent(self.ele_eta_pass.FindBin(x_value))
			all_events = self.ele_eta_all.GetBinContent(self.ele_eta_all.FindBin(x_value))
			y_value = 0.; y_err = 0.
			if all_events > 0. :
				y_value = passing_events/all_events
				if passing_events>0. :
					y_err = y_value*sqrt(1./passing_events+1./all_events)
			self.ele_eta_gr.SetPoint(i,x_value,y_value)
			self.ele_eta_gr.SetPointError(i,x_err,y_err)
		for i in range(n_pu_bins) :
			x_value = (pu_bins[i+1]+pu_bins[i])/2
			x_err   = (pu_bins[i+1]-pu_bins[i])/2
			passing_events = self.evt_pu_pass.GetBinContent(self.evt_pu_pass.FindBin(x_value))
			all_events = self.evt_pu_all.GetBinContent(self.evt_pu_all.FindBin(x_value))
			y_value = 0.; y_err = 0.
			if all_events > 0. :
				y_value = passing_events/all_events
				if passing_events>0. :
					y_err = y_value*sqrt(1./passing_events+1./all_events)
			self.pu_gr.SetPoint(i,x_value,y_value)
			self.pu_gr.SetPointError(i,x_err,y_err)
		#Fit the graphs with constants, then save the fit functions as graphs
		self.ele_pt_const = TF1('ele_pt_const','[0]',ele_pt_bins[0],ele_pt_bins[n_ele_pt_bins])
		self.ele_eta_const = TF1('ele_eta_const','[0]',eta_bins[0],eta_bins[n_eta_bins])
		self.pu_const = TF1('pu_const','[0]',pu_bins[0],pu_bins[n_pu_bins])
		self.ele_pt_gr.Fit('ele_pt_const')
		self.ele_eta_gr.Fit('ele_eta_const')
		self.pu_gr.Fit('pu_const')
		self.ele_pt_fit_value  = self.ele_pt_const.GetParameter(0)
		self.ele_pt_fit_err    = self.ele_pt_const.GetParError(0)
		self.ele_eta_fit_value = self.ele_eta_const.GetParameter(0)
		self.ele_eta_fit_err   = self.ele_eta_const.GetParError(0)
		self.pu_fit_value 	   = self.pu_const.GetParameter(0)
		self.pu_fit_err 	   = self.pu_const.GetParError(0)
		self.ele_pt_fit_gr 	= TGraphErrors(n_ele_pt_bins,ele_pt_x,ele_pt_y,ele_pt_xe,ele_pt_ye)
		self.ele_eta_fit_gr = TGraphErrors(n_eta_bins,ele_eta_x,ele_eta_y,ele_eta_xe,ele_eta_ye)
		self.pu_fit_gr 		= TGraphErrors(n_pu_bins,pu_x,pu_y,pu_xe,pu_ye)
		for i in range(n_ele_pt_bins) :
			x_value = (ele_pt_bins[i+1]+ele_pt_bins[i])/2
			x_err   = (ele_pt_bins[i+1]-ele_pt_bins[i])/2
			self.ele_pt_fit_gr.SetPoint(i,x_value,self.ele_pt_fit_value)
			self.ele_pt_fit_gr.SetPointError(i,x_err,self.ele_pt_fit_err)
		for i in range(n_eta_bins) :
			x_value = (eta_bins[i+1]+eta_bins[i])/2
			x_err   = (eta_bins[i+1]-eta_bins[i])/2
			self.ele_eta_fit_gr.SetPoint(i,x_value,self.ele_eta_fit_value)
			self.ele_eta_fit_gr.SetPointError(i,x_err,self.ele_eta_fit_err)
		for i in range(n_pu_bins) :
			x_value = (pu_bins[i+1]+pu_bins[i])/2
			x_err   = (pu_bins[i+1]-pu_bins[i])/2
			self.pu_fit_gr.SetPoint(i,x_value,self.pu_fit_value)
			self.pu_fit_gr.SetPointError(i,x_err,self.pu_fit_err)
		#Write the tree, histograms, and graphs
		outputfile.cd()
		self.tree.Write()
		for thing in self.histos_and_graphs :
			thing.Write()	

	def __getEvtWeight__(self) :
		#print '---------- getting event weight ----------' #DEBUG
		#if it's data just return "1"
		if self.name.lower().find('data')!=-1 :
		#	print 'data event' #DEBUG
			return 1.0
		#start with the stuff from the file
		#print 'weight=%.8f'%(self.rwbs['weight'][0][0]) #DEBUG
		#if self.rwbs['sf_pileup'][0][0]<0.5 : #DEBUG
		#	print 'pileup=%d, sf_pileup=%.8f'%(self.phyobjbs['evt_top'][0][0],self.rwbs['sf_pileup'][0][0]) #DEBUG
		#print 'sf_btag_eff=%.8f'%(self.rwbs['sf_btag_eff'][0][0]) #DEBUG
		#print 'sf_mu_R=%.8f, sf_mu_F=%.8f, sf_scale_comb=%.8f, sf_pdf_alphas=%.8f'%(self.rwbs['sf_mu_R'][0][0],self.rwbs['sf_mu_F'][0][0],self.rwbs['sf_scale_comb'][0][0],self.rwbs['sf_pdf_alphas'][0][0]) #DEBUG
		initialweight = self.rwbs['weight'][0][0]*self.rwbs['sf_pileup'][0][0]*self.rwbs['sf_btag_eff'][0][0]*self.rwbs['sf_mu_R'][0][0]
		initialweight*= self.rwbs['sf_mu_F'][0][0]*self.rwbs['sf_scale_comb'][0][0]*self.rwbs['sf_pdf_alphas'][0][0]
		#print 'initial weight = %.4f'%(initialweight) #DEBUG
		#Have to use the correctors to get other stuff though because they may use different objects than in the reconstructor
		if self.mode=='t' :
			#muon trigger efficiency
			mutrigeffbtof, scrap1, scrap2, mutrigeffgh, scrap3, scrap4 = self.corrector.getTrigEff(self.phyobjbs['evt_top'][0][0],mu',self.phyobjbs['themuon_pt'][0][0],self.phyobjbs['themuon_eta'][0][0])
		#	print 'trigeffbtof=%.4f, trigeffgh=%.4f'%(mutrigeffbtof,mutrigeffgh) #DEBUG
			#muon ID efficiency
			muideffbtof, scrap1, scrap2, muideffgh, scrap3, scrap4 = self.corrector.getIDEff(self.phyobjbs['evt_pu'][0][0],'mu',self.phyobjbs['themuon_pt'][0][0],self.phyobjbs['themuon_eta'][0][0])
		#	#print 'ideffbtof=%.4f, ideffgh=%.4f'%(muideffbtof,muideffgh) #DEBUG
			#muon isolation efficiency
			muisoeffbtof, scrap1, scrap2, muisoeffgh, scrap3, scrap4 = self.corrector.getIsoEff(self.phyobjbs['evt_pu'][0][0],self.phyobjbs['evt_top'][0][0],'mu',self.phyobjbs['themuon_pt'][0][0],self.phyobjbs['themuon_eta'][0][0])
		#	print 'isoeffbtof=%.4f, isoeffgh=%.4f'%(muisoeffbtof,muisoeffgh) #DEBUG
			#muon tracking efficiency
			mutrkeff, scrap1, scrap2 = self.corrector.getTrkEff(self.phyobjbs['evt_pu'][0][0],'mu',self.phyobjbs['themuon_pt'][0][0],self.phyobjbs['themuon_eta'][0][0])
		#	print 'trkeff=%.4f'%(mutrkeff) #DEBUG
			#get the lepton flavor
			ismu = self.phyobjbs['lepflavor'][0][0]==1; isel = self.phyobjbs['lepflavor'][0][0]==2
			#combine and return
			finalweight = initialweight*mutrkeff*(((LUMINOSITYBTOFMU*ismu+LUMINOSITYBTOFELE*isel)*mutrigeffbtof*muideffbtof*muisoeffbtof)+((LUMINOSITYGHMU*ismu+LUMINOSITYGHELE*isel)*mutrigeffgh*muideffgh*muisoeffgh))
		#	if self.rwbs['sf_pileup'][0][0]<0.5 : #DEBUG
		#		print 'finalweight=%.4f'%(finalweight) #DEBUG
			return finalweight
		elif self.mode=='id' :
			#tag ID efficiency
			tagideffbtof, scrap1, scrap2, tagideffgh, scrap3, scrap4 = self.corrector.getIDEff(self.phyobjbs['evt_pu'][0][0],'el',self.phyobjbs['tag_pt'][0][0],self.phyobjbs['tag_eta'][0][0])
			#tag isolation efficiency
			tagisoeffbtof, scrap1, scrap2, tagisoeffgh, scrap3, scrap4 = self.corrector.getIsoEff(self.phyobjbs['evt_pu'][0][0],self.phyobjbs['evt_top'][0][0],'el',self.phyobjbs['tag_pt'][0][0],self.phyobjbs['tag_eta'][0][0])
			#combine and return
			return initialweight*((LUMINOSITYBTOFELE*tagideffbtof*tagisoeffbtof)+(LUMINOSITYGHELE*tagideffgh*tagisoeffgh))
		elif self.mode=='iso' : #fix this later
			return None

	def getEventsSelected(self) :
		return self.ele_pt_all.Integral()
	def getEventsPassing(self) :
		return self.ele_pt_pass.Integral()
	def getGraphs(self) :
		return self.ele_pt_gr, self.ele_eta_gr, self.pu_gr
	def getFitGraphs(self) :
		return self.ele_pt_fit_gr, self.ele_eta_fit_gr, self.pu_fit_gr
	def getHisto2DAll(self) :
		return self.histo_2d_all
	def getHisto2DPass(self) :
		return self.histo_2d_pass

##########									Set up 									##########

print 'Opening files . . . '
#Set up the chains and add files
data_chain = TChain('tree')
mc_chain = TChain('tree')
filenamelist = []
for sample in samplenames :
	if options.usepreskim == 'no' :
		filenamelist+=glob('../'+sample+'/aggregated_'+sample+'_*.root')
	elif options.usepreskim == 'yes' :
		filenamelist+=glob('../total_ttree_files/'+sample+'_skim_all.root')
print 'Chaining files... '
for filename in filenamelist :
	if filename.find('SingleMu')!=-1 or filename.find('SingleEl')!=-1 :
		#print 'Adding file '+filename+' to data chain' #DEBUG
		data_chain.Add(filename)
	else :
		#print 'Adding file '+filename+' to MC chain' #DEBUG
		mc_chain.Add(filename)
#Set up output file
filename = 'electron'
if mode=='t' :
	filename+='_trigger'
elif mode=='id' :
	filename+='_ID'
if mode=='iso' :
	filename+='_isolation'
filename+='_efficiency_'+options.topology+'_'+options.ttbar_generator+'.root'
print 'Output file will be %s'%(filename)
outfile  = TFile(filename,'recreate')
#Set up file group objects and make histos/graphs
mc_group = filegroup('mc',mc_chain,outfile,mode)
data_group = filegroup('data',data_chain,outfile,mode)
#calculate numbers
data_events_selected = data_group.getEventsSelected()
data_events_passing  = data_group.getEventsPassing()
mc_events_selected 	 = mc_group.getEventsSelected()
mc_events_passing  	 = mc_group.getEventsPassing()
data_eff = data_events_passing/data_events_selected
data_eff_err = data_eff*sqrt(1./data_events_passing+1./data_events_selected)
mc_eff = mc_events_passing/mc_events_selected
mc_eff_err = mc_eff*sqrt(1./mc_events_passing+1./mc_events_selected)
data_mc_sf = data_eff/mc_eff
data_mc_sf_err = data_mc_sf*sqrt((data_eff_err/data_eff)**2+(mc_eff_err/mc_eff)**2)
print 'DATA EVENTS SELECTED = '+str(data_events_selected)
print 'DATA EVENTS PASSING = '+str(data_events_passing) 
print 'DATA EFFICIENCY = '+str(data_eff)+' +/- '+str(data_eff_err)
print 'MC EVENTS SELECTED = '+str(mc_events_selected)
print 'MC EVENTS PASSING = '+str(mc_events_passing)
print 'MC EFFICIENCY = '+str(mc_eff)+' +/- '+str(mc_eff_err)
print 'DATA/MC SCALEFACTOR = '+str(data_mc_sf)+' +/- '+str(data_mc_sf_err)
#make pretty plots
canvs = []
el_pt_canv = TCanvas('el_pt_canv','el_pt_canv',1100,900); canvs.append(el_pt_canv)
el_eta_canv = TCanvas('el_eta_canv','el_eta_canv',1100,900); canvs.append(el_eta_canv)
el_pu_canv = TCanvas('el_pu_canv','el_pu_canv',1100,900); canvs.append(el_pu_canv)
data_graphs = data_group.getGraphs()
mc_graphs = mc_group.getGraphs()
data_eff_lines = []
data_eff_lines.append(TLine(ele_pt_bins[0],data_eff,ele_pt_bins[n_ele_pt_bins],data_eff))
data_eff_lines.append(TLine(eta_bins[0],data_eff,eta_bins[n_eta_bins],data_eff))
data_eff_lines.append(TLine(pu_bins[0],data_eff,pu_bins[n_pu_bins],data_eff))
mc_eff_lines = []
mc_eff_lines.append(TLine(ele_pt_bins[0],mc_eff,ele_pt_bins[n_ele_pt_bins],mc_eff))
mc_eff_lines.append(TLine(eta_bins[0],mc_eff,eta_bins[n_eta_bins],mc_eff))
mc_eff_lines.append(TLine(pu_bins[0],mc_eff,pu_bins[n_pu_bins],mc_eff))
data_fit_graphs = data_group.getFitGraphs()
mc_fit_graphs = mc_group.getFitGraphs()
for i in range(len(canvs)) :
	canvs[i].cd()
	data_graphs[i].SetMarkerStyle(21)
	data_graphs[i].SetMarkerColor(kRed)
	data_graphs[i].SetFillStyle(3004)
	data_graphs[i].SetFillColor(kRed)
	mc_graphs[i].SetMarkerStyle(21)
	mc_graphs[i].SetMarkerColor(kBlue)
	mc_graphs[i].SetFillStyle(3005)
	mc_graphs[i].SetFillColor(kBlue)
	data_fit_graphs[i].SetLineWidth(3)
	data_fit_graphs[i].SetLineStyle(2)
	data_fit_graphs[i].SetLineColor(kRed)
	data_fit_graphs[i].SetFillColor(kRed)
	data_fit_graphs[i].SetFillStyle(3003)
	mc_fit_graphs[i].SetLineWidth(3)
	mc_fit_graphs[i].SetLineStyle(2)
	mc_fit_graphs[i].SetLineColor(kBlue)
	mc_fit_graphs[i].SetFillColor(kBlue)
	mc_fit_graphs[i].SetFillStyle(3003)
	data_eff_lines[i].SetLineWidth(3)
	data_eff_lines[i].SetLineColor(kRed)
	mc_eff_lines[i].SetLineWidth(3)
	mc_eff_lines[i].SetLineColor(kBlue)
	data_graphs[i].GetYaxis().SetRangeUser(0.5,1.3)
	data_graphs[i].SetTitle(data_graphs[i].GetTitle().lstrip(data_group.name))
	leg = TLegend(0.2,0.7,0.4,0.85)
	leg.AddEntry(data_graphs[i],'data','PEF')
	leg.AddEntry(data_fit_graphs[i],'data fit','LF')
	leg.AddEntry(data_eff_lines[i],'overall data','L')
	leg.AddEntry(mc_graphs[i],'MC','PEF')
	leg.AddEntry(mc_fit_graphs[i],'MC fit','LF')
	leg.AddEntry(mc_eff_lines[i],'overall MC','L')
	data_graphs[i].Draw('APE3')
	mc_graphs[i].Draw('PE3 SAME')
	data_eff_lines[i].Draw('SAME')
	mc_eff_lines[i].Draw('SAME')
	data_fit_graphs[i].Draw('LE3 SAME')
	mc_fit_graphs[i].Draw('LE3 SAME')
	leg.Draw('SAME')
	canvs[i].Update()
	outfile.cd()
	canvs[i].Write()
#make the 2D scalefactor plot
sf_plot = TH2D('sf_plot','Data/MC efficiency scalefactors over electron phase space; p_{T} (GeV); |#eta|',n_ele_pt_bins_2D,ele_pt_bins_2D,n_eta_bins_2D,eta_bins_2D)
data_2D_histo_all = data_group.getHisto2DAll(); data_2D_histo_pass = data_group.getHisto2DPass()
mc_2D_histo_all = mc_group.getHisto2DAll(); mc_2D_histo_pass = mc_group.getHisto2DPass()
for i in range(n_ele_pt_bins_2D) :
	for j in range(n_eta_bins_2D) :
		xvalue = (ele_pt_bins_2D[i+1]+ele_pt_bins_2D[i])/2.
		yvalue = (eta_bins_2D[j+1]+eta_bins_2D[j])/2.
		thisbin = mc_2D_histo_all.FindFixBin(xvalue,yvalue)
		mcall   = mc_2D_histo_all.GetBinContent(thisbin)
		mcpass  = mc_2D_histo_pass.GetBinContent(thisbin)
		dataall   = data_2D_histo_all.GetBinContent(thisbin)
		datapass  = data_2D_histo_pass.GetBinContent(thisbin)
		content = 1.0; err = 1.0
		if datapass!=0. and dataall!=0. and mcpass!=0. and mcall!=0. :
			content = (datapass/dataall)/(mcpass/mcall)
			err = content*sqrt(1./datapass+1./dataall+1./mcpass+1./mcall)
		sf_plot.SetBinContent(thisbin,content)
		sf_plot.SetBinError(thisbin,err)
sf_plot_canv = TCanvas('sf_plot_canv','sf_plot_canv',1200,900);
sf_plot_canv.cd()
sf_plot.Draw("COLZ")
outfile.cd()
sf_plot_canv.Write()
sf_plot.Write()

outfile.Close()