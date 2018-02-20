from ROOT import *
import os

#build a chain of the unskimmed files
#unskimmed_dir = '../powheg_TT/'
unskimmed_dir = '../SingleMu_Run2016Bv2/'
unskimmed_chain = TChain('tree')
filelist = []
for filename in os.listdir(unskimmed_dir) :
	if filename.find('JES')==-1 and filename.find('JER')==-1 and filename.find('.root')!=-1 :
		filelist.append(filename)
		unskimmed_chain.Add(os.path.join(unskimmed_dir,filename))
print 'unskimmed files = %s'%(filelist) #DEBUG

#open the skimmed file and get its tree
#skimmed_file_path = '../total_ttree_files/powheg_TT_skim_all.root'
skimmed_file_path = '../total_ttree_files/SingleMu_Run2016Bv2_skim_all.root'
skimmed_file = TFile.Open(skimmed_file_path)
skimmed_tree = skimmed_file.Get('tree')
print 'skimmed tree (from %s) = %s '%(skimmed_file_path,skimmed_tree) #DEBUG

#open output file
output_file_name = 'skim_check_plots.root'
output_file = TFile(output_file_name,'recreate')

#plot class
class Plot(object):
	"""docstring for Plot"""
	def __init__(self,canvas_name,title,unskimmed_cutstring,skimmed_cutstring):
		self._c_name = canvas_name
		self._title = title+'; c*; Events/0.1'
		self._u_cut = unskimmed_cutstring
		self._s_cut = skimmed_cutstring
	def savePlot(self) :
		print 'saving plot "%s"'%(self._c_name)
		#define the canvas
		self._canv = TCanvas(self._c_name,self._c_name,1100,900)
		self._canv.cd()
		#plot cstar for the unskimmed and skimmed events
		unskimmed_chain.Draw('cstar>>%s_u(20,-1.,1.)'%(self._c_name),'(%s)'%(self._u_cut))
		self._u_h = gROOT.FindObject('%s_u'%(self._c_name))
		self._u_h.SetLineColor(kRed); self._u_h.SetLineWidth(3); self._u_h.SetTitle(self._title)
		skimmed_tree.Draw('cstar>>%s_s(20,-1.,1.)'%(self._c_name),'(%s)'%(self._s_cut))
		self._s_h = gROOT.FindObject('%s_s'%(self._c_name))
		self._s_h.SetLineColor(kBlue); self._s_h.SetLineWidth(3)
		self._canv.cd()
		self._u_h.Draw('HIST'); self._s_h.Draw('HIST SAME')
		output_file.cd()
		self._canv.Write()

#defining cuts
skimcut = 'fullselection==1 || wjets_cr_selection==1 || qcd_A_SR_selection==1 || qcd_B_SR_selection==1 || qcd_C_SR_selection==1 || qcd_A_CR_selection==1 || qcd_B_CR_selection==1 || qcd_C_CR_selection==1'
fullselection = 'fullselection==1'

#define plots
allplots = []
allplots.append(Plot('all','all skimmed events',skimcut,'weight!=0.'))
allplots.append(Plot('all_v2','all skimmed events (v2)',skimcut,skimcut))
allplots.append(Plot('selected','fully selected events',fullselection,fullselection))

#plot plots
for p in allplots :
	p.savePlot()

#close the output file
output_file.Close()

print 'Done!'