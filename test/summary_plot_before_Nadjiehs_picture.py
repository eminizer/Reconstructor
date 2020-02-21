# This script makes the summary plot for the Tevatron A_FB measurements, 
# with mine at the bottom, compared to MC predictions

# Imports
from ROOT import *
from math import *
import CMS_lumi, tdrstyle
from array import array

# Style stuff
gROOT.SetBatch()
tdrstyle.setTDRStyle()
iPeriod = 0 # iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.cmsTextSize = 2.5
CMS_lumi.relPosX  = 0.008
CMS_lumi.relPosY = 0.008
CMS_lumi.relExtraDY = 0.85
CMS_lumi.extraOverCmsTextSize = 0.67

# output filename and plot canvas name
outfilename = 'Afb_summary_plot.root'
plotname = 'summary_plot_A_FB'
# x axis limits and number of bins to define background histogram
x_low = -0.18
x_hi = 0.32
nhistbinsx = 50
#where the middle line is between the entry boxes and the number boxes
xmid = 0.75
#length of the lines on the ends of the graph error bars
errbarlinelength = 0.0175
#graph colors
color_d0 = kRed+2
color_cdf = kBlue+2
color_cdfd0 = kGreen+2

# list of Tevatron measurements by category
tevatron_afbs = {
	'lepton+jets':{
		#'D0 (0.9 fb^{-1})':{
		#	'ref':'PRL 100 (2008) 142002',
		#	'Afb':0.12,
		#	'statErr':0.08,
		#	'systErr':0.01,
		#},
		#'CDF (1.9 fb^{-1})':{
		#	'ref':'PRL 101 (2008) 202001',
		#	'Afb':0.24,
		#	'statErr':0.13,
		#	'systErr':0.04,
		#},
		#'CDF (5.3 fb^{-1})':{
		#	'ref':'PRD 83 (2011) 112003',
		#	'Afb':0.158,
		#	'statErr':0.072,
		#	'systErr':0.017,
		#},
		#'D0 (5.4 fb^{-1})':{
		#	'ref':'PRD 84 (2011) 112005',
		#	'Afb':0.196,
		#	'statErr':0.060,
		#	'systErr':{'up':0.018,'dn':0.026},
		#},
		'CDF (9.4 fb^{-1})':{
			'ref':'PRD 87 (2013) 092002, #sqrt{s}=1.96 TeV',
			'Afb':0.164,
			'statErr':0.039,
			'systErr':0.026,
		},
		'D0 (9.7 fb^{-1})':{
			'ref':'PRD 90 (2014) 072011, #sqrt{s}=1.96 TeV',
			'Afb':0.106,
			'statErr':0.027,
			'systErr':0.013,
		},
	},
	'dilepton':{
		'D0 (9.7 fb^{-1})':{
			'ref':'PRD 92 (2015) 052007, #sqrt{s}=1.96 TeV',
			'Afb':0.175,
			'statErr':0.056,
			'systErr':0.031,
		},
		'CDF (9.1 fb^{-1})':{
			'ref':'PRD 93 (2016) 112005, #sqrt{s}=1.96 TeV',
			'Afb':0.12,
			'statErr':0.11,
			'systErr':0.07,
		},
	},
	'combination':{
		'D0 (9.7 fb^{-1})':{
			'ref':'PRD 92 (2015) 052007, #sqrt{s}=1.96 TeV',
			'Afb':0.118,
			'statErr':0.025,
			'systErr':0.013,
		},
		'CDF (9.1 fb^{-1})':{
			'ref':'PRD 93 (2016) 112005, #sqrt{s}=1.96 TeV',
			'Afb':0.160,
			'totErr':0.045,
		},
		'CDF+D0':{
			'ref':'PRL 120 (2018) 042001, #sqrt{s}=1.96 TeV',
			'Afb':0.128,
			'statErr':0.021,
			'systErr':0.014,
		},
	},
}

# Numbers for my measurement
my_afb = {
	'CMS (35.9 fb^{-1})':{
		'ref':'TOP-15-018 (2019), #sqrt{s}=13 TeV',
		'Afb':0.048,
		'statErr':{'up':0.088,'dn':0.084},
		'systErr':0.028,
	},
}

# Tevatron theory comparison
tevatron_theory = {
	'NNLO QCD (+ NLO EW)':{
		'ref':'Czakon et. al. PRL 115 (2015) 052001, #sqrt{s}=1.96 TeV',
		'Afb':0.095,
		'totErr':0.007,
	}
}

# Powheg theory comparison (for my measurement)
powheg_theory = {
	'POWHEGv2 NLO':{
		'ref':'q#bar{q}, event counting, #sqrt{s}=13 TeV',
		'Afb':0.0512,
		'totErr':0.0004
	}
}

# Helper function to parse error values in case they're asymmetric
def parseErrorValues(err) :
	errUp = 0.; errDn = 0.
	if isinstance(err,float)  :
		errUp = err
		errDn = err
	elif 'up' in err.keys() and 'dn' in err.keys() :
		errUp = err['up']
		errDn = err['dn']
	else :
		print 'WARNING: err = %s and cannot be parsed by parseErrorValues....'%(err)
	return errUp, errDn

# Entry class
class Entry(object) :

	def __init__(self,name,details,is_theory=False) :
		self._name = name
		self._number = None
		self._is_theory = is_theory
		self._collaboration = name.split()[0]
		self._refstring = details['ref']
		self._year = int(self._refstring.split('(')[1][:4]) if self._refstring.find('(')!=-1 else 3000
		self._top = 0.0
		self._bot = 0.0
		self._lbox = None
		self._rbox = None
		self._nbox = None
		self._Afb = details['Afb']
		self._statErrUp = 0.0
		self._statErrDn = 0.0
		self._systErrUp = 0.0
		self._systErrDn = 0.0
		self._totErrUp = 0.0
		self._totErrDn = 0.0
		keys = details.keys()
		if 'totErr' in keys :
			self._totErrUp, self._totErrDn = parseErrorValues(details['totErr'])			
		elif 'statErr' in keys and 'systErr' in keys :
			self._statErrUp, self._statErrDn = parseErrorValues(details['statErr'])
			self._systErrUp, self._systErrDn = parseErrorValues(details['systErr'])
			self._totErrUp = sqrt(self._statErrUp**2+self._systErrUp**2)
			self._totErrDn = sqrt(self._statErrDn**2+self._systErrDn**2)
		else :
			print 'WARNING: error listing scheme unidentifiable for entry %s with details:\n%s\n'%(name, details)

	def getYear(self) :
		return self._year
	def getIssue(self) :
		return int(self._refstring.split()[1])
	def getCollabName(self) :
		return self._collaboration
	def setNumber(self,n) :
		self._number = n
	def makeBoxes(self,lmargin,rmargin,top,h) :
		x1 = lmargin+0.005
		x2 = 1. - (rmargin+0.005)
		self._top = top - (self._number * h)
		self._bot = self._top - h
		lboxy2 = self._top
		lboxy1 = self._top - (2./3.)*h
		rboxy2 = lboxy1
		rboxy1 = lboxy2 - h
		self._lbox = TPaveText(x1,lboxy1,xmid,lboxy2,'NDC')
		self._rbox = TPaveText(x1,rboxy1,xmid,rboxy2,'NDC')
		self._nbox = TPaveText(xmid,self._bot,x2,self._top,'NDC')
		self._lbox.SetTextSize(0.034)
		self._rbox.SetTextSize(0.019)
		self._nbox.SetTextSize(0.025)
		allboxes = [self._lbox,self._rbox,self._nbox]
		for b in allboxes :
			b.SetBorderSize(0)
			b.SetFillColor(0)
			b.SetTextAlign(12)
			b.SetTextFont(42)
			b.SetMargin(0.01)
			#b.SetTextFont(4)
			if self._is_theory :
				b.SetTextColor(14)
		self._nbox.SetTextAlign(13)
		self._rbox.SetTextAlign(11)
		self._lbox.AddText(self._name)
		self._rbox.AddText(self._refstring)
		ndigits = len(str(self._Afb).split('.')[1].rstrip('0'))
		nboxstring = ''
		nboxstring+='{0:.{1}f}'.format(self._Afb,ndigits)
		if self._statErrUp!=0. :
			if self._statErrUp!=self._statErrDn :
				nboxstring+='^{{+{0:.{2}f}}}_{{-{1:.{2}f}}}'.format(self._statErrUp,self._statErrDn,ndigits)
			else :
				nboxstring+=' #pm {0:.{1}f}'.format(self._statErrUp,ndigits)
			if self._systErrUp!=0. :
				if self._systErrUp!=self._systErrDn :
					nboxstring+='^{{+{0:.{2}f}}}_{{-{1:.{2}f}}}'.format(self._systErrUp,self._systErrDn,ndigits)
				else :
					nboxstring+=' #pm {0:.{1}f}'.format(self._systErrUp,ndigits)
		else :
			if self._totErrUp!=self._totErrDn :
				nboxstring+='^{{+{0:.{2}f}}}_{{-{1:.{2}f}}}'.format(self._totErrUp,self._totErrDn,ndigits)
			else :
				nboxstring+=' #pm {0:.{1}f}'.format(self._totErrUp,ndigits)
		self._nbox.AddText(nboxstring)
	def getLbox(self) :
		return self._lbox
	def getRbox(self) :
		return self._rbox
	def getNbox(self) :
		return self._nbox
	def getBottom(self) :
		return self._bot
	def getTop(self) :
		return self._top
	def getAfb(self) :
		return self._Afb
	def getTotErr(self) :
		return 0.5*(self._totErrUp+self._totErrDn)
	def getStatErrHi(self) :
		return self._statErrUp
	def getStatErrLo(self) :
		return self._statErrDn
	def getTotErrHi(self) :
		return self._totErrUp
	def getTotErrLo(self) :
		return self._totErrDn

# Category class
class Category(object) :

	def __init__(self,name,entrydict) :
		self._name = name
		self._entries = []
		for ename,entry in entrydict.iteritems() :
			self._entries.append(Entry(ename,entry))
		self._entries.sort(key=lambda x: (x.getYear(),x.getIssue()))
		self._lboxes = []
		self._rboxes = []
		self._nboxes = []

	def getName(self) :
		return self._name
	def getCollabNames(self) :
		return [entry.getCollabName() for entry in self._entries]
	def getNpoints(self) :
		return len(self._entries)
	def numberEntriesFrom(self,start) :
		for ientry in range(len(self._entries)) :
			self._entries[ientry].setNumber(start+ientry)
	def makeEntryBoxes(self,lm,rm,top,h) :
		for entry in self._entries :
			entry.makeBoxes(lm,rm,top,h)
			self._lboxes.append(entry.getLbox())
			self._rboxes.append(entry.getRbox())
			self._nboxes.append(entry.getNbox())
	def getEntryBoxes(self) :
		return self._lboxes+self._rboxes+self._nboxes
	def getSectionBottom(self) :
		return self._entries[-1].getBottom()
	def getSectionTop(self) :
		return self._entries[0].getTop()
	def getEntries(self) :
		return self._entries

# Plot class
class Plot(object) :

	def __init__(self) :
		self._collaboration_names = []
		self._categories = []
		self._n_points = 0
		self._tevatron_theory = None
		self._my_measurement = None
		self._powheg_theory = None
		self._bg_histo = None
		self._divider_lines = []
		self._cat_labels = []
		self._legend_height = 0.10

	def addCategory(self,catname,entrydict) :
		self._categories.append(Category(catname,entrydict))
		for collab_name in self._categories[-1].getCollabNames() :
			if not collab_name in self._collaboration_names :
				self._collaboration_names.append(collab_name)

	def addTheory(self,ids,name,details) :
		newEntry = Entry(name,details,True)
		if ids == 'tevatron' :
			self._tevatron_theory = newEntry
		elif ids == 'powheg' :
			self._powheg_theory = newEntry
		else :
			print 'WARNING: theory ID string %s for theory dictionary %s with details %s not recognized by addTheory....'%(ids,name,details)

	def addMyAfb(self,afbd) :
		self._my_measurement = Entry(afbd.keys()[0],afbd[afbd.keys()[0]])

	def __countNpoints__(self) :
		for cat in self._categories :
			cat.numberEntriesFrom(self._n_points)
			self._n_points+=cat.getNpoints()
		if self._tevatron_theory != None :
			self._tevatron_theory.setNumber(self._n_points)
			self._n_points+=1
		if self._my_measurement != None :
			self._my_measurement.setNumber(self._n_points)
			self._n_points+=1
		if self._powheg_theory != None :
			self._powheg_theory.setNumber(self._n_points)
			self._n_points+=1

	def __makeAndPlotHistogram__(self) :
		self._bg_histo = TH2F('bg','',nhistbinsx,x_low,x_hi,100,self._canvas.GetBottomMargin(),1.-self._canvas.GetTopMargin())
		self._bg_histo.SetStats(0)
		self._bg_histo.GetXaxis().CenterTitle(True)
		self._bg_histo.GetXaxis().SetLabelSize(0.03)
		self._bg_histo.GetXaxis().SetTitleSize(0.05)
		self._bg_histo.GetXaxis().SetTitle('A_{FB}')
		self._bg_histo.GetYaxis().SetLabelSize(0.0)
		self._bg_histo.GetYaxis().SetTickLength(0.0)
		#self._bg_histo.GetYaxis().LabelsOption('v')
		self._bg_histo.Draw()

	def __makeAndPlotBoxes__(self) :
		lm = self._canvas.GetLeftMargin()
		rm = self._canvas.GetRightMargin()
		y1 = self._canvas.GetBottomMargin() + self._bg_histo.GetXaxis().GetTickLength()
		y2 = 1. - self._canvas.GetTopMargin() - self._legend_height #- self._bg_histo.GetXaxis().GetTickLength() 
		entryboxheight = (y2-y1)/float(self._n_points)
		allboxes = []
		for cat in self._categories :
			cat.makeEntryBoxes(lm,rm,y2,entryboxheight)
			allboxes+=cat.getEntryBoxes()
		if self._tevatron_theory != None :
			self._tevatron_theory.makeBoxes(lm,rm,y2,entryboxheight)
			allboxes+=[self._tevatron_theory.getLbox(),self._tevatron_theory.getRbox(),self._tevatron_theory.getNbox()]
		if self._my_measurement != None :
			self._my_measurement.makeBoxes(lm,rm,y2,entryboxheight)
			allboxes+=[self._my_measurement.getLbox(),self._my_measurement.getRbox(),self._my_measurement.getNbox()]
		if self._powheg_theory != None :
			self._powheg_theory.makeBoxes(lm,rm,y2,entryboxheight)
			allboxes+=[self._powheg_theory.getLbox(),self._powheg_theory.getRbox(),self._powheg_theory.getNbox()]
		for box in allboxes :
			box.Draw('SAME')

	def __makeAndPlotTheoryLines__(self) :
		tvly1 = self._tevatron_theory.getBottom()
		tvly2 = self._categories[0].getSectionTop()
		tvafb = self._tevatron_theory.getAfb()
		tvtoterr = self._tevatron_theory.getTotErr()
		self._tevatron_box = TBox(tvafb-tvtoterr,tvly1,tvafb+tvtoterr,tvly2)
		self._tevatron_line = TLine(tvafb,tvly1,tvafb,tvly2)
		ply1 = self._powheg_theory.getBottom()
		ply2 = self._my_measurement.getTop()
		pafb = self._powheg_theory.getAfb()
		ptoterr = self._powheg_theory.getTotErr()
		self._powheg_box = TBox(pafb-ptoterr,ply1,pafb+ptoterr,ply2)
		self._powheg_line = TLine(pafb,ply1,pafb,ply2)
		self._tevatron_box.SetFillColor(kBlue-9)
		self._powheg_box.SetFillColor(kMagenta+1)
		self._tevatron_line.SetLineColor(kBlue-7)
		self._powheg_line.SetLineColor(kMagenta+1)
		for line in [self._tevatron_line,self._powheg_line] :
			line.SetLineWidth(3)
		for obj in [self._tevatron_box,self._powheg_box,self._tevatron_line,self._powheg_line] :
			obj.Draw('SAME')

	def __makeAndPlotGraphs__(self) :
		#the list of the error bar lines
		self._all_err_bar_lines = []
		# find the points from each collaboration
		points_D0_x = []; points_D0_y = []; points_D0_xel_stat = []; points_D0_xeh_stat = []; points_D0_xel_tot = []; points_D0_xeh_tot = []
		points_CDF_x = []; points_CDF_y = []; points_CDF_xel_stat = []; points_CDF_xeh_stat = []; points_CDF_xel_tot = []; points_CDF_xeh_tot = []
		points_comb_x = []; points_comb_y = []; points_comb_xel_stat = []; points_comb_xeh_stat = []; points_comb_xel_tot = []; points_comb_xeh_tot = []
		for cat in self._categories :
			for entry in cat.getEntries() :
				collab = entry.getCollabName()
				px = entry.getAfb()
				etop = entry.getTop()
				ebot = entry.getBottom()
				py = etop-(1./3.)*(etop-ebot)
				pexl_stat = entry.getStatErrLo()
				pexh_stat = entry.getStatErrHi()
				pexl_tot = entry.getTotErrLo()
				pexh_tot = entry.getTotErrHi()
				newlines = []
				if pexl_stat!=0. :
					new_line_stat_low = TLine(px-pexl_stat,py-0.5*errbarlinelength,px-pexl_stat,py+0.5*errbarlinelength); newlines.append(new_line_stat_low)
				if pexh_stat!=0. :
					new_line_stat_high = TLine(px+pexh_stat,py-0.5*errbarlinelength,px+pexh_stat,py+0.5*errbarlinelength); newlines.append(new_line_stat_high)
				if pexl_tot!=0. :
					new_line_tot_low = TLine(px-pexl_tot,py-0.5*errbarlinelength,px-pexl_tot,py+0.5*errbarlinelength); newlines.append(new_line_tot_low)
				if pexh_tot!=0. :
					new_line_tot_high = TLine(px+pexh_tot,py-0.5*errbarlinelength,px+pexh_tot,py+0.5*errbarlinelength); newlines.append(new_line_tot_high)
				for line in newlines  :
					line.SetLineWidth(3)
					self._all_err_bar_lines.append(line)
				if collab=='D0' :
					points_D0_x.append(px); points_D0_y.append(py) 
					points_D0_xel_stat.append(pexl_stat); points_D0_xeh_stat.append(pexh_stat) 
					points_D0_xel_tot.append(pexl_tot); points_D0_xeh_tot.append(pexh_tot)
					for line in newlines :
						line.SetLineColor(color_d0)
				elif collab=='CDF' :
					points_CDF_x.append(px); points_CDF_y.append(py) 
					points_CDF_xel_stat.append(pexl_stat); points_CDF_xeh_stat.append(pexh_stat) 
					points_CDF_xel_tot.append(pexl_tot); points_CDF_xeh_tot.append(pexh_tot)
					for line in newlines :
						line.SetLineColor(color_cdf)
				elif collab=='CDF+D0' :
					points_comb_x.append(px); points_comb_y.append(py) 
					points_comb_xel_stat.append(pexl_stat); points_comb_xeh_stat.append(pexh_stat) 
					points_comb_xel_tot.append(pexl_tot); points_comb_xeh_tot.append(pexh_tot)
					for line in newlines :
						line.SetLineColor(color_cdfd0)
				else :
					print 'WARNING: collaboration name %s not recognized by __makeAndPlotGraphs__'%(collab)
		npoints_D0 = len(points_D0_x)
		self._gr_D0_stat = TGraphAsymmErrors(npoints_D0,array('d',points_D0_x),array('d',points_D0_y),
												array('d',points_D0_xel_stat),array('d',points_D0_xeh_stat),
												array('d',npoints_D0*[0.0]),array('d',npoints_D0*[0.0]))
		self._gr_D0_tot = TGraphAsymmErrors(npoints_D0,array('d',points_D0_x),array('d',points_D0_y),
												array('d',points_D0_xel_tot),array('d',points_D0_xeh_tot),
												array('d',npoints_D0*[0.0]),array('d',npoints_D0*[0.0]))
		npoints_CDF = len(points_CDF_x)
		self._gr_CDF_stat = TGraphAsymmErrors(npoints_CDF,array('d',points_CDF_x),array('d',points_CDF_y),
												array('d',points_CDF_xel_stat),array('d',points_CDF_xeh_stat),
												array('d',npoints_CDF*[0.0]),array('d',npoints_CDF*[0.0]))
		self._gr_CDF_tot = TGraphAsymmErrors(npoints_CDF,array('d',points_CDF_x),array('d',points_CDF_y),
												array('d',points_CDF_xel_tot),array('d',points_CDF_xeh_tot),
												array('d',npoints_CDF*[0.0]),array('d',npoints_CDF*[0.0]))
		npoints_comb = len(points_comb_x)
		self._gr_comb_stat = TGraphAsymmErrors(npoints_comb,array('d',points_comb_x),array('d',points_comb_y),
												array('d',points_comb_xel_stat),array('d',points_comb_xeh_stat),
												array('d',npoints_comb*[0.0]),array('d',npoints_comb*[0.0]))
		self._gr_comb_tot = TGraphAsymmErrors(npoints_comb,array('d',points_comb_x),array('d',points_comb_y),
												array('d',points_comb_xel_tot),array('d',points_comb_xeh_tot),
												array('d',npoints_comb*[0.0]),array('d',npoints_comb*[0.0]))
		for gr in [self._gr_D0_stat,self._gr_D0_tot] :
			gr.SetLineColor(color_d0)
			gr.SetMarkerColor(color_d0)
			gr.SetMarkerStyle(22)
		for gr in [self._gr_CDF_stat,self._gr_CDF_tot] :
			gr.SetLineColor(color_cdf)
			gr.SetMarkerColor(color_cdf)
			gr.SetMarkerStyle(23)
		for gr in [self._gr_comb_stat,self._gr_comb_tot] :
			gr.SetLineColor(color_cdfd0)
			gr.SetMarkerColor(color_cdfd0)
			gr.SetMarkerStyle(21)
		# add the graph of my data point
		points_mine_x = []; points_mine_y = []; points_mine_xel_stat = []; points_mine_xeh_stat = []; points_mine_xel_tot = []; points_mine_xeh_tot = []
		points_mine_x.append(self._my_measurement.getAfb()); points_mine_y.append(self._my_measurement.getTop() - (1./3.)*(self._my_measurement.getTop()-self._my_measurement.getBottom())) 
		points_mine_xel_stat.append(self._my_measurement.getStatErrLo()); points_mine_xeh_stat.append(self._my_measurement.getStatErrHi()) 
		points_mine_xel_tot.append(self._my_measurement.getTotErrLo()); points_mine_xeh_tot.append(self._my_measurement.getTotErrHi())
		npoints_mine = len(points_mine_x)
		newlines = []
		for i in range(npoints_mine) :
			if points_mine_xel_stat[i]!=0. :
				new_line_stat_low = TLine(points_mine_x[i]-points_mine_xel_stat[i],points_mine_y[i]-0.5*errbarlinelength,points_mine_x[i]-points_mine_xel_stat[i],points_mine_y[i]+0.5*errbarlinelength); newlines.append(new_line_stat_low)
			if points_mine_xeh_stat[i]!=0. :
				new_line_stat_high = TLine(points_mine_x[i]+points_mine_xeh_stat[i],points_mine_y[i]-0.5*errbarlinelength,points_mine_x[i]+points_mine_xeh_stat[i],points_mine_y[i]+0.5*errbarlinelength); newlines.append(new_line_stat_high)
			if points_mine_xel_tot[i]!=0. :
				new_line_tot_low = TLine(points_mine_x[i]-points_mine_xel_tot[i],points_mine_y[i]-0.5*errbarlinelength,points_mine_x[i]-points_mine_xel_tot[i],points_mine_y[i]+0.5*errbarlinelength); newlines.append(new_line_tot_low)
			if points_mine_xeh_tot[i]!=0. :
				new_line_tot_high = TLine(points_mine_x[i]+points_mine_xeh_tot[i],points_mine_y[i]-0.5*errbarlinelength,points_mine_x[i]+points_mine_xeh_tot[i],points_mine_y[i]+0.5*errbarlinelength); newlines.append(new_line_tot_high)
		for line in newlines  :
			line.SetLineWidth(3)
			self._all_err_bar_lines.append(line)
		self._gr_mine_stat = TGraphAsymmErrors(npoints_mine,array('d',points_mine_x),array('d',points_mine_y),
												array('d',points_mine_xel_stat),array('d',points_mine_xeh_stat),
												array('d',npoints_mine*[0.0]),array('d',npoints_mine*[0.0]))
		self._gr_mine_tot = TGraphAsymmErrors(npoints_mine,array('d',points_mine_x),array('d',points_mine_y),
												array('d',points_mine_xel_tot),array('d',points_mine_xeh_tot),
												array('d',npoints_mine*[0.0]),array('d',npoints_mine*[0.0]))
		for gr in [self._gr_mine_stat,self._gr_mine_tot] :
			gr.SetLineColor(kBlack)
			gr.SetMarkerColor(kBlack)
			gr.SetMarkerStyle(20)
		# set the generic graph attributes and draw them
		for gr in [self._gr_D0_stat,self._gr_D0_tot,self._gr_CDF_stat,self._gr_CDF_tot,self._gr_comb_stat,self._gr_comb_tot,self._gr_mine_stat,self._gr_mine_tot] :
			gr.SetMarkerSize(2)
			gr.SetLineWidth(3)
			gr.Draw('PE1 SAME')
		#draw the error bar lines
		for line in self._all_err_bar_lines :
			line.Draw('SAME')

	def __makeAndPlotDividers__(self) :
		#make dividers for the different types of Tevatron measurements
		for cat in self._categories :
			topy = cat.getSectionTop()
			newline = TLine(self._bg_histo.GetXaxis().GetXmin(),topy,self._bg_histo.GetXaxis().GetXmax(),topy)
			newline.SetLineWidth(2)
			newline.SetLineStyle(2)
			self._divider_lines.append(newline)
			newtex = TLatex(-0.05,topy-0.01,'p#bar{p} '+cat.getName())
			newtex.SetTextSize(0.025)
			newtex.SetTextAlign(13)
			self._cat_labels.append(newtex)
		#make divider for our measurement
		topy = self._my_measurement.getTop()
		boty = self._canvas.GetBottomMargin()+self._bg_histo.GetXaxis().GetTickLength()
		newline = TLine(self._bg_histo.GetXaxis().GetXmin(),topy,self._bg_histo.GetXaxis().GetXmax(),topy)
		newline.SetLineWidth(2)
		self._divider_lines.append(newline)
		newtex = TLatex(0.075,boty+0.02,'pp lepton+jets')
		newtex.SetTextSize(0.025)
		newtex.SetTextAlign(11)
		self._cat_labels.append(newtex)
		#draw all of the objects
		for line in self._divider_lines :
			line.Draw('SAME')
		for tex in self._cat_labels :
			tex.Draw('SAME')

	def __addLegend__(self) :
		exp_x = x_hi - 0.5*(x_hi - x_low) + 0.05
		exp_xer_stat = 0.03
		exp_xer_tot = 0.06
		#add the graph of the point in the legend
		leggry = 1. - self._canvas.GetTopMargin() - (0.5*self._legend_height) #- self._bg_histo.GetXaxis().GetTickLength()
		self._gr_leg_stat = TGraphAsymmErrors(1,array('d',[exp_x]),array('d',[leggry-0.01]),array('d',[exp_xer_stat]),array('d',[exp_xer_stat]),array('d',[0.00]),array('d',[0.00]))
		self._gr_leg_tot  = TGraphAsymmErrors(1,array('d',[exp_x]),array('d',[leggry-0.01]),array('d',[exp_xer_tot]),array('d',[exp_xer_tot]),array('d',[0.00]),array('d',[0.00]))
		for gr in [self._gr_leg_stat,self._gr_leg_tot] :
			gr.SetLineColor(14)
			gr.SetLineWidth(3)
			gr.SetMarkerColor(14)
			gr.SetMarkerStyle(20)
			gr.SetMarkerSize(2)
		#add its error bar lines
		newlines = []
		newlines.append(TLine(exp_x-exp_xer_stat,(leggry-0.01)-0.5*errbarlinelength,exp_x-exp_xer_stat,(leggry-0.01)+0.5*errbarlinelength))
		newlines.append(TLine(exp_x-exp_xer_tot,(leggry-0.01)-0.5*errbarlinelength,exp_x-exp_xer_tot,(leggry-0.01)+0.5*errbarlinelength))
		newlines.append(TLine(exp_x+exp_xer_stat,(leggry-0.01)-0.5*errbarlinelength,exp_x+exp_xer_stat,(leggry-0.01)+0.5*errbarlinelength))
		newlines.append(TLine(exp_x+exp_xer_tot,(leggry-0.01)-0.5*errbarlinelength,exp_x+exp_xer_tot,(leggry-0.01)+0.5*errbarlinelength))
		for line in newlines :
			line.SetLineWidth(3)
			line.SetLineColor(14)
			self._all_err_bar_lines.append(line)
		#add the labels to the legend graph
		self._stat_tex = TLatex(exp_x-exp_xer_stat,leggry-0.028,'stat')
		self._tot_tex = TLatex(exp_x-exp_xer_tot,leggry-0.028,'total')
		for tex in [self._stat_tex,self._tot_tex] :
			tex.SetTextSize(0.018)
			tex.SetTextColor(14)
			tex.SetTextAlign(22)
		#add the boxes for the legend also
		self._leg_tbox = TPaveText(self._canvas.GetLeftMargin()+0.005,leggry-0.3*self._legend_height,xmid,leggry+0.49*self._legend_height,'NDC')
		self._leg_bbox = TPaveText(self._canvas.GetLeftMargin()+0.005,leggry-0.49*self._legend_height,xmid,leggry-0.3*self._legend_height,'NDC')
		self._leg_nbox = TPaveText(xmid,leggry-0.49*self._legend_height,1.-(self._canvas.GetRightMargin()+0.005),leggry+0.495*self._legend_height,'NDC')
		#self._leg_tbox.AddText('CMS Preliminary')
		#self._leg_bbox.AddText('Tevatron  #sqrt{s} = 1.96 TeV; 2016 LHC #sqrt{s} = 13.0 TeV')
		self._leg_nbox.AddText('A_{FB} #pm (stat.) #pm (syst.)')
		for b in [self._leg_tbox,self._leg_bbox,self._leg_nbox] :
			b.SetBorderSize(0)
			b.SetFillColor(0)
			b.SetMargin(0.01)
			#b.SetTextFont(4)
			b.SetTextAlign(12)
		self._leg_tbox.SetTextSize(0.046)
		self._leg_tbox.SetTextAlign(13)
		self._leg_bbox.SetTextSize(0.017)
		self._leg_nbox.SetTextSize(0.025)
		self._leg_nbox.SetTextColor(14)
		self._leg_nbox.SetTextAlign(13)
		#add the boxes and labels for the theory  predictions
		boxx1 = -0.06
		boxsidelength = 0.0125
		tboxy2 = 1. - self._canvas.GetTopMargin() - 0.045
		pboxy2 = tboxy2 - boxsidelength - 0.01
		self._tevatron_theory_leg_box = TBox(boxx1,tboxy2-boxsidelength,boxx1+boxsidelength,tboxy2)
		self._powheg_theory_leg_box   = TBox(boxx1,pboxy2-boxsidelength,boxx1+boxsidelength,pboxy2)
		self._tevatron_theory_leg_box.SetFillColor(kBlue-7)
		self._powheg_theory_leg_box.SetFillColor(kMagenta+1)
		self._tevatron_theory_leg_tex = TLatex(boxx1+boxsidelength+0.001,tboxy2-0.5*boxsidelength,'Tevatron NNLO QCD')
		self._powheg_theory_leg_tex = TLatex(boxx1+boxsidelength+0.001,pboxy2-0.5*boxsidelength,'CMS POWHEGv2 NLO')
		for tex in [self._tevatron_theory_leg_tex,self._powheg_theory_leg_tex] :
			tex.SetTextSize(0.018)
			tex.SetTextAlign(12)
			tex.SetTextFont(42)
		#add the title
		self._title = TLatex(exp_x,1. - self._canvas.GetTopMargin() - 0.01,'top quark forward-backward asymmetry (parton-level)')
		self._title.SetTextSize(0.025)
		self._title.SetTextAlign(23)
		#draw the objects
		self._leg_tbox.Draw('SAME')
		self._leg_bbox.Draw('SAME')
		self._leg_nbox.Draw('SAME')
		self._gr_leg_stat.Draw('PE1 SAME')
		self._gr_leg_tot.Draw('PE1 SAME')
		for line in newlines :
			line.Draw('SAME')
		self._stat_tex.Draw('SAME')
		self._tot_tex.Draw('SAME')
		self._tevatron_theory_leg_box.Draw('SAME')
		self._powheg_theory_leg_box.Draw('SAME')
		self._tevatron_theory_leg_tex.Draw('SAME')
		self._powheg_theory_leg_tex.Draw('SAME')
		self._title.Draw('SAME')

	def plotOnCanv(self) :
		#start by counting up how many points we have total and numbering the entries
		self.__countNpoints__()
		#make the canvas
		self._canvas = TCanvas(plotname,plotname,900,900)
		self._canvas.SetTicky(0)
		self._canvas.SetTopMargin(0.025)
		self._canvas.SetBottomMargin(0.10)
		self._canvas.SetLeftMargin(0.025)
		self._canvas.SetRightMargin(self._canvas.GetLeftMargin())
		#make and plot the background histogram
		self.__makeAndPlotHistogram__()
		#make and plot the label and reference boxes for all the entries
		self.__makeAndPlotBoxes__()
		#make and plot the theory lines
		self.__makeAndPlotTheoryLines__()
		#make and plot the graphs
		self.__makeAndPlotGraphs__()
		#add the legend
		self.__addLegend__()
		#make and plot the dividers
		self.__makeAndPlotDividers__()
		#add the legend
		CMS_lumi.CMS_lumi(self._canvas,iPeriod,11)
		#Write the final canvas to the file
		self._canvas.Write()

# open the output file
outfile = TFile(outfilename,'recreate')

# initialize the plot
plotobj = Plot()
# add all the tevatron measurement categories to the plot
for cat in sorted(tevatron_afbs.keys(),reverse=True) : #make sure they get listed as 'lepton+jets', 'dilepton', 'combination'
	plotobj.addCategory(cat,tevatron_afbs[cat])
# add the theory value for comparisons with the tevatron measurements
plotobj.addTheory('tevatron',tevatron_theory.keys()[0],tevatron_theory[tevatron_theory.keys()[0]])
# add our measurement to the plot
plotobj.addMyAfb(my_afb)
# add the powheg generator value for comparison with our measurement
plotobj.addTheory('powheg',powheg_theory.keys()[0],powheg_theory[powheg_theory.keys()[0]])
# make the plot
plotobj.plotOnCanv()




# close the output file
outfile.Close()