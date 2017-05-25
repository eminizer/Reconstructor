#imports
from ROOT import *
from glob import glob
from datetime import date
import os
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--sample',  metavar='F', type='string', action='store', 
                              dest='sample') 
parser.add_option('--outtag',  metavar='F', type='string', action='store', 
                              default='', 
                              dest='outtag') ## name for output file
(options, args) = parser.parse_args()
sample = options.sample
outtag = options.outtag
if sample==None :
	optstring = raw_input('samplename (outfilename_append) >  ')
	optlist = optstring.split()
	if len(optlist)<1 :
		print 'must input sample name'
		exit()
	sample = optlist[0]
	if len(optlist)>1 :
		outtag = optlist[1]

if not os.path.isdir('../'+sample) :
	print 'sample name '+sample+' is not valid!'
	exit()

#open the input files
infiles = glob('../'+sample+'/aggregated_*.root')

#start up the output file
outfilename = 'kinfit_tuning_plots_'+sample+'_'+str(date.today())
if outtag!='' :
	outfilename+='_'+outtag
outfilename+='.root'
outfile = TFile(outfilename,'recreate')

#build the chain
fullchain = TChain('tree')
for f in infiles :
	print 'adding file '+str(f)+'...'
	fullchain.Add(f)

#declare the plots
#tmass_low=100.; tmass_high=350.; ntbins=75
tmass_low=100.; tmass_high=300.; ntbins=100
allhistos = []
t1leptM = TH1F('t1leptM','type-1 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t1leptM)
t1hadtM = TH1F('t1hadtM','type-1 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t1hadtM)
t2leptM = TH1F('t2leptM','type-2 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t2leptM)
t2hadtM = TH1F('t2hadtM','type-2 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t2hadtM)
t3leptM = TH1F('t3leptM','type-3 leptonic top mass, correct hypothesis; M_{top}^{lep} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t3leptM)
t3hadtM = TH1F('t3hadtM','type-3 hadronic top mass, correct hypothesis; M_{top}^{had} (GeV); fraction',ntbins,tmass_low,tmass_high); allhistos.append(t3hadtM)
for histo in allhistos :
	histo.SetDirectory(0)

#set up skimmed chain
com_cuts = 'fullselection==1 && ismatchable==1'
chain=fullchain.CopyTree(com_cuts)

#plot plots
weights = '35867.*weight*sf_pileup*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas'
print 'Drawing type-1 leptonic top mass...'
chain.Draw('leptcorprefitM>>t1leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==1)')
print 'Drawing type-1 hadronic top mass...'
chain.Draw('hadtcorprefitM>>t1hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==1)')
print 'Drawing type-2 leptonic top mass...'
chain.Draw('leptcorprefitM>>t2leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==2)')
print 'Drawing type-2 hadronic top mass...'
chain.Draw('hadtcorprefitM>>t2hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==2)')
print 'Drawing type-3 leptonic top mass...'
chain.Draw('leptcorprefitM>>t3leptM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3)')
print 'Drawing type-3 hadronic top mass...'
chain.Draw('hadtcorprefitM>>t3hadtM('+str(ntbins)+','+str(tmass_low)+','+str(tmass_high)+')','('+weights+')*('+com_cuts+' && eventTopology==3)')
print 'Done.'
for histo in allhistos :
	thisName = histo.GetName()
	#print 'histo name = %s'%(thisName) #DEBUG
	toAdd = gROOT.FindObject(thisName)
	#print '	object = %s'%(toAdd) #DEBUG
	histo.Add(toAdd)
	if histo.Integral()==0 :
		print thisName+' has an integral of zero'
	histo.Scale(1./histo.Integral())
	histo.SetMarkerStyle(20)
	histo.GetYaxis().SetTitleOffset(1.2)

#declare canvases
allcanvs = []
t1leptM_canv = TCanvas('t1leptM_canv','t1leptM_canv',1100,900); allcanvs.append(t1leptM_canv)
t1hadtM_canv = TCanvas('t1hadtM_canv','t1hadtM_canv',1100,900); allcanvs.append(t1hadtM_canv)
t2leptM_canv = TCanvas('t2leptM_canv','t2leptM_canv',1100,900); allcanvs.append(t2leptM_canv)
t2hadtM_canv = TCanvas('t2hadtM_canv','t2hadtM_canv',1100,900); allcanvs.append(t2hadtM_canv)
t3leptM_canv = TCanvas('t3leptM_canv','t3leptM_canv',1100,900); allcanvs.append(t3leptM_canv)
t3hadtM_canv = TCanvas('t3hadtM_canv','t3hadtM_canv',1100,900); allcanvs.append(t3hadtM_canv)
total_canv = TCanvas('total_canv','total_canv',1100,900)

colors = [kRed,kBlue,kGreen,kOrange,kMagenta,kCyan+2]
legEntries = ['type-1, leptonic','type-1, hadronic','type-2, leptonic','type-2, hadronic','type-3, leptonic','type-3, hadronic']

#draw plots on canvases and fit with gaussians, print results
lines = []
total_leg = TLegend(0.1,0.7,0.48,0.9)
for i in range(len(allhistos)) :
	histo = allhistos[i]
	print 'Doing '+str(allhistos[i].GetName())
	allcanvs[i].cd()
	histo.Draw()
	gausfunc1 = TF1('gausfunc1','gaus',150.,195.)
	histo.Fit(gausfunc1,"","",150.,195.)
	gmean = gausfunc1.GetParameter(1)
	gausfunc = TF1('gausfunc','gaus',gmean-20.,gmean+20.)
	#gausfunc = TF1('gausfunc','[0]*ROOT::Math::lognormal_pdf(x,[1],[2])',tmass_low,tmass_high)
	#gausfunc.SetParameters(10000000000.,172.5,0.5)
	gausfunc.SetLineColor(kRed)
	gausfunc.SetLineWidth(4)
	histo.Fit(gausfunc,"","",gmean-20.,gmean+20.)
	gausfunc.Draw('SAME')
	gmean = gausfunc.GetParameter(1); gwidth = gausfunc.GetParameter(2)
	lines.append('%s, Gaussian: %.1f +/- %.1f'%(allhistos[i].GetName(),gmean,gwidth))
	total_canv.cd()
	histo.SetMarkerColor(colors[i])
	histo.SetFillColor(colors[i])
	gausfunc.SetLineColor(colors[i])
	histo.Draw() if i==0 else histo.Draw('SAME')
	gausfunc.Draw('SAME')
	total_leg.AddEntry(histo,'#splitline{%s}{M_{t}=%.1f #pm %.1f}'%(legEntries[i],gmean,gwidth),'F')
total_canv.cd()
total_leg.Draw('SAME')

print '\n\n\nFit Results:'	
for line in lines :
	print line

#save plots and canvases
outfile.cd()
for histo in allhistos :
	histo.Write()
for canv in allcanvs :
	canv.Write()
total_canv.Write()
outfile.Write()
outfile.Close()