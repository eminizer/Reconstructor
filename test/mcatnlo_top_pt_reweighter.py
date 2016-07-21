import os, sys, glob
from ROOT import *
import math
from array import array
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--lodir',  metavar='F', type='string', action='store', 
                              default='/eos/uscms/store/user/eminizer/TT_CT10_TuneZ2star_8TeV-powheg-tauola/TTBar_Powheg_v1/151116_191806/0000', 
                              dest='lodir',  help='') ## Sets directory of leading-order files
parser.add_option('--nlodir', metavar='F', type='string', action='store', 
                              default='/eos/uscms/store/user/eminizer/TT_8TeV-mcatnlo/TTBar_MCatNLO_v3/160608_184637/0000', 
                              dest='nlodir', help='') ## Sets directory of next-to-leading-order files
parser.add_option('--outname', metavar='F', type='string', action='store', 
                              default='pt_distribution_comparisons', 
                              dest='outname', help='') ## Name of output file
parser.add_option('--total_jobs', metavar='F', type='int', action='store', 
                              default=1, 
                              dest='total_jobs', help='') ## Total number of jobs
parser.add_option('--n_job', metavar='F', type='int', action='store', 
                              default=0, 
                              dest='n_job', help='') ## Number index of this job
(options, args) = parser.parse_args()

#Set up handles and labels for use
#MC GenParticle variables
GenHandle = Handle( "vector<reco::GenParticle>" )
GenLabel = ( "prunedGenParticles", "" )
#MC GenEvent variables
GenEventHandle = Handle("GenEventInfoProduct")
GenEventLabel  = ("generator","")

#Function that checks if a particle is a top, and adds its info to the appropriate histograms
def add_to_histos(event,is_nlo) :
    #find the top quarks
    event.getByLabel(GenLabel,GenHandle)
    if not GenHandle.isValid() :
        print 'INVALID HANDLE!!'
        return -1.
    evtweight = 1.0
    if is_nlo :
        event.getByLabel(GenEventLabel,GenEventHandle) 
        if not GenEventHandle.isValid() :
            print 'INVALID GEN EVENT HANDLE!!'
            return -1.
        GenEvent = GenEventHandle.product()
        if GenEvent.weight()<0 :
            evtweight*=-1
        nlo_weight[0] = evtweight
    GenParticles = GenHandle.product()
    foundtop = False; foundatop = False
    for ig in GenParticles :
        if foundtop and foundatop :
            if is_nlo :
                nlo_tree.Fill()
            else :
                lo_tree.Fill()
            return 0.
        istop = ig.pdgId() == 6 and ig.status() == 3
        isatop = ig.pdgId() == -6 and ig.status() == 3
        if istop :
            foundtop = True
            tpt = ig.pt()
            if is_nlo :
                nlo_tpt[0] = tpt
                nlo_pt_hist.Fill(tpt,evtweight)
            else :
                lo_tpt[0] = tpt
                lo_pt_hist.Fill(tpt,evtweight)
        elif isatop :
            foundatop = True
            atpt = ig.pt()
            if is_nlo :
                nlo_atpt[0] = atpt
                nlo_pt_hist.Fill(atpt,evtweight)
            else :
                lo_atpt[0] = atpt
                lo_pt_hist.Fill(atpt,evtweight)
    return -1.

#Make output file
outfilename = options.outname
if options.total_jobs > 1 :
    outfilename+='_'+str(options.n_job)
if not outfilename.endswith('.root') :
    outfilename+='.root'
outfile = TFile(outfilename,'recreate')

#Set up the TTrees
lo_tree = TTree('lo_tree','lo_tree') 
lo_tpt = array('d',[-1.]); lo_tree.Branch('lo_tpt',lo_tpt,'lo_tpt/D')
lo_atpt = array('d',[-1.]); lo_tree.Branch('lo_atpt',lo_atpt,'lo_atpt/D')
nlo_tree = TTree('nlo_tree','nlo_tree')
nlo_tpt = array('d',[-1.]); nlo_tree.Branch('nlo_tpt',nlo_tpt,'nlo_tpt/D')
nlo_atpt = array('d',[-1.]); nlo_tree.Branch('nlo_atpt',nlo_atpt,'nlo_atpt/D')
nlo_weight = array('d',[0]); nlo_tree.Branch('nlo_weight',nlo_weight,'nlo_weight/D')

#Set up histograms
PTBINS = array('d',[0.,30.,45.,60.,75.,90.,105.,120.,135.,150.,165.,180.,200.,225.,250.,285.,310.,360.,410.,485.,560.,650.,1000.])
NBINS = len(PTBINS)-1
lo_pt_hist  = TH1D('lo_pt_hist', 'Generated top quark p_{T} distributions; p_{T} (GeV); ',NBINS,PTBINS)
lo_pt_hist.SetMarkerStyle(25); lo_pt_hist.SetLineWidth(4); lo_pt_hist.SetLineColor(kRed)
nlo_pt_hist = TH1D('nlo_pt_hist','',NBINS,PTBINS)
nlo_pt_hist.SetMarkerStyle(25); nlo_pt_hist.SetLineWidth(4); nlo_pt_hist.SetLineColor(kBlue)
sf_plot = TH1D('sf_plot','NLO/LO cross section scalefactor; top quark p_{T} (GeV); #frac{NLO events}{LO events}',NBINS,PTBINS)
sf_plot.SetLineWidth(3)

#Read in input files
LO_FILES  = glob.glob(options.lodir+'/*.root')
NLO_FILES = glob.glob(options.nlodir+'/*.root')
lo_events  = Events(LO_FILES)
nlo_events = Events(NLO_FILES)

print 'Set up complete, looping. . .'

#for every leading order event
print ' Doing LO events'
realcount = 0
count = 0
totalevents = lo_events.size()
for event in lo_events :
    #check the grid split
    realcount+=1
    if ((realcount-1)-options.n_job) % options.total_jobs != 0 :
        continue
    check = add_to_histos(event,False)
    if check<0 :
        print 'COULDN\'T FIND A TOP AND ANTITOP IN THIS EVENT'
    count+=1
    percent_done = 100.*count/(totalevents/options.total_jobs)
    if count%10000 == 0 :
        print '     count at %d of %d (%.2f%% done)'%(count,totalevents/options.total_jobs,percent_done)
print ' Done'

#for every NLO event
print ' Doing NLO events'
realcount = 0
count = 0
totalevents = nlo_events.size()
for event in nlo_events :
    #check the grid split
    realcount+=1
    if ((realcount-1)-options.n_job) % options.total_jobs != 0 :
        continue
    check = add_to_histos(event,True)
    if check<0 :
        print 'COULDN\'T FIND A TOP AND ANTITOP IN THIS EVENT'
    count+=1
    percent_done = 100.*count/(totalevents/options.total_jobs)
    if count%10000 == 0 :
        print '     count at %d of %d (%.2f%% done)'%(count,totalevents/options.total_jobs,percent_done)
print ' Done'

print 'Done looping.'

#Write the unscaled histograms
outfile.cd()
lo_pt_hist.Write()
nlo_pt_hist.Write()

#Normalize the p_T plots
lo_pt_hist.Sumw2()
nlo_pt_hist.Sumw2()
lo_pt_hist.Scale(19748.*245.8/lo_pt_hist.Integral())
nlo_pt_hist.Scale(19748.*245.8/nlo_pt_hist.Integral())

#Fill the scalefactor plot
for i in range(1,NBINS+1) :
    nlo = nlo_pt_hist.GetBinContent(i); lo = lo_pt_hist.GetBinContent(i)
    nlo_err = nlo_pt_hist.GetBinError(i); lo_err = lo_pt_hist.GetBinError(i)
    if nlo>0 and lo>0 :
        content = nlo/lo
        err = content*sqrt((nlo_err/nlo)**2 + (lo_err/lo)**2)
        sf_plot.SetBinContent(i,content)
        sf_plot.SetBinError(i,err)

#Write the scalefactor plot
outfile.cd()
sf_plot.Write()

allcanvases = []
#Make canvases and legends and plot plots
generated_top_pt_canv = TCanvas('generated_top_pt_canv','generated_top_pt_canv',1200,900); allcanvases.append(generated_top_pt_canv)
generated_top_pt_leg = TLegend(0.62,0.67,0.9,0.9)
generated_top_pt_leg.AddEntry(lo_pt_hist,'Leading Order (Powheg)','L')
generated_top_pt_leg.AddEntry(nlo_pt_hist,'Next to Leading Order (MC@NLO)','L')
generated_top_pt_canv.cd()
lo_pt_hist.Draw('HIST'); nlo_pt_hist.Draw("HIST SAME"); generated_top_pt_leg.Draw()
sf_plot_canv = TCanvas('sf_plot_canv','sf_plot_canv',1200,900); allcanvases.append(sf_plot_canv)
sf_plot_canv.cd()
sf_plot.Draw()

#Save all the canvases
outfile.cd()
for canv in allcanvases :
    canv.Write()

#Write the ttrees
outfile.cd()
lo_tree.Write()
nlo_tree.Write()
