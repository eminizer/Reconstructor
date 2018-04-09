from ROOT import *
from glob import glob
from datetime import date
import os
from math import sqrt
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--sample',  metavar='F', type='string', action='store', 
                              dest='sample',  help='')
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
	print 'sample name %s is not valid!'%(sample)
	exit()

#open the input ttbar files in a chain
allfiles = glob('../'+sample+'/aggregated_'+sample+'*.root')
filelist = [filename for filename in allfiles if filename.find('JES')==-1 and filename.find('JER')==-1]
#print 'filelist (%d files) = %s'%(len(filelist),filelist) #DEBUG
chain = TChain('tree')
for filename in filelist :
	chain.Add(filename)

#draw into histograms to get numbers of events
fully_selected_cut = 'eventType<2 && fullselection==1'
matchable_cut = 'ismatchable==1'
not_matchable_cut = 'ismatchable==0'
correct_cut = 'iscorrect==1'
matched_postfit_cut = 'ismatchedpostfit==1'
#type-1
chain.Draw('cstar>>t1_selected(40,-1.,1.)','('+fully_selected_cut+' && eventTopology==1)')
t1_selected = (gROOT.FindObject('t1_selected')).Integral()
#print 't1_selected = %d'%(t1_selected) #DEBUG
chain.Draw('cstar>>t1_matchable(40,-1.,1.)','('+fully_selected_cut+' && '+matchable_cut+' && eventTopology==1)')
t1_matchable = (gROOT.FindObject('t1_matchable')).Integral()
chain.Draw('cstar>>t1_correct(40,-1.,1.)','('+fully_selected_cut+' && '+matchable_cut+' && '+correct_cut+' && eventTopology==1)')
t1_correct = (gROOT.FindObject('t1_correct')).Integral()
chain.Draw('cstar>>t1_postfit(40,-1.,1.)','('+fully_selected_cut+' && '+not_matchable_cut+' && '+matched_postfit_cut+' && eventTopology==1)')
t1_postfit = (gROOT.FindObject('t1_postfit')).Integral()
#type-2
chain.Draw('cstar>>t2_selected(40,-1.,1.)','('+fully_selected_cut+' && eventTopology==2)')
t2_selected = (gROOT.FindObject('t2_selected')).Integral()
#print 't2_selected = %d'%(t2_selected) #DEBUG
chain.Draw('cstar>>t2_matchable(40,-1.,1.)','('+fully_selected_cut+' && '+matchable_cut+' && eventTopology==2)')
t2_matchable = (gROOT.FindObject('t2_matchable')).Integral()
chain.Draw('cstar>>t2_correct(40,-1.,1.)','('+fully_selected_cut+' && '+matchable_cut+' && '+correct_cut+' && eventTopology==2)')
t2_correct = (gROOT.FindObject('t2_correct')).Integral()
chain.Draw('cstar>>t2_postfit(40,-1.,1.)','('+fully_selected_cut+' && '+not_matchable_cut+' && '+matched_postfit_cut+' && eventTopology==2)')
t2_postfit = (gROOT.FindObject('t2_postfit')).Integral()
#type-3
chain.Draw('cstar>>t3_selected(40,-1.,1.)','('+fully_selected_cut+' && eventTopology==3)')
t3_selected = (gROOT.FindObject('t3_selected')).Integral()
#print 't3_selected = %d'%(t3_selected) #DEBUG
chain.Draw('cstar>>t3_matchable(40,-1.,1.)','('+fully_selected_cut+' && '+matchable_cut+' && eventTopology==3)')
t3_matchable = (gROOT.FindObject('t3_matchable')).Integral()
chain.Draw('cstar>>t3_correct(40,-1.,1.)','('+fully_selected_cut+' && '+matchable_cut+' && '+correct_cut+' && eventTopology==3)')
t3_correct = (gROOT.FindObject('t3_correct')).Integral()
chain.Draw('cstar>>t3_postfit(40,-1.,1.)','('+fully_selected_cut+' && '+not_matchable_cut+' && '+matched_postfit_cut+' && eventTopology==3)')
t3_postfit = (gROOT.FindObject('t3_postfit')).Integral()

#calculate numbers and percentages and put 'em in a dictionary
values = {}
t1_matchable_percent = 100.*t1_matchable/t1_selected
t1_matchable_error = t1_matchable_percent*sqrt(1./t1_matchable+1./t1_selected)
t1_correct_percent = 100.*t1_correct/t1_matchable
t1_correct_error = t1_correct_percent*sqrt(1./t1_correct+1./t1_matchable)
values['t1'] = {'matchable':'%s(%s)'%(t1_matchable_percent,t1_matchable_error),'correct':'%s(%s)'%(t1_correct_percent,t1_correct_error)}
t2_matchable_percent = 100.*t2_matchable/t2_selected
t2_matchable_error = t2_matchable_percent*sqrt(1./t2_matchable+1./t2_selected)
t2_correct_percent = 100.*t2_correct/t2_matchable
t2_correct_error = t2_correct_percent*sqrt(1./t2_correct+1./t2_matchable)
values['t2'] = {'matchable':'%s(%s)'%(t2_matchable_percent,t2_matchable_error),'correct':'%s(%s)'%(t2_correct_percent,t2_correct_error)}
t3_matchable_percent = 100.*t3_matchable/t3_selected
t3_matchable_error = t3_matchable_percent*sqrt(1./t3_matchable+1./t3_selected)
t3_correct_percent = 100.*t3_correct/t3_matchable
t3_correct_error = t3_correct_percent*sqrt(1./t3_correct+1./t3_matchable)
values['t3'] = {'matchable':'%s(%s)'%(t3_matchable_percent,t3_matchable_error),'correct':'%s(%s)'%(t3_correct_percent,t3_correct_error)}

#save lines of .tex to a file 
print 'Table lines:'
outname = 'kinfit_performance_table_lines_'+sample+'_'+str(date.today())
if outtag!='' :
	outname+='_'+outtag
outname+='.txt'
newline = 'Type-1 & %d & %d & %s & %d & %s & %d	\\\\'%(t1_selected,t1_matchable,values['t1']['matchable'],t1_correct,values['t1']['correct'],t1_postfit)
print newline
os.system('echo "'+newline+'" >> '+outname+'')
newline = 'Type-2 & %d & %d & %s & %d & %s & %d	\\\\'%(t2_selected,t2_matchable,values['t2']['matchable'],t2_correct,values['t2']['correct'],t2_postfit)
print newline
os.system('echo "'+newline+'" >> '+outname+'')
newline = 'Type-3 & %d & %d & %s & %d & %s & %d	\\\\'%(t3_selected,t3_matchable,values['t3']['matchable'],t3_correct,values['t3']['correct'],t3_postfit)
print newline
os.system('echo "'+newline+'" >> '+outname+'')
print 'Done :3c'