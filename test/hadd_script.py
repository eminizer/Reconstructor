import os, glob
from ROOT import *

sample_names = []
sample_names.append('mcatnlo_qq_semilep_TT')
#sample_names.append('mcatnlo_qq_semilep_TT_just_MC')
sample_names.append('mcatnlo_gg_semilep_TT')
#sample_names.append('mcatnlo_gg_semilep_TT_just_MC')
sample_names.append('mcatnlo_dilep_TT')
sample_names.append('mcatnlo_had_TT')
sample_names.append('WJets_HT-200to400')
sample_names.append('WJets_HT-400to600')
sample_names.append('WJets_HT-600to800')
sample_names.append('WJets_HT-800to1200')
#sample_names.append('WJets_HT-2500toInf')
sample_names.append('DYJets_M-50_HT-100to200')
sample_names.append('DYJets_M-50_HT-200to400')
sample_names.append('DYJets_M-50_HT-400to600')
#sample_names.append('DYJets_M-50_HT-600toInf')
sample_names.append('ST_s-c')
sample_names.append('ST_t-c_top')
sample_names.append('ST_tW-c_top')
sample_names.append('ST_t-c_antitop')
sample_names.append('ST_tW-c_antitop')
sample_names.append('WW_to_2L_2Nu')
sample_names.append('WW_to_L_Nu_2Q')
sample_names.append('WZ_to_L_Nu_2Q')
sample_names.append('WZ_to_L_3Nu')
sample_names.append('WZ_to_2L_2Q')
sample_names.append('WZ_to_3L_Nu')
sample_names.append('ZZ_to_2L_2Nu')
sample_names.append('ZZ_to_2L_2Q')
sample_names.append('ZZ_to_4L')
sample_names.append('SingleEl_Run2016Bv2')
sample_names.append('SingleEl_Run2016C')
sample_names.append('SingleEl_Run2016D')
sample_names.append('SingleMu_Run2016Bv2')
sample_names.append('SingleMu_Run2016C')
sample_names.append('SingleMu_Run2016D')

for name in sample_names :
    print 'doing '+name
    
#    #make directories
#    os.system('mkdir '+name)
    
    os.chdir(name)

#    #get rid of old files
#    os.system('bash cleanup.bash')
#    os.system('rm -rf output *.root')
#    os.system('mv ana.listOfJobs_all ana.listOfJobs')

#    #clean out the input and ana.listOfJobs, even
#    os.system('rm -rf input.txt')
#    os.system('rm -rf ana.listOfJobs')
    
#    #copy scripts
#    os.system('cp ../grid_sub.csh .; cp ../cleanup.bash .')
    
#    #make input file
#    os.system('rm -rf input.txt')
#    directory = raw_input('nTuple directory for '+name+': ')
#    os.system('python ../make_ttree_input_file.py --directory '+directory)
    
#    #make new ana.listOfJobs
#    os.system('rm -rf ana.listOfJobs')
#    nJobs = raw_input('number of jobs for '+name+': ')
#    xSec  = raw_input('cross section for '+name+': ')
#    generator = None
#    if name.lower().find('powheg')!=-1 :
#    	generator='powheg'
#    elif name.lower().find('mcatnlo')!=-1 :
#    	generator='mcatnlo'
#    else :
#	    generator = raw_input('MC generator for '+name+': ')
#    cmd = 'python ../make_list_of_jobs.py --n_jobs '+nJobs+' --name '+name+' --on_grid yes'
#    cmd+= ' --generator '+generator+' --xSec '+xSec
#    os.system(cmd)

#    #make list of failed jobs
#    os.system('python ../make_failed_job_list.py')
    
#    #submit jobs
#    os.system('tcsh grid_sub.csh')

    #skim files
    os.system('rm -rf *_skim_tree.root')
    filelist = glob.glob('*_tree.root')
    for i in range(len(filelist)) :
        print ' '+str(i)+': '+str(filelist[i])
        f = TFile(filelist[i]); t = f.Get('tree')
        newname = filelist[i].replace('_tree.root','')+'_skim_tree.root'
        newFile = TFile(newname,'recreate')
        newTree = t.CopyTree('weight!=0.')
        #newTree = t.CopyTree('fullselection==1')
        newTree.Write()
        newFile.Close()
    i = 0
    while i<len(filelist) :
        if filelist[i].find('JES')!=-1 or filelist[i].find('JER')!=-1 or filelist[i].find('skim')!=-1 :
            filelist.pop(i)
        else :
            i+=1
    cmd = 'hadd -f '+name+'_skim_all.root '+name+'_?_skim_tree.root'
    if len(filelist) > 10 :
        cmd += ' '+name+'_??_skim_tree.root'
        if len(filelist) > 100 :
            cmd += ' '+name+'_???_skim_tree.root'
            if len(filelist) > 1000 :
                cmd += ' '+name+'_????_skim_tree.root'
    os.system(cmd)
    if name.find('Run2012')==-1 :
        cmd = 'hadd -f '+name+'_JES_up_skim_all.root '+name+'_JES_up_*_skim_tree.root'
        os.system(cmd)
        cmd = 'hadd -f '+name+'_JES_down_skim_all.root '+name+'_JES_down_*_skim_tree.root'
        os.system(cmd)
        cmd = 'hadd -f '+name+'_JER_up_skim_all.root '+name+'_JER_up_*_skim_tree.root'
        os.system(cmd)
        cmd = 'hadd -f '+name+'_JER_down_skim_all.root '+name+'_JER_down_*_skim_tree.root'
        os.system(cmd)
    os.system('mv *_all.root ../total_ttree_files')
    os.system('rm -rf *_skim_tree.root')

#    #skim files (e/mu selection)
#    filelist = glob.glob('*_tree.root')
#    i = 0
#    while i<len(filelist) :
#        if filelist[i].find('JES')!=-1 or filelist[i].find('JER')!=-1 or filelist[i].find('skim')!=-1 :
#            filelist.pop(i)
#        else :
#            i+=1
#    for i in range(len(filelist)) :
#        print ' '+str(i)+''
#        f = TFile(filelist[i]); t = f.Get('tree')
#        cutstring = 'mu_trigger==1 && muon1_isLoose==1 && muon1_pt>40. && abs(muon1_eta)<2.4'
#        cutstring+= ' && (muon1_relPt>25. || muon1_dR>0.5)'
#        cutstring+= ' && ele1_pt>40. && abs(ele1_eta)<2.4'
#        cutstring+= ' && (ele1_relPt>25. || ele1_dR>0.5)'
#        newTree = t.CopyTree(cutstring)
#        newname = filelist[i].replace('_tree.root','')+'_emu_skim_tree.root'
#        newFile = TFile(newname,'recreate')
#        newTree.Write()
#        newFile.Close()
#    cmd = 'hadd -f '+name+'_emu_skim_all.root '+name+'_?_emu_skim_tree.root'
#    if len(filelist) > 10 :
#        cmd += ' '+name+'_??_emu_skim_tree.root'
#        if len(filelist) > 100 :
#            cmd += ' '+name+'_???_emu_skim_tree.root'
#            if len(filelist) > 1000 :
#                cmd += ' '+name+'_????_emu_skim_tree.root'
#    os.system(cmd)
#    os.system('mv *_all.root ../total_ttree_files')
#    os.system('rm -rf *emu_skim*.root')

#    #skim files (e+/e- selection)
#    filelist = glob.glob('*_tree.root')
#    i = 0
#    while i<len(filelist) :
#        if filelist[i].find('JES')!=-1 or filelist[i].find('JER')!=-1 or filelist[i].find('skim')!=-1 :
#            filelist.pop(i)
#        else :
#            i+=1
#    for i in range(len(filelist)) :
#        print ' '+str(i)+''
#        f = TFile(filelist[i]); t = f.Get('tree')
#        newname = filelist[i].replace('_tree.root','')+'_ee_skim_tree.root'
#        newFile = TFile(newname,'recreate')
#        ele1_tag_selec = '(ele1_pt>40. && abs(ele1_eta)<2.4 && ele1_isLoose==1 && (ele1_relPt>25. || ele1_dR>0.5))'
#        ele1_probe_selec = '(ele1_pt>40. && abs(ele1_eta)<2.4 && (ele1_relPt>25. || ele1_dR>0.5))'
#        ele2_tag_selec = '(ele2_pt>40. && abs(ele2_eta)<2.4 && ele2_isLoose==1 && (ele2_relPt>25. || ele2_dR>0.5))'
#        ele2_probe_selec = '(ele2_pt>40. && abs(ele2_eta)<2.4 && (ele2_relPt>25. || ele2_dR>0.5))'
#        cutstring = '(('+ele1_tag_selec+' && '+ele2_probe_selec+') || ('+ele2_tag_selec+' && '+ele1_probe_selec+'))'
#        #muon 1 and 2 rejection
#        cutstring += ' && (muon1_pt<40. || abs(muon1_eta)>2.4 || muon1_isLoose!=1 || (muon1_relPt<25. && muon1_dR<0.5))'
#        cutstring += ' && (muon2_pt<40. || abs(muon2_eta)>2.4 || muon2_isLoose!=1 || (muon2_relPt<25. && muon2_dR<0.5))'
#        #require two electrons to have opposite charges
#        cutstring += ' && (ele1_Q+ele2_Q==0)'
#        #other leptonic side cuts
#        cutstring+=' && lepW_pt[0]>50.'
#        #other jet requirements
#        cutstring+=' && lepb_pt>50. && hadt_pt>50. && max(lepb_pt,hadt_pt)>150.'
#        cutstring+=' && (lepb_csv>0.244 || hadt_csv>0.244)'
#        newTree = t.CopyTree(cutstring)
#        newTree.Write()
#        newFile.Close()
#    cmd = 'hadd -f '+name+'_ee_skim_all.root '+name+'_?_ee_skim_tree.root'
#    if len(filelist) > 10 :
#        cmd += ' '+name+'_??_ee_skim_tree.root'
#        if len(filelist) > 100 :
#            cmd += ' '+name+'_???_ee_skim_tree.root'
#            if len(filelist) > 1000 :
#                cmd += ' '+name+'_????_ee_skim_tree.root'
#    os.system(cmd)
#    os.system('mv *_all.root ../total_ttree_files')
#    os.system('rm -rf *_ee_skim*.root')

    os.chdir('..')
#os.chdir('total_ttree_files')
#os.system('hadd -f SingleMu_Run2012_all.root SingleMu_Run2012*_all.root')
#cmd = 'hadd -f Powheg_semilep_TT_all.root Powheg_qq_semilep_TT_all.root'
#cmd += ' Powheg_qq_semilep_TT_SC_all.root Powheg_gg_semilep_TT_all.root Powheg_gg_semilep_TT_SC_all.root'
#os.system(cmd)
#os.chdir('..')
