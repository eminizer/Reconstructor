#! /bin/sh
tar czvf tarball.tgz ./input.txt -C ../../python/ . -C ../test/other_input_files/ .
voms-proxy-init --voms cms

/uscms_data/d3/eminizer/runManySections/runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball.tgz \ana.listOfJobs commands.cmd
/uscms_data/d3/eminizer/runManySections/runManySections.py --submitCondor commands.cmd
condor_q eminizer
