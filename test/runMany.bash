#!/bin/bash

## Used to run executables.
if [ -z "$1" ] ; then
  echo Usage: $0 RunManyPath commandfile section# cluster#
  exit;
fi

#some debugging
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node

# environment setup
source /cvmfs/cms.cern.ch/cmsset_default.sh

# save current directory
export CURRENTDIR=`pwd -P`

export RUNMANYPATH=$1
shift
export COMMANDFILE=$1
shift
export CMSSW_LOC=$1
shift
export CMSSWVER=$1
shift
export SECTION=$1
shift
export Cluster=$1
export Process=$SECTION
export JID=JID_$Cluster\_$Process
shift
export CopyDirectory=$1

# figure out which directory runManySections.py is in
export RUNMANYSCRIPT=$RUNMANYPATH"/runManySections.py"

echo script, $RUNMANYSCRIPT

# Before doing anything else, untar the tarball if there is one
export TARBALL=`${RUNMANYSCRIPT} --tarball  ${COMMANDFILE}`
if [ -n "$TARBALL" ] ; then
    echo Tarball :$TARBALL:
    # Do we want to untar this in a subdirectory?
    export UNTARDIR=`${RUNMANYSCRIPT} --untardir  ${COMMANDFILE}`
    if [ -n $UNTARDIR ] ; then
        mkdir -p $UNTARDIR
        cd $UNTARDIR
    fi # untarred
    tar xzf $CURRENTDIR/$TARBALL
    cd $CURRENTDIR
fi # tarball

#set scram architecture
export SCRAM_ARCH=$CMSSW_LOC
echo "Scram architecture: "$SCRAM_ARCH

#make a new release area
echo "Pulling CMSSW version: "$CMSSWVER
scram project $CMSSWVER
echo "CD-ing to CMSSW release directory "
cd $CMSSWVER
echo "current directory: `pwd`"

# cmsenv
echo "doing cmsenv"
eval `scramv1 runtime -sh`
echo "CMSSW_BASE: "$CMSSW_BASE
cd $CURRENTDIR
echo "current directory: `pwd`"

# If we are on a system that has the 'storage.xml' file, tell
# runManySections.py use it.
export REPLACE=' '
if [ -f "$CMS_PATH/SITECONF/local/PhEDEx/storage.xml" ]; then
   export REPLACE="--replaceFilelist";
fi

# get the command and log 
export  COMMAND=`${RUNMANYSCRIPT} --command ${REPLACE} ${COMMANDFILE} ${SECTION}`
export  LOG=$CURRENTDIR"/"`${RUNMANYSCRIPT} --log  ${COMMANDFILE} ${SECTION}`


# setup environment
hostname > $LOG
pwd >> $LOG
echo $COMMAND  >> $LOG 
(time $COMMAND) >> $LOG 2>&1

if [ -n "$CopyDirectory" ] ; then
   cd $CURRENTDIR
   export COPYCOMMAND=`${RUNMANYSCRIPT} --copyCommand ${CopyDirectory} ${COMMANDFILE}`
   echo copy command ${COPYCOMMAND} >> $LOG 2>&1
   eval $COPYCOMMAND
fi
