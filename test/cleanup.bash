#!/bin/bash

rm -rf notneeded
rm -rf commands.cmd
rm -rf tarball.tgz
rm -rf garbage*root
mkdir output
mv output*.log output
mv condor-* output
