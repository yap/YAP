#!/bin/bash

# Download ROOT archive
wget https://root.cern.ch/download/root_v6.06.08.Linux-ubuntu14-x86_64-gcc4.8.tar.gz -O root.tar.gz

# Extract ROOT archive in this directory
tar xzf root.tar.gz 

# !! Remember to `. ./root/bin/thisroot.sh` in the .travis.yml file !!
