#!/usr/bin/bash
# 
# download genbank record
#
accession=NC_006273
if [ ! -z "$1" ]; then acccession=$1; fi

wget -O ${accession}.gbk "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=genbank&id=${accession}&conwithfeat=on&withparts=on&hide-cdd=on"
