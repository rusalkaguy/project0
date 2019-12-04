#!/usr/bin/bash
# 
# download genbank record
#
accession=NC_006273
if [ ! -z "$1" ]; then echo "accession=$1"; accession=$1; fi
wget -nv -O ${accession}.gbk "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=genbank&id=${accession}&conwithfeat=on&withparts=on&hide-cdd=on"
wget -nv -O ${accession}.fa "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=${accession}&conwithfeat=on&withparts=on&hide-cdd=on"

