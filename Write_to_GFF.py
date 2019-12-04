#!/usr/bin/env python
# Import adds modules to your script from the python argv
# argv = argument variable that holds the arguments you pass to your script when you run it
from sys import argv 

# Defining Accession Number, File_name, and Database
accession_number= 'NC_006273'
if len(argv) > 1:
    accession_number=argv[1]
file_name = accession_number + '.gbk'
db = 'nucleotide' 

# Write Genbank File to GFF file format: http://biopython.org/wiki/GFF_Parsing
from BCBio import GFF
from Bio import SeqIO
 
in_file = file_name
out_file = accession_number + '.gff'
in_handle = open(in_file)
out_handle = open(out_file, "w")
 
GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
 
in_handle.close()
out_handle.close()
