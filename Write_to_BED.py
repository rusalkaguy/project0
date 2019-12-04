#!/usr/bin/env python
# Pull, parse, and write genbank data in bed6 format.

# Import adds modules to your script from the python argv
# argv = argument variable that holds the arguments you pass to your script when you run it
from sys import argv 

# Defining Accession Number, File_name, and Database
accession_number = argv[1] # Upacks argv so that it gets assigned to 1 variable you can work with
file_name = accession_number + '.gbk'
db = 'nucleotide' 

# Write Genbank File to BED file format: https://gist.github.com/brantfaircloth/893580
from Bio import SeqIO

import pdb

def main():

	for record in SeqIO.parse(open(file_name, "rU"), "genbank") : #record is genome
		# add chrom name using split and join('v')
		ucsc_chrom='v'.join(accession_number.split('.'))
		# get accession_number from record
		outf = open( ucsc_chrom +'.bed', 'w')
		for feature in record.features:
			if feature.type == 'gene':
				start = feature.location.start.position
				stop = feature.location.end.position
				try:
					name = feature.qualifiers['gene'][0]
				except:
					# some features only have a locus tag
					name = feature.qualifiers['locus_tag'][0]
				if feature.strand < 0:
					strand = "-"
				else:
					strand = "+"
				bed_line = ucsc_chrom + "\t{0}\t{1}\t{2}\t1000\t{3}\n".format(start, stop, name, strand)
				outf.write(bed_line)
		outf.close()

if __name__ == '__main__':
	main()

# Run with $ python Write_to_BED.py NC_006273.2
