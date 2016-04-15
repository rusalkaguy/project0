# This program pulls and reads a genbank file by command line input of accession number,
# then formats and writes the relevant information in bed format.

# Import adds modules to your script from the python argv
# argv = argument variable that holds the arguments you pass to your script when you run it
from sys import argv 

# Defining Accession Number, File_name, and Database
accession_number = argv[1] # Upacks argv-> assigned to 1 variable you can work with
file_name = accession_number + '.gbk'
db = 'nucleotide' 

# Write Genbank File to BED file format: https://gist.github.com/brantfaircloth/893580
from Bio import SeqIO

import pdb

# Import regular expressions
import re

# for debugging output
import pprint
pp = pprint.PrettyPrinter(indent=0)

# Compile patterns for regular expressions
p = re.compile(r"""
    complement\(( # Start of a location
    .*)      # 0 or more(*) of anything(.)
    \)       # Trailing end parenthesis
    """, re.VERBOSE)

q = re.compile(r"""
    [0-9]+       # 1 or more(+) of any number
    \.\.         # .. \ denotes not anything (.)
    [0-9]+       # 1 or more(+) of any number
    """, re.VERBOSE)

q2 = re.compile(r"""
    ([0-9]+)       # 1 or more(+) of any number
    \.\.         # .. \ denotes not anything (.)
    ([0-9]+)       # 1 or more(+) of any number
    """, re.VERBOSE)


def main():
	for record in SeqIO.parse(open(file_name, "rU"), "genbank") : #record is genome
		# add chrom name using split and join('v')
		ucsc_chrom='v'.join(accession_number.split('.'))
		# get accession_number from record
		outf = open( ucsc_chrom +'.bed', 'w')
		gene_dict={} # Initialize outermost dictionary
		for feature in record.features: # record is genome
			if feature.type == 'gene':
				start = feature.location.start.position
			 	stop = feature.location.end.position
				try:
					gene_name=feature.qualifiers['gene'][0]
				except:
					# some features only have a locus tag
					gene_name = feature.qualifiers['locus_tag'][0]
				# Check to see if gene exits in dictionary
				if gene_name not in gene_dict: # if gene not in dict,
					gene_dict[gene_name]={'gene_name': gene_name} # set dictionary and store gene_name to process later
					#gene_dict[gene_name]['gene_name']=gene_name
				# Set gene_dict to key-value pair {gene_name:gene_def_dict} 
				gene_def_dict = gene_dict[gene_name]
				if feature.strand < 0:
					strand = "-"
				else:
					strand = "+"
				bed_line = ucsc_chrom + "\t{0}\t{1}\t{2}\t1000\t{3}\n".format(start, stop, gene_name, strand)
				outf.write(bed_line)
		outf.close()
	pp.pprint(gene_dict)

if __name__ == '__main__':
	main()
# Run with $ python write_to_bed.py NC_006273.2