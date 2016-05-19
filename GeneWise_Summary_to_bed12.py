# This program parses GeneWise Summary .txt files, and formats and outputs in 12 column .bed file.

# Import regular expressions
import re

import pprint
pp = pprint.PrettyPrinter(indent=0)
p=re.compile(r"""
	.+ 			# 1 or more(+) of anything
	_ex([0-9]+) # 1 or more(+) of any number after ex
	""", re.VERBOSE)

file_open = open ("HCMV-Merlin_annotation_details.txt", "r+")
file_read = file_open.readlines()
'''
Use Columns 1-4 to construct a 12 column bed. 
When there are exons, you get one line for the 'gene' where aa_seq = seq for whole gene 
and nt_seq = seq for whole gene, PLUS a line for each exon, where aa_seq is empty 
and nt_seq is the seq for just that exon. 
So, construct a single .bed file line from each set of GENE + GENE_ex*. 

The 4th  column, 'status' can be used for the score  and filtering
    Missing  = do not include gene in output bed file
    Partial = set BED score to 500
    Full = set BED score to 1000
'''
# Variables for column locations
gene_name_col=0
start_col=1
stop_col=2
status_col=3
aa_seq_col=4
strand='.'
gene_dict={}
cds_array = [] # make empty array to store strand, start, and stop info
# neccessary when more than one exon or CDS region exists for a gene
for line in file_read: # 1 gene per line
	line = line.split('#')[0]
	tab = line.split()  # 'tab' could be any any other variable
	# Split at tabs is default: i.e. empty brackets.
	# Insert space or comma in parenthesis () if that denotes separation.
	if len(tab) == 0:
		continue
	# Assign variables for items of interest
	gene_name = tab[gene_name_col]
	start_loc = tab[start_col]
	stop_loc = tab[stop_col]
	status = tab[status_col]
	aa_seq = tab[aa_seq_col]
	if status=='Missing':
		continue
	elif status=='Partial':
		score=500
	elif status=='Full':
		score = 1000

	s=p.search(gene_name)
	cds_def = {} 
	cds_array=[]# make empty array to store exon strand, start, and stop info
	# neccessary when more than one exon or CDS region exists for a gene
	if not s: # if no exons
		print 'No search hit'
		gene_dict[gene_name]={'gene_name':gene_name, 'start': start_loc, 'stop': stop_loc, 'score': score, 'CDS': [] }
		gene_def_dict=gene_dict[gene_name]

	else: # if exons present
		print 'search found. Span=', s.span(), ' Group=', s.group()
		gene_exon_name=str(s.group())
		gene_name_str_list = gene_exon_name.split("_")
		gene_name=gene_name_str_list[0]
		if gene_name in gene_dict:
			# fetch existing gene record from dictionary
			gene_def_dict = gene_dict[gene_name]
		else:
			# create gene in dictionary
			gene_dict[gene_name]=gene_def_dict
		# extract and create exon_def
		exon_num= gene_name_str_list[1]
		exon_def= {'exon_num': exon_num, 'strand':strand, 'start': start_loc, 'stop': stop_loc}
		cds_def[exon_num]=exon_def
		#cds_array.append(exon_def)
		#print cds_array
		gene_def_dict['CDS'].append(exon_def)
		gene_dict[gene_name]=gene_def_dict
'''
for key in gene_dict:


def format_bed12_line(gene_def_dict):
item_rgb='0,0,0'
    # columns 7-8 thickStart and thickStop
    cds_def=gene_def_dict['CDS']
    cds_exon_count= len(cds_def)
'''
pp.pprint(gene_dict)


# Close the original opened text file
file_open.close;