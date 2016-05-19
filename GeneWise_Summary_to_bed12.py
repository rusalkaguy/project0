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

# Parse GeneWise Summary Text File
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
		#print 'No search hit'
		gene_dict[gene_name]={'gene_name':gene_name, 'start': start_loc, 'stop': stop_loc, 'score': score, 'CDS': [], 'strand':strand, }
		gene_def_dict=gene_dict[gene_name]

	else: # if exons present
		#print 'search found. Span=', s.span(), ' Group=', s.group()
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
		exon_def= {'exon_num': exon_num, 'start': start_loc, 'stop': stop_loc}
		cds_def[exon_num]=exon_def
		#cds_array.append(exon_def)
		#print cds_array
		gene_def_dict['CDS'].append(exon_def)
		gene_dict[gene_name]=gene_def_dict
# Debugging of parseing
pp.pprint(gene_dict)
'''
When there are exons, you get one line for the 'gene' where aa_seq = seq for whole gene 
and nt_seq = seq for whole gene, PLUS a line for each exon, where aa_seq is empty 
and nt_seq is the seq for just that exon. 
So, construct a single .bed file line from each set of GENE + GENE_ex*. 
'''
# Format Output
chrom='NC_006273'

def format_bed12_line(gene_dict):
	bed6_str=chrom+'\t'+gene_dict['start']+'\t'+gene_dict['stop']+'\t'+ gene_dict['gene_name']+'\t'+str(score)+'\t'+'strand'
	print bed6_str
	'''
	# columns 7-8 thickStart and thickStop
	cds_def=gene_dict['CDS']
	print cds_def
	cds_exon_count= len(cds_def)
	if cds_exon_count>0 :
		# compute blockCount    blockSizes  blockStarts
		exon_list=[cds_def]# list of CDS, where CDS is a list of exon_def's
		for exon_def in exon_list:
			block_count = len(exon_def) # number of blocks for each mRNA
			thick_start = cds_def[0]['start']
			thick_stop = cds_def[cds_exon_count-1]['stop']
			bed6_str = chrom+'\t'+ exon_def[0]['start']+'\t'+\
				exon_def[block_count-1]['stop']+'\t'+gene_def_dict['gene_name']+'\t'+str(score)+'\t'+\
				[0]['strand']
			# create block size and block start lists
			block_sizes=[] # initialize array of integer sizes
			block_starts=[] # initialize array of integer starts
			for exon in mrna_def:
				block_size= str(int(exon_def['stop'])-int(exon_def['start']))
				block_sizes.append(block_size)
				block_start= str(int(exon_def['start'])-int(mrna_def[0]['start']))                
				block_starts.append(block_start)

			block_sizes_str=','.join(block_sizes)
			block_starts_str=','.join(block_starts)
			# format thickStart and thickEnd columns 7 and 8 and blocks
			mrna_ouput='\t'.join([bed6_str, thick_start, thick_stop, rgb, str(block_count), block_sizes_str, block_starts_str])
			mrna_ouputs.append(mrna_ouput)
		return '\n'.join(mrna_ouputs) 
	'''
# Write the formatted information to a new bed file.
bed_file= open('bed_file.bed','w') 

for key in gene_dict:
    #print '-----------key of gene_dict-------------'
    #print key
    print format_bed12_line(gene_dict[key])
    bed_file.write(format_bed12_line(gene_dict[key])+'\n')
		
# Close the original opened text file
file_open.close;

# Close the newly written bed file
bed_file.close