# Regular expressions find and parse particular pieces of information for bed formatting

# Variables for column location
feature_type_col = 0
gene_name_col = 1
location_col = 2

# Import regular expressions
import re

# for debugging output
import pprint
pp = pprint.PrettyPrinter(indent=0)

# Open the file of interest
file_open = open ("example_genbank_redux.txt", "r+")
file_read = file_open.readlines()

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

gene_dict={}
for line in file_read: # for feature in record=genome, instead of line
    line = line.split('#')[0]
    tab = line.split()  # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    if len(tab) == 0:
        continue
    # Assign variables for items of interest
    gene_name = tab[gene_name_col]
    feature_type = tab[feature_type_col]
    location = tab[location_col]

    # Check to see if gene exits
    if gene_name not in gene_dict: # if gene not in dict,
        gene_dict[gene_name]={'gene_name': gene_name} # set dictionary and store gene_name to process later
        #gene_dict[gene_name]['gene_name']=gene_name
    # Set gene_dict to key-value pair {gene_name:gene_def_dict} 
    gene_def_dict = gene_dict[gene_name]

    # Regular expressions identify the complement pattern and join pattern
    strand = '+'
    comp_hit = p.match(location) # finds complement at begining of row
    # did you find complement?
    if comp_hit:
        strand = '-'
        #print 'Match found. Span=', ' Group(0)=', comp_hit.group(), 'Group(1)=', comp_hit.group(1)
        spans= comp_hit.group(1)
    else:   
        spans = location
    #print "spans="+ spans
    '''
    # if comma--> join use findall
        # Identify the join pattern for exons
    # Use findall or split on comma operator ",".split(location)-->return array with 1+ 
    # Use findall because it ignores the joins
    fa = q.findall( spans )
    if fa !=[]: # if fa is NOT equal to an empty array
        print " findall=", fa
        for item in fa:
            start_stop = item.split("..")
            print start_stop
    '''
    # use group function --> more efficient because only parse loop once
    loc_array = [] # make empty array to store strand, start, and stop info
    # neccessary when more than one exon or CDS region exists for a gene
    hits = q2.finditer(spans) # returns objects for matches per group
    for match in hits:
        #print 'finditer found. Span=', ' Group=', match.group()
        loc_dict={ 'strand':strand, 'start':match.group(1), 'stop':match.group(2)}
        loc_array.append(loc_dict)

    # Add key-value pair {feature_type:loc_array} to gene_def dictionary
    # Check to see if 'mRNA' is already a key in the dictionary
    if feature_type == 'mRNA' and 'mRNA' in gene_def_dict:
        # If value associated with'mRNA' is already in dict and the value is a dictionary
        if len(gene_def_dict['mRNA']) > 0 and type(gene_def_dict['mRNA'][0]) is dict:
            # Add loc_array as the second element to the array containing the value for dict 
            gene_def_dict['mRNA'] = [gene_def_dict['mRNA'], loc_array]
        # If len(gene_def_dict['mRNA']) is an empty array
        else: 
            gene_def_dict['mRNA'].append(loc_array)
    # If value associated with'mRNA' is not already in dict and the value is a dictionary
    else:
        # Append loc_array to empty loc_array
        gene_def_dict[feature_type]=loc_array
    # Add key-value pair {'mRNA':[]} to gene_def_dict to make format conistent and aid formatting later.
    if 'mRNA' not in gene_def_dict:
        gene_def_dict['mRNA']=[]

#print "#--------- gene_dict { gene_def_dict { exon_dict } } ------------"
#pp.pprint(gene_dict)

'''
print "#--------- formatted gene_dict with map and function ------------"
def format_bed(x): # x is input argument
    return str(x['start']),'\t',str(x['stop']),'\t','\t',str(score), str(x['strand'])
for key in gene_dict:
    gene_list = gene_dict[key]
    print chrom+'\t'.join(map(format_bed,gene_list))
'''
# create the first 6 columns (easy) and columns 7/8 (more interesting)
def format_bed6_line(gene_def_dict): 
    # first 6 columns: 
    # how do i fix 'chrom','name', and 'score' without using variables in the function?
    # either make def more specific or call info outside of def and for loop, to later incoorporate
    
    return 'chrom'+'\t'+ gene_def_dict['gene'][0]['start']+'\t'+\
        gene_def_dict['gene'][0]['stop']+'\t'+gene_def_dict['gene_name']+'\t'+'score'+'\t'+\
        gene_def_dict['gene'][0]['strand']
    # how do i differentiate between CDS and genes or mRNA without a for loop or if statement?
    # nested calls to get values of dictionaries by dict[key]
    # str.join() is a similar function to split that will separate items by tab if i indicate column and delimeter
chrom='NC_006273v2'
rgb='0,0,0'
def format_bed12_line(gene_def_dict):

    #pp.pprint(gene_def_dict)
    # column 9:itemRgb
    score=1000
    item_rgb='0,0,0'
    # columns 7-8 thickStart and thickStop
    cds_def=gene_def_dict['CDS']
    cds_exon_count= len(cds_def)
    mrna_count=len(gene_def_dict['mRNA'])
    #print 'CDS_count = '+str(cds_count)
    if cds_exon_count>0 :
        # compute blockCount    blockSizes  blockStarts
        mrna_ouput=[]
        mrna_ouputs=[]

        virtual_mrna_list = [cds_def] # list of CDS, where CDS is a list of exon_def's

        #print('----------cds_def------------------')
        #pp.pprint(cds_def)
        # when no mRNA is present (mrna_count=0), cds defines block size and start
        if mrna_count>0:
            virtual_mrna_list = gene_def_dict['mRNA']
        for mrna_def in virtual_mrna_list:
            #print('----------mRNA_def-----------------')
            #pp.pprint(mrna_def)
            block_count = len(mrna_def) # number of blocks for each mRNA
            # only apply the CDS to the mRNA's that completely contain it.
            # Does the CDS start and stop fall within the mRNA start and stop?
            # If not, then set thick start and stop to be mRNA start and stop.
            if cds_def[0]['start']<mrna_def[0]['start'] or cds_def[0]['stop']>mrna_def[block_count-1]['stop']:
                thick_start=mrna_def[0]['start']
                thick_stop=mrna_def[block_count-1]['stop']
            else:          
                thick_start = cds_def[0]['start']
                thick_stop = cds_def[cds_exon_count-1]['stop']

            bed6_str = chrom+'\t'+ mrna_def[0]['start']+'\t'+\
                mrna_def[block_count-1]['stop']+'\t'+gene_def_dict['gene_name']+'\t'+str(score)+'\t'+\
                mrna_def[0]['strand']

            # If there are 2 mRNAs, we need to output 2 lines in the bed file
            # for m in range(0,mrna_count):
            #print "mrna_def"
            #pp.pprint(mrna_def)
            # create block size and block start lists
            block_sizes=[] # initialize array of integer sizes
            block_starts=[] # initialize array of integer starts
            for exon_def in mrna_def:
                #print "exon_def"
                #pp.pprint(exon_def)
                #print mrna_def[0]['start']
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


# Write the formatted information to a new bed file.
bed_file= open('bed_file.bed','w') 
# 'a' creates new file if it does not exist and does not truncate the file if it does exist
# 'w' creates the file if the file does not exist, but it will truncate the existing file


for key in gene_dict:
    #print '-----------key of gene_dict-------------'
    #print key
    print format_bed12_line(gene_dict[key])
    bed_file.write(format_bed12_line(gene_dict[key])+'\n')

# Close the original opened text file
file_open.close;
# Close the newly written bed file
bed_file.close