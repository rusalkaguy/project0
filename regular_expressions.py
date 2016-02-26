
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
for line in file_read:
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
        gene_dict[gene_name]={}    # create empty dict

    # Set gene_dict to key-value pair {gene_name:gene_def} 
    gene_def = gene_dict[gene_name]

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
    loc_array = [] # make empty array
    hits = q2.finditer(spans) # returns objects for matches per group
    for match in hits:
        #print 'finditer found. Span=', ' Group=', match.group()
        loc_dict={  'strand':strand, 'start':match.group(1), 'stop':match.group(2)}
        #loc_dict={}
        #loc_dict['start'] = match.group(1)
        #stop = match.group(2)
        #print loc_dict
        loc_array.append(loc_dict)
    #print loc_array

    # Add key-value pair {feature_type:loc_array} to gene_def dictionary
    gene_def[feature_type]=loc_array

print "#--------- gene_dict ------------"
pp.pprint(gene_dict)

''' 
Variables writing in bed format
chrom = NC_006273 # The name of the chromosome
chromStart =  # The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chomEnd =  # The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. 
# For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
name = # Defines the name of the BED line. 
score = # A score between 0 and 1000.
strand = # Defines the strand - either '+' or '-'complement.
thickStart = # The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). 
# When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
thickEnd = # The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
itemRgb = # An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. 
blockCount = # The number of blocks (exons) in the BED line.
blockSizes = # A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
blockStart = # A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
'''

# Variable for bed format
chrom = 'NC_006273'
score = 0


print "#--------- formatted gene_dict with map and function ------------"
def format_bed(x): # x is input argument
    return str(chrom),'\t',str(x['start']),'-',str(x['stop']),str(gene_dict[key]),str(score)
for key in gene_dict:
    gene_list = gene_dict[key]
    print key+'\t'+",".join(map(format_bed,gene_list))
    # error--> sol'n: http://stackoverflow.com/questions/18931315/typeerror-string-indices-must-be-integers-not-str-working-with-dict
    # as you go through each dictionary, add another nest
