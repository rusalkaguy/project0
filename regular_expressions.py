
feature_type_col = 0
gene_name_col = 1
location_col = 2

# Import regular expressions
import re

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
    join\(( # Start of a location
    .*)      # 0 or more(*) of anything(.)
    \)       # Trailing end parenthesis
    """, re.VERBOSE)

for line in file_read:
    line = line.split('#')[0]
    tab = line.split()  # 'tab' could be any any other variable
    # Split at tabs is default: i.e. empty brackets.
    # Insert space or comma in parenthesis () if that denotes separation.
    if len(tab) == 0:
        continue
    # print tab
    gene_name = tab[gene_name_col]
    feature_type = tab[feature_type_col]
    location = tab[location_col]
    # Regular expressions identify the complement pattern and join pattern
    hits = p.finditer(location)
    for match in hits:
        # Print start and stop location of complement w/o including 'complement'
        print 'finditer found.',' complement=', match.group(1) 
        # debug start_stop= re.split(r""".*\.\..*""",match)
        print start_stop
    # Identify the join pattern for exons
    # Use findall or split on comma operator ",".split(location)-->return array with 1+ 
    # Use findall because it ignores the joins
    fa = q.findall( location )
    if fa !=[]: # if fa is NOT equal to an empty array
        strand = '+'
        print " findall=", fa
        for item in fa:
            "..".split(item)
            print fa
    