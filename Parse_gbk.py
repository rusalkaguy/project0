#Defining Accession Number and Database
assession_number='U00096.gbk'
db='nucleotide'
#Parsing Genbank Files

from Bio import SeqIO #Use SeqIO.read if there is only one genome (or sequence) in the file, and SeqIO.parse if there are multiple sequences. Since we're using genbank files, there typically (I think) only be a single giant sequence of the genome.
genome=SeqIO.read(assession_number,'genbank') #you MUST tell SeqIO what format is being read

print genome.features[:10] #prints a short description of the first ten features

#Feature Types
feats=set()
for feat in genome.features:
	feats.add(feat.type)

feats

#Feature Qualifiers
#Typical information will be 'product' (for genes), 'gene' (name) , and 'note' for misc. crap. Her's the qualifier dictionary for the first coding sequence (feature.type=='CDS'):

feat=genome.features[4]

feat.qualifiers

#How to use information from above in practice: 'product' and 'function' provide the current knowledge of what the gene (is thought to) make and what it (is thought to) do. 
predicted_genes=[]
for i,feature in enumerate(genome.features):
    if feature.type=='CDS':
        if 'product' in feature.qualifiers: #verify it codes for a real protein (not pseudogene)
            tot_genes+=1
            #data in qualifiers are all lists, even if only 1 string, so [0] converts to string
            #use lower() to make sure weirdly capitilized strings get selected as well
            product=feature.qualifiers['product'][0].lower()

            if product=='predicted protein':
                num_predicted+=1
                predicted_genes.append(feature)

#Grabbing the coding sequence
seq=feature.extract(genome.seq)
protein=seq.translate() #can also use optional switches like cds=True to raise errors if not a proper coding sequence (lacks start and stop codons)

sequence=feature.extract(genome.seq) # simple if you just want the gene sequence; works on antisense strand genes
#things get trickier if you want information outside of gene
location=feature.location
#position gives the actual python slice location of the start and end. need to reverse if on negative strand
#could use to get upstream or downstream regions
start=location.start.position
end=location.end.position

#Grabbing genomes from Genbank
from Bio import Entrez
#Replace with your real email 
Entrez.email = 'bheater@uab.edu'
handle=Entrez.efetch(db='nucleotide',id=assession_number,rettype='gb') # Accession id works, returns genbank format, looks in the 'nucleotide' database
#store locally
local_file=open(assession_number,'w')
local_file.write(handle.read())
handle.close()
local_file.close()

#Use NCBI Eutils to Fetch Records: http://biopython.org/DIST/docs/api/Bio.Entrez-pysrc.html
#Fetches Entrez results, which are returned as a handle.
def efetch (db, **keywords): #EFetch retreives records in the requested formate from a list of UIs or frmo the environment of the user.
	"""
	Example:
		>>> from Bio import Entrez
		>>>Entrez.email = "Your.Name.Here@example.org" 
		>>>handle = Entrez.efetch(db="nucleotide", id="57240072", rettype="gb", retmode="text")
		>>>print(handle.readline().strip())
		LOCUS       AY851612                 892 bp    DNA     linear   PLN 10-APR-2007
		>>> handle.close() 
	**Warning:** NCBI changed the default retmode in Feb 2002, so databased now give XML instead of text output.
	"""
	cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
	variables = {'db': db} 
	variables.update(keywords)
	post = False 
	try:
		ids = variables["id"]
	except KeyError:
		pass #Raises IOError exception if a network error occurs.
	else:
		if isinstance(ids,list):
			ids = ",".join(ids)
			variables["id"] = ids
		if ids.count(",") >=200:
			# NCBI prefers an HTTP POST instead of an HTTP GET if there are >~200 IDs
			post = True
	return _open(cgi, variables, post) #Returns a handle to results.
	#http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc114
		