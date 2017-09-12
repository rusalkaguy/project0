table Virus
"Chromosomes and Genes of Viral Strain SR10-01 for GenomeAccession"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  name;		"Name of gene"
uint    score;		"Score"
char[1] strand;		"+ or - for strand"
uint	thickStart;	"Starting position of thick feature"
uint	thickEnd;	"Ending position of thick feature"
uint	reserved;	"RGB value determines the display color"
int	blockCount;	"Number of blocks or exons"
int[blockCount]	blockSizes;	"List of block sizes"
int[blockCount]	chromStarts;	"List of block starting positions"
string  genomeAccession;	"Genome accession"
)