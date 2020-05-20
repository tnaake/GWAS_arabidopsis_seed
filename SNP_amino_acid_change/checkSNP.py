## check AA change...
## Note: use TAIR 9 XML files

from collections import namedtuple
from parse_xml import TAIRparser
import string
import gzip

########################################################################
# input data

HMM_FILES   = "349_t9_mafSNP_C%d.hmp.txt.gz"
TAIR_FILES  = "ch%d.xml.gz"  
INPUT_FILE  = "GeneList.txt"
OUTPUT_FILE = "Result.txt"

########################################################################

def getSNP():
	snp = namedtuple("snp", ['rs', 'allele','chrom', 'pos', 'strand', 'col0']) 
	hmm  = HMM_FILES
	ret  = []

	def parseHMM(i):
		head = True
		with gzip.open(hmm % i) as fp:
			for line in fp:
				x = line.strip().split("\t")
				if head:
					col0 = x.index("ecotype.6909")
					head = False
					continue
				ret.append( snp(x[0], x[1], x[2], int(x[3]), x[4], x[col0]) )

	for i in [1,2,3,4,5]:
		parseHMM(i)
	return ret

def translate(s):
	codons = {
		'AAA': 'K',  'AAC': 'N',  'AAG': 'K',  'AAT': 'N',
		'ACA': 'T',  'ACC': 'T',  'ACG': 'T',  'ACT': 'T',
		'AGA': 'R',  'AGC': 'S',  'AGG': 'R',  'AGT': 'S',
		'ATA': 'I',  'ATC': 'I',  'ATG': 'M',  'ATT': 'I',

		'CAA': 'Q',  'CAC': 'H',  'CAG': 'Q',  'CAT': 'H',
		'CCA': 'P',  'CCC': 'P',  'CCG': 'P',  'CCT': 'P',
		'CGA': 'R',  'CGC': 'R',  'CGG': 'R',  'CGT': 'R',
		'CTA': 'L',  'CTC': 'L',  'CTG': 'L',  'CTT': 'L',

		'GAA': 'E',  'GAC': 'D',  'GAG': 'E',  'GAT': 'D',
		'GCA': 'A',  'GCC': 'A',  'GCG': 'A',  'GCT': 'A',
		'GGA': 'G',  'GGC': 'G',  'GGG': 'G',  'GGT': 'G',
		'GTA': 'V',  'GTC': 'V',  'GTG': 'V',  'GTT': 'V',

		'TAA': '*',  'TAC': 'Y',  'TAG': '*',  'TAT': 'Y',
		'TCA': 'S',  'TCC': 'S',  'TCG': 'S',  'TCT': 'S',
		'TGA': '*',  'TGC': 'C',  'TGG': 'W',  'TGT': 'C',
		'TTA': 'L',  'TTC': 'F',  'TTG': 'L',  'TTT': 'F'
	}
	i = 0
	p = []
	while i < len(s):
		p.append( codons[ s[i:(i+3)] ]  )
		i += 3

	return "".join(p)

def revseq(s):
	s = s.translate(string.maketrans("ATGC","TACG"))
	return s[::-1]

def buildCDS(x, model, inputseq):
	if len(x.seq) != len(inputseq):
		return ""

	seq = []
	for e in model.exons:
		if e.cds is not None:
			if x.strand == 1:
				seq.append( inputseq[ (e.cds.start - x.start):(e.cds.end - x.start + 1): ] )
			else:
				seq.append( inputseq[ (e.cds.end - x.end):(e.cds.start - x.end + 1): ] )
		
	seq = "".join(seq)
	if x.strand == 1:
		return seq
	return revseq(seq)

def replaceStr(s, t, i):
	return s[:i] + t + s[(i+1):]

#####################################################################
#####################################################################

SNP = getSNP()

TU = []
for i in [1,2,3,4,5]:
	TU.extend( TAIRparser(TAIR_FILES % i ) )

loci = map(lambda x: x.locus, TU)

Gene = namedtuple('Gene', ['code','name'])
genes = []

with open(INPUT_FILE) as fp:
	for line in fp:
		if line[:2] == "AT":
			x = line.strip().split("\t")
			genes.append( Gene(x[0], x[1] ) )

fp = open(OUTPUT_FILE, "wt")
fp.write("AGI_code Gene_Name Gene_Model SNP Allele Col-0 Change From To\n".replace(" ", "\t"))

for gene in genes:
	x = TU[ loci.index( gene.code ) ]

	## find SNP
	chrom = gene.code[2]
	sel = []
	stt, end = (x.start, x.end) if x.strand == 1 else (x.end, x.start)

	for z in SNP:
		if z.chrom == chrom and z.pos >= stt and z.pos <= end:
			sel.append(z)

	if len(sel) == 0:
		continue

	## replace nucleotides 1 per 1...
	for m in x.model:
		for z in sel:
			if x.seq[ z.pos - stt ] != z.col0:
				print "error"
				print z
				next
			nt = filter(lambda w: w != z.col0, z.allele.split("/"))[0]
			S = replaceStr(x.seq, nt, z.pos - stt)
			cds = buildCDS(x, m, S)
			prot = translate(cds)
			if prot == m.prot:
				change, frm, to = False, '.', '.'
			else:
				change = True
				tmp = filter(lambda x: x[0] != x[1], zip(m.prot, prot))
				if len(tmp) != 1:
					print "Error =>"
				frm, to = tmp[0]
			fp.write("\t".join(( gene.code, gene.name, m.locus, z.rs, z.allele, z.col0, str(change), frm, to)) + "\n")

fp.close()

