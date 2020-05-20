#!/usr/bin/env python2

## parses the XML files from TAIR 9 or 10
## ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/Tair10_XML/
## ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/Tair9_XML/

from lxml import etree
import sys

class _tu:
	def __init__(self, locus="", seq="", start=0, end=0, model=[], gtype=''):
		self.locus = locus
		self.seq   = seq
		self.start = start
		self.end   = end
		self.model = model
		self.gtype = gtype

class _model:
	def __init__(self, locus="", exons=[], seq="", prot = ""):
		self.locus  = locus
		self.exons  = exons
		self.seq    = seq
		self.prot   = prot

class _exon:
	def __init__(self, name="", start=0, end=0, cds=None):
		self.name  = name
		self.start = start
		self.end   = end
		self.cds   = cds

class _cds:
	def __init__(self, start=0, end=0):
		self.start = start
		self.end   = end

def printGene(t):
	print "%s: start:%d end:%d strand:%d type:%s" % \
		(t.locus, t.start, t.end, t.strand, t.gtype)
	for m in t.model:
		print "  %s:" % m.locus
		for e in m.exons:
			print "    %s: start:%s end:%s" % (e.name, e.start, e.end)
			if e.cds is not None:
				print "      start:%s end:%s" % (e.cds.start, e.cds.end)

def getCoordset(x):
	start, end = 0, 0
	for y in x:
		if y.tag == "END5":
			start = int(y.text)
		elif y.tag == "END3":
			end = int(y.text)
	return (start, end)

def getLocusName(x):
	for i in x:
		if i.tag == "PUB_LOCUS":
			return i.text
	return 'NA'

def parseCDS(x):
	start,end = 0,0
	for i in x:
		if i.tag=="COORDSET":
			start, end = getCoordset(i)
	return _cds(start, end)

def parseModel(x):
	m = _model(exons=[])
	for i in x:
		if i.tag == "PUB_LOCUS":
			m.locus = i.text
		elif i.tag == "EXON":
			ex = _exon()
			for j in i:
				if j.tag == "FEAT_NAME":
					ex.name = j.text
				if j.tag == "COORDSET":
					ex.start, ex.end = getCoordset(j)
				if j.tag == "CDS":
					ex.cds = parseCDS(j)
			m.exons.append(ex)
		elif i.tag == "CDS_SEQUENCE":
			m.seq = i.text
		elif i.tag == "PROTEIN_SEQUENCE":
			m.prot = i.text

	return m

def parseGene(z, gtype):
	t = _tu(model=[], gtype=gtype)
	for x in z:
		## print x.tag
		if x.tag == "GENE_INFO":
			t.locus = getLocusName(x)
		elif x.tag == "COORDSET":
			t.start, t.end = getCoordset(x)
		elif x.tag == "MODEL":
			m = parseModel(x)
			t.model.append(m)
		elif x.tag == "TRANSCRIPT_SEQUENCE":
			t.seq = x.text
	return t

def TAIRparser(infile):
	parser = etree.XMLParser(huge_tree = True)
	tree   = etree.parse(infile, parser)
	root   = tree.getroot()
	pseudo = root[0]

	for asm in pseudo:
		if asm.tag == "ASSEMBLY":
			break

	for x in asm:
		if x.tag == "GENE_LIST":
			glist = x
			break

	for x in glist:
		if x.tag == "PROTEIN_CODING":
			prot = x
		if x.tag == "RNA_GENES":
			rna  = x
		if x.tag == "TRANSPOSABLE_ELEMENT_GENES":
			trna = x

	genes = []
	for tu in prot:
		t = parseGene(tu, "protein coding")
		genes.append(t)

	for r in rna:
		t = parseGene(r, "RNA gene")
		genes.append(t)

	for tr in trna:
		t = parseGene(tr, "transposable element")
		genes.append(t)

	for t in genes:
		if t.start > t.end:
			t.strand = -1
		else:
			t.strand = 1

	## sort the list by start position
	genes.sort(key=lambda t: t.start)
	return genes

## a simple binary search to search
def checkPos(t, pos):
	start, end = t.start, t.end
	if t.strand == -1:
		start, end = t.end, t.start
	if pos >= start and pos <= end:
		return 0
	elif pos > end:
		return 1
	else:
		return -1

def searchLocus(loci, pos):
	n = len(loci)
	i = int(n / 2)
	imin, imax = 0, n-1
	while imax - imin != 1:
		t = loci[i]
		d = checkPos(t, pos)
		if d == 0:
			return t
		if d == 1:
			imin = i
		else:
			imax = i
		i = imin + int((imax - imin) / 2)
	return None

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "Usage: %s <input.xml>" % (sys.argv[0])
		sys.exit()
	tus = TAIRparser(sys.argv[1])

