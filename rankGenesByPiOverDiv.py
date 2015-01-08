#!/usr/bin/env python

import sys,os

def usage():
    return '''Calculates pi/divergence for the longest transcript of each gene in the user-specified gff file, ranks genes based on these values, and prints the results to stdout.

usage: python rankGenesByPiOverDiv.py afDir ancDir refDir gffFileName maskFileName

afDir: the path to a directory of tab-delimited files containing the derived allele frequency at each site within a single chromosome (see README for format); files must be named <chrname>.daf
ancDir: the path to a directory of fasta-formatted files each containing the inferred ancestral sequence of a single chromosome. Each file must be naemd <chrname>.fa
refDir: the path to a directory of fasta-formatted files each containing the reference sequence of a single chromosome. The header line of each file contains the name of the chromosome.
gffFileName: a gffFile containing transcript annotations. Following FlyBase, transcript names are <geneName>-R<transcriptNumber>.
maskFileName: the path to a tab-delimited file with coordinates of regions masked by RepeatMasker (see README for format); these sites are omitted from calculations
'''

try:
    afDir, ancDir, refDir, gffFileName, maskFileName = sys.argv[1:]
except Exception:
    sys.exit(usage())

def processRefChr(c,seq,fixationH,ancH,maskedH):
    maskedH[c] = []
    for i in range(len(seq)):
        maskedH[c].append(0)
        if ancH[c][i] in "ACGT" and seq[i] in "ACGT" and seq[i] != ancH[c][i]:
            pos = i+1
            fixationH[(c,pos)] = [ancH[c][i],seq[i]]

ancH = {}
for ancFaFileName in os.listdir(ancDir):
    if ancFaFileName.split(".")[-1] in ["fa","fasta"]:
        c = ancFaFileName.split(".")[0]
        sys.stderr.write("reading %s\n" %(ancDir+ancFaFileName))
        ancFile = open(ancDir+ancFaFileName)
        first = 1
        for line in ancFile.xreadlines():
            if line.startswith(">"):
                assert first
                seq = ""
                first = 0
            else:
                seq += line.strip().upper()
        ancFile.close()
        ancH[c] = seq

fixationH = {}
maskedH = {}
for refFaFileName in os.listdir(refDir):
    refFaFile = open(refDir +"/"+ refFaFileName)
    for line in refFaFile.xreadlines():
        if line.startswith(">"):
            seq = ""
            c = line.strip().split()[0][1:]
        else:
            seq += line.strip()
    refFaFile.close()
    processRefChr(c,seq,fixationH,ancH,maskedH)

gffFile = open(gffFileName)
gene2transcripts = {}
for line in gffFile.xreadlines():
    #A01     blat    CDS     4354155 4354677 100     -       .       Parent=Bra011027
    #C1      JCVI    CDS     977890  978273  .       +       0       ID=Bo1g004770.1_cds_2;Parent=Bo1g004770.1
    #C1      JCVI    mRNA    998080  998833  .       +       .       ID=Bo1g004820.1;Parent=Bo1g004820;Note=Iron-sulfur cluster insertion protein erpA
    if line.startswith(">"):
        break
    elif not line.startswith("#"):
        line = line.strip().split("\t")
        c,source,annotType,s,e,score,strand,frame,info = line
        if not c.startswith("chr"):
            c = "chr" + c
        s,e = int(s),int(e)
        if strand in "+-" and annotType == 'mRNA' and maskedH.has_key(c):
            transcript = info.split("ID=")[1].split(";")[0]
            geneName = "-R".join(info.split("Name=")[1].split(";")[0].split("-R")[:-1])
#            sys.stderr.write("geneName: %s\n" %(geneName))
            if not gene2transcripts.has_key(geneName):
                gene2transcripts[geneName] = []
            gene2transcripts[geneName].append((transcript,c,s,e))
gffFile.close()
#sys.stderr.write("gene2transcripts len: %s\n" %(len(gene2transcripts.keys())))

snpH = {}
for filename in os.listdir(afDir):
    c = filename.split(".")[0]
    if not "nonsyn" in filename:
        sys.stderr.write("reading %s\n" %(filename))
        afFile = open(afDir+filename)
        first = 1
        for line in afFile.xreadlines():
            if not first:
                site,der,anc,daf,sampleSize = line.strip().split("\t")
                site = int(site)+1
                daf = int(daf)
                if (site % 1000000) == 0:
                    sys.stderr.write("checked %s sites\r" %(site))
                sampleSize = int(sampleSize)
                if anc in "ACGTacgt" and sampleSize > 1 and der in "ACGTacgt" and daf < sampleSize:
                    p = daf/float(sampleSize)
                    snpH[(c,site)] = (anc,der,2*p*(1-p)*(sampleSize/(sampleSize-1.0)))
            first = 0
        afFile.close()
        sys.stderr.write("\n")

maskFile = open(maskFileName)
for line in maskFile.xreadlines():
    c,s,e = line.strip().split()[:3]
    if maskedH.has_key(c):
        s = int(s)+1
        e = int(e)
        for pos in range(s,e+1):
            maskedH[c][pos-1]=1
maskFile.close()
sys.stderr.write("maskedH keys: %s\n" %(maskedH.keys()))

genes = gene2transcripts.keys()
genes.sort()
outLs = []
outH = {}
for gene in genes:
    transcripts = gene2transcripts[gene]
    maxTrans,maxTransLen = False,0
    for transcript,c,s,e in transcripts:
        translen = e-s+1
        if translen > maxTransLen:
            maxTransLen = translen
            maxTrans = (transcript,c,s,e)
    assert maxTrans

    transcript,c,s,e = maxTrans
    if maskedH.has_key(c):
        pi = 0
        div = 0
        for pos in range(s,e+1):
            if not maskedH[c][pos-1]:
                if snpH.has_key((c,pos)):
                    pi += snpH[(c,pos)][2]
                if fixationH.has_key((c,pos)):
                    div += 1
        if not 0 in [pi,div]:
            polyDiv = pi/float(div)
            outH[gene] = 1
            outLs.append((gene,pi,div,polyDiv))

print "gene\tpoly\tdiv\tpolyDiv\trank"
denom = float(len(outLs))
cumD = 1/denom
outLs.sort(lambda x,y: cmp(x[-1],y[-1]))
currOutLs = ["%s\t%s\t%s\t%s" %outLs[0]]
prevRatio = polyDiv
outH = {}
for gene,pi,div,polyDiv in outLs[1:]:
    outH[gene] = 1
    if polyDiv == prevRatio:
        currOutLs.append("%s\t%s\t%s\t%s" %(gene,pi,div,polyDiv))
    else:
        for outline in currOutLs:
            print outline + "\t" + str(cumD)
        currOutLs = ["%s\t%s\t%s\t%s" %(gene,pi,div,polyDiv)]
    cumD += 1/denom
    prevRatio = polyDiv
for outline in currOutLs:
    print outline + "\t" + str(cumD)
assert abs(cumD)-1.0 < 1e-5
