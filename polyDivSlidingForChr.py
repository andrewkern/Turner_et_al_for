#!/usr/bin/env python

import sys

def usage():
    return '''Calculates the ratio of pi:divergence in sliding windows within a chromosome, and prints the results to stdout.

usage: python polyDivSlidingForChr.py refFileName dafFileName maskFileName winSize stepSize

refFileName: the path to a fasta-formatted file containing the reference sequence of a single chromosome. The header line contains the name of the chromosome.
dafFileName: the path to a tab-delimited file containing the derived allele frequency at each site within a single chromosome (see README for format)
maskFileName: the path to a tab-delimited file with coordinates of regions masked by RepeatMasker (see README for format); these sites are omitted fr
om calculations
winSize: the length of windows in which the ratio of pi:divergence is calculated
stepSize: the distance between starting locations of windows to be included in the output
'''

try:
    refFileName, dafFileName, maskFileName, winSize, stepSize = sys.argv[1:]
except Exception:
    sys.exit(usage())

winSize = int(winSize)
stepSize = int(stepSize)

refFile = open(refFileName)
first = 1
for line in refFile.xreadlines():
    if line.startswith(">"):
        assert first
        refSeq = ""
        c = line.strip().split()[0][1:]
        first = 0
    else:
        refSeq += line.strip()
refFile.close()
refSeq = list(refSeq)

maskFile = open(maskFileName)
for line in maskFile.xreadlines():
    currC,s,e = line.strip().split()[:3]
    if currC == c:
        s = int(s)+1
        e = int(e)
        for pos in range(s,e+1):
            refSeq[pos-1]="N"
maskFile.close()
refSeq = "".join(refSeq)

chrLen = len(refSeq)
piDenom = [0]*chrLen
piNum = [0]*chrLen

#sys.stderr.write("reading %s\n" %(dafFileName))
dafFile = open(dafFileName)
first = True
ancSeq = ""
for line in dafFile.xreadlines():
    if not first:
        site,der,anc,daf,sampleSize = line.strip().split("\t")
        site = int(site)
        ancSeq += anc
        sampleSize = int(sampleSize)
        if anc in "ACGTacgt" and sampleSize > 1:
            piDenom[site] = 1
            if der in "ACGTacgt":
                p = float(daf)/sampleSize
                piNum[site] = 2*p*(1-p)*(sampleSize/(sampleSize-1.0))
    first = False
dafFile.close()

refSeq = refSeq.upper()
ancSeq = ancSeq.upper()
assert len(refSeq) == len(ancSeq)

divNum,divDenom = [],[]
for i in range(chrLen):
    if refSeq[i] != "N" and ancSeq[i] != "N":
        divDenom.append(1)
        if refSeq[i] == ancSeq[i]:
            divNum.append(0)
        else:
            divNum.append(1)
    else:
        divNum.append(0)
        divDenom.append(0)

divDenomSum = 0
divNumSum = 0
piDenomSum = 0
piNumSum = 0
for i in range(winSize):

    divDenomSum += divDenom[i]
    divNumSum += divNum[i]

    piDenomSum += piDenom[i]
    piNumSum += piNum[i]

badDiv,badPi = False,False
if divDenomSum:
    div = divNumSum/divDenomSum
else:
    badDiv = True
if piDenomSum:
    pi = piNumSum/piDenomSum
else:
    badPi = True
if not badPi and not badDiv and div != 0:
    print "%s\t%s\t%s\t%s" %(c,1,winSize,pi/div)

def pickSandE(winS,winE,stepSize):
    midp = (winE+winS)/2
    s = midp - ((stepSize/2)-1)
    e = midp + (stepSize/2)
    return s,e

for winStart in range(1,chrLen-winSize+1):
    divDenomSum = divDenomSum - divDenom[winStart-1] + divDenom[winStart+winSize-1]
    divNumSum = divNumSum - divNum[winStart-1] + divNum[winStart+winSize-1]
    piDenomSum = piDenomSum - piDenom[winStart-1] + piDenom[winStart+winSize-1]
    piNumSum = piNumSum - piNum[winStart-1] + piNum[winStart+winSize-1]

    badDiv,badPi = False,False
    if divDenomSum:
        div = divNumSum/float(divDenomSum)
    else:
        badDiv = True
    if piDenomSum:
        pi = piNumSum/piDenomSum
    else:
        badPi = True
    if not badPi and not badDiv and winStart % stepSize == 0 and div != 0:
        s,e = pickSandE(winStart,winStart+1+winSize,stepSize)
        print "%s\t%s\t%s\t%s" %(c,s,e,pi/div)
