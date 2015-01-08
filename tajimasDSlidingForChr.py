#!/usr/bin/env python

import sys

def usage():
    return '''Calculates Tajima's D statistic in sliding windows within a chromosome, and prints the results to stdout.

usage: python tajimasDSlidingForChr.py refFileName dafFileName maskFileName winSize stepSize

refFileName: the path to a fasta-formatted file containing the reference sequence of a single chromosome. The header line contains the name of the chromosome.
dafFileName: the path to a tab-delimited file containing the derived allele frequency at each site within a single chromosome (see README for format)
maskFileName: the path to a tab-delimited file with coordinates of regions masked by RepeatMasker (see README for format); these sites are omitted fr
om calculations
winSize: the length of windows in which Tajima's D is calculated
stepSize: the distance between starting locations of windows to be included in the output
'''

try:
    refFileName, dafFileName, maskFileName, winSize, stepSize = sys.argv[1:]
except Exception:
    sys.exit(usage())
winSize = int(winSize)
stepSize = int(stepSize)

def calcThetaW(n):
    thetaW = 0
    if len(n) > 1:
        for j in range(len(n)):
            if n[j] != "NA":
                denom = 0
                for i in range(1,n[j]):
                    denom += 1/float(i)
                thetaW += 1 / denom
    return thetaW

def a1f(n):
    a1 = 0.0
    for i in range(1,n):
        a1 += 1.0/i
    return a1

def a2f(n):
    a2 = 0.0
    for i in range(1,n):
        a2 += 1.0/(i*i)
    return a2

def b1f(n):
    b1 = (n + 1.0)/(3.0*(n-1.0))
    return b1

def b2f(n):
    b2 = (2*(n*n + n + 3.0))/(9*n*(n - 1))
    return b2

def e1f(a1, c1):
    e1 = c1/a1
    return e1

def e2f(a1, a2, c2):
    e2 = c2/((a1*a1)+a2)
    return e2

def c1f(a1, b1):
    c1 = b1 - (1/a1)
    return c1

def c2f(n, a1, a2, b2):
    c2 = b2 - ((n+2)/(a1*n)) + (a2/(a1 * a1))
    return c2

def calcTajimasD(segSites,thetaW,pi,n):
    if thetaW == 0:
        return 0.0
    else:
        nNum, nDenom = 0, 0
        narray = n
        for snpCov in n:
            if snpCov != "NA":
                nNum += snpCov
                nDenom += 1
        n = int(nNum/float(nDenom))
        if n > 3:
            a1 = a1f(n)
            a2 = a2f(n)
            b1 = b1f(n)
            b2 = b2f(n)
            c1 = c1f(a1, b1)
            c2 = c2f(n, a1, a2, b2)
            e1 = e1f(a1, c1)
            e2 = e2f(a1, a2, c2)
            #try:
            ans = ((pi - thetaW)/((e1*segSites) + ((e2*segSites)*(segSites-1)))**0.5)
            #except Exception:
            #    print a1, a2, b1, b2, c1, c2, e1, e2
            #    print n, narray
            #    print e1, e2, segSites, thetaW, pi
            #    sys.exit()
            return ans
        else:
            return "NA"

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

refSeq = refSeq.upper()
chrLen = len(refSeq)

#sys.stderr.write("reading %s\n" %(dafFileName))
dafFile = open(dafFileName)
first = 1
piNum = [0]*chrLen
piDenom = [0]*chrLen
sNum = [0]*chrLen
sampleSizes = [0]*chrLen
coverages = ["NA"]*chrLen#array of SNP coverages ("NA" if not a snp)
for line in dafFile.xreadlines():
    if not first:
        site,der,anc,daf,sampleSize = line.strip().split("\t")
        site = int(site)
        #if (site % 1000000) == 0:
            #sys.stderr.write("checked %s sites\r" %(site))
        sampleSize = int(sampleSize)
        if anc in "ACGTacgt" and sampleSize > 1:
            piDenom[site] = 1
            if der in "ACGTacgt":
                coverages[site] = sampleSize
                p = float(daf)/sampleSize
                piNum[site] = 2*p*(1-p)*(sampleSize/(sampleSize-1.0))
                sNum[site] = 1
    first = 0
dafFile.close()
#sys.stderr.write("\n")
sSum = 0
piDenomSum = 0
piSum = 0
currCoverages = coverages[:winSize]
for i in range(winSize):
    piDenomSum += piDenom[i]
    piSum += piNum[i]
    sSum += sNum[i]

if piDenomSum:
    badPi = False
else:
    badPi = True

if not badPi:
    tajD =  calcTajimasD(sSum, calcThetaW(currCoverages), piSum, currCoverages)
    if tajD != "NA":
        print "%s\t%s\t%s\t%s" %(c, 1, winSize, tajD)

def pickSandE(winS,winE,stepSize):
    midp = (winE+winS)/2
    s = midp - ((stepSize/2)-1)
    e = midp + (stepSize/2)
    return s,e

for winStart in range(1,chrLen-winSize+1):
    piDenomSum = piDenomSum - piDenom[winStart-1] + piDenom[winStart+winSize-1]
    piSum = piSum - piNum[winStart-1] + piNum[winStart+winSize-1]
    sSum = sSum - sNum[winStart-1] + sNum[winStart+winSize-1]
    ditchedCoverages = currCoverages[0]
    currCoverages = currCoverages[1:] + [coverages[winStart+winSize-1]]
    snpCount = 0
    for snpCov in currCoverages:
        if snpCov != "NA":
            snpCount += 1
    if snpCount != sSum:
        print sSum, currCoverages
        print sNum[winStart-1],sNum[winStart+winSize-1]
        print ditchedCoverages,currCoverages[-1]
        print coverages[winStart-1], coverages[winStart+winSize-1]
        print winStart
        raise Exception

    if piDenomSum:
        badPi = False
    else:
        badPi = True

    if not badPi and winStart % stepSize == 0:
        #print (s, e, winStart, winStart+winSize)
        #print currCoverages, sSum
        s,e = pickSandE(winStart,winStart+1+winSize,stepSize)
        tajD =  calcTajimasD(sSum, calcThetaW(currCoverages), piSum, currCoverages)
        if tajD != "NA":
            print "%s\t%s\t%s\t%s" %(c, s, e, tajD)
