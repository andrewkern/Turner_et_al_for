//pgStats for a bedFile rather than windows

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "stringWrap.h"
#include "sequenceMatrix.h"
#include "pgSummaryStats.h"
#include "bedFile.h"


void usage();

int main(int argc, char *argv[]){
	struct sequenceMatrix *aSeqMat, *ancSeqMat;
	int i, bedElNumber, h, ss, j, start, stop, step;
	struct bedEl data[100000];
	double pi, theta_h, z, w, wins[50], max, winSize, taj, H;
	int nwins = 9;
	

	
	if(argc < 2){
		usage();
		exit(1);
	}

	//open fastaFile and bedFile
	aSeqMat = sequenceMatrix_importFasta(argv[1]);
	ancSeqMat = sequenceMatrix_importFasta(argv[2]);
	bedElNumber = bedFileImport3(argv[3], data);

	//normal header
	printf("chrom\tchromStart\tchromEnd\tpi\tsegSites\tthetaH\ttajD\tfayWuH\tHapCount\tZnS\tOmega");
	for(i=0;i<nwins;i++)printf("\tpiWin%d",i);
	printf("\n");

	//loop through beds; the adjustments to end are to honor the zero indexed half open bed convention
	for(i=0;i<bedElNumber;i++){
		printf("%s\t%ld\t%ld\t",data[i].chrom,data[i].chromStart,data[i].chromEnd);
		pi = nucdivFromTo(aSeqMat, data[i].chromStart, data[i].chromEnd+1);
		ss = segSiteCountFromTo(aSeqMat, data[i].chromStart, data[i].chromEnd+1);
		taj = tajd(aSeqMat->sampleSize,ss,pi);
		theta_h = thetaHFromTo(aSeqMat,ancSeqMat,data[i].chromStart, data[i].chromEnd+1);
		H = theta_h - pi;
		h = nHaplotypes(aSeqMat,data[i].chromStart, data[i].chromEnd+1);
		z = ZnSFromTo(aSeqMat,data[i].chromStart, data[i].chromEnd+1);
	//	w = omegaMaxFromTo(aSeqMat,data[i].chromStart,data[i].chromEnd+1);
		w = omegaAtCenter(aSeqMat, data[i].chromStart, data[i].chromEnd+1, data[i].chromStart + 5000.0);
		
		printf("%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f",pi, ss, theta_h, taj, H, h, z, w);

		//window stats
	
		winSize = floor((data[i].chromEnd+1 - data[i].chromStart) / 9);
		step = winSize;
		start = data[i].chromStart;
		stop = start + step;
		/*pi_itr = 0;*/

		for (j = 0; j < nwins; ++j){
			wins[j] = nucdivFromTo(aSeqMat, start, stop);
			start += step;
			stop += step;
		}
		//print normalized windows
		max=0.0;
		for( j=0; j<nwins ; j++){
			if (wins[j] > max){
				max = wins[j];
			}
		} 
		if (max > 0) {
			for( j=0; j<nwins ; j++) printf("\t%f",wins[j]/max) ;    
		}
		else{
			for( j=0; j<nwins ; j++) printf("\t%f",0.0) ;    
		}
		
		
		printf("\n");
	}
	sequenceMatrix_free(aSeqMat);
	return(0);
}	

void usage(){
        printf("Computes various population genetic summary statistics in regions of a chromosome specified in a .bed file.\n\n");
	printf("usage: pgStatsBed ingroupFastaFile ancestorFastaFile bedFile\n\n");
        printf("ingroupFastaFile: a fasta file with a population sample of chromosomes\n");
        printf("ancestorFastaFile: a fasta file with the chromosomes inferred ancestral state\n");
        printf("bedFile: a file with the coordinates in which summary stats are computed, in bed format\n");
}
