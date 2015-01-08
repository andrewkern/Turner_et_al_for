CC = gcc
CFLAGS = -O3 -Wall -lm -lgsl -lgslcblas -lc

BASICS = rhoCalc.c stringWrap.c sequenceMatrix.c pgSummaryStats.c vector.c snpFile.c numerical.c nrutil.c bedFile.c

all: pgStatsBed

pgStatsBed:	pgStatsBed.c $(BASICS)
	$(CC) pgStatsBed.c $(BASICS) -o pgStatsBed  $(CFLAGS) -g

clean:
	rm -f pgStatsBed

