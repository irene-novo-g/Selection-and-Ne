// SUMM_REP_BOpin.c

#include "libhdr"
#define BB 1001

int x, a, c, b, nchrom, NGEN, NBIN, PHASE, SAM;
unsigned long long int sumi[BB];
double w, sumcc[BB], sumd2[BB], sumD2[BB], sumVV[BB];
double cc[BB], d2[BB], D2[BB], VV[BB], F;

FILE *fin, *fout, *fnsnp;

main()
{
	getintandskip("NGEN :",&NGEN, 1, 5000);
	getintandskip("NBIN :",&NBIN, 1, 1000);
	getintandskip("n.chrom :",&nchrom, 1, 1000);

	/* ***************** readfile ******************** */

	fout = fopen ("outfileLD","w");
	fin = fopen ("CHROM","r");

	for (c=1; c<=nchrom; c++)
	{
		fscanf(fin,"%d", &x);
		PHASE = x;
		fscanf(fin,"%d", &x);
		SAM = x;
		fscanf(fin,"%lf", &w);
		if (w != -99)	F += w;

		for (b=1; b<=NBIN; b++)
		{
			fscanf(fin,"%d", &a);
			sumi[b] += a;
			fscanf(fin,"%lf", &w);
			sumcc[b] += w*a;
			fscanf(fin,"%lf", &w);
			sumd2[b] += w*a;
			fscanf(fin,"%d", &x);
			fscanf(fin,"%lf", &w);
			sumD2[b] += w*a;
			fscanf(fin,"%lf", &w);
			sumVV[b] += w*a;
		}
	}
	for (b=1; b<=NBIN; b++)
	{
		cc[b] = sumcc[b] / sumi[b];
		d2[b] = sumd2[b] / sumi[b];
		D2[b] = sumD2[b] / sumi[b];
		VV[b] = sumVV[b] / sumi[b];
	}

	fclose(fin);

	/* ***** outfileLD ********************************************************* */

	fprintf(fout,"%d\n", PHASE);
	fprintf(fout,"%d\n", SAM);
	fprintf(fout,"%f\n", F/nchrom);
	for (b=1; b<=NBIN; b++)
	if (sumi[b] != 0)
	{
			if (b <= 5)	fprintf(fout,"%lld %f %f %d\n", sumi[b], cc[b], d2[b], 2*b);
			else	fprintf(fout,"%lld %f %f %d\n", sumi[b], cc[b], d2[b], ((NGEN/NBIN)*(b-5))+10);
	}
}
