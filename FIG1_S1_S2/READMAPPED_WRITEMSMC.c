// READSLIMOUT2.c

#include "libhdr"

#define NN 150000 // Maximum 150000 SNPs segregating
#define CC 4001   // Maximum N=2000 individuals

int i, j, s, sb, m, x, l, nind, ind, NCHR;
int pos[NN], posb[NN], homs[8];

double w;
int numSNP, numSNPb;
char ch, crom[CC][NN], cromb[8][NN];

FILE *fmap, *fped, *fmsmc;

main()
{
	//getintandskip("NIND :",&nind,2,2000);
	//getintandskip("Genome length :",&genomelen,1,infinity);
	getintandskip("NCHR :",&NCHR,1,infinity);
	//getrealandskip("cMMb :",&cMMb,0.0,(double)infinity);

	readfiles();
	MSMC_file();

	return(0);
}

/* **************************************************************************** */

readfiles()
{
	// MSMC FILE
	fmsmc = fopen ("datamsmc","w");

	// ********** read .map **********

	fmap = fopen ("dataBP.map","r");

	while (!feof(fmap))
	{
		s ++;
		fscanf(fmap,"%d", &x);		//chr
		fscanf(fmap,"%lf", &w);		//rec rate cM
		fscanf(fmap,"%d", &x);		//pos bp
		pos[s] = x;
		if (pos[s] == pos[s-1])
		{
			pos[s] = pos[s-1] + 1;
		}
	}
	numSNP = s - 1;

	fclose(fmap);

	// ********** read .ped **********

	fped = fopen ("dataBP.ped","r");
	
	i=1;
	ind = 1;
	while (!feof(fped))
	{
		fscanf(fped,"%d", &x);		//sex
		fscanf(fped,"%d", &x);		//-9
		for(s=1; s<=numSNP; s++)
		{
			fscanf(fped,"%c", &ch); //space
			fscanf(fped,"%c", &ch); //1st allele
			crom[i][s] = ch;	
			fscanf(fped,"%c", &ch);	//space
			fscanf(fped,"%c", &ch);	//2nd allele
			crom[i+1][s] = ch;
		}
		fscanf(fped,"%c", &ch); 	//\n
		i+=2;
		ind++;

	}

	fclose(fped);

	return(0);
}

/* **************************************************************************** */


MSMC_file()
{
	sb = 0;
	for (s=1; s<=numSNP; s++)
	{
		if ((crom[1][s] == crom[2][s]) && (crom[2][s] == crom[3][s]) && (crom[3][s] == crom[4][s]) && (crom[4][s] == crom[5][s]) && (crom[5][s] == crom[6][s]) && (crom[6][s] == crom[7][s]) && (crom[7][s] == crom[8][s]))
		{
		} else 
		{
			sb++;
			posb[sb] = pos[s];
			for (i=1; i<=8; i++)
			{
				cromb[i][sb] = crom[i][s];
			}
		}
	}
	numSNPb = sb;

	for (sb=1; sb<=numSNPb; sb++)
	{
		fprintf(fmsmc,"%d\t%d\t%d\t", NCHR, posb[sb], posb[sb]-posb[sb-1]);
		for (i=1; i<=8; i++)
		{
			fprintf(fmsmc, "%c", cromb[i][sb]);
		}
		fprintf(fmsmc,"\n");
	}
	
	fclose(fmsmc);
	return(0);
}

