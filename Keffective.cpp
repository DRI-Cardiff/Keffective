#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


int z(int i, int j, int r);
int c(int i, int j, int k);
int EM(double h[4], int g[9]);
double r(int g[9]);


int main(int argc, char *argv[])
{
FILE *fo, *fi;
int j, N, win, w;
char *infile, MarName[50], outKeff[50];
double Keff, alpha;
long Nmar, m;

fpos_t pos;
int ii, k, gg[9], *gen1, gen, a1, a2;
double d, rr, P_I, eps;
double 	*maxR;


if(argc>5)
{
	infile = argv[1];
	Nmar = atol(argv[2]);
	N = atoi(argv[3]);
	P_I = atof(argv[4]);
	win = atoi(argv[5]);
}
else
{
  fprintf(stderr,"Input parameters: FileName Nmarkers Nindividuals SigLevel WindowSize\n");
  exit(1);
}

if (win>Nmar) win=Nmar;

gen1=new int[N];
sprintf(outKeff,"%s.res", infile);

printf("Estimation of the number of effective tests...\n");
printf("The authors: V Moskvina and KM Schmidt.\n");
printf("N of markers:\t\t%i\n",Nmar);
printf("N of individuals:\t%i\n",N);
printf("Overall significance level:\t%f\n",P_I);
printf("Window Size:\t%i\n\n",win);


	maxR=new double[Nmar];
	for (m=0;m<Nmar;m++) maxR[m]=-1;

	fi = fopen(infile, "r");
	fgetpos(fi,&pos);
	fclose(fi);

	w=win;
	for (m=0;m<Nmar;m++)
	{
		fi = fopen(infile, "r");
		fsetpos(fi,&pos);

		a1=0; a2=0; m--;
		while (!feof(fi) && (a1==0 || a2==0))
		{
			fscanf(fi, "%s", &MarName);		
		
			a1=0; a2=0;
			for (j=0;j<N;j++)
			{
				fscanf(fi, "%i", &gen1[j]);
				if (gen1[j]==1) a1+=2;
				if (gen1[j]==2) {a1++; a2++;} 
				if (gen1[j]==3) a2+=2;
			}
			m++;
			if (m>=Nmar-1) break; 		
		}		
		maxR[m]=0;		
		fgetpos(fi,&pos);

		if (m==Nmar-1)
			int a = 0;

		if (m+w > Nmar) w=Nmar-m;
		ii=m+1;
		k =m+1;
		while (ii<m+w)
		{
			if (k>=Nmar) break;
			a1=0; a2=0;
			while (!feof(fi) && (a1==0 || a2==0))
			{
				fscanf(fi, "%s", &MarName);
				for (j=0;j<9;j++) gg[j]=0;
				a1=0; a2=0;
				for (j=0;j<N;j++)
				{
					fscanf(fi, "%i", &gen);
					if (gen==1 && gen1[j]==1) gg[0]++; // 11 11
					if (gen==1 && gen1[j]==2) gg[1]++; // 11 12
					if (gen==1 && gen1[j]==3) gg[2]++; // 11 22
					if (gen==2 && gen1[j]==1) gg[3]++; // 12 11
					if (gen==2 && gen1[j]==2) gg[4]++; // 12 12
					if (gen==2 && gen1[j]==3) gg[5]++; // 12 22
					if (gen==3 && gen1[j]==1) gg[6]++; // 22 11
					if (gen==3 && gen1[j]==2) gg[7]++; // 22 12
					if (gen==3 && gen1[j]==3) gg[8]++; // 22 22
				
					if (gen==1) a1+=2;
					if (gen==2) {a1++; a2++;} 
					if (gen==3) a2+=2;
				}

				k++;
				if (k>=Nmar) break; 
			}

			rr=r(gg); 
			if (rr<0) rr=-rr;
			
			if (rr<=1)
			{
				ii++;				
				if (maxR[m]<rr) maxR[m] = rr;
			}
			k++;
		}
		fclose(fi);
	}


	alpha = 1/(double)Nmar;
	Keff  = Nmar;
	d = 1-pow((1-alpha),Keff);

	eps = alpha/100; //precision

	while (d>P_I)
	{
		alpha -= eps;
		if (alpha > 0)
		{
			Keff = 0;
			for (k=0; k<Nmar; k++) 
				if (maxR[k]>=0)
					Keff += sqrt(1-pow(maxR[k],-1.31*log10(alpha)));

			d = 1-pow((1-alpha),Keff);
		}
		else 
			break;
	}

	fo = fopen(outKeff, "w");
	fprintf(fo,"Estimation of the number of effective tests...\n");
	fprintf(fo,"The authors: V Moskvina and KM Schmidt. (Version: September 2007).\n");
	fprintf(fo,"N of markers:\t\t%i\n",Nmar);
	fprintf(fo,"N of individuals:\t%i\n",N);
	fprintf(fo,"Overall significance level:\t%f\n",P_I);
	fprintf(fo,"Window Size:\t%i\n\n",win);
	fprintf(fo,"Keff=\t%f\n",Keff);
	fclose(fo);
	delete maxR;

delete gen1;

printf("The number of effective tests is:\t%.2f\n",Keff);
printf("(the results are written to %s)\n",outKeff);

return 0;
}





















double r(int g[9])
{
double dd, p1, p2, q1, q2, f[4];

	if (EM(f, g)==0)
	{
		dd = f[0]*f[3]-f[1]*f[2];
		p1 = f[0]+f[1];
		p2 = f[2]+f[3];
		q1 = f[0]+f[2];
		q2 = f[1]+f[3];
		if (p1>0 && p2 >0 && q1>0 && q2>0)
		{
			return dd/sqrt(p1*p2*q1*q2);
		}
		else return -2.0;
	}
	else return -2.0;
}


int c(int i,int j,int k)
{ //0; 1 - if haplotype pair (i,j) is compatible to genotype k
	if      (i==0 && j==0 && k==0) return 1;
	else if (i==0 && j==1 && k==1) return 1;
	else if (i==0 && j==2 && k==3) return 1;
	else if (i==0 && j==3 && k==4) return 1;
	
	else if (i==1 && j==0 && k==1) return 1;
	else if (i==1 && j==1 && k==2) return 1;
	else if (i==1 && j==2 && k==4) return 1;
	else if (i==1 && j==3 && k==5) return 1;

	else if (i==2 && j==0 && k==3) return 1;
	else if (i==2 && j==1 && k==4) return 1;
	else if (i==2 && j==2 && k==6) return 1;
	else if (i==2 && j==3 && k==7) return 1;

	else if (i==3 && j==0 && k==4) return 1;
	else if (i==3 && j==1 && k==5) return 1;
	else if (i==3 && j==2 && k==7) return 1;
	else if (i==3 && j==3 && k==8) return 1;

	else return 0;
}


int z(int i, int j, int r)
{ //0,1,2 how often haplotype r occurs in the haplotype pair (i,j)
	int n = 0;
	if (i==r) n++;
	if (j==r) n++;
	return n;
}


int EM(double h[4], int g[9])
{
int i, j, k, r, t, iter, N, Nhap, Ngen;
double *ht, d, d1, d2, eps;

Nhap = 4;
Ngen = 9;

N=0; for (i=0;i<Ngen;i++) N+=g[i];

iter=10000;
eps=0.000001;

ht = new double[Nhap];

for (r=0; r<Nhap; r++) 	
	h[r]=1/(double)(Nhap);

for (t=0;t<iter;t++)
{
	for (r=0; r<Nhap; r++)
	{
		d = 0;
		for (k=0; k<Ngen; k++)
		{
			d1=0;d2=0;
			if (g[k]>0)
			{
				for (i=0; i<Nhap; i++)
					for (j=0; j<Nhap; j++)
					{
						d1 += c(i,j,k)*z(i,j,r)*h[i]*h[j];
						d2 += c(i,j,k)*h[i]*h[j];
					}
				d+= g[k]*d1/d2;
			}
		}
		ht[r] = d/2/N;
	}

	k=0;
	d=0;
	for (r=0; r<Nhap; r++)
	{
		d1=ht[r]-h[r];
		if (d1<=eps) k++;
		if (d1>d) d=d1;
		h[r]=ht[r];
	}
	if (k==Nhap)
	{
		//printf("Number of iterations %i, eps=%f\n", t, eps);
		break;
	}
	
}
if (t>iter)
{
	printf("EM did not converge.\nNumber of iterations %i, final precision=%f\n", iter, d);
	return -1;
}

delete ht;
return 0;
}


