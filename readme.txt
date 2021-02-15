Please cite: Moskvina V, Schmidt KM, On multiple-testing correction in genome-wide association studies. Genet Epidemiol. 2008 Sep; 32(6):567-73.

The input file format is tab delimited with
the number of rows equals to the Number of markers (Nmar)
and the number of columns equals to the number of individuals plus 1 (Nind+1).
First column - marker names;
the other columns are genotypes for each individual.

THERE ARE NO LIMITATIONS FOR THE NUMBERS OF MARKERS AND INDIVIDUALS

For example:
Marker1 1 2 0 3 ...

The genotypes codes as:
1 for 11 genotype
2 for 12 genotype
3 for 22 genotype
0 - missing value


WINDOWS:
To run the program type in the command window:
Keffective.exe FileName Nmar Nind SigLevel WindowSize

For example:
Keffective.exe example.dat 500 1000 0.05 20


LINUX:
To compile the program type: g++ -o Keff Keffective.cpp
To pun the program type: ./Keff example.dat 500 1000 0.05 10


The presision of the results can be changed in the text, line 142: 
eps = alpha/100; //precision
