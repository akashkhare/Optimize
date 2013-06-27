#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <stdio.h>
#include <sstream> 

using namespace std;

int main (int argc, char* argv[]) {

	//set path to NGILA
	char cmd_ngila[100] = "ngila";

	vector <char> sequenceA;
	vector <char> sequenceC;
	vector <double> answers;
	char cmd [100];
	char remove [50];
	FILE * pFile;
	string word, output, file; 
	char letter;
	int flag = 0; 
	int mismatch = 0;
	int match = 0;
	int indels = 0;
	int indel_length = 0;
	double nuisance_a, distance_t, rate_r, mean_q, parameter_l;

	sprintf (cmd, "%s %s -o Output.fasta --pairs each -m geo"
			, cmd_ngila, argv[1]);
	system ((char *)cmd);

	file = "Output.fasta";
	ifstream myfile(file.c_str());	
	
	if(!myfile) 
    {
                cout<<"Error opening output file"<<endl;
                system("pause");
                return -1;
    }
	
	sequenceA.push_back('Z');
	while(myfile.get(letter))
    {
		if (letter == '>') {
			flag = flag + 1;
			if ((flag % 2) == 0)
				sequenceA.push_back('Z');
			if ((flag % 2) == 1)
				sequenceC.push_back('Z');
			getline (myfile, word);
		}

		else if ((flag % 2) == 1 && letter != '\n')  
			sequenceA.push_back(letter);
		
		else if ((flag % 2) == 0 && letter != '\n')
			sequenceC.push_back(letter);
    }
	sequenceC.push_back('Z');
	
	for (int i = 0; i < sequenceC.size(); i++) {
	
		if (sequenceC[i] == '-') {
			indel_length = indel_length + 1;
			if (sequenceC[i-1] == 'A' || sequenceC[i-1] == 'T' || sequenceC[i-1] == 'C' 
			    || sequenceC[i-1] == 'G' || sequenceC[i-1] == 'Z') 
						indels = indels + 1;
		}		
		if (sequenceA[i] == '-') {
			indel_length = indel_length + 1;
			if (sequenceA[i-1] == 'A' || sequenceA[i-1] == 'T' || sequenceA[i-1] == 'C' 
			    || sequenceA[i-1] == 'G' || sequenceA[i-1] == 'Z') 
						indels = indels + 1;
		}
			
		if ((sequenceA[i] == 'A' || sequenceA[i] == 'T' || sequenceA[i] == 'C' || sequenceA[i] == 'G') && 
		   (sequenceC[i] == 'A' || sequenceC[i] == 'T' || sequenceC[i] == 'C' || sequenceC[i] == 'G')) {
		   
			if (sequenceA[i] == sequenceC[i])
				match = match + 1;
			else 
				mismatch = mismatch + 1;
		}
	}
	
	nuisance_a = ((double)(match + mismatch + indels) / 100.0);
	mean_q = (double) indel_length / (double) indels;
	distance_t = -0.75 * log (1.0 - ((double)(4 * (mismatch)) / (double)(3 * (match + mismatch))));
	rate_r = -1 * (log10 ((1.0 - ((double) indels / (double) (indels + mismatch + match))))) / (2 * distance_t);
	parameter_l = rate_r * distance_t;
	answers.push_back(distance_t);
	answers.push_back(mean_q);
	answers.push_back(parameter_l);
	answers.push_back(nuisance_a);

  	pFile = fopen ("sim_values.txt", "a");
	for (int i = 0; i < 4; i++) {
		fprintf (pFile, "%-8.6f\t", answers[i]);	
	}
	fprintf (pFile, "\n"); 
	fclose (pFile);

	sprintf (remove, "rm %s", argv[1]);
	system ((char *)remove);

}
