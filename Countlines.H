#include <iostream>
#include <fstream>

int Countlines(char *filename)
{
	ifstream ReadFile;
	int n=0;
	string tmp;

	ReadFile.open(filename,std::ios::in);

	while(getline(ReadFile,tmp,'\n'))
	{
		n++;
	}
	ReadFile.close();

	return n;
}

scalar analyticalRegion::B = 0.1;
scalar analyticalRegion::rho_l = 684;
