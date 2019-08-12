ifstream file;
int lines;
char filename[512] = "data";

file.open(filename, std::ios::in);

lines = Countlines(filename);

float data[lines][2] = { 0 };

Info<<lines<<endl;


for (int i = 0; i<lines; i++)
{
	for (int j = 0; j < 2; j++)
	{
		file >> data[i][j];
	}
}

//Info<<data[9][1]<<endl;

file.close();

int assimilationTimeIndex = 0;
scalar assimilationObservation = data[assimilationTimeIndex][1];
