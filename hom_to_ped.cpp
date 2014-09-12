#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int main (int argc, char* argv[]) 
{

	string line , meta[6];
	stringstream ss;
	vector< char > seq[2];
	char snp;
	
	while( getline( cin , line ) )
	{
		for ( int i = 0 ; i < 2 ; i++ )
		{
			if ( i == 1 ) getline( cin , line );
			ss.clear(); ss.str( line );
			for ( int m = 0 ; m < 6 ; m++ ) ss >> meta[m];
			seq[i].clear();
			while ( !ss.eof() )
			{
				snp = '\0';
				ss >> snp >> snp;
				if ( snp == '\0' ) continue; 
				
				seq[i].push_back( snp );
			}
		}
		
		cout << meta[0] << ' ' << meta[1].substr(0,meta[1].length()-2);
		for ( int m = 2 ; m < 6 ; m++ ) cout << ' ' << meta[m];
		for ( int s = 0 ; s < seq[0].size() ; s++ ) cout << ' ' << seq[0][s] << ' ' << seq[1][s];
		cout << endl;
	}

	return 1;
}
