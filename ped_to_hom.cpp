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
		ss.clear(); ss.str( line );
		for ( int i = 0 ; i < 6 ; i++ ) ss >> meta[i];
		
		seq[0].clear(); seq[1].clear();
		
		while ( !ss.eof() )
		{
			for ( int i = 0 ; i < 2 ; i++ )
			{
				snp = '\0';
				ss >> snp;
				if ( snp == '\0' ) continue;
				seq[i].push_back( snp );
			}
		}
		
		for ( int i = 0 ; i < 2 ; i++ )
		{
			cout << meta[0] << ' ' << meta[1] << '.' << i;
			for ( int m = 2 ; m < 6 ; m++ ) cout << ' ' << meta[m];
			
			for ( int s = 0 ; s < seq[i].size() ; s++ ) cout << ' ' << seq[i][s] << ' ' << seq[i][s];
			cout << endl;
		}
	}

	return 1;
}
