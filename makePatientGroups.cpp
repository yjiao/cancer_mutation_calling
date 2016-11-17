// 2016 11 08
// Yunxin Joy Jiao
// Simple script that create bam lists for samtools mpileup

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <set>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <memory.h>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <string>
#include <unistd.h>
#include <cassert>

using namespace std;

typedef long long LL; 
typedef pair<int,int> PII;
typedef pair<PII,char> PIIC;

#define MP make_pair
#define PB push_back
#define FF first
#define SS second

#define FORN(i, n) for (int i = 0; i <  (int)(n); i++)
#define ALL(c) (c).begin(), (c).end()
#define FOR1(i, n) for (int i = 1; i <= (int)(n); i++)
#define FORD(i, n) for (int i = (int)(n) - 1; i >= 0; i--)
#define FOREACH(i, c) for (typeof((c).begin()) i = (c).begin(); i != (c).end(); i++)

void split(const string &input, char delim, vector<string> &output) {
    // split a string "input" based on delimiter, return vector "output"
    stringstream input_stream;
    input_stream.str(input);
    string token;
    while (getline(input_stream, token, delim)) {
	output.PB(token);
    }
}

static inline void parse_sample_info(string &path_samples, map<string, vector<string> > &patmap) {
    // parse the sample info file to get a list of patients, store in set "patmap"
    ifstream sample_info(path_samples);
    string line;

    // this assumes that the patient ID is the second field, and that the file is csv
    // Edit accordingly if this is not the case
    char delim = ',';
    size_t idx_pid = 1;
    size_t idx_bam = 2;
    getline(sample_info, line);  // get rid of header line
    while (!sample_info.eof()) {
	getline(sample_info, line);
	//printf("%s\n", line.c_str());

	vector <string> tokens;
	split(line, delim, tokens);

	// for now, ignore patients with MMD1_etc, not sure what these are
	// also ignore whole blood files
	if (tokens.size() >=idx_pid && 
		tokens[idx_pid].find('_') > tokens[idx_pid].length() &&
		tokens[idx_bam].find("WB") > tokens[idx_bam].length()) {
	    patmap[tokens[idx_pid]].PB( tokens[idx_bam] );
	}
    }
}

static inline void write_bam_list(map<string, vector<string> > &patmap, string &prefix_output) {
    for (auto it = patmap.begin(); it != patmap.end(); it++) {
	string pid = it->FF;
	string filename = prefix_output + "/" + pid + ".bamlist";
	ofstream fh(filename);
	vector<string> samps = it->SS;
	
	for (auto samp = samps.begin(); samp != samps.end(); samp++) {
	    fh << *samp << "\n";
	}
	fh.close();
    }
}

int main(int argc, char* argv[]) {
    string path_samples = argv[1];
    string prefix_output = argv[2];
    
    if (access(path_samples.c_str(), F_OK) == -1) {
	printf("Input file does not exist: %s\n", path_samples.c_str());
	return 1;
    }

    map <string, vector<string> > pat_sets;
    // get a set of all patients
    parse_sample_info(path_samples, pat_sets);
    write_bam_list(pat_sets, prefix_output);
}