// 2016 11 08
// Yunxin Joy Jiao
// This file implements a simple voting algorithm for mutect, strelka, and vardict
// Outputs maflist format for marginal and passed (>= 2 callers), as well as BED format file
// needed for downstream force calling using samtools pileup
// Note that assert statements can be turned off for more speed

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
#include <queue>
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

typedef struct {
    uint8_t count;
    bool mutect;
    bool strelka;
    bool vardict;
} callstats_t;

typedef struct {
    uint8_t chr;  // for now, ignore X, Y, contigs
    uint64_t start;
    uint64_t end;
    string ref;
    string alt;
} mutation_t;

struct mutation_compare {
    // custom compare function for mutation_t objects
    bool operator() (const mutation_t &lhs, const mutation_t &rhs) const {
	// returns true if lhs is before rhs
	if (lhs.chr < rhs.chr) return true;
	if (lhs.chr > rhs.chr) return false;

	// now they have the same chr
	if (lhs.start < rhs.start) return true;
	if (lhs.start > rhs.start) return false;

	// now they have the same chr and same start
	if (lhs.end < rhs.end) return true;
	if (lhs.end > rhs.end) return false;

	// same chr, start, end, so the ref seq will also be the same
	if (lhs.alt < rhs.alt) return true;
	return false;
    }
};

void split(const string &input, char delim, vector<string> &output) {
    // split a string "input" based on delimiter, return vector "output"
    stringstream input_stream;
    input_stream.str(input);
    string token;
    while (getline(input_stream, token, delim)) {
	output.PB(token);
    }
}

static inline void parse_sample_info(string &path_samples, set<string> &all_patients) {
    // parse the sample info file to get a list of patients, store in set "all_patients"
    ifstream sample_info(path_samples);
    string line;

    // this assumes that the patient ID is the second field, and that the file is csv
    // Edit accordingly if this is not the case
    char delim = ',';
    size_t field = 1;
    getline(sample_info, line);  // get rid of header line
    while (!sample_info.eof()) {
	getline(sample_info, line);

	vector <string> split_line;
	split(line, delim, split_line);

	// for now, ignore patients with MMD1_etc, not sure what these are
	if (split_line.size() >=field && split_line[field].find('_') > split_line[field].length()) {
	    all_patients.insert(split_line[field]);
	}
    }
}

callstats_t new_callstats() {
    callstats_t callstats;
    callstats.mutect = 0;
    callstats.strelka = 0;
    callstats.vardict = 0;
    callstats.count = 1;
    return callstats;
}

mutation_t new_mutation(uint8_t chr, uint64_t start, uint64_t end, string ref, string alt) {
    mutation_t mut;
    mut.chr = chr;
    mut.start = start;
    mut.end = end;
    mut.ref = ref;
    mut.alt = alt;
    return mut;
}

static void parse_vardict(string &path_vardict, map<mutation_t, callstats_t, mutation_compare> &mutations) {
    string line;
    ifstream fh(path_vardict);
    while (!fh.eof()) {
	getline(fh, line);
	while (line[0] == '#') getline(fh, line);

	vector<string> tokens;
	split(line, '\t', tokens);
	// chrom pos id ref alt qual filter etc
	uint8_t chr;
	if (tokens.size() >= 5) {  // line is not eof
	    //printf("%s\n", line.c_str());
	    try {
		chr = stoi(tokens[0]);
	    } catch (invalid_argument& e) {
		break;
		// non-numeric chromosome, pass
	    }
	    uint64_t start = stoi(tokens[1]);
	    string ref = tokens[3];
	    string alt = tokens[4];
	    uint64_t end = start + ref.length() - 1;
	    
	    // create new mutation_t object
	    mutation_t mut = new_mutation(chr, start, end, ref, alt);

	    if (mutations.find(mut) == mutations.end()) {  // the current mut is not in the set of mutations found
		//printf("new in vardict: %d, %d, %s, %s\n", chr, start, ref.c_str(), alt.c_str());
		callstats_t callstats = new_callstats();
		callstats.vardict = 1;
		mutations[mut] = callstats;
	    } else {
		mutations[mut].vardict = 1;
		mutations[mut].count++;
	    }
	}
    }
    printf("  parsed vardict, mutation set size = %zu\n", mutations.size());
}

static void parse_strelka(string &path_strelka, map<mutation_t, callstats_t, mutation_compare> &mutations) {
    string line;
    ifstream fh(path_strelka);
    while (!fh.eof()) {
	getline(fh, line);
	while (line[0] == '#') getline(fh, line);

	vector<string> tokens;
	split(line, '\t', tokens);
	// chrom pos id ref alt qual filter etc
	uint8_t chr;
	if (tokens.size() >= 5) {  // line is not eof
	    //printf("%s\n", line.c_str());
	    try {
		chr = stoi(tokens[0]);
	    } catch (invalid_argument& e) {
		break;
		// non-numeric chromosome, pass
	    }
	    uint64_t start = stoi(tokens[1]);
	    string ref = tokens[3];
	    string alt = tokens[4];
	    uint64_t end = start + ref.length() - 1;
	    
	    // create new mutation_t object
	    mutation_t mut = new_mutation(chr, start, end, ref, alt);

	    if (mutations.find(mut) == mutations.end()) {  // the current mut is not in the set of mutations found
		//printf("new in strelka: %d, %d, %s, %s\n", chr, start, ref.c_str(), alt.c_str());
		callstats_t callstats = new_callstats();
		callstats.strelka = 1;
		mutations[mut] = callstats;
	    } else {
		mutations[mut].strelka = 1;
		mutations[mut].count++;
	    }
	}
    }
    printf("  parsed strelka, mutation set size = %zu\n", mutations.size());
}

static void parse_mutect(string &path_mutect, map<mutation_t, callstats_t, mutation_compare> &mutations) {
    string line;
    ifstream fh(path_mutect);
    while (!fh.eof()) {
	getline(fh, line);

	vector<string> tokens;
	split(line, ' ', tokens);
	if (tokens.size() >= 4) {  // line is not eof
	    uint8_t chr;
	    try {
		chr = stoi(tokens[0]);
	    } catch (invalid_argument& e) {
		break;
		// non-numeric chromosome, pass
	    }
	    uint64_t start = stoi(tokens[1]);
	    uint64_t end = start;
	    string ref = tokens[2];
	    string alt = tokens[3];
	    
	    // create new mutation_t object
	    mutation_t mut = new_mutation(chr, start, end, ref, alt);

	    if (mutations.find(mut) == mutations.end()) {  // the current mut is not in the set of mutations found
		callstats_t callstats = new_callstats();
		callstats.mutect = 1;
		mutations[mut] = callstats;
	    } else {
		mutations[mut].mutect = 1;
		mutations[mut].count++;
	    }
	}
    }
    printf("  parsed mutect, mutation set size = %zu\n", mutations.size());
}

static inline void write_bed(ofstream &fh, mutation_t mut) {
    // Bed format is 0-based, half-open
    fh << "chr" << to_string(mut.chr) << '\t';
    fh << mut.start - 1 << '\t';
    fh << mut.end + 1;
    fh << "\n";
}

static inline void write_region(ofstream &fh, mutation_t mut) {
    // Bed format is 0-based, half-open
    fh << to_string(mut.chr) << ':';
    fh << mut.start << '-' << mut.end;
    fh << "\n";
}
static inline void write_line(ofstream &fh, mutation_t mut, callstats_t cs) {
    fh << to_string(mut.chr) << '\t';
    fh << mut.start << '\t';
    fh << mut.end << '\t';
    fh << mut.ref << '\t';
    fh << mut.alt << '\t';
//    printf("%d %s %d %s %d %s\n", cs.mutect, to_string(cs.mutect).c_str(), cs.strelka, to_string(cs.strelka).c_str(), cs.vardict, to_string(cs.vardict).c_str());
    fh << to_string(cs.mutect) << '\t';
    fh << to_string(cs.strelka) << '\t';
    fh << to_string(cs.vardict) << '\t';
    fh << to_string(cs.count);
    fh << "\n";
}

static inline void write_header(ofstream &fh) {
    fh << "chr\tstart\tend\tref_allele\talt_allele\tmutect\tstrelka\tvardict\tcount\n";
}

static inline void process_patient(string &pid, string &prefix_mutect, string &prefix_strelka, string &prefix_vardict, string &prefix_output) {

    map<mutation_t, callstats_t, mutation_compare> mutations;
    printf("pid: %s---------------\n", pid.c_str());
    string path_mutect = prefix_mutect + pid + "_mutect_combined.maf";
    string path_strelka = prefix_strelka + pid + "_strelka_combined.vcf";
    string path_vardict = prefix_vardict + pid + "_vardict_combined.vcf";

    if (access(path_mutect.c_str(), F_OK) != -1) parse_mutect(path_mutect, mutations);
    else printf("!!!!!!!!!!WARNING mutect union does not exist %s\n", path_mutect.c_str());
    if (access(path_strelka.c_str(), F_OK) != -1) parse_strelka(path_strelka, mutations);
    else printf("!!!!!!!!!!WARNING strelka union does not exist %s\n", path_strelka.c_str());
    if (access(path_vardict.c_str(), F_OK) != -1) parse_vardict(path_vardict, mutations);
    else printf("!!!!!!!!!!WARNING vardict union does not exist %s\n", path_vardict.c_str());


    ofstream fh_pass(prefix_output + pid + "_pass.maf");
    ofstream fh_margin(prefix_output + pid + "_marginal.maf");
    ofstream fh_bed(prefix_output + pid + "_pass.bed");
    ofstream fh_region(prefix_output + pid + "_pass.region");
    write_header(fh_pass);
    write_header(fh_margin);
    for (auto it = mutations.begin(); it != mutations.end(); it++) {
	mutation_t mut = it->FF;
	callstats_t callstats = it->SS;
	assert(callstats.mutect == 0 || callstats.mutect == 1);
	assert(callstats.strelka == 0 || callstats.strelka == 1);
	assert(callstats.vardict == 0 || callstats.vardict == 1);
	assert(callstats.mutect + callstats.strelka + callstats.vardict == callstats.count);
	if (callstats.count >= 2){
	    write_line(fh_pass, mut, callstats);
	    write_bed(fh_bed, mut);
	    write_region(fh_region, mut);
	}
	else write_line(fh_margin, mut, callstats);
    }
    fh_pass.close();
    fh_margin.close();
}

int main(int argc, char* argv[]) {
    string path_samples = argv[1];
    string prefix_mutect  = argv[2];
    string prefix_strelka = argv[3];
    string prefix_vardict = argv[4];
    string prefix_output  = argv[5];

    // get a set of all patients
    set<string> all_patients;
    parse_sample_info(path_samples, all_patients);

    // add a dash at the end just in case
    prefix_mutect += "/";
    prefix_strelka += "/";
    prefix_vardict += "/";
    prefix_output += "/";

    // print out the parameters passed into the program, for debugging
    printf("path_samples: %s\n", path_samples.c_str());
    printf("prefix_mutect: %s\n", prefix_mutect.c_str());
    printf("prefix_strelka: %s\n", prefix_strelka.c_str());
    printf("prefix_vardict: %s\n", prefix_vardict.c_str());
    printf("prefix_output: %s\n", prefix_output.c_str());

    for (auto it = all_patients.begin(); it != all_patients.end(); it++) {
	// iterate through list of patients
	// get the paths to the union of all mutations found for mutect, strelka, and vardict
	string pid(*it);
	process_patient(pid, prefix_mutect, prefix_strelka, prefix_vardict, prefix_output);
    }
}