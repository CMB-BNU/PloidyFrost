#ifndef CDBG_HPP_
#define CDBG_HPP_

#include "CompactedDBG.hpp"
#include "MyUnitig.hpp"
#include "kmc_file.h"

class CDBG
{

private:
	CompactedDBG<MyUnitig> &cdbg;
	CKMCFile kmer_databas;
	size_t complex_size;
	double match;
	double mismatch;
	double gap;

public:
	CDBG(CompactedDBG<MyUnitig> &cdbg, const size_t &complexsize, double &m, double &d, double &g, string kmc_db = "");
	~CDBG();
	pair<double, bool> readCov(const string &u, const uint &low = 0, const uint &up = 10001);
	pair<double, uint> readCov(const UnitigMap<MyUnitig> &u);
	void sortSeq_simple(vector<double> &cov, vector<UnitigMap<MyUnitig>> &unitig_vec, int low, int high);
	void sortSeq_branching(vector<string> &str_vec, int low, int high);
	void setUnitigId(const string &outpre, const string &, const size_t &thr);
	void printInfo(const bool &verbose, const string &outpre);
	void clearMarking(const vector<UnitigMap<MyUnitig>> &set_km_seen);
	void clearMarking();
	void findSuperBubble_multithread_ptr(const string &outpre, const size_t &thread);
	void extractSuperBubble_multithread_ptr(const UnitigMap<MyUnitig> &s, mutex &m);
	void setNoBubble_multithread_ptr(const vector<UnitigMap<MyUnitig>> &set_km_seen, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &, mutex &x);
	void ploidyEstimation_multithread_ptr(const string &outpre, const int &lower, const int &upper, const size_t &thr);
	void setNoBubble_multithread_ptr(const vector<UnitigMap<MyUnitig>> &set_km_seen, mutex &x, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &);
	void setNoBubble_multithread_ptr_cycle(const vector<UnitigMap<MyUnitig>> &set_km_seen, mutex &x, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &);
	void setNoBubble_ptr_cycle(const vector<UnitigMap<MyUnitig>> &set_km_seen, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &);
	void setNoBubble_ptr(const vector<UnitigMap<MyUnitig>> &set_km_seen, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &);
	void setNoBubble_ptr(pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &, const vector<UnitigMap<MyUnitig>> &set_km_seen);
	void ploidyEstimation_ptr(const string &outpre, const int &lower, const int &upper);
	void findSuperBubble_ptr(const string &outpre);
	void extractSuperBubble_ptr(const UnitigMap<MyUnitig> &s);
};

#endif
