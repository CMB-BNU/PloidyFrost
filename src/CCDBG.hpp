#ifndef CCDBG_HPP_
#define CCDBG_HPP_
#include "ColoredCDBG.hpp"
#include "MyUnitig.hpp"
#include "kmc_file.h"
class CCDBG
{
private:
	ColoredCDBG<MyUnitig> &cdbg;
	CKMCFile kmer_databas;
	size_t complex_size;
	vector<CKMCFile *> kmer_databas_vector;
	double match;
	double mismatch;
	double gap;

public:
	CCDBG(ColoredCDBG<MyUnitig> &cdbg, const size_t &complexsize, double &m, double &d, double &g, string kmc_db = "", const size_t &thread = 1);
	~CCDBG();
	void setUnitigId(const string &outpre, const string &graphfile, const size_t &thr);
	void findSuperBubble_multithread_ptr(const string &outpre, const size_t &thread, const size_t &cov = 10000);
	void sortSeq_simple(vector<size_t> &path_color, vector<UnitigColorMap<MyUnitig>> &unitig_vec, vector<std::vector<double>> &cov_vec, int low, int high);
	void sortSeq_branching(vector<string> &str_vec, int low, int high);
	void extractSuperBubble_multithread_ptr(const UnitigColorMap<MyUnitig> &s, mutex &m);
	void ploidyEstimation_multithread_ptr(const string &outpre, const vector<pair<int, int>> &cutoff, const size_t &thr);
	void setNoBubble_multithread_ptr(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, mutex &x, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p);
	void setNoBubble_multithread_ptr(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p, mutex &x);
	double computeCramerVCoefficient(const vector<double> &A, const vector<double> &B);
	pair<double, bool> readCov(const string &s, const uint &low, const uint &up, const size_t &color_id);
	pair<double, bool> readCovUni(const UnitigColorMap<MyUnitig> &u, const uint &low, const uint &up, const size_t &color_id);
	void setNoBubble_multithread_ptr_cycle(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, mutex &x, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p);
	void setUnitigId(const string &outpre, const string &graphfile);
	void setNoBubble_ptr_cycle(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p);
	void printInfo(const bool &verbose, const string &outpre);
	void clearMarking(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen);
	void ploidyEstimation_ptr(const string &outpre, const vector<pair<int, int>> &cutoff);
	void findSuperBubble_ptr(const string &outpre);
	void extractSuperBubble_ptr(const UnitigColorMap<MyUnitig> &s);
	void setNoBubble_ptr(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &);
	void setNoBubble_ptr(pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &, const vector<UnitigColorMap<MyUnitig>> &vec_km_seen);
};
#endif