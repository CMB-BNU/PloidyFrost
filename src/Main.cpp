#include <string>
#include <fstream>
#include "ColoredCDBG.hpp"
#include "CDBG.hpp"
#include "CCDBG.hpp"
#include "MyUnitig.hpp"
#include "GmmModel.hpp"
#include "kmc_file.h"
#include <float.h>
using namespace std;
void PrintUsage()
{
	cout << "Usage: PloidyFrost -g <BifrostGraph> -d <KMCDatabase> -o <outfile_prefix> ..." << endl
		 << endl;
	cout << "parameters with required argument:" << endl
		 << endl
		 << "  -g,             Input Bifrost Graph file (GFA format)"
		 << endl
		 << "  -f,             Input Bifrost color file (BFG_COLORS format,colored CDBG)"
		 << endl
		 << "  -o,             Prefix for Output files (default : 'output')"
		 << endl
		 << "  -t,             Number of Threads (default is 1)"
		 << endl
		 << "  -d,             Load KMC Database"
		 << endl
		 << "  -h,             kmer Histogram file"
		 << endl
		 << "  -l,             Lower coverage threshold (default : 10 )"
		 << endl
		 << "  -u,             Upper coverage threshold (default : 1000 )"
		 << endl
		 << "  -C,             Input Coverage thresholds file (colored CDBG)"
		 << endl
		 << "  -z,             Maximum number of unitigs in superbubble (default : 8 )"
		 << endl
		 << "  -M,             Match score (default : 2 )"
		 << endl
		 << "  -D,             Mismatch penalty (default : -1 )"
		 << endl
		 << "  -G,             Gap penalty (default : -3 )"
		 << endl
		 << endl
		 << "parameters with no argument:" << endl
		 << endl
		 << "  -v,             Print information messages during construction"
		 << endl
		 << "  -i,             Output Information about Bifrost graph"
		 << endl
		 << endl;

	cout << "Usage: PloidyFrost model" << endl
		 << "GMM model"
		 << endl
		 << "  -f,             Prefix of coverage files"
		 << endl
		 << "  -g,             Allele frequency file"
		 << endl
		 << "  -l,             Minimum ploidy level"
		 << endl
		 << "  -u,             Maximum ploidy level"
		 << endl
		 << "  -q,             Minimum allele frequency"
		 << endl
		 << "  -m,             Weigth minimum threshold ( 1/(p-1)/m , default value of m is 5)"
		 << endl
		 << "  -n,             Weigth minimum threshold ( first_weight > maximum_weight/(p-1)/n , default value of n is 2)"
		 << endl
		 << "  -k,             Maximum iterations"
		 << endl
		 << "  -a,             Maximum delta"
		 << endl
		 << "  -o,             Output prefix"
		 << endl
		 << endl;
	cout << "Usage: PloidyFrost cutoffL  " << endl
		 << "Compute Lower coverage threshold" << endl
		 << "PloidyFrost cutoffL kmer_histogram_file " << endl
		 << endl;
	cout << "Usage: PloidyFrost cutoffU " << endl
		 << "Compute Upper coverage threshold" << endl
		 << "PloidyFrost cutoffU kmer_histogram_file [quantile[<1 ,default:0.998]]" << endl
		 << endl;
}
/**
 * @param graphfile: GFA format graph file
 * @param colorfile: BFG_COLORS format file of colored CDBG
 * @param nb_threads: number of threads (default: 1)
 * @param coverage_lower: [CDBG] a cutoff value to filter the "error" structures with low-coverage unitigs
 * @param coverage_upper: [CDBG] a cutoff value to filter the "repeat" structures with high-coverage unitigs
 **/
struct Options
{
	string graphfile;
	string colorfile;
	size_t nb_threads;
	bool verbose;
	int coverage_lower;
	int coverage_upper;
	size_t complex_size;
	string coveragefile;
	double frequency;
	string outprefix;
	int k;
	bool info;
	string db;
	bool bubble;
	double delta;
	vector<pair<int, int>> coverage_vec;
	string hist;
	bool p;
	double mthreshold;
	double nthreshold;
	double match;
	double mismatch;
	double gap;
	Options() : delta(0.01), nb_threads(1), verbose(false), outprefix("output"),
				k(25), info(false), bubble(false), p(true), coverage_lower(10), frequency(0.998),
				coverage_upper(1000), complex_size(8), mthreshold(5.0), nthreshold(2.0), match(2), mismatch(-1), gap(-3){};
};
void parseOptions(int argc, char **argv, Options &opt)
{
	int oc;
	while ((oc = getopt(argc, argv, "M:D:G:z:a:l:q:u:e:C:R:o:t:g:f:k:d:m:n:h:ibvpNSc")) != -1)
	{
		switch (oc)
		{
		case 'z':
			opt.complex_size = atoi(optarg);
			break;
		case 'q':
			opt.frequency = atof(optarg);
			break;
		case 'm':
			opt.mthreshold = atof(optarg);
			break;
		case 'n':
			opt.nthreshold = atof(optarg);
			break;
		case 'M':
			opt.match = atof(optarg);
			break;
		case 'D':
			opt.mismatch = atof(optarg);
			break;
		case 'G':
			opt.gap = atof(optarg);
			break;
		case 'u':
			opt.coverage_upper = atoi(optarg);
		case 'C':
			opt.coveragefile = optarg;
			break;
		case 'a':
			opt.delta = atof(optarg);
			break;
		case 'h':
			opt.hist = optarg;
			break;
		case 'g':
			opt.graphfile = optarg;
			break;
		case 'f':
			opt.colorfile = optarg;
			break;
		case 'o':
			opt.outprefix = optarg;
			break;
		case 'l':
			opt.coverage_lower = atoi(optarg);
			break;
		case 't':
			opt.nb_threads = atoi(optarg);
			break;
		case 'k':
			opt.k = atoi(optarg);
			break;
		case 'v':
			opt.verbose = true;
			break;
		case 'd':
			opt.db = optarg;
			break;
		case 'i':
			opt.info = true;
			break;
		case 'b':
			opt.bubble = true;
			break;
		case 'p':
			opt.p = true;
			break;
		default:
			cout << "Invalid option" << endl;
			PrintUsage();
			exit(EXIT_FAILURE);
		}
	}
}
int cutoffL(const string &file)
{
	ifstream f;
	f.open(file, ios::in);
	if (!f.is_open())
	{
		cout << "ERROR:Open Histogram File " << file << " error!" << endl;
		exit(EXIT_FAILURE);
	}
	vector<size_t> v;

	string s;
	while (getline(f, s, '\n'))
	{
		int pos1 = string::npos;
		pos1 = s.find("\t");
		if (pos1 != string::npos)
		{
			v.emplace_back((atoll(s.substr(pos1 + 1).c_str())));
		}
		else
		{
			cerr << "Error: Histogram File is badly Formatted." << endl;
			exit(EXIT_FAILURE);
		}
	}
	int peakPositions;
	for (peakPositions = 1; peakPositions < v.size(); peakPositions++)
	{
		if (v[peakPositions - 1] < v[peakPositions])
		{
			break;
		}
	}
	return int(round(1.25 * (peakPositions - 1)));
}
int cutoffH(const string &file, double frequency = 0.998)
{
	ifstream f;
	f.open(file, ios::in);
	if (!f.is_open())
	{
		cout << "ERROR:Open Histogram File " << file << " error!" << endl;
		exit(EXIT_FAILURE);
	}
	vector<size_t> v;
	v.emplace_back(0);
	string s;
	while (getline(f, s, '\n'))
	{
		int pos1 = string::npos;
		pos1 = s.find("\t");
		if (pos1 != string::npos)
		{
			v.emplace_back((atoll(s.substr(pos1 + 1).c_str())) + v.back());
		}
		else
		{
			cerr << "Error: Histogram File is badly Formatted." << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (v.size() <= 2)
	{
		cerr << "Error: Histogram File is badly Formatted." << endl;
		exit(EXIT_FAILURE);
	}
	size_t cf = frequency * (v.back() - v[1]) + v[1];
	size_t peakPositions;
	for (peakPositions = 2; peakPositions < v.size(); peakPositions++)
	{
		if (v[peakPositions] > cf)
		{
			break;
		}
	}
	return peakPositions;
}
bool check_ProgramOptions(Options &opt)
{
	bool ret = true;
	size_t max_threads = std::thread::hardware_concurrency();
	if (opt.nb_threads <= 0)
	{
		cerr << "Error: Number of threads cannot be less than or equal to 0." << endl;
		ret = false;
	}
	if (opt.nb_threads > max_threads)
	{
		cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
		ret = false;
	}
	if (opt.p == true)
	{
		if (!opt.db.empty())
		{
			uint8_t kmc_db_num = 0;
			auto checkKmerDb = [&](const string &name) -> bool
			{
				FILE *fp = fopen(name.c_str(), "rb");
				if (fp == NULL)
				{
					cerr << "Error: Could not read the input kmc database " << name << "." << endl;
					ret = false;
					return false;
				}
				else
					fclose(fp);
				return true;
			};
			if (opt.colorfile.empty())
			{
				if (!check_file_exists(opt.db + ".kmc_pre") || !check_file_exists(opt.db + ".kmc_suf"))
				{
					ret = false;
				}
				else
				{
					checkKmerDb(opt.db + ".kmc_pre");
					checkKmerDb(opt.db + ".kmc_suf");
				}
			}
			else
			{
				if (!check_file_exists(opt.db))
				{
					ret = false;
				}
				else
				{
					ifstream infile;
					infile.open(opt.db);
					if (infile.fail())
					{
						ret = false;
					}
					else
					{
						string name;
						vector<string> db;
						while (getline(infile, name, '\n'))
						{
							db.push_back(name);
							if (!checkKmerDb(name + ".kmc_pre"))
								break;
							if (!checkKmerDb(name + ".kmc_suf"))
								break;
						}
						kmc_db_num = db.size();
					}
				}
			}
			if (!opt.hist.empty())
			{
				if (opt.colorfile.empty())
				{
					opt.coverage_lower = max(10, cutoffL(opt.hist));
					opt.coverage_upper = cutoffH(opt.hist, opt.frequency);
				}
				else
				{
					if (!check_file_exists(opt.hist))
					{
						ret = false;
					}
					else
					{
						ifstream infile;
						infile.open(opt.hist);
						if (infile.fail())
						{
							ret = false;
						}
						else
						{
							string name;
							uint8_t i = 0;
							while (getline(infile, name, '\n'))
							{
								opt.coverage_vec.push_back(pair<int, int>(max(10, cutoffL(name)), cutoffH(name, opt.frequency)));
								if (opt.coverage_vec[i].first > opt.coverage_vec[i].second)
								{
									cerr << "Error: lower cutoff need be smaller than upper cutoff " << endl;
									ret = false;
								}
								i++;
							}
							if (i != kmc_db_num)
							{
								cerr << "ERROR: the numbers of kmc databases and hist files are not equal! " << endl;
								exit(EXIT_FAILURE);
							}
						}
					}
				}
			}
			else
			{
				if (!opt.colorfile.empty())
				{
					if (!opt.coveragefile.empty())
					{
						if (!check_file_exists(opt.coveragefile))
						{
							ret = false;
						}
						else
						{
							ifstream infile;
							infile.open(opt.coveragefile);
							if (infile.fail())
							{
								ret = false;
							}
							else
							{
								string name;
								uint8_t i = 0;
								while (getline(infile, name, '\n'))
								{
									int pos1 = string::npos;
									pos1 = name.find("\t");
									if (pos1 != string::npos)
									{
										opt.coverage_vec.push_back(pair<int, int>(atoll(name.substr(0, pos1).c_str()), atoll(name.substr(pos1 + 1).c_str())));
									}
									else
									{
										cerr << "Error: Coverage File is badly Formatted." << endl;
										exit(EXIT_FAILURE);
									}
									if (opt.coverage_vec[i].first < 0 || opt.coverage_vec[i].second < 0)
									{
										cerr << "Error: Filter coverage need a positive number." << endl;
										ret = false;
									}
									if (opt.coverage_vec[i].first > opt.coverage_vec[i].second)
									{
										cerr << "Error: lower cutoff need be smaller than upper cutoff " << endl;
										ret = false;
									}
									i++;
								}
								if (i != kmc_db_num)
								{
									cerr << "ERROR: the numbers of kmc databases and coverages are not equal! " << endl;
									exit(EXIT_FAILURE);
								}
							}
						}
					}
					else
					{
						opt.coverage_vec.insert(opt.coverage_vec.end(), kmc_db_num, pair<int, int>(10, 1000));
					}
				}
			}
		}
		else
		{
			cerr << "Error: Need input a kmc database prefix!\n";
			ret = false;
		}
		opt.bubble = true;
	}
	if (opt.complex_size < 4)
	{
		cerr << "Error: Maximum number of unitigs in superbubble is at least 4 !" << endl;
		ret = false;
	}
	if (opt.mismatch > opt.match)
	{
		cerr << "Error: Mismatch penalty should be smaller than match score !" << endl;
		ret = false;
	}
	if (opt.gap > opt.match)
	{
		cerr << "Error: Gap penalty should be smaller than match score !" << endl;
		ret = false;
	}
	if (opt.outprefix.length() == 0)
	{
		cerr << "Error: No output filename prefix given." << endl;
		ret = false;
	}
	if (opt.coverage_lower < 0 || opt.coverage_upper < 0)
	{
		cerr << "Error: Filter coverage need a positive number." << endl;
		ret = false;
	}
	if (opt.coverage_lower > opt.coverage_upper)
	{
		cerr << "Error: lower cutoff need be smaller than upper cutoff " << endl;
		ret = false;
	}
	if (opt.frequency < 0 || opt.frequency > 1)
	{
		cerr << "Error: Filter coverage need a positive number." << endl;
		ret = false;
	}
	if (opt.graphfile.empty())
	{
		cerr << "Error: No graph file was provided in input." << endl;
		ret = false;
	}
	else if (!check_file_exists(opt.graphfile))
	{
		cerr << "Error: The graph file does not exist." << endl;
		ret = false;
	}
	else
	{
		FILE *fp = fopen(opt.graphfile.c_str(), "r");
		if (fp == NULL)
		{
			cerr << "Error: Could not read input graph file " << opt.graphfile << "." << endl;
			ret = false;
		}
		else
			fclose(fp);
	}
	if (!opt.colorfile.empty())
	{
		if (!check_file_exists(opt.colorfile))
		{
			cerr << "Error: The input color file does not exist." << endl;
			ret = false;
		}
		else
		{
			FILE *fp = fopen(opt.colorfile.c_str(), "rb");
			if (fp == NULL)
			{
				cerr << "Error: Could not read input color file " << opt.colorfile << "." << endl;
				ret = false;
			}
			else
				fclose(fp);
		}
	}
	return ret;
}
bool check_Model_ProgramOptions(Options &opt)
{
	bool ret = true;
	if (opt.coverage_lower > opt.coverage_upper)
	{
		cerr << "Error:  min gauss <= max gauss  " << endl;
		ret = false;
	}
	if (opt.coverage_lower < 1 || opt.coverage_upper < 1)
	{
		cerr << "Error: gauss > 0  " << endl;
		ret = false;
	}
	if (opt.frequency >= 0.5)
	{
		cerr << "Error: frequency cutoff value should < 0.5  " << endl;
		ret = false;
	}
	if (opt.k < 0)
	{
		cerr << "Error: iterate count should > 0 " << endl;
		ret = false;
	}
	if (opt.delta < 0)
	{
		cerr << "Error: iterate delta should > 0 " << endl;
		ret = false;
	}
	if (opt.mthreshold < 0)
	{
		cerr << "Error: minimum threshold should > 0 " << endl;
		ret = false;
	}
	if (opt.nthreshold < 0)
	{
		cerr << "Error: minimum threshold should > 0 " << endl;
		ret = false;
	}
	if (opt.colorfile.empty() && opt.graphfile.empty())
	{
		cout << "ERROR: input a frequency or coverage file " << endl;
		ret = false;
	}
	if (!opt.graphfile.empty())
	{
		ifstream f;
		f.open(opt.graphfile, ios::in);
		if (!f.is_open())
		{
			cout << "ERROR: open frequency file " << opt.graphfile << " error!" << endl;
			ret = false;
		}
		f.close();
	}
	if (!opt.colorfile.empty())
	{
		ifstream f;
		f.open(opt.colorfile + "_bicov.txt", ios::in);
		if (!f.is_open())
		{
			cout << "ERROR: open coverage file " << opt.colorfile + "_bicov.txt"
				 << " error!" << endl;
			ret = false;
		}
		f.close();
		f.open(opt.colorfile + "_tricov.txt", ios::in);
		if (!f.is_open())
		{
			cout << "ERROR: open coverage file " << opt.colorfile + "_tricov.txt"
				 << " error!" << endl;
			ret = false;
		}
		f.close();
		f.open(opt.colorfile + "_tetracov.txt", ios::in);
		if (!f.is_open())
		{
			cout << "ERROR: open coverage file " << opt.colorfile + "_tetracov.txt"
				 << " error!" << endl;
			ret = false;
		}
	}
	return ret;
}

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		PrintUsage();
	}
	else
	{
		Options opt;

		if (strcmp(argv[1], "model") == 0)
		{
			opt.coverage_lower = 1;
			opt.coverage_upper = 9;
			opt.frequency = 0;
			opt.k = 1000;
			opt.delta = 0.01;
			parseOptions(argc, argv, opt);
			if (check_Model_ProgramOptions(opt))
			{
				GmmModel model;
				model.setMThreshold(opt.mthreshold);
				model.setNThreshold(opt.nthreshold);
				model.setMaxIterNum(opt.k);
				model.setMaxDeltaNum(opt.delta);
				if (opt.colorfile.empty() == false)
				{
					model.readCovFile(opt.colorfile, opt.frequency);
				}
				else
				{
					model.readFreFile(opt.graphfile, opt.frequency);
				}
				ofstream outfile;
				outfile.open(opt.outprefix + "_model_result.txt", ios::out | ios::trunc);
				if (!outfile.is_open())
				{
					cout << "ERROR: open output file " << opt.outprefix << "_model_result.txt error!" << endl;
					exit(EXIT_FAILURE);
				}
				double maxll = DBL_MIN;
				double minaic = DBL_MAX;
				double ll_p = 0;
				double aic_p = 0;
				for (int i = opt.coverage_lower; i <= opt.coverage_upper; i++)
				{
					model.resize(i);
					model.emIterate();

					model.output(outfile);
					if (model.getLogLikelihood() > maxll)
					{
						maxll = model.getLogLikelihood();
						ll_p = i + 1;
					}
					if (model.getAIC() < minaic)
					{
						minaic = model.getAIC();
						aic_p = i + 1;
					}
				}
				outfile << "max loglikelihood : " << maxll << "\tploidy : " << ll_p << endl;
				outfile << "min AIC : " << minaic << "\tploidy : " << aic_p << endl;
				outfile << "estimated ploidy level is : " << aic_p << endl;
				outfile.close();
			}
			else
			{
				cout << "Usage: PloidyFrost model" << endl
					 << "GMM model"
					 << endl
					 << "  -f,             Prefix of coverage files"
					 << endl
					 << "  -g,             Allele frequency file"
					 << endl
					 << "  -l,             Minimum ploidy level"
					 << endl
					 << "  -u,             Maximum ploidy level"
					 << endl
					 << "  -q,             Minimum allele frequency"
					 << endl
					 << "  -m,             Weigth minimum threshold ( each_weight > 1/(p-1)/m , default value of m is 5)"
					 << endl
					 << "  -n,             Weigth minimum threshold ( first_weight > maximum_weight/(p-1)/n , default value of n is 2)"
					 << endl
					 << "  -k,             Maximum iterations"
					 << endl
					 << "  -a,             Maximum delta"
					 << endl
					 << "  -o,             Output prefix"
					 << endl
					 << endl;
			}
			return 0;
		}
		if (strcmp(argv[1], "cutoffL") == 0)
		{
			if (argc != 3)
			{
				cout << "Usage:PloidyFrost cutoffL kmer_histogram_file" << endl;
				exit(EXIT_FAILURE);
			}
			cout << max(10, cutoffL(argv[2])) << endl;
			return 0;
		}
		if (strcmp(argv[1], "cutoffU") == 0)
		{
			if (argc == 3)
			{
				cout << cutoffH(argv[2]) << endl;
			}
			else if (argc == 4)
			{
				double y;
				try
				{
					y = stod(argv[3]);
				}
				catch (const exception &)
				{
					cout << "Usage:PloidyFrost cutoffU kmer_histogram_file (quantile[<1 ,default:0.998]) " << endl;
					exit(EXIT_FAILURE);
				}
				if (y >= 1)
				{
					cout << "Usage:PloidyFrost cutoffU kmer_histogram_file (quantile[<1 ,default:0.998]) " << endl;
					exit(EXIT_FAILURE);
				}
				cout << cutoffH(argv[2], y);
			}
			else
			{
				cout << "Usage:PloidyFrost cutoffU kmer_histogram_file (quantile[<1 ,default:0.998]) " << endl;
				exit(EXIT_FAILURE);
			}
			return 0;
		}

		parseOptions(argc, argv, opt);
		if (check_ProgramOptions(opt) == false)
		{
			PrintUsage();
			return 0;
		}
		if (opt.graphfile.empty())
		{
			cout << "No input file given to load Bifrost graph!" << endl;
			exit(EXIT_FAILURE);
		}
		if (!opt.colorfile.empty())
		{
			ColoredCDBG<MyUnitig> cdbg;
			time_t start_time = time(NULL);
			if (cdbg.read(opt.graphfile, opt.colorfile, opt.nb_threads, opt.verbose))
			{
				cout << "ColoredCDBG::read(): Graph loading successful" << endl;
			}
			else
			{
				cout << "ColoredCDBG::read(): Graph could not be loaded! Exit." << endl;
				exit(EXIT_FAILURE);
			}
			time_t end_time = time(NULL);
			cout << "CCDBG: Graph loading Real time : "
				 << (double)difftime(end_time, start_time) << "s" << endl;
			CCDBG g(cdbg, opt.complex_size, opt.match, opt.mismatch, opt.gap, opt.db, opt.nb_threads);
			g.setUnitigId(opt.outprefix, opt.graphfile, opt.nb_threads);
			if (opt.info)
			{
				g.printInfo(opt.verbose, opt.outprefix);
			}

			if (opt.bubble)
			{
				g.findSuperBubble_multithread_ptr(opt.outprefix, opt.nb_threads);
			}
			if (opt.p)
			{
				for (size_t i = 0; i < opt.coverage_vec.size(); i++)
				{
					cout << "CCDBG:: Database " << i << " Minimum Coverage:" << opt.coverage_vec[i].first << endl;
					cout << "CCDBG:: Maximum Coverage:" << opt.coverage_vec[i].second << endl;
				}
				if (opt.bubble)
					g.ploidyEstimation_multithread_ptr(opt.outprefix, opt.coverage_vec, opt.nb_threads);
			}
		}
		else
		{
			CompactedDBG<MyUnitig> cdbg;
			time_t start_time = time(NULL);
			if (cdbg.read(opt.graphfile, opt.nb_threads, opt.verbose))
			{
				cout << "CompactedDBG::read(): Graph loading successful" << endl;
			}
			else
			{
				cout << "CompactedDBG::read(): Graph could not be loaded! Exit." << endl;
				exit(EXIT_FAILURE);
			}
			time_t end_time = time(NULL);
			cout << "CDBG: Graph loading Real time : "
				 << (double)difftime(end_time, start_time) << "s" << endl;
			CDBG g(cdbg, opt.complex_size, opt.match, opt.mismatch, opt.gap, opt.db);
			g.setUnitigId(opt.outprefix, opt.graphfile, opt.nb_threads);
			if (opt.info)
			{
				g.printInfo(opt.verbose, opt.outprefix);
			}

			if (opt.bubble)
			{
				if (opt.p)
				{
					g.findSuperBubble_multithread_ptr(opt.outprefix, opt.nb_threads);
				}
				else
					g.findSuperBubble_multithread_ptr(opt.outprefix, opt.nb_threads);
			}
			if (opt.p)
			{
				cout << "CDBG:: Minimum Coverage:" << opt.coverage_lower << endl;
				cout << "CDBG:: Maximum Coverage:" << opt.coverage_upper << endl;
				g.ploidyEstimation_multithread_ptr(opt.outprefix, opt.coverage_lower, opt.coverage_upper, opt.nb_threads);
			}
		}
	}
}
