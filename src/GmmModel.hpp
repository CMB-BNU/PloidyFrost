#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
class GmmModel
{
public:
    void resize(size_t g);
    void setMThreshold(double m) { m_thre = m; }
    void setNThreshold(double n) { n_thre = n; }

    void setMaxIterNum(int i) { emMaxIter = i; }
    void setMaxDeltaNum(double i) { emMaxDelta = i; }
    inline double getProbability(double mean, double var, double x)
    {
        return 1 / (sqrt(2 * M_PI * var)) * exp(-(pow(x - mean, 2) / (2 * var)));
    }
    double computeLogLikelihood();
    void emIterate();
    void emStep();
    inline double getLogLikelihood() { return logLikelihood; }
    void readData(vector<double>);
    void readFreFile(const string &filename, const double &freq);
    void readCovFile(const string &filename, const double &freq);
    void output(ofstream &);
    void print();
    inline double computeAIC()
    {
        aic = (2 * (gauss * 2 - 1) - 2 * logLikelihood) / allele_fre.size();
        return aic;
    }
    inline double getAIC()
    {
        return aic;
    }

private:
    vector<double> allele_fre;
    size_t gauss;
    vector<double> weights;
    vector<double> means;
    vector<double> vars;
    double m_thre = 5.0;
    double n_thre = 2.0;
    int emMaxIter = 1000;
    double emMaxDelta = 0.01;
    double logLikelihood;
    double aic = 0;
    vector<vector<int>> cov_vec;
};