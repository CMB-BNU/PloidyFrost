#include <iostream>
#include "GmmModel.hpp"
#include <float.h>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>

void GmmModel::resize(size_t g)
{
    gauss = g;
    vars.resize(g);
    means.resize(g);
    weights.resize(g);
    for (int i = 1; i <= g; i++)
    {
        means[i - 1] = ((double)i / (g + 1));
        weights[i - 1] = ((double)1 / g);
        vars[i - 1] = 0.01;
    }
}
void GmmModel::readCovFile(const string &name, const double &frequency)
{
    allele_fre.clear();
    ifstream bicov;
    ifstream tricov;
    ifstream tetracov;
    ifstream pentacov;

    vector<int> count(10000, 0);
    bicov.open(name + "_bicov.txt", ios::in);
    tricov.open(name + "_tricov.txt", ios::in);
    tetracov.open(name + "_tetracov.txt", ios::in);
    pentacov.open(name + "_pentacov.txt", ios::in);

    if (!bicov.is_open() || !tricov.is_open() || !tetracov.is_open() || !pentacov.is_open())
    {
        cout << "Model::readCovFile() : Open cov file error" << endl;
        exit(EXIT_FAILURE);
    }
    string s = "";
    while (getline(bicov, s, '\n'))
    {
        int pos1 = string::npos;
        pos1 = s.find("\t");
        int pos2 = string::npos;
        pos2 = s.find("\t", pos1 + 1);
        if (pos2 != string::npos)
        {
            vector<int> cov;
            cov.push_back(atoi(s.substr(0, pos1).c_str()));
            cov.push_back(atoi(s.substr(pos1 + 1, pos2).c_str()));
            int cov_sum = cov[0] + cov[1];
            if (cov_sum < 10000)
            {
                if (cov[0] / cov_sum >= frequency && cov[0] / cov_sum <= 1 - frequency)
                {
                    count[cov_sum]++;
                    cov_vec.push_back(cov);
                    allele_fre.push_back(double(cov[0]) / cov_sum);
                    allele_fre.push_back(double(cov[1]) / cov_sum);
                }
            }
        }
        else
        {
            continue;
            cout << "Model::readCovFile() : Open cov file error" << endl;
            bicov.close();
            tricov.close();
            tetracov.close();
            exit(EXIT_FAILURE);
        }
    }
    bicov.close();
    while (getline(tricov, s, '\n'))
    {
        int pos1 = string::npos;
        int pos2 = string::npos;
        int pos3 = string::npos;
        pos1 = s.find("\t");
        pos2 = s.find("\t", pos1 + 1);
        pos3 = s.find("\t", pos2 + 1);
        if (pos3 != string::npos)
        {
            vector<int> cov;
            cov.push_back(atoi(s.substr(0, pos1).c_str()));
            cov.push_back(atoi(s.substr(pos1 + 1, pos2).c_str()));
            cov.push_back(atoi(s.substr(pos2 + 1, pos3).c_str()));
            int cov_sum = cov[0] + cov[1] + cov[2];
            if (cov_sum < 10000)
            {
                int min = cov[0];
                if (cov[1] < cov[0])
                {
                    min = cov[1];
                }
                if (cov[2] < cov[1])
                {
                    min = cov[2];
                }
                if (min / cov_sum >= frequency && min / cov_sum <= 1 - frequency)
                {
                    count[cov_sum]++;
                    cov_vec.push_back(cov);
                    allele_fre.push_back(double(cov[0]) / cov_sum);
                    allele_fre.push_back(double(cov[1]) / cov_sum);
                    allele_fre.push_back(double(cov[2]) / cov_sum);
                }
            }
        }
        else
        {
            continue;
            cout << "Model::readCovFile() : Open cov file error" << endl;
            tricov.close();
            tetracov.close();
            exit(EXIT_FAILURE);
        }
    }
    tricov.close();
    while (getline(tetracov, s, '\n'))
    {
        int pos1 = string::npos;
        int pos2 = string::npos;
        int pos3 = string::npos;
        int pos4 = string::npos;
        pos1 = s.find("\t");
        pos2 = s.find("\t", pos1 + 1);
        pos3 = s.find("\t", pos2 + 1);
        pos4 = s.find("\t", pos3 + 1);
        if (pos4 != string::npos)
        {
            vector<int> cov;
            cov.push_back(atoi(s.substr(0, pos1).c_str()));
            cov.push_back(atoi(s.substr(pos1 + 1, pos2).c_str()));
            cov.push_back(atoi(s.substr(pos2 + 1, pos3).c_str()));
            cov.push_back(atoi(s.substr(pos3 + 1, pos4).c_str()));
            int cov_sum = cov[0] + cov[1] + cov[2] + cov[3];
            if (cov_sum < 10000)
            {
                int min = cov[0];
                if (cov[1] < cov[0])
                {
                    min = cov[1];
                }
                if (cov[2] < cov[1])
                {
                    min = cov[2];
                }
                if (cov[3] < cov[2])
                {
                    min = cov[3];
                }
                if (min / cov_sum >= frequency && min / cov_sum <= 1 - frequency)
                {
                    count[cov_sum]++;
                    cov_vec.push_back(cov);
                    allele_fre.push_back(double(cov[0]) / cov_sum);
                    allele_fre.push_back(double(cov[1]) / cov_sum);
                    allele_fre.push_back(double(cov[2]) / cov_sum);
                    allele_fre.push_back(double(cov[3]) / cov_sum);
                }
            }
        }
        else
        {
            continue;
            cout << "Model::readCovFile() : Open cov file error" << endl;
            tetracov.close();
            exit(EXIT_FAILURE);
        }
    }
    pentacov.close();

    while (getline(pentacov, s, '\n'))
    {
        int pos1 = string::npos;
        int pos2 = string::npos;
        int pos3 = string::npos;
        int pos4 = string::npos;
        int pos5 = string::npos;

        pos1 = s.find("\t");
        pos2 = s.find("\t", pos1 + 1);
        pos3 = s.find("\t", pos2 + 1);
        pos4 = s.find("\t", pos3 + 1);
        pos5 = s.find("\t", pos4 + 1);

        if (pos5 != string::npos)
        {
            vector<int> cov;
            cov.push_back(atoi(s.substr(0, pos1).c_str()));
            cov.push_back(atoi(s.substr(pos1 + 1, pos2).c_str()));
            cov.push_back(atoi(s.substr(pos2 + 1, pos3).c_str()));
            cov.push_back(atoi(s.substr(pos3 + 1, pos4).c_str()));
            cov.push_back(atoi(s.substr(pos4 + 1, pos5).c_str()));

            int cov_sum = cov[0] + cov[1] + cov[2] + cov[3] + cov[4];
            if (cov_sum < 10000)
            {
                int min = cov[0];
                if (cov[1] < cov[0])
                {
                    min = cov[1];
                }
                if (cov[2] < cov[1])
                {
                    min = cov[2];
                }
                if (cov[3] < cov[2])
                {
                    min = cov[3];
                }
                if (cov[4] < cov[3])
                {
                    min = cov[4];
                }
                if (min / cov_sum >= frequency && min / cov_sum <= 1 - frequency)
                {
                    count[cov_sum]++;
                    cov_vec.push_back(cov);
                    allele_fre.push_back(double(cov[0]) / cov_sum);
                    allele_fre.push_back(double(cov[1]) / cov_sum);
                    allele_fre.push_back(double(cov[2]) / cov_sum);
                    allele_fre.push_back(double(cov[3]) / cov_sum);
                    allele_fre.push_back(double(cov[4]) / cov_sum);
                }
            }
        }
        else
        {
            continue;
            cout << "Model::readCovFile() : Open cov file error" << endl;
            pentacov.close();
            exit(EXIT_FAILURE);
        }
    }
    pentacov.close();
}

void GmmModel::readFreFile(const string &filename, const double &frequency)
{
    ifstream fileReader;
    fileReader.open(filename, ios::in);
    if (!fileReader.is_open())
    {
        cout << "ERROR: open frequency file error!" << endl;
        exit(EXIT_FAILURE);
    }
    double a;
    while (!fileReader.eof())
    {
        fileReader >> a;
        if (a >= frequency && a <= 1 - frequency)
            allele_fre.push_back(a);
    }
}
double GmmModel::computeLogLikelihood()
{
    double sum = 0;
    for (auto af : allele_fre)
    {
        double af_sum = 0;
        for (int i = 0; i < gauss; i++)
        {
            af_sum += (weights[i] * getProbability(means[i], vars[i], af));
        }
        if (af_sum == 0.0)
            af_sum = DBL_MIN;
        sum += log(af_sum);
    }
    return sum;
}
void GmmModel::emStep()
{
    vector<double> gaussPart(gauss, 0);
    vector<double> gaussSum(gauss, 0);
    vector<double> varSum(gauss, 0);
    vector<double> meanSum(gauss, 0);
    double sum = 0;
    for (auto af : allele_fre)
    {
        double rowsum = 0;
        for (int i = 0; i < gauss; i++)
        {
            gaussPart[i] = weights[i] * getProbability(means[i], vars[i], af);
            if (gaussPart[i] == 0.0)
                gaussPart[i] = DBL_MIN;
            rowsum += gaussPart[i];
        }
        for (int i = 0; i < gauss; i++)
        {
            gaussPart[i] /= rowsum;
            gaussSum[i] += gaussPart[i];
            varSum[i] += gaussPart[i] * pow(af - means[i], 2);
            meanSum[i] += gaussSum[i] * af;
            sum += gaussPart[i];
        }
    }
    vector<double> new_means = means;
    vector<double> new_vars(gauss, 0);
    vector<double> new_weights(gauss, 0);
    for (int i = 0; i < gauss; i++)
    {

        double var = 1 / gaussSum[i] * varSum[i];
        double weight = gaussSum[i] / sum;
        double mean = 1 / gaussSum[i] * meanSum[i];
        if (var == 0.0)
        {
            var = DBL_MIN;
        }
        new_vars[i] = var;
        new_weights[i] = weight;
    }

    double max_weight = *max_element(new_weights.begin(), new_weights.end());
    if (max_weight != new_weights[0] && max_weight != new_weights[gauss - 1])
    {
        double min_weight = *min_element(new_weights.begin(), new_weights.end());
        if (min_weight < (double)1 / gauss / m_thre)
        {
            return;
        }
        if (min_weight < max_weight / gauss / n_thre)
        {
            return;
        }
    }
    means = new_means;
    vars = new_vars;
    weights = new_weights;
}
void GmmModel::print()
{
    cout << "ploidy:\t" << gauss + 1 << "\tgauss:\t" << gauss << endl;
    cout << "avg loglikelihood:\t" << getLogLikelihood() / allele_fre.size() << endl;
    cout << "AIC:\t" << computeAIC() << endl;
    cout << "means:\t" << endl
         << "\t";
    for (int i = 0; i < gauss; i++)
        cout << means[i] << "\t";
    cout << endl;
    cout << "weights:\t" << endl
         << "\t";
    for (int i = 0; i < gauss; i++)
        cout << weights[i] << "\t";
    cout << endl;
    cout << "variances:\t" << endl
         << "\t";
    for (int i = 0; i < gauss; i++)
        cout << vars[i] << "\t";
    cout << endl
         << endl;
}
void GmmModel::output(ofstream &os)
{
    os << "ploidy : " << gauss + 1 << "\tgauss : " << gauss << endl;
    os << "avg loglikelihood : " << getLogLikelihood() / allele_fre.size() << endl;
    os << "AIC : " << getAIC() << endl;
    os << "means :\t" << endl
       << "\t";
    for (int i = 0; i < gauss; i++)
        os << means[i] << "\t";
    os << endl;
    os << "weights :\t" << endl
       << "\t";
    for (int i = 0; i < gauss; i++)
        os << weights[i] << "\t";
    os << endl;
    os << "variances :\t" << endl
       << "\t";
    for (int i = 0; i < gauss; i++)
        os << vars[i] << "\t";
    os << endl
       << "-----------------------------------" << endl;
}
void GmmModel::emIterate()
{
    logLikelihood = computeLogLikelihood();
    double last = logLikelihood;
    double deltaLogl = DBL_MAX;
    size_t count = 0;
    while ((deltaLogl > emMaxDelta) && (count < emMaxIter))
    {
        emStep();
        last = logLikelihood;
        logLikelihood = computeLogLikelihood();
        deltaLogl = logLikelihood - last;
        ++count;
    }
    computeAIC();
}
void GmmModel::readData(vector<double> data)
{
    allele_fre = data;
}
