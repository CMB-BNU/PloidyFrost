#include "CCDBG.hpp"
#include <string>
#include "kmc_file.h"
#include <stack>
#include <map>
#include <iostream>
#include "SeqAlign.hpp"
#include <algorithm>
#include <float.h>
using namespace std;
CCDBG::CCDBG(ColoredCDBG<MyUnitig> &graph, const size_t &complexsize, double &m, double &d, double &g, string db, const size_t &thread) : cdbg(graph), complex_size(complexsize), match(m), mismatch(d), gap(g)
{
    if (!db.empty())
    {
        kmer_databas_vector.resize(cdbg.getNbColors());
        ifstream infile;
        infile.open(db);
        if (infile.fail())
        {
            cerr << "CCDBG::CCDBG():Error: Open kmc database name file error" << endl;
            exit(EXIT_FAILURE);
        }
        else
        {
            if (thread == 1)
            {
                for (uint8_t i = 0; i < cdbg.getNbColors(); i++)
                {
                    string name;
                    getline(infile, name, '\n');
                    CKMCFile *kmer_db = new CKMCFile;
                    if (!kmer_db->OpenForRA(name))
                    {
                        cerr << "CCDBG::CCDBG():Error: Open kmc database " << name << " error" << endl;
                        exit(EXIT_FAILURE);
                    }
                    else
                    {
                        cout << "CCDBG::CCDBG(): kmc database " << name << " initialized" << endl;
                        kmer_databas_vector[i] = (kmer_db);
                    }
                }
            }
            else
            {
                vector<std::thread> workers;
                mutex m;
                int i = 0;
                for (size_t t = 0; t < thread; ++t)
                {
                    workers.emplace_back([&]
                                         {
                                             string name;
                                             int threadid = 0;
                                             while (1)
                                             {
                                                 {
                                                     lock_guard<mutex> lock(m);
                                                     if (getline(infile, name, '\n'))
                                                     {
                                                         threadid = i;
                                                         i++;
                                                     }
                                                     else
                                                     {
                                                         break;
                                                     }
                                                 }
                                                 CKMCFile *kmer_db = new CKMCFile;
                                                 if (!kmer_db->OpenForRA(name))
                                                 {
                                                     cerr << "CCDBG::CCDBG():Error: Open kmc database " << name << " error" << endl;
                                                     exit(EXIT_FAILURE);
                                                 }
                                                 else
                                                 {
                                                     cout << "CCDBG::CCDBG(): kmc database " << name << " initialized" << endl;
                                                     kmer_databas_vector[threadid] = (kmer_db);
                                                 }
                                             } });
                }
                for (auto &t : workers)
                    t.join();
            }
        }
    }
    cout << "CCDBG::CCDBG():CCDBG initialized!" << endl;
}
pair<double, bool> CCDBG::readCov(const string &s, const uint &low, const uint &up, const size_t &color_id)
{
    double sum = 0;
    CKmerAPI kmer_object(cdbg.getK());
    uint count = 0;
    if (kmer_databas_vector[color_id]->GetBothStrands())
    {
        for (int i = 0; i <= s.length() - cdbg.getK(); i++)
        {
            kmer_object.from_string(s.substr(i, cdbg.getK()));
            if (!kmer_databas_vector[color_id]->IsKmer(kmer_object))
            {
                kmer_object.reverse();
            }
            if (kmer_databas_vector[color_id]->CheckKmer(kmer_object, count))
            {
                if (count > low && count < up)
                    sum += count;
                else
                {

                    return pair<double, bool>(0, false);
                }
            }
            else
            {
                cout << "CCDBG::readCov():" << kmer_object.to_string() << " kmer can not find" << endl;
                return pair<double, bool>(0, false);
                exit(EXIT_FAILURE);
            }
        }
    }
    return pair<double, bool>(sum / (s.length() - cdbg.getK() + 1), true);
}
pair<double, bool> CCDBG::readCovUni(const UnitigColorMap<MyUnitig> &u, const uint &low, const uint &up, const size_t &color_id)
{
    double sum = 0;
    CKmerAPI kmer_object(cdbg.getK());
    uint count = 0;
    if (kmer_databas_vector[color_id]->GetBothStrands())
    {
        for (int i = 0; i < u.len; i++)
        {
            kmer_object.from_string(u.referenceUnitigToString().substr(i, cdbg.getK()));
            if (!kmer_databas_vector[color_id]->IsKmer(kmer_object))
            {
                kmer_object.reverse();
            }
            if (kmer_databas_vector[color_id]->CheckKmer(kmer_object, count))
            {
                if (count > low && count < up)
                    sum += count;
                else
                {

                    return pair<double, bool>(0, false);
                }
            }
            else
            {
                cout << "CCDBG::readCov():" << kmer_object.to_string() << " kmer can not find" << endl;
                return pair<double, bool>(0, false);
                exit(EXIT_FAILURE);
            }
        }
    }
    return pair<double, bool>(sum / u.len, true);
}
CCDBG::~CCDBG()
{
    for (const auto &db : kmer_databas_vector)
        db->Close();
}

void CCDBG::findSuperBubble_multithread_ptr(const string &outpre, const size_t &thr, const size_t &cov)
{
    if (thr > 1)
    {
        cout << "CCDBG::findSuperBubble(): Finding superbubbles" << endl;
        ofstream outfile;
        if (access("PloidyFrost_output", 0))
        {
            system("mkdir ./PloidyFrost_output");
        }
        outfile.open("PloidyFrost_output/" + outpre + "_super_bubble.txt", ios::out | ios::trunc);
        if (outfile.fail())
        {
            cout << "CCDBG:: Open super_bubble file error";
            exit(EXIT_FAILURE);
        }
        std::mutex mtx;
        std::mutex mtx_uni;
        vector<thread> workers;
        unitigIterator<DataAccessor<MyUnitig>, DataStorage<MyUnitig>> all_iter = cdbg.begin();
        clock_t start_clock = clock();
        double start_time = time(NULL);
        size_t nb_unitig_processed = 0;
        for (size_t t = 0; t < thr; ++t)
        {
            workers.emplace_back([&]
                                 {
                unitigIterator<DataAccessor<MyUnitig>, DataStorage<MyUnitig>> iter= cdbg.end();
                {
                    lock_guard<std::mutex> lock(mtx_uni);
                    if (all_iter != cdbg.end())
                    {
                        iter = all_iter;
                        all_iter++;
                    }
                    else
                    {
                        return;
                    }
                }
                while (iter != cdbg.end())
                {
                    UnitigMap<DataAccessor<MyUnitig>, DataStorage<MyUnitig>> u(*iter);
                    ++nb_unitig_processed;
                    if (nb_unitig_processed % 100000 == 0)
                    {
                        cout << "CCDBG::findSuperBubble(): Processed " << nb_unitig_processed << " unitigs " << endl;
                    }
                            u.strand = true;
                            if (u.getSuccessors().cardinality()>1&&u.getData()->getData(u)->get_plus() == NULL)
                            {
                                extractSuperBubble_multithread_ptr(u, mtx);
                            }                                
                            u.strand = false;
                            if (u.getSuccessors().cardinality()>1&&u.getData()->getData(u)->get_minus() == NULL)
                            {
                                extractSuperBubble_multithread_ptr(u, mtx);
                            }
                        
                        
                    lock_guard<std::mutex> lock(mtx_uni);
                    if (all_iter != cdbg.end())
                    {
                        iter = all_iter;
                        all_iter++;
                    }
                    else
                    {
                        return;
                    }
                } });
        }
        cout << "CCDBG::findSuperBubble(): There are " << cdbg.size() << " unitigs " << endl;
        for (auto &t : workers)
            t.join();
        time_t end_time = time(NULL);
        cout << "CCDBG::findSuperBubble(): Finding superbubbles Cpu time : "
             << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
        cout << "CCDBG::findSuperBubble(): Finding superbubbles Real time : "
             << (double)difftime(end_time, start_time) << "s" << endl;
        outfile << "BubbleId\tEntrance\tStrand\tExit\tisSimple\tisComplex" << endl;
        workers.clear();
        std::mutex mtx_file;
        all_iter = cdbg.begin();
        atomic<size_t> nb_super_bubble(0);
        for (size_t t = 0; t < thr; ++t)
        {
            workers.emplace_back([&]
                                 {
                unitigIterator<DataAccessor<MyUnitig>, DataStorage<MyUnitig>> iter;
                size_t bid=0;
                {
                    lock_guard<std::mutex> lock(mtx_uni);
                    if (all_iter != cdbg.end())
                    {
                        iter = all_iter;
                        all_iter++;
                    }
                    else
                    {
                        return;
                    }
                }
                while (iter != cdbg.end())
                {
                    MyUnitig *data = iter->getData()->getData(*iter);
                    if (data==NULL||!data->is_super())
                    {
                        lock_guard<std::mutex> lock(mtx_uni);
                        if (all_iter != cdbg.end())
                        {
                            iter = all_iter;
                            all_iter++;
                        }
                        else
                        {
                            return;
                        }
                        continue;
                    }
                    if (data->get_plus() != NULL)
                    {
                        bid=nb_super_bubble.fetch_add(1);
                        lock_guard<std::mutex> lock(mtx_file);                    
                        outfile << bid << "\t" << data->get_id() << "\t"
                                << "+"
                                << "\t" << data->get_plus()->get_id() 
                                << "\t" <<(data->isStrict(true)?"1":"0")
                                << "\t" <<(data->isComplex(true)?"1":"0")//<< "\t" <<(data->iscut()?"1":"0")
                                <<endl;
                      
                    }
                    if (data->get_minus() != NULL)
                    {
                        bid=nb_super_bubble.fetch_add(1);
                        lock_guard<std::mutex> lock(mtx_file);
                        outfile << bid << "\t" << data->get_id() << "\t"
                                << "-"
                                << "\t" << data->get_minus()->get_id()
                                << "\t" <<(data->isStrict(false)?"1":"0")
                                << "\t" <<(data->isComplex(false)?"1":"0")//<< "\t" <<(data->iscut()?"1":"0")
                                <<endl;
         
                    }
                
                    lock_guard<std::mutex> lock(mtx_uni);
                    if (all_iter != cdbg.end())
                    {
                        iter = all_iter;
                        all_iter++;
                    }
                    else
                    {
                        return;
                    }
                } });
        }
        for (auto &t : workers)
            t.join();
        outfile.close();
        cout << "CCDBG::findSuperBubble(): " << nb_super_bubble << "  SuperBubbles Found" << endl;
    }
    else
    {
        findSuperBubble_ptr(outpre);
    }
}
double CCDBG::computeCramerVCoefficient(const vector<double> &A, const vector<double> &B)
{
    double n = 0;
    double nA = 0;
    double nB = 0;
    uint8_t count = 0;
    vector<double> p(A.size(), 0);
    double chi = 0;
    for (int i = 0; i < A.size(); i++)
    {
        nA += A[i];
        nB += B[i];
        p[i] = A[i] + B[i];
        n = n + p[i];
        if (p[i] != 0)
        {
            ++count;
        }
    }
    if (count < 2)
    {
        return 0;
    }
    for (int i = 0; i < A.size(); i++)
    {
        if (p[i] == 0)
        {
            continue;
        }
        double exA = nA * p[i] / n;
        double exB = nB * p[i] / n;
        chi += pow(A[i] - exA, 2) / exA;
        chi += pow(B[i] - exB, 2) / exB;
    }

    return sqrt(chi / n);
}

void CCDBG::sortSeq_simple(vector<size_t> &path_color, vector<UnitigColorMap<MyUnitig>> &unitig_vec, vector<std::vector<double>> &cov_vec, int low, int high)
{
    if (high <= low)
        return;
    int i = low;
    int j = high;
    while (true)
    {
        while (path_color[i] >= path_color[low])
        {
            if (path_color[i] > path_color[low])
            {
                i++;
            }
            else
            {
                if (unitig_vec[i].referenceUnitigToString().length() > unitig_vec[low].referenceUnitigToString().length())
                {
                    i++;
                }
                else if (unitig_vec[i].referenceUnitigToString().length() == unitig_vec[low].referenceUnitigToString().length())
                {
                    if (strcmp(unitig_vec[i].referenceUnitigToString().c_str(), unitig_vec[low].referenceUnitigToString().c_str()) > 0)
                    {
                        i++;
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
            if (i == high)
            {
                break;
            }
        }
        while (path_color[j] <= path_color[low])
        {
            if (path_color[j] < path_color[low])
            {
                j--;
            }
            else
            {

                if (unitig_vec[j].referenceUnitigToString().length() < unitig_vec[low].referenceUnitigToString().length())
                {
                    j--;
                }
                else if (unitig_vec[j].referenceUnitigToString().length() == unitig_vec[low].referenceUnitigToString().length())
                {
                    if (strcmp(unitig_vec[j].referenceUnitigToString().c_str(), unitig_vec[low].referenceUnitigToString().c_str()) < 0)
                    {
                        j--;
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
            if (j == low)
            {
                break;
            }
        }
        if (i >= j)
            break;
        double t = path_color[i];
        path_color[i] = path_color[j];
        path_color[j] = t;
        UnitigColorMap<MyUnitig> u = unitig_vec[i];
        unitig_vec[i] = unitig_vec[j];
        unitig_vec[j] = u;
        for (int k = 0; k < cov_vec.size(); k++)
        {
            t = cov_vec[k][i];
            cov_vec[k][i] = cov_vec[k][j];
            cov_vec[k][j] = t;
        }
    }
    double t = path_color[low];
    path_color[low] = path_color[j];
    path_color[j] = t;
    UnitigColorMap<MyUnitig> u = unitig_vec[low];
    unitig_vec[low] = unitig_vec[j];
    unitig_vec[j] = u;
    for (int k = 0; k < cov_vec.size(); k++)
    {
        t = cov_vec[k][low];
        cov_vec[k][low] = cov_vec[k][j];
        cov_vec[k][j] = t;
    }
    sortSeq_simple(path_color, unitig_vec, cov_vec, low, j - 1);
    sortSeq_simple(path_color, unitig_vec, cov_vec, j + 1, high);
}
void CCDBG::sortSeq_branching(vector<string> &str_vec, int low, int high)
{
    if (high <= low)
        return;
    int i = low;
    int j = high;
    while (true)
    {
        while (str_vec[i].length() >= str_vec[low].length())
        {
            if (str_vec[i].length() > str_vec[low].length())
            {
                i++;
            }
            else
            {
                if (strcmp(str_vec[i].c_str(), str_vec[low].c_str()) > 0)
                {
                    i++;
                }
                else
                {
                    break;
                }
            }
            if (i == high)
            {
                break;
            }
        }
        while (str_vec[j].length() <= str_vec[low].length())
        {
            if (str_vec[j].length() < str_vec[low].length())
            {
                j--;
            }
            else
            {
                if (strcmp(str_vec[j].c_str(), str_vec[low].c_str()) < 0)
                {
                    j--;
                }
                else
                {
                    break;
                }
            }
            if (j == low)
            {
                break;
            }
        }
        if (i >= j)
            break;
        string s = str_vec[j];
        str_vec[j] = str_vec[i];
        str_vec[i] = s;
    }
    string s = str_vec[j];
    str_vec[j] = str_vec[low];
    str_vec[low] = s;
    sortSeq_branching(str_vec, low, j - 1);
    sortSeq_branching(str_vec, j + 1, high);
}

void CCDBG::ploidyEstimation_multithread_ptr(const string &outpre, const vector<pair<int, int>> &cutoff, const size_t &thr)
{
    if (thr > 1)
    {
        clock_t start_clock = clock();
        double start_time = time(NULL);
        cout << "CCDBG::PloidyEstimation():  Analyzing superbubbles to generate sites' information" << endl;
        atomic<size_t> var_count_all(0);
        atomic<size_t> coreNum(0);
        atomic<size_t> coreCov(0);
        atomic<size_t> allele4(0);
        atomic<size_t> allele2(0);
        atomic<size_t> allele3(0);
        atomic<size_t> allele5(0);
        if (access("PloidyFrost_output", 0))
        {
            system("mkdir ./PloidyFrost_output");
        }
        ofstream bifre;
        ofstream trifre;
        ofstream tetrafre;
        ofstream bicov;
        ofstream tricov;
        ofstream tetracov;
        ofstream allfre;
        ofstream pentafre;
        ofstream pentacov;
        pentafre.open("PloidyFrost_output/" + outpre + "_pentafre.txt", ios::out | ios::trunc);
        pentacov.open("PloidyFrost_output/" + outpre + "_pentacov.txt", ios::out | ios::trunc);
        allfre.open("PloidyFrost_output/" + outpre + "_allele_frequency.txt", ios::out | ios::trunc);
        bifre.open("PloidyFrost_output/" + outpre + "_bifre.txt", ios::out | ios::trunc);
        trifre.open("PloidyFrost_output/" + outpre + "_trifre.txt", ios::out | ios::trunc);
        tetrafre.open("PloidyFrost_output/" + outpre + "_tetrafre.txt", ios::out | ios::trunc);
        bicov.open("PloidyFrost_output/" + outpre + "_bicov.txt", ios::out | ios::trunc);
        tricov.open("PloidyFrost_output/" + outpre + "_tricov.txt", ios::out | ios::trunc);
        tetracov.open("PloidyFrost_output/" + outpre + "_tetracov.txt", ios::out | ios::trunc);
        ofstream s_var;
        s_var.open("PloidyFrost_output/" + outpre + "_alignseq.txt", ios::out | ios::trunc);
        if (!allfre.is_open() || !bifre.is_open() || !trifre.is_open() ||
            !tetrafre.is_open() || !bicov.is_open() || !tricov.is_open() || !tetracov.is_open() || !s_var.is_open())
        {
            cout << "CCDBG:: PloidyEstimation():Open file error" << endl;
            exit(EXIT_FAILURE);
        }

        vector<thread> workers;
        unitigIterator<DataAccessor<MyUnitig>, DataStorage<MyUnitig>> all_iter = cdbg.begin();
        mutex mtx_uni;
        mutex mtx_file;
        mutex mtx_bicov;
        mutex mtx_tricov;
        mutex mtx_tetracov;
        mutex mtx_pentacov;
        mutex mtx_allfre;
        mutex mtx_var;
        size_t nb_unitig_processed = 0;
        for (size_t t = 0; t < thr; ++t)
        {
            workers.emplace_back([&]
                                 {
                unitigIterator<DataAccessor<MyUnitig>, DataStorage<MyUnitig>> iter;
                {
                    lock_guard<std::mutex> lock(mtx_uni);
                    if (all_iter != cdbg.end())
                    {
                        iter = all_iter;
                        all_iter++;
                    }
                    else
                    {
                        return;
                    }
                }
                while (iter != cdbg.end())
                {
                    ++nb_unitig_processed;
                    if (nb_unitig_processed % 100000 == 0)
                    {
                        cout << "CCDBG::PloidyEstimation(): Processed " << nb_unitig_processed << " unitigs " << endl;
                    }
                    UnitigColorMap<MyUnitig> u(*iter);
                    if ( u.getData()->getData(u)->is_both_visited())
                    
                    {
                        lock_guard<std::mutex> lock(mtx_uni);
                        if (all_iter != cdbg.end())
                        {
                            iter = all_iter;
                            ++all_iter;
                            continue;
                        }
                        else
                        {
                            return;
                        }
                    }
                    while (!u.getData()->getData(u)->is_both_visited())
                    {
                        if(u.getData()->getData(u)->is_plus_visited() != true)
                        {
                            u.strand = true;
                             if(u.getData()->getData(u)->isComplex(true))
                             {
                                 u.getData()->getData(u)->set_plus_visited();
                                 continue;
                             }
                        }
                        else if(u.getData()->getData(u)->is_minus_visited() != true)
                        {
                            u.strand = false;
                             if(u.getData()->getData(u)->isComplex(false))
                             {
                                 u.getData()->getData(u)->set_minus_visited();
                                 break;
                             }
                        }
                        else 
                        {
                            break ;
                        }
                        bool isStrict = u.getData()->getData(u)->isStrict(u.strand);
                        UnitigColorMap<MyUnitig> exit_uni=NULL;
                        bool flag=true;
                        pair<double, uint> core =pair<double,uint>(0,0);
                        for(size_t i=0 ;i<cdbg.getNbColors();i++){
                            pair<double, bool> temp = readCovUni(u, cutoff[i].first, cutoff[i].second,i);
                            if(temp.second==true)
                            {core.first+=temp.first;
                            core.second+=temp.second;}
                            else
                            {
                                break;
                            }
                        }
                        if(flag)
                        {                
                            if(isStrict)
                            {                               
                                exit_uni = (*u.getSuccessors().begin()->getSuccessors().begin());
                
                                if(u.referenceUnitigToString().compare(exit_uni.referenceUnitigToString())<0){
                                    if(u.strand)
                                        u.getData()->getData(u)->set_plus_visited();
                                    else
                                        u.getData()->getData(u)->set_minus_visited();
                                    continue;
                                }                            
                                
                                vector< vector< double > > cov_vec(cdbg.getNbColors());
                                for(int i=0;i<cov_vec.size();i++){
                                    cov_vec[i].resize(u.getSuccessors().cardinality(),0);
                                }            
                                vector< size_t > path_color;
                                uint8_t path=0;
                                vector<UnitigColorMap<MyUnitig>> unitig_vec;
                                
                                for (const auto &uu : u.getSuccessors())
                                {
                                    unitig_vec.push_back(uu);
                                    size_t j=0;
                                    for(size_t i=0;i<cdbg.getNbColors();i++){
                                        if(uu.getData()->getUnitigColors(uu)->contains(uu,i)==true){
                                            j++;
                                            pair<double, bool> inside = readCovUni(uu, cutoff[i].first, cutoff[i].second,i);
                                            if(inside.second==true){
                                                cov_vec[i][path]=inside.first;
                                            }
                                            else{
                                                flag=false;
                                                break;
                                            }
                                        }
                                    }                    
                                    if(flag==false){
                                        break;
                                    }
                                    if(uu.getData()->getUnitigColors(uu)->size(uu)!=j*uu.len){
                                        flag=false;
                                        break;
                                    }
                                    ++path;
                                    path_color.push_back(j);
                                }
                                
                                if (flag == true)
                                {
                                    flag=false;
                                    for(auto v:cov_vec){
                                        int c=0;
                                        for(auto d:v){
                                            if(d!=0.0){
                                                ++c;
                                            }
                                        }
                                        if(c>1){
                                            flag=true;
                                            break;              
                                        }
                                    }
                                }
                                if (flag == true)
                                {

                                    sortSeq_simple(path_color,unitig_vec,cov_vec,0,path_color.size()-1);

                                    SeqAlign seqalign(match,mismatch,gap);
                                    vector<uint> var_site;
                                    vector<uint> snp_pos;
                                    vector<uint> indel_pos;
                                    vector<uint> indel_len;
                                    vector<vector<unsigned short>> partition;
                                    vector<string> str_vec;
                                    for (const auto  & u : unitig_vec)
                                    {
                                        str_vec.push_back(u.mappedSequenceToString());
                                    }
                                    seqalign.SequenceAlignment(str_vec, snp_pos, indel_pos, partition, indel_len);
                                    if (str_vec.size() != 0)
                                    {        
                                        stringstream var_info_stream; 
                                        size_t var_count=var_count_all.fetch_add(1);
                                        for (auto s : str_vec)
                                        {
                                            var_info_stream << var_count << "\t" <<1<<"\t"<<u.getData()->getData(u)->get_id()<<"\t"<<exit_uni.getData()->getData(exit_uni)->get_id()<< "\t" << s << endl;
                                        }
                                        vector<size_t> allele_info(4,0);
                                        stringstream bifre_info; 
                                        stringstream bicov_info; 
                                        stringstream tricov_info; 
                                        stringstream trifre_info; 
                                        stringstream tetracov_info; 
                                        stringstream tetrafre_info; 
                                        stringstream pentafre_info; 
                                        stringstream pentacov_info;  
                                        coreCov.fetch_add((size_t)core.first);
                                        coreNum.fetch_add(1);
                                        for (int i = 0; i < partition.size(); i++)
                                        {
                                            if (partition[i].back() > 0)
                                            {
                                                var_site.push_back(i);
                                            } 
                                        }
                                        
                                        
                                        uint var_distance = 0;
                                        
                                        uint indel = 0;
                                        double coefficient=0;
                                        
                                        for(size_t ci=0;ci<cov_vec.size()-1;ci++){                                        
                                            for(size_t cj=ci+1;cj<cov_vec.size();cj++){ 
                                                coefficient=max(coefficient,computeCramerVCoefficient(cov_vec[ci],cov_vec[cj]));
                                            }
                                        }
                                        for (uint i = 0; i < var_site.size(); i++)
                                        {
                                            stringstream cov_info_all; 
                                            vector<unsigned short> part = partition[var_site[i]];
                                            unsigned short maxnum = *max_element(part.begin(), part.end());
                                            vector<double> temp_cov(maxnum, 0);
                                            if (i == 0)
                                            {
                                                if (i != var_site.size() - 1)
                                                {
                                                    var_distance = min((size_t)(var_site[i + 1] - var_site[i] - 1), u.size);
                                                }
                                                else
                                                {
                                                    var_distance = min(u.size, exit_uni.size);
                                                }
                                            }
                                            else if (i == var_site.size() - 1)
                                            {
                                                var_distance = min((size_t)(var_site[i] - var_site[i - 1] - 1), exit_uni.size);
                                            }
                                            else
                                            {
                                                var_distance = min(var_site[i] - var_site[i - 1] - 1, var_site[i + 1] - var_site[i] - 1);
                                            }
                                            if (find(indel_pos.begin(), indel_pos.end(), var_site[i]) != indel_pos.end())
                                            {
                                                ++indel;
                                                cov_info_all << 1 << "\t"  << indel_len[indel - 1] << "\t" << var_count << "\t" << var_site.size() << "\t" <<coefficient<< "\t"<< var_distance<<"\t" << endl;
                                            }
                                            else
                                            {
                                                cov_info_all << 1 << "\t" << "0\t" << var_count << "\t" << var_site.size() <<"\t"<<coefficient<< "\t" << var_distance <<"\t" << endl;
                                            }
                                            for(int s =0 ;s <cdbg.getNbColors();s++){
                                                vector<double> temp_cov(maxnum, 0.0);
                                                vector<double> res_cov;
                                                double sum = 0;
                                                for (int j = 0; j < part.size(); j++)
                                                {
                                                    temp_cov[part[j] - 1] += cov_vec[s][j];
                                                }
                                                for(auto c:temp_cov){
                                                    if(c>0.0){
                                                        res_cov.push_back(c);
                                                        sum+=c;
                                                    }
                                                }
                                                if(res_cov.size()<2){
                                                    continue;
                                                }                                    
                                                stringstream cov_info; 
                                                stringstream fre_info; 
                                                for (auto c : res_cov)
                                                {
                                                    cov_info << c << "\t";
                                                    fre_info << (double)c / sum << "\n";
                                                }
                                                cov_info<<s<<"\t"<<cov_info_all.str();
                                                switch (res_cov.size())
                                                {
                                                    case 2:
                                                        ++allele_info[0];
                                                        bifre_info << fre_info.str();
                                                        bicov_info << cov_info.str();
                                                        
                                                        
                                                        break;
                                                    case 3:
                                                        ++allele_info[1];
                                                        trifre_info << fre_info.str();
                                                        tricov_info << cov_info.str();
                                                        
                                                        
                                                        break;
                                                    case 4:
                                                        ++allele_info[2];
                                                        tetrafre_info << fre_info.str();
                                                        tetracov_info << cov_info.str();
                                                        break;
                                                    case 5:
                                                        ++allele_info[3];
                                                        pentafre_info << fre_info.str();
                                                        pentacov_info << cov_info.str();
                                                        
                                                }
                                               
                                            }
                                        }
                                        allele2.fetch_add(allele_info[0]);
                                        allele3.fetch_add(allele_info[1]);
                                        allele4.fetch_add(allele_info[2]);
                                        allele5.fetch_add(allele_info[3]);
                                        {
                                            lock_guard<mutex> lock(mtx_var);    
                                            s_var<<var_info_stream.str();
                                        }
                                        {
                                            lock_guard<mutex> lock(mtx_allfre);    
                                            allfre<<bifre_info.str()<<trifre_info.str()<<tetrafre_info.str();              
                                        }    
                                        {
                                            lock_guard<mutex> lock(mtx_bicov);    
                                            bicov<<bicov_info.str();                                                  
                                            bifre<<bifre_info.str();   
    
                                        }           
                                                    
                                        {
                                            lock_guard<mutex> lock(mtx_tricov);    
                                            tricov<<tricov_info.str(); 
                                            trifre<<trifre_info.str(); 
                                        }              
                                              
                                        {
                                            lock_guard<mutex> lock(mtx_tetracov);    
                                            tetracov<<tetracov_info.str();   
                                            tetrafre<<tetrafre_info.str();      
                                        }
                                        {
                                            lock_guard<mutex> lock(mtx_pentacov);    
                                            pentacov<<pentacov_info.str();                                     
                                            pentafre<<pentafre_info.str();  

                                        }
                                    }
                                }
                            }
                            else
                        {
                                exit_uni=*u.getSuccessors().begin();
                            while(exit_uni.getData()->getData(exit_uni)!=u.getData()->getData(u)->get_plus_or_minus(u.strand))
                            {
                                exit_uni=*exit_uni.getSuccessors().begin();                    
                            } 
                            if(u.referenceUnitigToString().compare(exit_uni.referenceUnitigToString())<0){
                                if(u.strand)
                                    u.getData()->getData(u)->set_plus_visited();
                                else
                                    u.getData()->getData(u)->set_minus_visited();
                                continue;
                            }
                            vector<string> str_vec;
                            stack<UnitigColorMap<MyUnitig>> major;
                            stack<UnitigColorMap<MyUnitig>> minor;
                            string bubbleStr = "";
                            
                            minor.push(u);
                            while (!minor.empty())
                            {
                                UnitigColorMap<MyUnitig> umi = minor.top();
                                minor.pop();
                                major.push(umi);
                                string str = umi.mappedSequenceToString();
                                bubbleStr += str.substr(0,umi.len);
                                if(umi.isSameReferenceUnitig(exit_uni))
                                {
                                    bubbleStr += str.substr(umi.len);
                                    str_vec.push_back(bubbleStr.substr(u.len - 1, bubbleStr.length()  - u.len + 1-umi.len+1));  
                                    
                                    
                                    
                                    bubbleStr = bubbleStr.substr(0, bubbleStr.length() - str.length());
                                    major.pop();
                                    while (!major.empty() && !minor.empty())
                                    {
                                        bool f = false;
                                        for (auto uma : major.top().getSuccessors())
                                        {
                                            if (uma == minor.top())
                                            {
                                                f = true;
                                                break;
                                            }
                                        }
                                        if (f == false)
                                        { 
                                            
                                            bubbleStr = bubbleStr.substr(0, bubbleStr.length() - major.top().len);
                                            major.pop();
                                        }
                                        else
                                        {
                                            break;
                                        }
                                    }
                                }
                                else
                                {
                                    for (auto u : umi.getSuccessors())
                                    {
                                        minor.push(u);
                                    }
                                }
                            }
                            
                            
                            
                            
                            
                            sortSeq_branching(str_vec, 0, str_vec.size() - 1);
                            SeqAlign seqalign(match,mismatch,gap);
                            
                            vector<uint> var_site;
                            vector<uint> snp_pos;
                            vector<uint> indel_pos;
                            vector<uint> indel_len;
                            vector<vector<unsigned short>> partition;
                            seqalign.SequenceAlignment(str_vec, snp_pos, indel_pos, partition, indel_len);
                            if (str_vec.size() != 0)
                            {
                                coreNum.fetch_add(1);
                                coreCov.fetch_add((size_t)core.first);                        
                                stringstream var_info_stream; 
                                size_t var_count=var_count_all.fetch_add(1);
                                for (auto s : str_vec)
                                {
                                    var_info_stream << var_count << "\t" <<0<<"\t"<<u.getData()->getData(u)->get_id()<<"\t"<<exit_uni.getData()->getData(exit_uni)->get_id()<< "\t" << s << endl;
                                }                                     
                                vector<size_t> allele_info(4,0);
                                stringstream bifre_info; 
                                stringstream bicov_info; 
                                stringstream tricov_info; 
                                stringstream trifre_info; 
                                stringstream tetracov_info; 
                                stringstream tetrafre_info; 
                                stringstream pentacov_info; 
                                stringstream pentafre_info;        
                                for (int i = 0; i < partition.size(); i++)
                                {
                                    if (partition[i].back() > 0)
                                    {
                                        var_site.push_back(i);
                                    }
                                }
                                uint indel = 0;
                                for (uint i = 0; i < var_site.size(); i++)
                                {
                            
                                    stringstream cov_info; 
                                    stringstream fre_info; 
                                    uint var_distance = 0;
                                    unsigned short maxnum = *max_element(partition[var_site[i]].begin(), partition[var_site[i]].end());
                                    vector<string> k_length_str(str_vec.size());
                                    vector<set<string>> string_set(maxnum); 
                                    if (i == 0)
                                    {
                                        if (i != var_site.size() - 1)
                                        {
                                            var_distance = min((size_t)(var_site[i + 1] - var_site[i] - 1), u.size);
                                        }
                                        else
                                        {
                                            var_distance = min(u.size, exit_uni.size);
                                        }
                                    }
                                    else if (i == var_site.size() - 1)
                                    {
                                        var_distance = min((size_t)(var_site[i] - var_site[i - 1] - 1), exit_uni.size);
                                    }
                                    else
                                    {
                                        var_distance = min(var_site[i] - var_site[i - 1] - 1, var_site[i + 1] - var_site[i] - 1);
                                    }
                                    if (find(indel_pos.begin(), indel_pos.end(), var_site[i]) != indel_pos.end())
                                    {
                                        vector<int > site_vec(str_vec.size(),var_site[i]);
                                        while(true){
                                            set<char> site_char;
                                            for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                            {
                                                string c=str_vec[s_size].substr(site_vec[s_size],1);
                                                while(c.compare("-")==0)
                                                {
                                                    site_vec[s_size]+=1;
                                                    c=str_vec[s_size].substr(site_vec[s_size],1);
                                                }
                                                site_vec[s_size]+=1;
                                                k_length_str[s_size]+=c;
                                                site_char.insert(c[0]);
                                            }
                                            if(site_char.size()>1){
                                                break;
                                            }
                                        }
                                        if(indel==0)
                                        {
                                            for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                            {
                                                int indel_i=k_length_str[s_size].length();
                                                k_length_str[s_size]=str_vec[s_size].substr(var_site[i]-cdbg.getK()+indel_i,cdbg.getK()-indel_i)+k_length_str[s_size];
                                            }
                                        }
                                        else 
                                        {
                                            for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                            {
                                                int indel_i=k_length_str[s_size].length();
                                                string temp_str;
                                                temp_str = str_vec[s_size].substr(0, var_site[i]);
                                                temp_str.erase(std::remove(temp_str.begin(), temp_str.end(), '-'), temp_str.end());
                                                if(temp_str.length()<cdbg.getK()-indel_i){
                                                    k_length_str[s_size] =temp_str+k_length_str[s_size];
                                                    for(int extend_site=site_vec[s_size];k_length_str[s_size].length()<cdbg.getK();extend_site++)
                                                    {
                                                        string c=str_vec[s_size].substr(extend_site,1);
                                                        if(c.compare("-")!=0)
                                                        {
                                                            k_length_str[s_size]+=c;
                                                        }
                                                    }
                                                }
                                                else
                                                    k_length_str[s_size] =temp_str.substr(temp_str.length()-cdbg.getK()+indel_i,cdbg.getK()-indel_i)+k_length_str[s_size];
                                            
                                            }
                                        }
                                        ++indel;
                                        
                                        
                                        for (int part_i = 0; part_i < partition[var_site[i]].size(); part_i++)
                                        {
                                            string_set[partition[var_site[i]][part_i] - 1].insert(k_length_str[part_i]);
                                        }
                                        vector< vector< double > > cov_vec(cdbg.getNbColors());
                                        for(int i=0;i<cov_vec.size();i++){
                                            cov_vec[i].resize(maxnum,0);
                                        }
                                        set<size_t> color_set;
                                        bool flag_site_cov = true;
                                        for (int string_i = 0; string_i < maxnum&& flag_site_cov; string_i++)
                                        {
                                            for (auto s : string_set[string_i])
                                            {
                                                UnitigColorMap<MyUnitig> path_uni=cdbg.findUnitig(s.c_str(),0,s.length());
                                                for(size_t ci=0;ci<cdbg.getNbColors();ci++){
                                                    if(path_uni.getData()->getUnitigColors(path_uni)->contains(path_uni,ci)){
                                                        pair<double, bool> cov_res = readCov(s, cutoff[ci].first, cutoff[ci].second,ci);
                                                        if (cov_res.second == false)
                                                        {                                                            
                                                            flag_site_cov = false;
                                                            break;
                                                        }                                                        
                                                        color_set.insert(ci);
                                                        cov_vec[ci][string_i]+=cov_res.first;
                                                    }
                                                }
                                                if(flag_site_cov==false){
                                                    break;
                                                }
                                            }
                                        }                                            
                                        if(color_set.size()!=cdbg.getNbColors()){
                                            continue;
                                        }
                                        if (flag_site_cov == false)
                                        {
                                            continue;
                                        }
                                        double coefficient=0;
                                
                                        for(size_t ci=0;ci<cov_vec.size()-1;ci++){                                        
                                            for(size_t cj=ci+1;cj<cov_vec.size();cj++){ 
                                                coefficient=max(coefficient,computeCramerVCoefficient(cov_vec[ci],cov_vec[cj]));
                                            }
                                        }   
                                        for(int ci=0;ci<cdbg.getNbColors();ci++){
                                            vector<double> temp_cov;      
                                            double sum=0;
                                            for(auto c:cov_vec[ci]){
                                                if(c>0.0){
                                                    sum+=c;
                                                    temp_cov.push_back(c);
                                                }
                                            }
                                            if(temp_cov.size()<2){
                                                continue;
                                            }
                                            stringstream cov_info; 
                                            stringstream fre_info;
                                            for (auto c : temp_cov)
                                            {
                                                cov_info << c << "\t";
                                                fre_info << (double)c / sum << "\n";
                                            }
                                            cov_info <<ci<<"\t"<< 0 << "\t" << indel_len[indel - 1] << "\t" << var_count << "\t" << var_site.size() << "\t" << coefficient<<"\t"<< var_distance << "\t" << endl;
                                 
                                            switch (temp_cov.size())
                                            {
                                                case 2:
                                                    ++allele_info[0];
                                                    bifre_info << fre_info.str();
                                                    bicov_info << cov_info.str();
                                                    
                                                    
                                                    break;
                                                case 3:
                                                    ++allele_info[1];
                                                    trifre_info << fre_info.str();
                                                    tricov_info << cov_info.str();
                                                    
                                                    
                                                    break;
                                                case 4:
                                                    ++allele_info[2];
                                                    tetrafre_info << fre_info.str();
                                                    tetracov_info << cov_info.str();
                                                    break;
                                                case 5:
                                                    ++allele_info[3];
                                                    pentafre_info << fre_info.str();
                                                    pentacov_info << cov_info.str(); 
                                                    
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (indel > 0)
                                        {
         
                                            for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                            {
                                              
                                                string temp_str;
                                                temp_str = str_vec[s_size].substr(0, var_site[i]+1);
                                                temp_str.erase(std::remove(temp_str.begin(), temp_str.end(), '-'), temp_str.end());                                    
                                                if(temp_str.length()<cdbg.getK()){
                                                    k_length_str[s_size] =temp_str;
                                                    for(int extend_site=var_site[i]+1;k_length_str[s_size].length()<cdbg.getK();extend_site++)
                                                    {
                                                        string c=str_vec[s_size].substr(extend_site,1);
                                                        if(c.compare("-")!=0)
                                                        {
                                                            k_length_str[s_size]+=c;
                                                        }
                                                    }
                                                }
                                                else
                                                    k_length_str[s_size] = (temp_str.substr(temp_str.length() - cdbg.getK(), cdbg.getK()));
                                            }
                                        }
                                        else
                                        {
                                            for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                            {
                                                k_length_str[s_size] = (str_vec[s_size].substr(var_site[i] - cdbg.getK() + 1, cdbg.getK()));
                                            }
                                        }
                                    
                                        
                                        for (int part_i = 0; part_i < partition[var_site[i]].size(); part_i++)
                                        {
                                            string_set[partition[var_site[i]][part_i] - 1].insert(k_length_str[part_i]);
                                        }
                                        vector< vector< double > > cov_vec(cdbg.getNbColors());
                                        for(int i=0;i<cov_vec.size();i++){
                                            cov_vec[i].resize(maxnum,0);
                                        }
                                        set<size_t> color_set;
                                        bool flag_site_cov = true;
                                        for (int string_i = 0; string_i < string_set.size() && flag_site_cov; string_i++)
                                        {
                                            for (auto s : string_set[string_i])
                                            {
                                                UnitigColorMap<MyUnitig> path_uni=cdbg.findUnitig(s.c_str(),0,s.length());
                                                for(size_t ci=0;ci<cdbg.getNbColors();ci++){
                                                    if(path_uni.getData()->getUnitigColors(path_uni)->contains(path_uni,ci)){
                                                        pair<double, bool> cov_res = readCov(s, cutoff[ci].first, cutoff[ci].second,ci);
                                                        if (cov_res.second == false)
                                                        {
                                                            flag_site_cov = false;
                                                            break;
                                                        }
                                                        cov_vec[ci][string_i]+=cov_res.first;     
                                                        color_set.insert(ci);
                                                    }
                                                }
                                                if(flag_site_cov==false){
                                                    break;
                                                }
                                            }
                                        }                                            
                                        if(color_set.size()!=cdbg.getNbColors()){
                                            continue;
                                        }
                                        if (flag_site_cov == false)
                                        {
                                            continue;
                                        }
                                        double coefficient=0;
                                
                                        for(size_t ci=0;ci<cov_vec.size()-1;ci++){                                        
                                            for(size_t cj=ci+1;cj<cov_vec.size();cj++){ 
                                                coefficient=max(coefficient,computeCramerVCoefficient(cov_vec[ci],cov_vec[cj]));
                                            }
                                        }                                           
                                        
                                        for(int ci=0;ci<cdbg.getNbColors();ci++){
                                            vector<double> temp_cov;      
                                            double sum=0;
                                            for(auto c:cov_vec[ci]){
                                                if(c>0.0){
                                                    sum+=c;
                                                    temp_cov.push_back(c);
                                                }
                                            }
                                            if(temp_cov.size()<2){
                                                continue;
                                            }
                                            stringstream cov_info; 
                                            stringstream fre_info;
                                            for (auto c : temp_cov)
                                            {
                                                cov_info << c << "\t";
                                                fre_info << (double)c / sum << "\n";
                                            }
                                            cov_info <<ci<<"\t"<< 0 << "\t" << "0\t" << var_count << "\t" << var_site.size() << "\t" <<coefficient<<"\t"<< var_distance << "\t" << endl;
                                           
                                            switch (temp_cov.size())
                                            {
                                                case 2:
                                                    ++allele_info[0];
                                                    bifre_info << fre_info.str();
                                                    bicov_info << cov_info.str();
                                                    
                                                    
                                                    break;
                                                case 3:
                                                    ++allele_info[1];
                                                    trifre_info << fre_info.str();
                                                    tricov_info << cov_info.str();
                                                    
                                                    
                                                    break;
                                                case 4:
                                                    ++allele_info[2];
                                                    tetrafre_info << fre_info.str();
                                                    tetracov_info << cov_info.str();
                                                    break;
                                                case 5:
                                                    ++allele_info[3];
                                                    pentafre_info << fre_info.str();
                                                    pentacov_info << cov_info.str(); 
                                                    
                                            }
                                        }
                                    }
                                }
                                    allele5.fetch_add(allele_info[3]);
                                    allele2.fetch_add(allele_info[0]);
                                    allele3.fetch_add(allele_info[1]);
                                    allele4.fetch_add(allele_info[2]);
                                    {
                                        lock_guard<mutex> lock(mtx_var);    
                                        s_var<<var_info_stream.str();
                                    }
                                    {
                                        lock_guard<mutex> lock(mtx_allfre);    
                                        allfre<<bifre_info.str()<<trifre_info.str()<<tetrafre_info.str();              
                                    }    
                                    {
                                        lock_guard<mutex> lock(mtx_bicov);    
                                        bicov<<bicov_info.str();        
                                                                        bifre<<bifre_info.str();   
          
                                    }           
                                    {
                                        lock_guard<mutex> lock(mtx_tricov);    
                                        tricov<<tricov_info.str();                                   
                                         trifre<<trifre_info.str();              
                                    }              
                                    {
                                        lock_guard<mutex> lock(mtx_tetracov);    
                                        tetracov<<tetracov_info.str();                                    
                                         tetrafre<<tetrafre_info.str();              
                                    }
                                   {
                                        lock_guard<mutex> lock(mtx_pentacov);    
                                        pentacov<<pentacov_info.str();                                     
                                        pentafre<<pentafre_info.str();              
                                    }
                            }
                        }
                        }
                        if(u.strand){
                            
                            u.getData()->getData(u)->set_plus_visited();
                            if(exit_uni.strand==true)
                            {
                                exit_uni.getData()->getData(exit_uni)->set_minus_visited();
                            }
                            else
                            {
                                exit_uni.getData()->getData(exit_uni)->set_plus_visited();
                            }
                        }else{
                            
                            u.getData()->getData(u)->set_minus_visited();
                            if(exit_uni.strand==true)
                            {
                                exit_uni.getData()->getData(exit_uni)->set_minus_visited();
                            }
                            else
                            {
                                exit_uni.getData()->getData(exit_uni)->set_plus_visited();
                            }
                        }
                    
                    }
                    {
                        lock_guard<std::mutex> lock(mtx_uni);
                        if (all_iter != cdbg.end())
                        {
                            iter = all_iter;
                            all_iter++;
                        }
                        else
                        {
                            return;
                        }
                    }
                } });
        }
        for (auto &w : workers)
            w.join();
        time_t end_time = time(NULL);
        cout << "CCDBG::PloidyEstimation(): Cpu time : "
             << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
        cout << "CCDBG::PloidyEstimation(): Real time : "
             << (double)difftime(end_time, start_time) << "s" << endl;
        cout << "CCDBG::PloidyEstimation(): Alleles in SuperBubbles  :\t"
             << "2 :" << allele2 << "\t"
             << "3 :" << allele3 << "\t"
             << "4 :" << allele4 << "\t"
             << "5 :" << allele5 << endl;

        bifre.close();
        trifre.close();
        tetrafre.close();
        bicov.close();
        tricov.close();
        tetracov.close();
        pentacov.close();
        pentafre.close();
        allfre.close();
        if (coreNum != 0)
        {
            const int avg = coreCov / coreNum;
            cout << "CCDBG::PloidyEstimation(): Sites' Average Coverage:" << avg << endl;
        }
    }
    else
    {
        ploidyEstimation_ptr(outpre, cutoff);
    }
}
void CCDBG::setNoBubble_multithread_ptr(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p, mutex &x)
{
    if (vec_km_seen.size() < 4)
    {
        return;
    }
    lock_guard<std::mutex> lock(x);
    if (p.second.getData()->getData(p.second)->get_plus_or_minus(!p.second.strand) == p.first.getData()->getData(p.first) && !p.first.getData()->getData(p.first)->is_non_super() && !p.second.getData()->getData(p.second)->is_non_super())
    {
        return;
    }
    if (p.second.getData()->getData(p.second)->get_plus_or_minus(!p.second.strand) != NULL || p.second.getData()->getData(p.second)->is_non_super())
    {
        for (const auto &ucm : vec_km_seen)
        {
            if (ucm == p.first)
            {
                if (p.first.strand == true)
                {
                    p.first.getData()->getData(p.first)->set_plus_self();
                }
                else
                {
                    p.first.getData()->getData(p.first)->set_minus_self();
                }
                continue;
            }
            if (ucm == p.second)
            {
                if (p.second.strand == true)
                {
                    p.second.getData()->getData(p.second)->set_minus_self();
                }
                else
                {
                    p.second.getData()->getData(p.second)->set_plus_self();
                }
                continue;
            }
            if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
            {
                MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
                if (ex->get_plus() == ucm.getData()->getData(ucm))
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->getData(ucm)->set_plus_self();
            if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
            {
                MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
                if (ex->get_plus() == ucm.getData()->getData(ucm))
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->getData(ucm)->set_minus_self();
            ucm.getData()->getData(ucm)->set_non_super();
        }
        return;
    }
    bool flag = true;
    if (vec_km_seen.size() <= 6)
    {
        for (const auto &ucm : vec_km_seen)
        {
            if (ucm == p.first || ucm == p.second)
            {
                continue;
            }
            if (ucm.getPredecessors().cardinality() == 1 && (ucm.getPredecessors().begin())->isSameReferenceUnitig(p.first) && ucm.getSuccessors().cardinality() == 1 && ucm.getSuccessors().begin()->isSameReferenceUnitig(p.second))
            {
                continue;
            }
            else
            {
                flag = false;
                break;
            }
        }
        if (flag == true)
        {
            p.first.getData()->getData(p.first)->setStrict(p.first.strand);
            p.second.getData()->getData(p.second)->setStrict(!p.second.strand);
        }
    }
    if (vec_km_seen.size() > complex_size)
    {
        p.first.getData()->getData(p.first)->setComplex(p.first.strand);
        p.second.getData()->getData(p.second)->setComplex(!p.second.strand);
    }
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm == p.first || ucm == p.second)
        {
            continue;
        }
        if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->getData(ucm)->set_plus_self();
        if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->getData(ucm)->set_minus_self();
        ucm.getData()->getData(ucm)->set_non_super();
    }
    bool f = true;
    if (p.first.getData()->getUnitigColors(p.first)->size(p.first) != p.first.len * cdbg.getNbColors())
    {
        f = false;
        p.first.getData()->getData(p.first)->set_non_super();
        if (p.first.strand == true)
        {
            p.first.getData()->getData(p.first)->set_plus_self();
        }
        else
        {
            p.first.getData()->getData(p.first)->set_minus_self();
        }
        if (p.second.strand == false)
        {
            p.second.getData()->getData(p.second)->set_plus_self();
        }
        else
        {
            p.second.getData()->getData(p.second)->set_minus_self();
        }
    }
    if (p.second.getData()->getUnitigColors(p.second)->size(p.first) != p.second.len * cdbg.getNbColors())
    {
        f = false;
        p.second.getData()->getData(p.second)->set_non_super();
        if (p.first.strand == true)
        {
            p.first.getData()->getData(p.first)->set_plus_self();
        }
        else
        {
            p.first.getData()->getData(p.first)->set_minus_self();
        }
        if (p.second.strand == false)
        {
            p.second.getData()->getData(p.second)->set_plus_self();
        }
        else
        {
            p.second.getData()->getData(p.second)->set_minus_self();
        }
    }
    if (f)
    {
        map<int, vector<int>> map;
        map[p.first.getData()->getData(p.first)->get_id()].resize(cdbg.getNbColors());
        map[p.second.getData()->getData(p.second)->get_id()].resize(cdbg.getNbColors());
        for (int i = 0; i < cdbg.getNbColors(); i++)
        {
            map[p.first.getData()->getData(p.first)->get_id()][i] = i;
            map[p.second.getData()->getData(p.second)->get_id()][i] = i;
        }
        for (const auto &ucm : vec_km_seen)
        {
            if (ucm == p.second)
            {
                continue;
            }
            if (map.find(ucm.getData()->getData(ucm)->get_id()) == map.end())
            {
                vector<int> color;
                for (int i = 0; i < cdbg.getNbColors(); i++)
                {
                    if (ucm.getData()->getUnitigColors(ucm)->contains(ucm, i))
                    {
                        color.push_back(i);
                    }
                }
                map[ucm.getData()->getData(ucm)->get_id()] = color;
            }
            set<int> suc_color;
            for (const auto &suc : ucm.getSuccessors())
            {
                for (int i = 0; i < map[ucm.getData()->getData(ucm)->get_id()].size(); i++)
                {
                    if (suc.getData()->getUnitigColors(suc)->contains(suc, map[ucm.getData()->getData(ucm)->get_id()][i]))
                    {
                        suc_color.insert(map[ucm.getData()->getData(ucm)->get_id()][i]);
                    }
                }
            }
            if (suc_color.size() == map[ucm.getData()->getData(ucm)->get_id()].size())
            {
                continue;
            }
            else
            {
                f = false;
                break;
            }
        }
        if (f)
        {
            if (p.first.strand == true)
            {
                p.first.getData()->getData(p.first)->set_plus(p.second.getData()->getData(p.second));
            }
            else
            {
                p.first.getData()->getData(p.first)->set_minus(p.second.getData()->getData(p.second));
            }
            if (p.second.strand == true)
            {
                p.second.getData()->getData(p.second)->set_minus(p.first.getData()->getData(p.first));
            }
            else
            {
                p.second.getData()->getData(p.second)->set_plus(p.first.getData()->getData(p.first));
            }
        }
        else
        {
            if (p.first.strand == true)
            {
                p.first.getData()->getData(p.first)->set_plus_self();
            }
            else
            {
                p.first.getData()->getData(p.first)->set_minus_self();
            }
            if (p.second.strand == false)
            {
                p.second.getData()->getData(p.second)->set_plus_self();
            }
            else
            {
                p.second.getData()->getData(p.second)->set_minus_self();
            }
        }
    }
}
void CCDBG::setNoBubble_multithread_ptr(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, mutex &x, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p)
{
    lock_guard<std::mutex> lock(x);
    if (p.first.strand == true)
    {
        if (p.first.getData()->getData(p.first)->get_plus() != NULL)
        {
            if (p.first.getData()->getData(p.first)->get_plus()->get_plus() == p.first.getData()->getData(p.first))
            {
                p.first.getData()->getData(p.first)->get_plus()->set_plus_self();
            }
            else
            {
                p.first.getData()->getData(p.first)->get_plus()->set_minus_self();
            }
        }
        p.first.getData()->getData(p.first)->set_plus_self();
    }
    else
    {
        if (p.first.getData()->getData(p.first)->get_minus() != NULL)
        {
            if (p.first.getData()->getData(p.first)->get_minus()->get_plus() == p.first.getData()->getData(p.first))
            {
                p.first.getData()->getData(p.first)->get_minus()->set_plus_self();
            }
            else
            {
                p.first.getData()->getData(p.first)->get_minus()->set_minus_self();
            }
        }
        p.first.getData()->getData(p.first)->set_minus_self();
    }
    if (p.second.strand == false)
    {
        if (p.second.getData()->getData(p.second)->get_plus() != NULL)
        {
            if (p.second.getData()->getData(p.second)->get_plus()->get_plus() == p.second.getData()->getData(p.second))
            {
                p.second.getData()->getData(p.second)->get_plus()->set_plus_self();
            }
            else
            {
                p.second.getData()->getData(p.second)->get_plus()->set_minus_self();
            }
        }
        p.second.getData()->getData(p.second)->set_plus_self();
    }
    else
    {
        if (p.second.getData()->getData(p.second)->get_minus() != NULL)
        {
            if (p.second.getData()->getData(p.second)->get_minus()->get_plus() == p.second.getData()->getData(p.second))
            {
                p.second.getData()->getData(p.second)->get_minus()->set_plus_self();
            }
            else
            {
                p.second.getData()->getData(p.second)->get_minus()->set_minus_self();
            }
        }
        p.second.getData()->getData(p.second)->set_minus_self();
    }
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm == p.first || ucm == p.second)
        {
            continue;
        }
        if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->getData(ucm)->set_plus_self();
        if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->getData(ucm)->set_minus_self();
        ucm.getData()->getData(ucm)->set_non_super();
    }
}
void CCDBG::extractSuperBubble_multithread_ptr(const UnitigColorMap<MyUnitig> &s, mutex &m)
{
    bool flag_cycle = false;
    bool flag_tip = false;
    bool flag_large_cycle = false;
    vector<UnitigColorMap<MyUnitig>> vertices_visit;
    pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> p;
    vector<UnitigColorMap<MyUnitig>> vec_km_seen;
    map<size_t, uint8_t> state_map;
    map<size_t, bool> strand_map;
    unordered_set<UnitigColorMap<MyUnitig>, UnitigMapHash<DataAccessor<MyUnitig>, DataStorage<MyUnitig>>> cycle_unitig_set;
    UnitigColorMap<MyUnitig> v(s);
    vertices_visit.push_back(v);
    vec_km_seen.push_back(v);
    while (!vertices_visit.empty())
    {
        v = vertices_visit.back();
        vertices_visit.pop_back();
        state_map[v.getData()->getData(v)->get_id()] = 0x01;
        strand_map[v.getData()->getData(v)->get_id()] = v.strand;
        if (!v.getSuccessors().hasSuccessors())
        {
            flag_tip = true;
        }
        else
        {
            for (auto &u : v.getSuccessors())
            {
                if (u == s)
                {
                    flag_cycle = true;
                    cycle_unitig_set.insert(s);
                    cycle_unitig_set.insert(v);
                    continue;
                }
                if (state_map.find(u.getData()->getData(u)->get_id()) == state_map.end())
                {
                    vec_km_seen.push_back(u);
                    strand_map[u.getData()->getData(u)->get_id()] = u.strand;
                    state_map[u.getData()->getData(u)->get_id()] = 0x02;
                    bool all_predecessor_visited = true;
                    for (const auto &predecessor : u.getPredecessors())
                    {
                        if (state_map.find(predecessor.getData()->getData(predecessor)->get_id()) != state_map.end())
                        {
                            if (state_map[predecessor.getData()->getData(predecessor)->get_id()] != 0x01)
                            {
                                all_predecessor_visited = false;
                            }
                            if (strand_map[predecessor.getData()->getData(predecessor)->get_id()] != predecessor.strand)
                            {
                                flag_cycle = true;
                                cycle_unitig_set.insert(u);
                                cycle_unitig_set.insert(predecessor);
                            }
                        }
                        else
                        {
                            all_predecessor_visited = false;
                        }
                    }
                    if (all_predecessor_visited)
                        vertices_visit.push_back(u);
                }
                else if (state_map[u.getData()->getData(u)->get_id()] != 0x01)
                {
                    if (strand_map[u.getData()->getData(u)->get_id()] != u.strand)
                    {
                        flag_cycle = true;
                        cycle_unitig_set.insert(u);
                        cycle_unitig_set.insert(v);
                    }
                    state_map[u.getData()->getData(u)->get_id()] = 0x02;
                    bool all_predecessor_visited = true;
                    for (const auto &predecessor : u.getPredecessors())
                    {
                        if (state_map.find(predecessor.getData()->getData(predecessor)->get_id()) != state_map.end())
                        {
                            if (state_map[predecessor.getData()->getData(predecessor)->get_id()] != 0x01)
                            {
                                all_predecessor_visited = false;
                            }
                            if (strand_map[predecessor.getData()->getData(predecessor)->get_id()] != predecessor.strand)
                            {
                                flag_cycle = true;
                                cycle_unitig_set.insert(u);
                                cycle_unitig_set.insert(predecessor);
                            }
                        }
                        else
                        {
                            all_predecessor_visited = false;
                        }
                    }
                    if (all_predecessor_visited)
                        vertices_visit.push_back(u);
                }
                else
                {
                    flag_cycle = true;
                    cycle_unitig_set.insert(v);
                    cycle_unitig_set.insert(u);
                }
            }
        }
        if (vertices_visit.size() == 1)
        {
            bool not_seen = true;
            for (const auto &cucm : vec_km_seen)
            {
                if (cucm != vertices_visit[0])
                {
                    if (state_map[cucm.getData()->getData(cucm)->get_id()] == 0x02)
                    {
                        not_seen = false;
                        break;
                    }
                }
            }
            if (not_seen)
            {
                p.first = s;
                p.second = vertices_visit[0];
                for (const auto &successor : vertices_visit[0].getSuccessors())
                {
                    if (successor == s)
                    {
                        setNoBubble_multithread_ptr_cycle(vec_km_seen, m, p);
                        return;
                    }
                }
                if (flag_cycle == true || flag_tip == true)
                {
                    setNoBubble_multithread_ptr(vec_km_seen, m, p);
                    return;
                }
                setNoBubble_multithread_ptr(vec_km_seen, p, m);
                return;
            }
        }
    }
    if (flag_cycle == true)
    {
        lock_guard<std::mutex> lock(m);
        for (const auto &ucm : cycle_unitig_set)
        {
            if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
            {
                MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
                if (ex->get_plus() == ucm.getData()->getData(ucm))
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->getData(ucm)->set_plus_self();
            if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
            {
                MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
                if (ex->get_plus() == ucm.getData()->getData(ucm))
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->getData(ucm)->set_minus_self();
            ucm.getData()->getData(ucm)->set_non_super();
        }
        if (s.strand == true)
        {
            s.getData()->getData(s)->set_plus_self();
        }
        else
        {
            s.getData()->getData(s)->set_minus_self();
        }
    }
    return;
}

void CCDBG::setUnitigId(const string &outpre, const string &graphfile, const size_t &thr)
{
    if (access("PloidyFrost_output", 0))
    {
        system("mkdir ./PloidyFrost_output");
    }
    ofstream output("PloidyFrost_output/" + outpre + "_Unitig_Id.txt", ofstream::out);
    cout << "CCDBG::setUnitigId(): Setting Unitig Id" << endl;
    double start_time = time(NULL);
    clock_t start_clock = clock();
    size_t id = 0;
    for (const auto &unitig : cdbg)
    {
        (unitig.getData()->getData(unitig))->set_id(++id);
        output << id << "\t" << unitig.referenceUnitigToString() << endl;
    }
    output.close();
    time_t end_time = time(NULL);
    cout << "CCDBG::setUnitigId(): Cpu time : "
         << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
    cout << "CCDBG::setUnitigId(): Real time : "
         << (double)difftime(end_time, start_time) << "s" << endl;
}
void CCDBG::printInfo(const bool &verbose, const string &outpre)
{
    if (verbose)
    {
        cout << ">>>>>>>>>Bifrost Graph Information>>>>>>>>>" << endl;
        cout << "k:" << cdbg.getK() << "\t";
        cout << "g:" << cdbg.getG() << "\t";
        cout << "NbColors:" << cdbg.getNbColors() << "\t";
        cout << "nbKmer:" << cdbg.nbKmers() << "\t";
        cout << "nbUnitig:" << cdbg.size() << "\t";
        cout << "length:" << cdbg.length() << "\n";
        vector<string> colors = cdbg.getColorNames();
        for (string &col : colors)
        {
            cout << col << endl;
        }
        cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    }
    ofstream output(outpre + "_graph_info.txt", ofstream::out);
    output << "k:" << cdbg.getK() << "\t";
    output << "g:" << cdbg.getG() << "\t";
    output << "NbColors:" << cdbg.getNbColors() << "\t";
    output << "nbKmer:" << cdbg.nbKmers() << "\t";
    output << "nbUnitig:" << cdbg.size() << "\t";
    output << "length:" << cdbg.length() << "\n";
    vector<string> colors = cdbg.getColorNames();
    for (string &col : colors)
    {
        output << col << endl;
    }
    output.close();
}
void CCDBG::clearMarking(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen)
{
    for (const auto &ucm : vec_km_seen)
    {
        MyUnitig *da_ucm = ucm.getData()->getData(ucm);
        da_ucm->clear(ucm);
    }
}
void CCDBG::findSuperBubble_ptr(const string &outpre)
{
    cout << "CCDBG::findSuperBubble(): Finding superbubbles" << endl;
    clock_t start_clock = clock();
    double start_time = time(NULL);
    size_t nb_super_bubble = 0;
    ofstream outfile;
    if (access("PloidyFrost_output", 0))
    {
        system("mkdir ./PloidyFrost_output");
    }
    outfile.open("PloidyFrost_output/" + outpre + "_super_bubble.txt", ios::out | ios::trunc);
    if (outfile.fail())
    {
        cout << "CCDBG:: Open super_bubble file error";
        exit(EXIT_FAILURE);
    }
    size_t nb_unitig_processed = 0;
    cout << "CCDBG::findSuperBubble(): There are " << cdbg.size() << " unitigs " << endl;
    for (const auto &unitig : cdbg)
    {
        ++nb_unitig_processed;
        if (nb_unitig_processed % 100000 == 0)
        {
            cout << "CCDBG::findSuperBubble(): Processed " << nb_unitig_processed << " unitigs " << endl;
        }
        UnitigColorMap<MyUnitig> u(unitig);
        u.strand = true;
        if (u.getSuccessors().cardinality() > 1 && u.getData()->getData(u)->get_plus() == NULL)
        {
            extractSuperBubble_ptr(u);
        }
        u.strand = false;
        if (u.getSuccessors().cardinality() > 1 && u.getData()->getData(u)->get_minus() == NULL)
        {
            extractSuperBubble_ptr(u);
        }
    }
    time_t end_time = time(NULL);
    cout << "CCDBG::findSuperBubble(): Finding superbubbles Cpu time : "
         << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
    cout << "CCDBG::findSuperBubble(): Finding superbubbles Real time : "
         << (double)difftime(end_time, start_time) << "s" << endl;
    outfile << "BubbleId\tEntrance\tStrand\tExit\tisSimple\tisComplex" << endl;
    for (const auto &unitig : cdbg)
    {
        if (!unitig.getData()->getData(unitig)->is_super())
        {
            continue;
        }
        if (unitig.getData()->getData(unitig)->get_plus() != NULL)
        {
            ++nb_super_bubble;
            outfile << nb_super_bubble << "\t" << unitig.getData()->getData(unitig)->get_id() << "\t"
                    << "+"
                    << "\t" << unitig.getData()->getData(unitig)->get_plus()->get_id()
                    << "\t" << (unitig.getData()->getData(unitig)->isStrict(true) ? "1" : "0")
                    << "\t" << (unitig.getData()->getData(unitig)->isComplex(true) ? "1" : "0") //<< "\t" << (unitig.getData()->getData(unitig)->iscut() ? "1" : "0")
                    << endl;
        }
        if (unitig.getData()->getData(unitig)->get_minus() != NULL)
        {
            ++nb_super_bubble;
            outfile << nb_super_bubble << "\t" << unitig.getData()->getData(unitig)->get_id() << "\t"
                    << "-"
                    << "\t" << unitig.getData()->getData(unitig)->get_minus()->get_id()
                    << "\t" << (unitig.getData()->getData(unitig)->isStrict(false) ? "1" : "0")
                    << "\t" << (unitig.getData()->getData(unitig)->isComplex(false) ? "1" : "0") //<< "\t" << (unitig.getData()->getData(unitig)->iscut() ? "1" : "0")
                    << endl;
        }
    }
    outfile.close();
    cout << "CCDBG::findSuperBubble(): " << nb_super_bubble << "  SuperBubbles Found" << endl;
}
void CCDBG::extractSuperBubble_ptr(const UnitigColorMap<MyUnitig> &s)
{
    bool flag_cycle = false;
    bool flag_tip = false;
    bool flag_large_cycle = false;
    vector<UnitigColorMap<MyUnitig>> vertices_visit;
    pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> p;
    vector<UnitigColorMap<MyUnitig>> vec_km_seen;
    map<size_t, uint8_t> state_map;
    map<size_t, bool> strand_map;
    unordered_set<UnitigColorMap<MyUnitig>, UnitigMapHash<DataAccessor<MyUnitig>, DataStorage<MyUnitig>>> cycle_unitig_set;
    UnitigColorMap<MyUnitig> v(s);
    vertices_visit.push_back(v);
    vec_km_seen.push_back(v);
    while (!vertices_visit.empty())
    {
        v = vertices_visit.back();
        vertices_visit.pop_back();
        state_map[v.getData()->getData(v)->get_id()] = 0x01;
        strand_map[v.getData()->getData(v)->get_id()] = v.strand;
        if (!v.getSuccessors().hasSuccessors())
        {
            flag_tip = true;
        }
        else
        {
            for (auto &u : v.getSuccessors())
            {
                if (u == s)
                {
                    flag_cycle = true;
                    cycle_unitig_set.insert(s);
                    cycle_unitig_set.insert(v);
                    continue;
                }
                if (state_map.find(u.getData()->getData(u)->get_id()) == state_map.end() || state_map[u.getData()->getData(u)->get_id()] != 0x01)
                {
                    if (state_map.find(u.getData()->getData(u)->get_id()) == state_map.end())
                    {
                        vec_km_seen.push_back(u);
                        strand_map[u.getData()->getData(u)->get_id()] = u.strand;
                    }
                    else
                    {
                        if (strand_map[u.getData()->getData(u)->get_id()] != u.strand)
                        {
                            flag_cycle = true;
                            cycle_unitig_set.insert(u);
                            cycle_unitig_set.insert(v);
                        }
                    }
                    state_map[u.getData()->getData(u)->get_id()] = 0x02;
                    bool all_predecessor_visited = true;
                    for (const auto &predecessor : u.getPredecessors())
                    {
                        if (state_map.find(predecessor.getData()->getData(predecessor)->get_id()) != state_map.end())
                        {
                            if (state_map[predecessor.getData()->getData(predecessor)->get_id()] != 0x01)
                            {
                                all_predecessor_visited = false;
                            }
                            if (strand_map[predecessor.getData()->getData(predecessor)->get_id()] != predecessor.strand)
                            {
                                flag_cycle = true;
                                cycle_unitig_set.insert(u);
                                cycle_unitig_set.insert(predecessor);
                            }
                        }
                        else
                        {
                            all_predecessor_visited = false;
                        }
                    }
                    if (all_predecessor_visited)
                        vertices_visit.push_back(u);
                }
                else
                {
                    flag_cycle = true;
                    cycle_unitig_set.insert(v);
                    cycle_unitig_set.insert(u);
                }
            }
        }
        if (vertices_visit.size() == 1)
        {
            bool not_seen = true;
            for (const auto &cucm : vec_km_seen)
            {
                if (cucm != vertices_visit[0])
                {
                    if (state_map[cucm.getData()->getData(cucm)->get_id()] == 0x02)
                    {
                        not_seen = false;
                        break;
                    }
                }
            }
            if (not_seen)
            {
                p.first = s;
                p.second = vertices_visit[0];
                for (const auto &successor : vertices_visit[0].getSuccessors())
                {
                    if (successor == s)
                    {
                        setNoBubble_ptr_cycle(vec_km_seen, p);
                        return;
                    }
                }
                if (flag_cycle == true || flag_tip == true)
                {
                    setNoBubble_ptr(vec_km_seen, p);
                    return;
                }
                setNoBubble_ptr(p, vec_km_seen);
                return;
            }
        }
    }
    if (flag_cycle == true)
    {
        for (const auto &ucm : cycle_unitig_set)
        {
            if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
            {
                MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
                if (ex->get_plus() == ucm.getData()->getData(ucm))
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->getData(ucm)->set_plus_self();
            if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
            {
                MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
                if (ex->get_plus() == ucm.getData()->getData(ucm))
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->getData(ucm)->set_minus_self();
            ucm.getData()->getData(ucm)->set_non_super();
        }
        if (s.strand == true)
        {
            s.getData()->getData(s)->set_plus_self();
        }
        else
        {
            s.getData()->getData(s)->set_minus_self();
        }
    }
    return;
}
void CCDBG::setNoBubble_multithread_ptr_cycle(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, mutex &x, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p)
{
    lock_guard<std::mutex> lock(x);
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
            ucm.getData()->getData(ucm)->set_plus_self();
        }
        if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
            ucm.getData()->getData(ucm)->set_minus_self();
        }
        ucm.getData()->getData(ucm)->set_non_super();
    }
    if (p.first.strand == true)
    {
        p.first.getData()->getData(p.first)->set_plus_self();
    }
    else
    {
        p.first.getData()->getData(p.first)->set_minus_self();
    }
    if (p.second.strand == false)
    {
        p.second.getData()->getData(p.second)->set_plus_self();
    }
    else
    {
        p.second.getData()->getData(p.second)->set_minus_self();
    }
}
void CCDBG::setNoBubble_ptr_cycle(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p)
{
    for (const auto &ucm : vec_km_seen)
    {
        // if (ucm == p.first || ucm == p.second)
        // {
        //     continue;
        // }
        if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
            ucm.getData()->getData(ucm)->set_plus_self();
        }
        if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
            ucm.getData()->getData(ucm)->set_minus_self();
        }
        ucm.getData()->getData(ucm)->set_non_super();
    }
    if (p.first.strand == true)
    {
        p.first.getData()->getData(p.first)->set_plus_self();
    }
    else
    {
        p.first.getData()->getData(p.first)->set_minus_self();
    }
    if (p.second.strand == false)
    {
        p.second.getData()->getData(p.second)->set_plus_self();
    }
    else
    {
        p.second.getData()->getData(p.second)->set_minus_self();
    }
}
void CCDBG::setNoBubble_ptr(pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p, const vector<UnitigColorMap<MyUnitig>> &vec_km_seen)
{
    if (vec_km_seen.size() < 4)
    {
        return;
    }
    if (p.second.getData()->getData(p.second)->is_non_super() || p.first.getData()->getData(p.first)->is_non_super())
    {
        for (const auto &ucm : vec_km_seen)
        {
            if (ucm == p.first)
            {
                if (p.first.strand == true)
                {
                    p.first.getData()->getData(p.first)->set_plus_self();
                }
                else
                {
                    p.first.getData()->getData(p.first)->set_minus_self();
                }
                continue;
            }
            if (ucm == p.second)
            {
                if (p.second.strand == true)
                {
                    p.second.getData()->getData(p.second)->set_minus_self();
                }
                else
                {
                    p.second.getData()->getData(p.second)->set_plus_self();
                }
                continue;
            }
            if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
            {
                MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
                if (ex->get_plus() == ucm.getData()->getData(ucm))
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->getData(ucm)->set_plus_self();
            if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
            {
                MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
                if (ex->get_plus() == ucm.getData()->getData(ucm))
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->getData(ucm)->set_minus_self();
            ucm.getData()->getData(ucm)->set_non_super();
        }
        return;
    }
    bool flag = true;
    if (vec_km_seen.size() <= 6)
    {
        for (const auto &ucm : vec_km_seen)
        {
            if (ucm == p.first || ucm == p.second)
            {
                continue;
            }
            if (ucm.getPredecessors().cardinality() == 1 && (ucm.getPredecessors().begin())->isSameReferenceUnitig(p.first) && ucm.getSuccessors().cardinality() == 1 && ucm.getSuccessors().begin()->isSameReferenceUnitig(p.second))
            {
                continue;
            }
            else
            {
                flag = false;
                break;
            }
        }
        if (flag == true)
        {
            p.first.getData()->getData(p.first)->setStrict(p.first.strand);
            p.second.getData()->getData(p.second)->setStrict(!p.second.strand);
        }
    }
    if (vec_km_seen.size() > complex_size)
    {
        p.first.getData()->getData(p.first)->setComplex(p.first.strand);
        p.second.getData()->getData(p.second)->setComplex(!p.second.strand);
    }
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm == p.first || ucm == p.second)
        {
            continue;
        }
        if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->getData(ucm)->set_plus_self();
        if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->getData(ucm)->set_minus_self();
        ucm.getData()->getData(ucm)->set_non_super();
    }
    bool f = true;
    if (p.first.getData()->getUnitigColors(p.first)->size(p.first) != p.first.len * cdbg.getNbColors())
    {
        f = false;
        p.first.getData()->getData(p.first)->set_non_super();
        if (p.first.strand == true)
        {
            p.first.getData()->getData(p.first)->set_plus_self();
        }
        else
        {
            p.first.getData()->getData(p.first)->set_minus_self();
        }
        if (p.second.strand == false)
        {
            p.second.getData()->getData(p.second)->set_plus_self();
        }
        else
        {
            p.second.getData()->getData(p.second)->set_minus_self();
        }
    }
    if (p.second.getData()->getUnitigColors(p.second)->size(p.first) != p.second.len * cdbg.getNbColors())
    {
        f = false;
        p.second.getData()->getData(p.second)->set_non_super();
        if (p.first.strand == true)
        {
            p.first.getData()->getData(p.first)->set_plus_self();
        }
        else
        {
            p.first.getData()->getData(p.first)->set_minus_self();
        }
        if (p.second.strand == false)
        {
            p.second.getData()->getData(p.second)->set_plus_self();
        }
        else
        {
            p.second.getData()->getData(p.second)->set_minus_self();
        }
    }
    if (f)
    {
        map<int, vector<int>> map;
        map[p.first.getData()->getData(p.first)->get_id()].resize(cdbg.getNbColors());
        map[p.second.getData()->getData(p.second)->get_id()].resize(cdbg.getNbColors());
        for (int i = 0; i < cdbg.getNbColors(); i++)
        {
            map[p.first.getData()->getData(p.first)->get_id()][i] = i;
            map[p.second.getData()->getData(p.second)->get_id()][i] = i;
        }
        for (const auto &ucm : vec_km_seen)
        {
            if (ucm == p.second)
            {
                continue;
            }
            if (map.find(ucm.getData()->getData(ucm)->get_id()) == map.end())
            {
                vector<int> color;
                for (int i = 0; i < cdbg.getNbColors(); i++)
                {
                    if (ucm.getData()->getUnitigColors(ucm)->contains(ucm, i))
                    {
                        color.push_back(i);
                    }
                }
                map[ucm.getData()->getData(ucm)->get_id()] = color;
            }
            set<int> suc_color;
            for (const auto &suc : ucm.getSuccessors())
            {
                for (int i = 0; i < map[ucm.getData()->getData(ucm)->get_id()].size(); i++)
                {
                    if (suc.getData()->getUnitigColors(suc)->contains(suc, map[ucm.getData()->getData(ucm)->get_id()][i]))
                    {
                        suc_color.insert(map[ucm.getData()->getData(ucm)->get_id()][i]);
                    }
                }
            }
            if (suc_color.size() == map[ucm.getData()->getData(ucm)->get_id()].size())
            {
                continue;
            }
            else
            {
                f = false;
                break;
            }
        }
        if (f)
        {
            if (p.first.strand == true)
            {
                p.first.getData()->getData(p.first)->set_plus(p.second.getData()->getData(p.second));
            }
            else
            {
                p.first.getData()->getData(p.first)->set_minus(p.second.getData()->getData(p.second));
            }
            if (p.second.strand == true)
            {
                p.second.getData()->getData(p.second)->set_minus(p.first.getData()->getData(p.first));
            }
            else
            {
                p.second.getData()->getData(p.second)->set_plus(p.first.getData()->getData(p.first));
            }
        }
        else
        {
            if (p.first.strand == true)
            {
                p.first.getData()->getData(p.first)->set_plus_self();
            }
            else
            {
                p.first.getData()->getData(p.first)->set_minus_self();
            }
            if (p.second.strand == false)
            {
                p.second.getData()->getData(p.second)->set_plus_self();
            }
            else
            {
                p.second.getData()->getData(p.second)->set_minus_self();
            }
        }
    }
}
void CCDBG::setNoBubble_ptr(const vector<UnitigColorMap<MyUnitig>> &vec_km_seen, pair<UnitigColorMap<MyUnitig>, UnitigColorMap<MyUnitig>> &p)
{
    if (p.first.strand == true)
    {
        if (p.first.getData()->getData(p.first)->get_plus() != NULL)
        {
            if (p.first.getData()->getData(p.first)->get_plus()->get_plus() == p.first.getData()->getData(p.first))
            {
                p.first.getData()->getData(p.first)->get_plus()->set_plus_self();
            }
            else
            {
                p.first.getData()->getData(p.first)->get_plus()->set_minus_self();
            }
        }
        p.first.getData()->getData(p.first)->set_plus_self();
    }
    else
    {
        if (p.first.getData()->getData(p.first)->get_minus() != NULL)
        {
            if (p.first.getData()->getData(p.first)->get_minus()->get_plus() == p.first.getData()->getData(p.first))
            {
                p.first.getData()->getData(p.first)->get_minus()->set_plus_self();
            }
            else
            {
                p.first.getData()->getData(p.first)->get_minus()->set_minus_self();
            }
        }
        p.first.getData()->getData(p.first)->set_minus_self();
    }
    if (p.second.strand == false)
    {
        if (p.second.getData()->getData(p.second)->get_plus() != NULL)
        {
            if (p.second.getData()->getData(p.second)->get_plus()->get_plus() == p.second.getData()->getData(p.second))
            {
                p.second.getData()->getData(p.second)->get_plus()->set_plus_self();
            }
            else
            {
                p.second.getData()->getData(p.second)->get_plus()->set_minus_self();
            }
        }
        p.second.getData()->getData(p.second)->set_plus_self();
    }
    else
    {
        if (p.second.getData()->getData(p.second)->get_minus() != NULL)
        {
            if (p.second.getData()->getData(p.second)->get_minus()->get_plus() == p.second.getData()->getData(p.second))
            {
                p.second.getData()->getData(p.second)->get_minus()->set_plus_self();
            }
            else
            {
                p.second.getData()->getData(p.second)->get_minus()->set_minus_self();
            }
        }
        p.second.getData()->getData(p.second)->set_minus_self();
    }
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm == p.first || ucm == p.second)
        {
            continue;
        }
        if (ucm.getData()->getData(ucm)->get_plus() != NULL && ucm.getData()->getData(ucm)->get_plus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_plus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->getData(ucm)->set_plus_self();
        if (ucm.getData()->getData(ucm)->get_minus() != NULL && ucm.getData()->getData(ucm)->get_minus() != ucm.getData()->getData(ucm))
        {
            MyUnitig *ex = ucm.getData()->getData(ucm)->get_minus();
            if (ex->get_plus() == ucm.getData()->getData(ucm))
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->getData(ucm)->set_minus_self();
        ucm.getData()->getData(ucm)->set_non_super();
    }
}
void CCDBG::ploidyEstimation_ptr(const string &outpre, const vector<pair<int, int>> &cutoff)
{
    clock_t start_clock = clock();
    double start_time = time(NULL);
    cout << "CCDBG::PloidyEstimation():  Analyzing superbubbles to generate sites' information" << endl;
    size_t coreNum = 0;
    size_t coreCov = 0;
    vector<size_t> allele(4, 0);
    if (access("PloidyFrost_output", 0))
    {
        system("mkdir ./PloidyFrost_output");
    }
    ofstream bifre;
    ofstream trifre;
    ofstream tetrafre;
    ofstream bicov;
    ofstream tricov;
    ofstream tetracov;
    ofstream allfre;
    ofstream pentafre;
    ofstream pentacov;
    allfre.open("PloidyFrost_output/" + outpre + "_allele_frequency.txt", ios::out | ios::trunc);
    bifre.open("PloidyFrost_output/" + outpre + "_bifre.txt", ios::out | ios::trunc);
    trifre.open("PloidyFrost_output/" + outpre + "_trifre.txt", ios::out | ios::trunc);
    tetrafre.open("PloidyFrost_output/" + outpre + "_tetrafre.txt", ios::out | ios::trunc);
    bicov.open("PloidyFrost_output/" + outpre + "_bicov.txt", ios::out | ios::trunc);
    tricov.open("PloidyFrost_output/" + outpre + "_tricov.txt", ios::out | ios::trunc);
    pentafre.open("PloidyFrost_output/" + outpre + "_pentafre.txt", ios::out | ios::trunc);
    pentacov.open("PloidyFrost_output/" + outpre + "_pentacov.txt", ios::out | ios::trunc);
    tetracov.open("PloidyFrost_output/" + outpre + "_tetracov.txt", ios::out | ios::trunc);
    ofstream s_var;
    s_var.open("PloidyFrost_output/" + outpre + "_alignseq.txt", ios::out | ios::trunc);
    if (!allfre.is_open() || !bifre.is_open() || !trifre.is_open() ||
        !tetrafre.is_open() || !bicov.is_open() || !tricov.is_open() || !tetracov.is_open() || !s_var.is_open() ||
        !pentafre.is_open() || !pentacov.is_open())
    {
        cout << "CCDBG:: PloidyEstimation():Open file error" << endl;
        exit(EXIT_FAILURE);
    }
    size_t var_count = 0;
    size_t nb_unitig_processed = 0;

    for (auto u : cdbg)
    {
        ++nb_unitig_processed;
        if (nb_unitig_processed % 100000 == 0)
        {
            cout << "CCDBG::PloidyEstimation(): Processed " << nb_unitig_processed << " unitigs " << endl;
        }
        if (u.getData()->getData(u)->is_both_visited())
        {
            continue;
        }
        while (!u.getData()->getData(u)->is_both_visited())
        {
            if (u.getData()->getData(u)->is_plus_visited() != true)
            {
                u.strand = true;
                if (u.getData()->getData(u)->isComplex(true))
                {
                    u.getData()->getData(u)->set_plus_visited();
                    continue;
                }
            }
            else if (u.getData()->getData(u)->is_minus_visited() != true)
            {
                u.strand = false;
                if (u.getData()->getData(u)->isComplex(false))
                {
                    u.getData()->getData(u)->set_minus_visited();
                    break;
                }
            }
            else
            {
                break;
            }
            bool isStrict = u.getData()->getData(u)->isStrict(u.strand);
            UnitigColorMap<MyUnitig> exit_uni = NULL;
            pair<double, uint> core = pair<double, uint>(0, 0);
            bool flag = true;
            for (size_t i = 0; i < cdbg.getNbColors(); i++)
            {
                pair<double, bool> temp = readCovUni(u, cutoff[i].first, cutoff[i].second, i);
                if (temp.second == true)
                {
                    core.first += temp.first;
                    core.second += temp.second;
                }
                else
                {
                    flag == false;
                    break;
                }
            }
            if (flag)
            {
                if (isStrict)
                {
                    exit_uni = (*u.getSuccessors().begin()->getSuccessors().begin());
                    if (u.referenceUnitigToString().compare(exit_uni.referenceUnitigToString()) < 0)
                    {
                        if (u.strand)
                            u.getData()->getData(u)->set_plus_visited();
                        else
                            u.getData()->getData(u)->set_minus_visited();
                        continue;
                    }
                    vector<vector<double>> cov_vec(cdbg.getNbColors());
                    for (int i = 0; i < cov_vec.size(); i++)
                    {
                        cov_vec[i].resize(u.getSuccessors().cardinality(), 0);
                    }
                    vector<size_t> path_color;
                    uint8_t path = 0;
                    vector<UnitigColorMap<MyUnitig>> unitig_vec;
                    for (const auto &uu : u.getSuccessors())
                    {
                        unitig_vec.push_back(uu);
                        size_t j = 0;
                        for (size_t i = 0; i < cdbg.getNbColors(); i++)
                        {
                            if (uu.getData()->getUnitigColors(uu)->contains(uu, i) == true)
                            {
                                j++;
                                pair<double, bool> inside = readCovUni(uu, cutoff[i].first, cutoff[i].second, i);
                                if (inside.second == true)
                                {
                                    cov_vec[i][path] = inside.first;
                                }
                                else
                                {
                                    flag = false;
                                    break;
                                }
                            }
                        }
                        if (flag == false)
                        {
                            break;
                        }
                        if (uu.getData()->getUnitigColors(uu)->size(uu) != j * uu.len)
                        {
                            flag = false;
                            break;
                        }
                        ++path;
                        path_color.push_back(j);
                    }
                    if (flag == true)
                    {
                        flag = false;
                        for (auto v : cov_vec)
                        {
                            int c = 0;
                            for (auto d : v)
                            {
                                if (d != 0.0)
                                {
                                    ++c;
                                }
                            }
                            if (c > 1)
                            {
                                flag = true;
                                break;
                            }
                        }
                    }
                    if (flag == true)
                    {

                        sortSeq_simple(path_color, unitig_vec, cov_vec, 0, path_color.size() - 1);

                        SeqAlign seqalign(match, mismatch, gap);
                        vector<uint> var_site;
                        vector<uint> snp_pos;
                        vector<uint> indel_pos;
                        vector<uint> indel_len;
                        vector<vector<unsigned short>> partition;
                        vector<string> str_vec;
                        for (auto u : unitig_vec)
                        {
                            str_vec.push_back(u.mappedSequenceToString());
                        }
                        seqalign.SequenceAlignment(str_vec, snp_pos, indel_pos, partition, indel_len);
                        if (str_vec.size() != 0)
                        {
                            ++var_count;
                            for (auto s : str_vec)
                            {
                                s_var << var_count << "\t" << 1 << "\t" << u.getData()->getData(u)->get_id() << "\t" << exit_uni.getData()->getData(exit_uni)->get_id() << "\t" << s << endl;
                            }
                            coreCov += (size_t)core.first;
                            coreNum++;
                            for (int i = 0; i < partition.size(); i++)
                            {
                                if (partition[i].back() > 0)
                                {
                                    var_site.push_back(i);
                                }
                            }
                            uint var_distance = 0;
                            uint indel = 0;
                            double coefficient = 0;
                            for (size_t ci = 0; ci < cov_vec.size() - 1; ci++)
                            {
                                for (size_t cj = ci + 1; cj < cov_vec.size(); cj++)
                                {
                                    coefficient = max(coefficient, computeCramerVCoefficient(cov_vec[ci], cov_vec[cj]));
                                }
                            }
                            for (uint i = 0; i < var_site.size(); i++)
                            {
                                stringstream cov_info_all;
                                vector<unsigned short> part = partition[var_site[i]];
                                unsigned short maxnum = *max_element(part.begin(), part.end());
                                vector<double> temp_cov(maxnum, 0);
                                if (i == 0)
                                {
                                    if (i != var_site.size() - 1)
                                    {
                                        var_distance = min((size_t)(var_site[i + 1] - var_site[i] - 1), u.size);
                                    }
                                    else
                                    {
                                        var_distance = min(u.size, exit_uni.size);
                                    }
                                }
                                else if (i == var_site.size() - 1)
                                {
                                    var_distance = min((size_t)(var_site[i] - var_site[i - 1] - 1), exit_uni.size);
                                }
                                else
                                {
                                    var_distance = min(var_site[i] - var_site[i - 1] - 1, var_site[i + 1] - var_site[i] - 1);
                                }
                                if (find(indel_pos.begin(), indel_pos.end(), var_site[i]) != indel_pos.end())
                                {
                                    ++indel;
                                    cov_info_all << 1 << "\t" << indel_len[indel - 1] << "\t" << var_count << "\t" << var_site.size() << "\t" << coefficient << "\t" << var_distance << "\t" << endl;
                                }
                                else
                                {
                                    cov_info_all << 1 << "\t"
                                                 << "0\t" << var_count << "\t" << var_site.size() << "\t" << coefficient << "\t" << var_distance << "\t" << endl;
                                }
                                for (int s = 0; s < cdbg.getNbColors(); s++)
                                {
                                    vector<double> temp_cov(maxnum, 0.0);
                                    vector<double> res_cov;
                                    double sum = 0;
                                    for (int j = 0; j < part.size(); j++)
                                    {
                                        temp_cov[part[j] - 1] += cov_vec[s][j];
                                    }
                                    for (auto c : temp_cov)
                                    {
                                        if (c > 0.0)
                                        {
                                            res_cov.push_back(c);
                                            sum += c;
                                        }
                                    }
                                    if (res_cov.size() < 2)
                                    {
                                        continue;
                                    }
                                    stringstream cov_info;
                                    stringstream fre_info;
                                    for (auto c : res_cov)
                                    {
                                        cov_info << c << "\t";
                                        fre_info << (double)c / sum << "\n";
                                    }
                                    cov_info << s << "\t" << cov_info_all.str();
                                    allfre << fre_info.str();
                                    switch (res_cov.size())
                                    {
                                    case 2:
                                        allele[0] = allele[0] + 1;
                                        bifre << fre_info.str();
                                        bicov << cov_info.str();
                                        break;
                                    case 3:
                                        allele[1] = allele[1] + 1;
                                        trifre << fre_info.str();
                                        tricov << cov_info.str();
                                        break;
                                    case 4:
                                        allele[2] = allele[2] + 1;
                                        tetrafre << fre_info.str();
                                        tetracov << cov_info.str();
                                        break;
                                    case 5:
                                        allele[3] = allele[3] + 1;
                                        pentafre << fre_info.str();
                                        pentacov << cov_info.str();
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    exit_uni = *u.getSuccessors().begin();
                    while (exit_uni.getData()->getData(exit_uni)->get_id() != u.getData()->getData(u)->get_bubble_id(u.strand))
                    {
                        exit_uni = *exit_uni.getSuccessors().begin();
                    }
                    if (u.referenceUnitigToString().compare(exit_uni.referenceUnitigToString()) < 0)
                    {
                        if (u.strand)
                            u.getData()->getData(u)->set_plus_visited();
                        else
                            u.getData()->getData(u)->set_minus_visited();
                        continue;
                    }

                    vector<string> str_vec;
                    stack<UnitigColorMap<MyUnitig>> major;
                    stack<UnitigColorMap<MyUnitig>> minor;
                    string bubbleStr = "";
                    int k = cdbg.getK();
                    minor.push(u);
                    while (!minor.empty())
                    {
                        UnitigColorMap<MyUnitig> umi = minor.top();
                        minor.pop();
                        major.push(umi);
                        string str = umi.mappedSequenceToString();
                        bubbleStr += str.substr(0, umi.len);
                        if (umi.isSameReferenceUnitig(exit_uni))
                        {
                            bubbleStr += str.substr(umi.len);
                            str_vec.push_back(bubbleStr.substr(u.len - 1, bubbleStr.length() - u.len + 1 - umi.len + 1));
                            bubbleStr = bubbleStr.substr(0, bubbleStr.length() - str.length());
                            major.pop();
                            while (!major.empty() && !minor.empty())
                            {
                                bool f = false;
                                for (auto uma : major.top().getSuccessors())
                                {
                                    if (uma == minor.top())
                                    {
                                        f = true;
                                        break;
                                    }
                                }
                                if (f == false)
                                {
                                    bubbleStr = bubbleStr.substr(0, bubbleStr.length() - major.top().len);
                                    major.pop();
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }
                        else
                        {
                            for (auto u : umi.getSuccessors())
                            {
                                minor.push(u);
                            }
                        }
                    }
                    sortSeq_branching(str_vec, 0, str_vec.size() - 1);
                    SeqAlign seqalign(match, mismatch, gap);
                    vector<uint> var_site;
                    vector<uint> snp_pos;
                    vector<uint> indel_pos;
                    vector<uint> indel_len;
                    vector<vector<unsigned short>> partition;
                    seqalign.SequenceAlignment(str_vec, snp_pos, indel_pos, partition, indel_len);
                    if (str_vec.size() != 0)
                    {
                        ++coreNum;
                        coreCov += (size_t)core.first;
                        ++var_count;
                        for (auto s : str_vec)
                        {
                            s_var << var_count << "\t" << 0 << "\t" << u.getData()->getData(u)->get_id() << "\t" << exit_uni.getData()->getData(exit_uni)->get_id() << "\t" << s << endl;
                        }
                        for (int i = 0; i < partition.size(); i++)
                        {
                            if (partition[i].back() > 0)
                            {
                                var_site.push_back(i);
                            }
                        }
                        uint indel = 0;
                        for (uint i = 0; i < var_site.size(); i++)
                        {
                            stringstream cov_info;
                            stringstream fre_info;
                            uint var_distance = 0;
                            unsigned short maxnum = *max_element(partition[var_site[i]].begin(), partition[var_site[i]].end());
                            vector<string> k_length_str(str_vec.size());
                            vector<set<string>> string_set(maxnum);
                            if (i == 0)
                            {
                                if (i != var_site.size() - 1)
                                {
                                    var_distance = min((size_t)(var_site[i + 1] - var_site[i] - 1), u.size);
                                }
                                else
                                {
                                    var_distance = min(u.size, exit_uni.size);
                                }
                            }
                            else if (i == var_site.size() - 1)
                            {
                                var_distance = min((size_t)(var_site[i] - var_site[i - 1] - 1), exit_uni.size);
                            }
                            else
                            {
                                var_distance = min(var_site[i] - var_site[i - 1] - 1, var_site[i + 1] - var_site[i] - 1);
                            }
                            if (find(indel_pos.begin(), indel_pos.end(), var_site[i]) != indel_pos.end())
                            {
                                vector<int> site_vec(str_vec.size(), var_site[i]);
                                while (true)
                                {
                                    set<char> site_char;
                                    for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                    {
                                        string c = str_vec[s_size].substr(site_vec[s_size], 1);
                                        while (c.compare("-") == 0)
                                        {
                                            site_vec[s_size] += 1;
                                            c = str_vec[s_size].substr(site_vec[s_size], 1);
                                        }
                                        site_vec[s_size] += 1;
                                        k_length_str[s_size] += c;
                                        site_char.insert(c[0]);
                                    }
                                    if (site_char.size() > 1)
                                    {
                                        break;
                                    }
                                }
                                if (indel == 0)
                                {
                                    for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                    {
                                        int indel_i = k_length_str[s_size].length();
                                        k_length_str[s_size] = str_vec[s_size].substr(var_site[i] - cdbg.getK() + indel_i, cdbg.getK() - indel_i) + k_length_str[s_size];
                                    }
                                }
                                else
                                {
                                    for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                    {
                                        int indel_i = k_length_str[s_size].length();
                                        string temp_str;
                                        temp_str = str_vec[s_size].substr(0, var_site[i]);
                                        temp_str.erase(std::remove(temp_str.begin(), temp_str.end(), '-'), temp_str.end());
                                        if (temp_str.length() < cdbg.getK() - indel_i)
                                        {
                                            k_length_str[s_size] = temp_str + k_length_str[s_size];
                                            for (int extend_site = site_vec[s_size]; k_length_str[s_size].length() < cdbg.getK(); extend_site++)
                                            {
                                                string c = str_vec[s_size].substr(extend_site, 1);
                                                if (c.compare("-") != 0)
                                                {
                                                    k_length_str[s_size] += c;
                                                }
                                            }
                                        }
                                        else
                                            k_length_str[s_size] = temp_str.substr(temp_str.length() - cdbg.getK() + indel_i, cdbg.getK() - indel_i) + k_length_str[s_size];
                                    }
                                }
                                ++indel;
                                for (int part_i = 0; part_i < partition[var_site[i]].size(); part_i++)
                                {
                                    string_set[partition[var_site[i]][part_i] - 1].insert(k_length_str[part_i]);
                                }
                                vector<vector<double>> cov_vec(cdbg.getNbColors());
                                for (int i = 0; i < cov_vec.size(); i++)
                                {
                                    cov_vec[i].resize(maxnum, 0.0);
                                }
                                set<size_t> color_set;
                                bool flag_site_cov = true;
                                for (int string_i = 0; string_i < maxnum && flag_site_cov; ++string_i)
                                {
                                    for (auto s : string_set[string_i])
                                    {
                                        UnitigColorMap<MyUnitig> path_uni = cdbg.findUnitig(s.c_str(), 0, s.length());
                                        for (size_t ci = 0; ci < cdbg.getNbColors(); ci++)
                                        {
                                            if (path_uni.getData()->getUnitigColors(path_uni)->contains(path_uni, ci))
                                            {
                                                color_set.insert(ci);
                                                pair<double, bool> cov_res = readCov(s, cutoff[ci].first, cutoff[ci].second, ci);
                                                if (cov_res.second == false)
                                                {
                                                    flag_site_cov = false;
                                                    break;
                                                }
                                                cov_vec[ci][string_i] = cov_res.first + cov_vec[ci][string_i];
                                            }
                                            else
                                            {
                                                cov_vec[ci][string_i] = cov_vec[ci][string_i] + 0;
                                            }
                                        }
                                        if (flag_site_cov == false)
                                        {
                                            break;
                                        }
                                    }
                                }
                                if (color_set.size() != cdbg.getNbColors())
                                {
                                    continue;
                                }
                                if (flag_site_cov == false)
                                {
                                    continue;
                                }
                                double coefficient = 0;
                                for (size_t ci = 0; ci < cov_vec.size() - 1; ci++)
                                {
                                    for (size_t cj = ci + 1; cj < cov_vec.size(); cj++)
                                    {
                                        coefficient = max(coefficient, computeCramerVCoefficient(cov_vec[ci], cov_vec[cj]));
                                    }
                                }
                                for (int ci = 0; ci < cdbg.getNbColors(); ci++)
                                {
                                    vector<double> temp_cov;
                                    double sum = 0;
                                    for (auto c : cov_vec[ci])
                                    {
                                        if (c > 0.0)
                                        {
                                            sum += c;
                                            temp_cov.push_back(c);
                                        }
                                    }
                                    if (temp_cov.size() < 2)
                                    {
                                        continue;
                                    }
                                    stringstream cov_info;
                                    stringstream fre_info;
                                    for (auto c : temp_cov)
                                    {
                                        cov_info << c << "\t";
                                        fre_info << (double)c / sum << "\n";
                                    }
                                    cov_info << ci << "\t" << 0 << "\t" << indel_len[indel - 1] << "\t" << var_count << "\t" << var_site.size() << "\t" << coefficient << "\t" << var_distance << "\t" << endl;
                                    allfre << fre_info.str();
                                    switch (temp_cov.size())
                                    {
                                    case 2:
                                        ++allele[0];
                                        bifre << fre_info.str();
                                        bicov << cov_info.str();
                                        break;
                                    case 3:
                                        ++allele[1];
                                        trifre << fre_info.str();
                                        tricov << cov_info.str();
                                        break;
                                    case 4:
                                        ++allele[2];
                                        tetrafre << fre_info.str();
                                        tetracov << cov_info.str();
                                        break;
                                    case 5:
                                        ++allele[3];
                                        pentafre << fre_info.str();
                                        pentacov << cov_info.str();
                                    }
                                }
                            }
                            else
                            {
                                if (indel > 0)
                                {

                                    for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                    {
                                        string temp_str;
                                        temp_str = str_vec[s_size].substr(0, var_site[i] + 1);
                                        temp_str.erase(std::remove(temp_str.begin(), temp_str.end(), '-'), temp_str.end());
                                        if (temp_str.length() < cdbg.getK())
                                        {
                                            k_length_str[s_size] = temp_str;
                                            for (int extend_site = var_site[i] + 1; k_length_str[s_size].length() < cdbg.getK(); extend_site++)
                                            {
                                                string c = str_vec[s_size].substr(extend_site, 1);
                                                if (c.compare("-") != 0)
                                                {
                                                    k_length_str[s_size] += c;
                                                }
                                            }
                                        }
                                        else
                                            k_length_str[s_size] = (temp_str.substr(temp_str.length() - cdbg.getK(), cdbg.getK()));
                                    }
                                }
                                else
                                {
                                    for (int s_size = 0; s_size < str_vec.size(); s_size++)
                                    {
                                        k_length_str[s_size] = (str_vec[s_size].substr(var_site[i] - cdbg.getK() + 1, cdbg.getK()));
                                    }
                                }
                                for (int part_i = 0; part_i < partition[var_site[i]].size(); part_i++)
                                {
                                    string_set[partition[var_site[i]][part_i] - 1].insert(k_length_str[part_i]);
                                }
                                vector<vector<double>> cov_vec(cdbg.getNbColors());
                                for (int i = 0; i < cov_vec.size(); i++)
                                {
                                    cov_vec[i].resize(maxnum, 0);
                                }
                                set<size_t> color_set;
                                bool flag_site_cov = true;
                                for (int string_i = 0; string_i < string_set.size() && flag_site_cov; string_i++)
                                {
                                    for (auto s : string_set[string_i])
                                    {

                                        UnitigColorMap<MyUnitig> path_uni = cdbg.findUnitig(s.c_str(), 0, s.length());
                                        for (size_t ci = 0; ci < cdbg.getNbColors(); ci++)
                                        {
                                            if (path_uni.getData()->getUnitigColors(path_uni)->contains(path_uni, ci))
                                            {
                                                color_set.insert(ci);
                                                pair<double, bool> cov_res = readCov(s, cutoff[ci].first, cutoff[ci].second, ci);
                                                if (cov_res.second == false)
                                                {
                                                    flag_site_cov = false;
                                                    break;
                                                }
                                                cov_vec[ci][string_i] += cov_res.first;
                                            }
                                        }
                                        if (flag_site_cov == false)
                                        {
                                            break;
                                        }
                                    }
                                }
                                if (color_set.size() != cdbg.getNbColors())
                                {
                                    continue;
                                }
                                if (flag_site_cov == false)
                                {
                                    continue;
                                }
                                double coefficient = 0;
                                for (size_t ci = 0; ci < cov_vec.size() - 1; ci++)
                                {
                                    for (size_t cj = ci + 1; cj < cov_vec.size(); cj++)
                                    {
                                        coefficient = max(coefficient, computeCramerVCoefficient(cov_vec[ci], cov_vec[cj]));
                                    }
                                }
                                for (int ci = 0; ci < cdbg.getNbColors(); ci++)
                                {
                                    vector<double> temp_cov;
                                    double sum = 0;
                                    for (auto c : cov_vec[ci])
                                    {
                                        if (c > 0.0)
                                        {
                                            sum += c;
                                            temp_cov.push_back(c);
                                        }
                                    }
                                    if (temp_cov.size() < 2)
                                    {
                                        continue;
                                    }
                                    stringstream cov_info;
                                    stringstream fre_info;
                                    for (auto c : temp_cov)
                                    {
                                        cov_info << c << "\t";
                                        fre_info << (double)c / sum << "\n";
                                    }
                                    cov_info << ci << "\t" << 0 << "\t"
                                             << "0\t" << var_count << "\t" << var_site.size() << "\t" << coefficient << "\t" << var_distance << "\t" << endl;
                                    allfre << fre_info.str();
                                    switch (temp_cov.size())
                                    {
                                    case 2:
                                        ++allele[0];
                                        bifre << fre_info.str();
                                        bicov << cov_info.str();
                                        break;
                                    case 3:
                                        ++allele[1];
                                        trifre << fre_info.str();
                                        tricov << cov_info.str();
                                        break;
                                    case 4:
                                        ++allele[2];
                                        tetrafre << fre_info.str();
                                        tetracov << cov_info.str();
                                        break;
                                    case 5:
                                        ++allele[3];
                                        pentafre << fre_info.str();
                                        pentacov << cov_info.str();
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (u.strand)
            {
                u.getData()->getData(u)->set_plus_visited();
                if (exit_uni.strand == true)
                {
                    exit_uni.getData()->getData(exit_uni)->set_minus_visited();
                }
                else
                {
                    exit_uni.getData()->getData(exit_uni)->set_plus_visited();
                }
            }
            else
            {
                u.getData()->getData(u)->set_minus_visited();
                if (exit_uni.strand == true)
                {
                    exit_uni.getData()->getData(exit_uni)->set_minus_visited();
                }
                else
                {
                    exit_uni.getData()->getData(exit_uni)->set_plus_visited();
                }
            }
        }
    }
    time_t end_time = time(NULL);
    cout << "CCDBG::PloidyEstimation(): Cpu time : "
         << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
    cout << "CCDBG::PloidyEstimation(): Real time : "
         << (double)difftime(end_time, start_time) << "s" << endl;
    cout << "CCDBG::PloidyEstimation(): Alleles in SuperBubbles  :\t"
         << "2 :" << allele[0] << "\t"
         << "3 :" << allele[1] << "\t"
         << "4 :" << allele[2] << "\t"
         << "5 :" << allele[3] << endl;
    bifre.close();
    trifre.close();
    tetrafre.close();
    bicov.close();
    tricov.close();
    tetracov.close();
    pentacov.close();
    pentafre.close();
    allfre.close();
    if (coreNum != 0)
    {
        const int avg = coreCov / coreNum;
        cout << "CCDBG::PloidyEstimation(): Sites' Average Coverage:" << avg << endl;
    }
}
