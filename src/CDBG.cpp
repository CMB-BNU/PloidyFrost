#include "CDBG.hpp"
#include <string>
#include "kmc_file.h"
#include <stack>
#include <thread>
#include <ctime>
#include <time.h>
#include <sstream>
#include "SeqAlign.hpp"
#include <queue>
#include <algorithm>
using namespace std;

CDBG::CDBG(CompactedDBG<MyUnitig> &graph, const size_t &complexsize, double &m, double &d, double &g, string db) : cdbg(graph), complex_size(complexsize), match(m), mismatch(d), gap(g)
{
    if (!db.empty())
        if (!kmer_databas.OpenForRA(db))
        {
            cerr << "CDBG::CDBG():Error: Open kmc database error ." << endl;
            exit(EXIT_FAILURE);
        }
    cout << "CDBG::CDBG():CDBG initialized!" << endl;
}
CDBG::~CDBG(void)
{
    kmer_databas.Close();
}

pair<double, bool> CDBG::readCov(const string &s, const uint &low, const uint &up)
{
    double sum = 0;
    CKmerAPI kmer_object(cdbg.getK());
    uint count = 0;
    if (kmer_databas.GetBothStrands())
    {
        for (int i = 0; i <= s.length() - cdbg.getK(); i++)
        {
            kmer_object.from_string(s.substr(i, cdbg.getK()));
            if (!kmer_databas.IsKmer(kmer_object))
            {
                kmer_object.reverse();
            }
            if (kmer_databas.CheckKmer(kmer_object, count))
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
                cout << "CDBG::readCov():" << kmer_object.to_string() << " kmer can not found ." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    return pair<double, bool>(sum / (s.length() - cdbg.getK() + 1), true);
}
/**
 * @brief get unitig average coverage and max kmer coverage difference from kmc database
 * @param u    unitig
 * @return  (unitig average coverage,max coverage difference)
 */
pair<double, uint> CDBG::readCov(const UnitigMap<MyUnitig> &u)
{
    double sum = 0;
    CKmerAPI kmer_object(cdbg.getK());
    uint count = 0;
    uint min = 10000;

    if (kmer_databas.GetBothStrands())
    {
        for (int i = 0; i < u.len; i++)
        {
            kmer_object.from_string(u.referenceUnitigToString().substr(i, cdbg.getK()));
            if (!kmer_databas.IsKmer(kmer_object))
            {
                kmer_object.reverse();
            }
            if (kmer_databas.CheckKmer(kmer_object, count))
            {

                sum += count;

                if (count < min)
                {
                    min = count;
                }
            }
            else
            {
                cout << "CDBG::readCov():" << kmer_object.to_string() << " kmer can not found ." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    else
    {
        for (int i = 0; i < u.len; i++)
        {
            kmer_object.from_string(u.mappedSequenceToString().substr(i, cdbg.getK()));
            if (!kmer_databas.CheckKmer(kmer_object, count))
            {
                cout << "CDBG::readCov():" << kmer_object.to_string() << " kmer can not found ." << endl;
                exit(EXIT_FAILURE);
            }
            else
            {
                sum += count;
                if (count < min)
                {
                    min = count;
                }
            }
        }
    }
    return pair<double, uint>(sum / u.len, min);
}
void CDBG::setUnitigId(const string &outpre, const string &graphfile, const size_t &thr)
{
    if (access("PloidyFrost_output", 0))
    {
        system("mkdir ./PloidyFrost_output");
    }
    ofstream output("PloidyFrost_output/" + outpre + "_Unitig_Id.txt", ofstream::out);
    cout << "CDBG::setUnitigId(): Setting Unitig Id" << endl;
    double start_time = time(NULL);
    clock_t start_clock = clock();
    size_t id = 0;
    for (const auto &unitig : cdbg)
    {
        unitig.getData()->set_id(++id);
        output << id << "\t" << unitig.referenceUnitigToString() << endl;
    }
    output.close();
    time_t end_time = time(NULL);
    cout << "CDBG::setUnitigId(): Cpu time : "
         << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
    cout << "CDBG::setUnitigId(): Real time : "
         << (double)difftime(end_time, start_time) << "s" << endl;
}
void CDBG::printInfo(const bool &verbose, const string &outpre)
{
    if (verbose)
    {
        cout << ">>>>>>>>>Bifrost Graph Information>>>>>>>>>" << endl;
        cout << "k:" << cdbg.getK() << "\t";
        cout << "g:" << cdbg.getG() << "\t";
        cout << "nbKmer:" << cdbg.nbKmers() << "\t";
        cout << "nbUnitig:" << cdbg.size() << "\t";
        cout << "length:" << cdbg.length() << "\t" << endl;
    }
    ofstream output(outpre + "_graph_info.txt", ofstream::out);
    output << "k:" << cdbg.getK() << "\t";
    output << "g:" << cdbg.getG() << "\t";
    output << "nbKmer:" << cdbg.nbKmers() << "\t";
    output << "nbUnitig:" << cdbg.size() << "\t";
    output << "length:" << cdbg.length() << "\t";
    output.close();
}
void CDBG::clearMarking(const vector<UnitigMap<MyUnitig>> &set_km_seen)
{
    for (const auto &ucm : set_km_seen)
    {
        ucm.getData()->clear(ucm);
    }
}
void CDBG::clearMarking()
{
    for (const auto &unitig : cdbg)
    {
        MyUnitig *data_ucm = unitig.getData();
        data_ucm->clear(unitig);
    }
}
void CDBG::findSuperBubble_ptr(const string &outpre)
{
    cout << "CDBG::findSuperBubble(): Finding superbubbles" << endl;
    size_t nb_super_bubble = 0;
    ofstream outfile;
    if (access("PloidyFrost_output", 0))
    {
        system("mkdir ./PloidyFrost_output");
    }
    outfile.open("PloidyFrost_output/" + outpre + "_super_bubble.txt", ios::out | ios::trunc);
    if (outfile.fail())
    {
        cout << "CDBG:: Open super_bubble file error";
        exit(EXIT_FAILURE);
    }
    clock_t start_clock = clock();
    double start_time = time(NULL);
    size_t nb_unitig_processed = 0;
    cout << "CDBG::findSuperBubble(): There are " << cdbg.size() << " unitigs " << endl;
    for (const auto &unitig : cdbg)
    {
        ++nb_unitig_processed;
        if (nb_unitig_processed % 100000 == 0)
        {
            cout << "CDBG::findSuperBubble(): Processed " << nb_unitig_processed << " unitigs " << endl;
        }
        UnitigMap<MyUnitig> u(unitig);
        u.strand = true;
        if (u.getSuccessors().cardinality() > 1 && u.getData()->get_plus() == NULL)
        {
            extractSuperBubble_ptr(u);
        }
        u.strand = false;
        if (u.getSuccessors().cardinality() > 1 && u.getData()->get_minus() == NULL)
        {
            extractSuperBubble_ptr(u);
        }
    }
    time_t end_time = time(NULL);
    cout << "CDBG::findSuperBubble():  Cpu time : "
         << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
    cout << "CDBG::findSuperBubble():  Real time : "
         << (double)difftime(end_time, start_time) << "s" << endl;
    outfile << "BubbleId\tEntrance\tStrand\tExit\tisSimple\tisComplex" << endl;
    for (const auto &unitig : cdbg)
    {
        MyUnitig *data = unitig.getData();
        if (data->is_both_visited())
        {
            continue;
        }
        if (!data->is_plus_visited())
        {
            ++nb_super_bubble;
            outfile << nb_super_bubble << "\t" << data->get_id() << "\t"
                    << "+"
                    << "\t" << data->get_plus()->get_id()
                    << "\t" << (data->isStrict(true) ? "1" : "0")
                    << "\t" << (data->isComplex(true) ? "1" : "0")
                    << endl;
        }
        if (!data->is_minus_visited())
        {
            ++nb_super_bubble;
            outfile << nb_super_bubble << "\t" << data->get_id() << "\t"
                    << "-"
                    << "\t" << data->get_minus()->get_id()
                    << "\t" << (data->isStrict(false) ? "1" : "0")
                    << "\t" << (data->isComplex(false) ? "1" : "0")
                    << endl;
        }
    }
    outfile.close();
    cout << "CDBG::findSuperBubble(): " << nb_super_bubble << "  SuperBubbles Found" << endl;
}
void CDBG::extractSuperBubble_ptr(const UnitigMap<MyUnitig> &s)
{
    bool flag_cycle = false;
    bool flag_tip = false;
    bool flag_large_cycle = false;
    vector<UnitigMap<MyUnitig>> vertices_visit;
    pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> p;
    vector<UnitigMap<MyUnitig>> vec_km_seen;
    map<size_t, uint8_t> state_map;
    map<size_t, bool> strand_map;
    unordered_set<UnitigMap<MyUnitig>, UnitigMapHash<MyUnitig>> cycle_unitig_set;
    UnitigMap<MyUnitig> v(s);
    vertices_visit.push_back(v);
    vec_km_seen.push_back(v);
    while (!vertices_visit.empty())
    {
        v = vertices_visit.back();
        vertices_visit.pop_back();
        state_map[v.getData()->get_id()] = 0x01;
        strand_map[v.getData()->get_id()] = v.strand;
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
                if (state_map.find(u.getData()->get_id()) == state_map.end() || state_map[u.getData()->get_id()] != 0x01)
                {
                    if (state_map.find(u.getData()->get_id()) == state_map.end())
                    {
                        vec_km_seen.push_back(u);
                        strand_map[u.getData()->get_id()] = u.strand;
                    }
                    else
                    {
                        if (strand_map[u.getData()->get_id()] != u.strand)
                        {
                            flag_cycle = true;
                            cycle_unitig_set.insert(u);
                            cycle_unitig_set.insert(v);
                        }
                    }
                    state_map[u.getData()->get_id()] = 0x02;
                    bool all_predecessor_visited = true;
                    for (const auto &predecessor : u.getPredecessors())
                    {
                        if (state_map.find(predecessor.getData()->get_id()) != state_map.end())
                        {
                            if (state_map[predecessor.getData()->get_id()] != 0x01)
                            {
                                all_predecessor_visited = false;
                            }
                            if (strand_map[predecessor.getData()->get_id()] != predecessor.strand)
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
                    if (state_map[cucm.getData()->get_id()] == 0x02)
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
            if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
            {
                MyUnitig *ex = ucm.getData()->get_plus();
                if (ex->get_plus() == ucm.getData())
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->set_plus_self();
            if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
            {
                MyUnitig *ex = ucm.getData()->get_minus();
                if (ex->get_plus() == ucm.getData())
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->set_minus_self();
            ucm.getData()->set_non_super();
        }
        if (s.strand == true)
        {
            s.getData()->set_plus_self();
        }
        else
        {
            s.getData()->set_minus_self();
        }
    }
    return;
}

void CDBG::sortSeq_branching(vector<string> &str_vec, int low, int high)
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

void CDBG::sortSeq_simple(vector<double> &cov, vector<UnitigMap<MyUnitig>> &unitig_vec, int low, int high)
{
    if (high <= low)
        return;
    int i = low;
    int j = high;
    while (true)
    {
        while (cov[i] >= cov[low])
        {
            if (cov[i] > cov[low])
            {
                i++;
            }
            else
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
            if (i == high)
            {
                break;
            }
        }
        while (cov[j] <= cov[low])
        {
            if (cov[j] < cov[low])
            {
                j--;
            }
            else
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
            if (j == low)
            {
                break;
            }
        }
        if (i >= j)
            break;
        double t = cov[i];
        cov[i] = cov[j];
        cov[j] = t;
        UnitigMap<MyUnitig> u = unitig_vec[j];
        unitig_vec[j] = unitig_vec[i];
        unitig_vec[i] = u;
    }
    double t = cov[low];
    cov[low] = cov[j];
    cov[j] = t;
    UnitigMap<MyUnitig> u = unitig_vec[j];
    unitig_vec[j] = unitig_vec[low];
    unitig_vec[low] = u;
    sortSeq_simple(cov, unitig_vec, low, j - 1);
    sortSeq_simple(cov, unitig_vec, j + 1, high);
}
void CDBG::setNoBubble_ptr_cycle(const vector<UnitigMap<MyUnitig>> &vec_km_seen, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &p)
{
    for (const auto &ucm : vec_km_seen)
    {

        if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_plus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_plus_self();
        if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_minus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_minus_self();

        ucm.getData()->set_non_super();
    }
    if (p.first.strand == true)
    {
        p.first.getData()->set_plus_self();
    }
    else
    {
        p.first.getData()->set_minus_self();
    }
    if (p.second.strand == false)
    {
        p.second.getData()->set_plus_self();
    }
    else
    {
        p.second.getData()->set_minus_self();
    }
}
void CDBG::setNoBubble_ptr(const vector<UnitigMap<MyUnitig>> &vec_km_seen, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &p)
{
    if (p.first.strand == true)
    {
        if (p.first.getData()->get_plus() != NULL)
        {
            if (p.first.getData()->get_plus()->get_plus() == p.first.getData())
            {
                p.first.getData()->get_plus()->set_plus_self();
            }
            else
            {
                p.first.getData()->get_plus()->set_minus_self();
            }
        }
        p.first.getData()->set_plus_self();
    }
    else
    {
        if (p.first.getData()->get_minus() != NULL)
        {
            if (p.first.getData()->get_minus()->get_plus() == p.first.getData())
            {
                p.first.getData()->get_minus()->set_plus_self();
            }
            else
            {
                p.first.getData()->get_minus()->set_minus_self();
            }
        }
        p.first.getData()->set_minus_self();
    }
    if (p.second.strand == false)
    {
        if (p.second.getData()->get_plus() != NULL)
        {
            if (p.second.getData()->get_plus()->get_plus() == p.second.getData())
            {
                p.second.getData()->get_plus()->set_plus_self();
            }
            else
            {
                p.second.getData()->get_plus()->set_minus_self();
            }
        }
        p.second.getData()->set_plus_self();
    }
    else
    {
        if (p.second.getData()->get_minus() != NULL)
        {
            if (p.second.getData()->get_minus()->get_plus() == p.second.getData())
            {
                p.second.getData()->get_minus()->set_plus_self();
            }
            else
            {
                p.second.getData()->get_minus()->set_minus_self();
            }
        }
        p.second.getData()->set_minus_self();
    }
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm == p.first || ucm == p.second)
        {
            continue;
        }
        if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_plus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_plus_self();
        if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_minus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_minus_self();
        ucm.getData()->set_non_super();
    }
}
void CDBG::setNoBubble_ptr(pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &p, const vector<UnitigMap<MyUnitig>> &vec_km_seen)
{
    if (vec_km_seen.size() < 4)
    {
        return;
    }
    if (p.second.getData()->is_non_super() || p.first.getData()->is_non_super())
    {
        for (const auto &ucm : vec_km_seen)
        {
            if (ucm == p.first)
            {
                if (p.first.strand == true)
                {
                    p.first.getData()->set_plus_self();
                }
                else
                {
                    p.first.getData()->set_minus_self();
                }
                continue;
            }
            if (ucm == p.second)
            {
                if (p.second.strand == true)
                {
                    p.second.getData()->set_minus_self();
                }
                else
                {
                    p.second.getData()->set_plus_self();
                }
                continue;
            }
            if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
            {
                MyUnitig *ex = ucm.getData()->get_plus();
                if (ex->get_plus() == ucm.getData())
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->set_plus_self();
            if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
            {
                MyUnitig *ex = ucm.getData()->get_minus();
                if (ex->get_plus() == ucm.getData())
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->set_minus_self();
            ucm.getData()->set_non_super();
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
            p.first.getData()->setStrict(p.first.strand);
            p.second.getData()->setStrict(!p.second.strand);
        }
    }
    if (vec_km_seen.size() > complex_size)
    {
        p.first.getData()->setComplex(p.first.strand);
        p.second.getData()->setComplex(!p.second.strand);
    }
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm == p.first || ucm == p.second)
        {
            continue;
        }
        if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_plus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_plus_self();
        if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_minus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_minus_self();
        ucm.getData()->set_non_super();
    }
    MyUnitig *data = p.first.getData();
    if (p.first.strand == true)
    {
        data->set_plus(p.second.getData());
    }
    else
    {
        data->set_minus(p.second.getData());
    }
    data = p.second.getData();
    if (p.second.strand == true)
    {
        data->set_minus(p.first.getData());
    }
    else
    {
        data->set_plus(p.first.getData());
    }
}
void CDBG::setNoBubble_multithread_ptr(const vector<UnitigMap<MyUnitig>> &vec_km_seen, mutex &x, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &p)
{
    lock_guard<std::mutex> lock(x);
    if (p.first.strand == true)
    {
        if (p.first.getData()->get_plus() != NULL)
        {
            if (p.first.getData()->get_plus()->get_plus() == p.first.getData())
            {
                p.first.getData()->get_plus()->set_plus_self();
            }
            else
            {
                p.first.getData()->get_plus()->set_minus_self();
            }
        }
        p.first.getData()->set_plus_self();
    }
    else
    {
        if (p.first.getData()->get_minus() != NULL)
        {
            if (p.first.getData()->get_minus()->get_plus() == p.first.getData())
            {
                p.first.getData()->get_minus()->set_plus_self();
            }
            else
            {
                p.first.getData()->get_minus()->set_minus_self();
            }
        }
        p.first.getData()->set_minus_self();
    }
    if (p.second.strand == false)
    {
        if (p.second.getData()->get_plus() != NULL)
        {
            if (p.second.getData()->get_plus()->get_plus() == p.second.getData())
            {
                p.second.getData()->get_plus()->set_plus_self();
            }
            else
            {
                p.second.getData()->get_plus()->set_minus_self();
            }
        }
        p.second.getData()->set_plus_self();
    }
    else
    {
        if (p.second.getData()->get_minus() != NULL)
        {
            if (p.second.getData()->get_minus()->get_plus() == p.second.getData())
            {
                p.second.getData()->get_minus()->set_plus_self();
            }
            else
            {
                p.second.getData()->get_minus()->set_minus_self();
            }
        }
        p.second.getData()->set_minus_self();
    }
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm == p.first || ucm == p.second)
        {
            continue;
        }
        if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_plus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_plus_self();

        if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_minus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_minus_self();

        ucm.getData()->set_non_super();
    }
}
void CDBG::setNoBubble_multithread_ptr(const vector<UnitigMap<MyUnitig>> &vec_km_seen, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &p, mutex &m)
{
    if (vec_km_seen.size() < 4)
    {
        return;
    }
    lock_guard<std::mutex> lock(m);

    if (p.second.getData()->get_plus_or_minus(!p.second.strand) == p.first.getData() && !p.second.getData()->is_non_super() && !p.first.getData()->is_non_super())
    {
        return;
    }
    if (p.second.getData()->get_plus_or_minus(!p.second.strand) != NULL || p.second.getData()->is_non_super() || p.first.getData()->is_non_super())
    {

        for (const auto &ucm : vec_km_seen)
        {
            if (ucm == p.first)
            {
                if (p.first.strand == true)
                {
                    p.first.getData()->set_plus_self();
                }
                else
                {
                    p.first.getData()->set_minus_self();
                }
                continue;
            }
            if (ucm == p.second)
            {
                if (p.second.strand == true)
                {
                    p.second.getData()->set_minus_self();
                }
                else
                {
                    p.second.getData()->set_plus_self();
                }
                continue;
            }
            if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
            {
                MyUnitig *ex = ucm.getData()->get_plus();
                if (ex->get_plus() == ucm.getData())
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->set_plus_self();
            if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
            {
                MyUnitig *ex = ucm.getData()->get_minus();
                if (ex->get_plus() == ucm.getData())
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->set_minus_self();
            ucm.getData()->set_non_super();
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
            p.first.getData()->setStrict(p.first.strand);
            p.second.getData()->setStrict(!p.second.strand);
        }
    }
    if (vec_km_seen.size() > complex_size)
    {
        p.first.getData()->setComplex(p.first.strand);
        p.second.getData()->setComplex(!p.second.strand);
    }
    for (const auto &ucm : vec_km_seen)
    {
        if (ucm == p.first || ucm == p.second)
        {
            continue;
        }
        if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_plus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_plus_self();
        if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_minus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_minus_self();
        ucm.getData()->set_non_super();
    }
    MyUnitig *data = p.first.getData();
    if (p.first.strand == true)
    {
        data->set_plus(p.second.getData());
    }
    else
    {
        data->set_minus(p.second.getData());
    }
    data = p.second.getData();
    if (p.second.strand == true)
    {
        data->set_minus(p.first.getData());
    }
    else
    {
        data->set_plus(p.first.getData());
    }
}
void CDBG::ploidyEstimation_ptr(const string &outpre, const int &lower, const int &upper)
{
    clock_t start_clock = clock();
    double start_time = time(NULL);
    size_t coreNum = 0;
    size_t coreCov = 0;
    size_t allele[4] = {};
    cout << "CDBG::PloidyEstimation():  Analyzing superbubbles to generate sites' information" << endl;
    if (access("PloidyFrost_output", 0))
    {
        system("mkdir ./PloidyFrost_output");
    }
    ofstream outfile;
    ofstream bifre;
    ofstream trifre;
    ofstream tetrafre;
    ofstream pentafre;
    ofstream bicov;
    ofstream tricov;
    ofstream tetracov;
    ofstream pentacov;

    ofstream allfre;
    ofstream s_var;
    allfre.open("PloidyFrost_output/" + outpre + "_allele_frequency.txt", ios::out | ios::trunc);
    bifre.open("PloidyFrost_output/" + outpre + "_bifre.txt", ios::out | ios::trunc);
    trifre.open("PloidyFrost_output/" + outpre + "_trifre.txt", ios::out | ios::trunc);
    tetrafre.open("PloidyFrost_output/" + outpre + "_tetrafre.txt", ios::out | ios::trunc);
    pentafre.open("PloidyFrost_output/" + outpre + "_pentafre.txt", ios::out | ios::trunc);
    pentacov.open("PloidyFrost_output/" + outpre + "_pentacov.txt", ios::out | ios::trunc);
    bicov.open("PloidyFrost_output/" + outpre + "_bicov.txt", ios::out | ios::trunc);
    tricov.open("PloidyFrost_output/" + outpre + "_tricov.txt", ios::out | ios::trunc);
    tetracov.open("PloidyFrost_output/" + outpre + "_tetracov.txt", ios::out | ios::trunc);
    outfile.open("PloidyFrost_output/" + outpre + "_allele_frequency.txt", ios::out | ios::trunc);
    s_var.open("PloidyFrost_output/" + outpre + "_alignseq.txt", ios::out | ios::trunc);
    if (!allfre.is_open() || !outfile.is_open() || !bifre.is_open() || !trifre.is_open() ||
        !tetrafre.is_open() || !bicov.is_open() || !tricov.is_open() || !tetracov.is_open() || !s_var.is_open() ||
        !pentafre.is_open() || !pentacov.is_open())
    {
        cout << "CDBG:: PloidyEstimation():Open file error" << endl;
        exit(EXIT_FAILURE);
    }
    size_t var_count = 0;
    size_t nb_unitig_processed = 0;

    for (const auto &unitig : cdbg)
    {
        ++nb_unitig_processed;
        if (nb_unitig_processed % 100000 == 0)
        {
            cout << "CDBG::PloidyEstimation(): Processed " << nb_unitig_processed << " unitigs " << endl;
        }
        UnitigMap<MyUnitig> u(unitig);
        MyUnitig *u_data = u.getData();
        if (u_data->is_both_visited())
        {
            continue;
        }
        while (!u_data->is_both_visited())
        {
            if (u_data->is_plus_visited() != true)
            {
                u.strand = true;
                if (u_data->isComplex(u.strand))
                {
                    u_data->set_plus_visited();
                    continue;
                }
            }
            else if (u_data->is_minus_visited() != true)
            {
                u.strand = false;
                if (u_data->isComplex(u.strand))
                {
                    u_data->set_minus_visited();
                    break;
                }
            }
            else
            {
                break;
            }
            double sum = 0;
            bool isStrict = u_data->isStrict(u.strand);
            UnitigMap<MyUnitig> exit_uni = NULL;
            pair<double, uint> core = readCov(u);
            if (isStrict)
            {
                exit_uni = (*u.getSuccessors().begin()->getSuccessors().begin());
                if (u.referenceUnitigToString().compare(exit_uni.referenceUnitigToString()) < 0)
                {
                    if (u.strand)
                    {
                        u.getData()->set_plus_visited();
                    }
                    else
                    {
                        u.getData()->set_minus_visited();
                    }
                    continue;
                }
                vector<UnitigMap<MyUnitig>> unitig_vec;
                bool flag = true;
                vector<double> cov;
                uint8_t allele_count = 0;
                for (const auto &uu : u.getSuccessors())
                {
                    unitig_vec.push_back(uu);
                    pair<double, uint> inside = readCov(uu);
                    if (inside.second > lower && inside.second < upper)
                    {
                        cov.push_back(inside.first);
                        sum += inside.first;
                        allele_count++;
                    }
                    else
                    {
                        flag = false;
                        break;
                    }
                }
                if (flag == true)
                {
                    vector<double> cov_pre;
                    if (u.getPredecessors().hasPredecessors())
                    {
                        for (const auto &uu : u.getPredecessors())
                        {
                            pair<double, uint> inside = readCov(uu);
                            cov_pre.push_back(inside.first);
                        }
                    }
                    if (cov_pre.size() < cov.size())
                    {
                        for (int i = cov_pre.size(); i < cov.size(); i++)
                        {
                            cov_pre.push_back(0);
                        }
                    }
                    sortSeq_simple(cov, unitig_vec, 0, cov.size() - 1);
                    SeqAlign seqalign(match, mismatch, gap);
                    uint variantNum = 0;
                    vector<uint> var_site;
                    vector<uint> snp_pos;
                    vector<uint> indel_pos;
                    vector<uint> indel_len;
                    vector<vector<unsigned short>> partition;
                    vector<string> str_vec;
                    for (auto uni : unitig_vec)
                    {
                        str_vec.push_back(uni.mappedSequenceToString());
                    }
                    seqalign.SequenceAlignment(str_vec, snp_pos, indel_pos, partition, indel_len);
                    if (str_vec.size() != 0)
                    {
                        ++var_count;
                        for (auto s : str_vec)
                        {
                            s_var << var_count << "\t" << 1 << "\t" << u.getData()->get_id() << "\t" << exit_uni.getData()->get_id() << "\t" << s << endl;
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
                        uint pre_var_pos = 0;
                        uint indel = 0;
                        for (uint i = 0; i < var_site.size(); i++)
                        {
                            stringstream cov_info;
                            stringstream fre_info;
                            vector<unsigned short> part = partition[var_site[i]];
                            unsigned short maxnum = *max_element(part.begin(), part.end());
                            vector<double> temp_cov(maxnum, 0);
                            uint var_distance = 0;
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
                            for (int j = 0; j < part.size(); j++)
                            {
                                temp_cov[part[j] - 1] += cov[j];
                            }
                            for (auto c : temp_cov)
                            {
                                cov_info << c << "\t";
                                fre_info << (double)c / sum << "\n";
                            }
                            if (find(indel_pos.begin(), indel_pos.end(), var_site[i]) != indel_pos.end())
                            {
                                ++indel;
                                cov_info << 1 << "\t" << indel_len[indel - 1] << "\t" << var_count << "\t" << var_site.size() << "\t" << var_distance << "\t" << endl;
                            }
                            else
                            {
                                cov_info << 1 << "\t"
                                         << "0\t" << var_count << "\t" << var_site.size() << "\t" << var_distance << "\t" << endl;
                            }
                            allfre << fre_info.str();
                            switch (maxnum)
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
            else
            {
                exit_uni = *u.getSuccessors().begin();
                while (exit_uni.getData()->get_id() != u.getData()->get_bubble_id(u.strand))
                {
                    exit_uni = *exit_uni.getSuccessors().begin();
                }
                if (u.referenceUnitigToString().compare(exit_uni.referenceUnitigToString()) < 0)
                {
                    if (u.strand)
                    {
                        u.getData()->set_plus_visited();
                    }
                    else
                    {
                        u.getData()->set_minus_visited();
                    }
                    continue;
                }
                vector<string> str_vec;
                stack<UnitigMap<MyUnitig>> major;
                stack<UnitigMap<MyUnitig>> minor;
                string bubbleStr = "";
                int k = cdbg.getK();
                minor.push(u);
                while (!minor.empty())
                {
                    UnitigMap<MyUnitig> umi = minor.top();
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
                uint variantNum = 0;
                vector<uint> var_site;
                vector<uint> snp_pos;
                vector<uint> indel_pos;
                vector<uint> indel_len;
                vector<vector<unsigned short>> partition;
                seqalign.SequenceAlignment(str_vec, snp_pos, indel_pos, partition, indel_len);
                if (str_vec.size() != 0)
                {
                    ++var_count;
                    for (auto s : str_vec)
                    {
                        s_var << var_count << "\t" << 0 << "\t" << u.getData()->get_id() << "\t" << exit_uni.getData()->get_id() << "\t" << s << endl;
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
                    uint pre_var_pos = 0;
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
                            vector<double> temp_cov(maxnum, 0);
                            double sum = 0;
                            bool flag_site_cov = true;
                            for (int string_i = 0; string_i < string_set.size() && flag_site_cov; string_i++)
                            {
                                for (auto s : string_set[string_i])
                                {
                                    pair<double, bool> cov_res = readCov(s, lower, upper);
                                    if (cov_res.second == false)
                                    {
                                        flag_site_cov = false;
                                        break;
                                    }
                                    temp_cov[string_i] += cov_res.first;
                                }
                                sum += temp_cov[string_i];
                            }
                            if (flag_site_cov == false)
                            {
                                continue;
                            }
                            for (auto tc : temp_cov)
                            {
                                cov_info << tc << "\t";
                                fre_info << tc / sum << endl;
                            }
                            cov_info << 0 << "\t" << indel_len[indel - 1] << "\t" << var_count << "\t" << var_site.size() << "\t" << var_distance << "\t" << endl;
                        }
                        else
                        {
                            if (indel > 0)
                            {
                                uint indel_len_sum = 0;
                                for (int len_i = 0; len_i < indel; len_i++)
                                {
                                    indel_len_sum += indel_len[len_i];
                                }
                                int start = var_site[i] - cdbg.getK() - indel_len_sum + 1;
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
                            vector<double> temp_cov(maxnum, 0);
                            double sum = 0;
                            bool flag_site_cov = true;
                            for (int string_i = 0; string_i < string_set.size() && flag_site_cov; string_i++)
                            {
                                for (auto s : string_set[string_i])
                                {
                                    pair<double, bool> cov_res = readCov(s, lower, upper);
                                    if (cov_res.second == false)
                                    {
                                        flag_site_cov = false;
                                        break;
                                    }
                                    temp_cov[string_i] += cov_res.first;
                                }
                                sum += temp_cov[string_i];
                            }
                            if (flag_site_cov == false)
                            {
                                continue;
                            }
                            for (auto tc : temp_cov)
                            {
                                cov_info << tc << "\t";
                                fre_info << tc / sum << endl;
                            }
                            cov_info << 0 << "\t"
                                     << "0\t" << var_count << "\t" << var_site.size() << "\t" << var_distance << "\t" << endl;
                        }
                        allfre << fre_info.str();
                        switch (maxnum)
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
            if (u.strand)
            {
                u.getData()->set_plus_visited();
                if (exit_uni.strand == true)
                {
                    exit_uni.getData()->set_minus_visited();
                }
                else
                {
                    exit_uni.getData()->set_plus_visited();
                }
            }
            else
            {
                u.getData()->set_minus_visited();
                if (exit_uni.strand == true)
                {
                    exit_uni.getData()->set_minus_visited();
                }
                else
                {
                    exit_uni.getData()->set_plus_visited();
                }
            }
        }
    }
    time_t end_time = time(NULL);
    cout << "CDBG::PloidyEstimation():  Cpu time : "
         << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
    cout << "CDBG::PloidyEstimation():  Real time : "
         << (double)difftime(end_time, start_time) << "s" << endl;
    cout << "CDBG::PloidyEstimation(): Alleles in SuperBubbles  :\t"
         << "2 :" << allele[0] << "\t"
         << "3 :" << allele[1] << "\t"
         << "4 :" << allele[2] << "\t"
         << "5 :" << allele[3] << endl;
    bifre.close();
    pentafre.close();
    pentacov.close();
    trifre.close();
    tetrafre.close();
    bicov.close();
    tricov.close();
    tetracov.close();
    outfile.close();
    allfre.close();
    s_var.close();
    const int avg = coreCov / coreNum;
    cout << "CDBG::PloidyEstimation(): Sites' Average Coverage:" << avg << endl;
}

void CDBG::findSuperBubble_multithread_ptr(const string &outpre, const size_t &thr)
{
    if (thr > 1)
    {
        cout << "CDBG::findSuperBubble(): Finding superbubbles" << endl;
        ofstream outfile;
        if (access("PloidyFrost_output", 0))
        {
            system("mkdir ./PloidyFrost_output");
        }
        outfile.open("PloidyFrost_output/" + outpre + "_super_bubble.txt", ios::out | ios::trunc);
        if (outfile.fail())
        {
            cout << "CDBG:: Open super_bubble file error";
            exit(EXIT_FAILURE);
        }
        std::mutex mtx;
        std::mutex mtx_uni;
        vector<thread> workers;
        unitigIterator<MyUnitig> all_iter = cdbg.begin();
        clock_t start_clock = clock();
        double start_time = time(NULL);
        size_t nb_unitig_processed = 0;
        for (size_t t = 0; t < thr; ++t)
        {
            workers.emplace_back([&]
                                 {
                unitigIterator<MyUnitig> iter;
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
                    UnitigMap<MyUnitig> u(*iter);
                    ++nb_unitig_processed;
                    if (nb_unitig_processed % 100000 == 0)
                    {
                        cout << "CDBG::findSuperBubble(): Processed " << nb_unitig_processed << " unitigs " << endl;
                    }
                    u.strand = true;
                    if (u.getSuccessors().cardinality()>1 && u.getData()->get_plus() == NULL)
                    {
                        extractSuperBubble_multithread_ptr(u, mtx);
                    }                        
                    u.strand=false;
                    if (u.getSuccessors().cardinality()>1 && u.getData()->get_minus() == NULL)
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
        cout << "CDBG::findSuperBubble(): There are " << cdbg.size() << " unitigs " << endl;
        for (auto &t : workers)
            t.join();
        time_t end_time = time(NULL);
        cout << "CDBG::findSuperBubble(): Finding superbubbles Cpu time : "
             << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
        cout << "CDBG::findSuperBubble(): Finding superbubbles Real time : "
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
                unitigIterator<MyUnitig> iter;
                size_t bid = 0;
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
                    
                    MyUnitig *data = iter->getData();
                    if (data->is_both_visited())
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
                    if (!data->is_plus_visited())
                    {
                        bid = nb_super_bubble.fetch_add(1);
                        lock_guard<std::mutex> lock(mtx_file);
                        outfile << bid << "\t" << data->get_id() << "\t"
                                << "+"
                                << "\t" << data->get_plus()->get_id() 
                                << "\t" <<(data->isStrict(true)?"1":"0")
                                << "\t" <<(data->isComplex(true)?"1":"0") 
                                <<endl;
                                
                    }
                    if (!data->is_minus_visited())
                    {
                        bid = nb_super_bubble.fetch_add(1);
                        lock_guard<std::mutex> lock(mtx_file);       
                        outfile << bid << "\t" << data->get_id() << "\t"
                            << "-"
                            << "\t" << data->get_minus()->get_id() 
                            << "\t" <<(data->isStrict(false)?"1":"0")
                            << "\t" <<(data->isComplex(false)?"1":"0")
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
        cout << "CDBG::findSuperBubble(): " << nb_super_bubble << "  SuperBubbles Found" << endl;
    }
    else
    {
        findSuperBubble_ptr(outpre);
    }
}
void CDBG::ploidyEstimation_multithread_ptr(const string &outpre, const int &lower, const int &upper, const size_t &thr)
{
    if (thr > 1)
    {
        clock_t start_clock = clock();
        double start_time = time(NULL);
        cout << "CDBG::PloidyEstimation():  Analyzing superbubbles to generate sites' information" << endl;
        atomic<size_t> coreNum(0);
        atomic<size_t> coreCov(0);
        atomic<size_t> allele5(0);
        atomic<size_t> allele4(0);
        atomic<size_t> allele2(0);
        atomic<size_t> allele3(0);
        atomic<size_t> var_count_all(0);
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
        pentafre.open("PloidyFrost_output/" + outpre + "_pentafre.txt", ios::out | ios::trunc);
        tetrafre.open("PloidyFrost_output/" + outpre + "_tetrafre.txt", ios::out | ios::trunc);
        bicov.open("PloidyFrost_output/" + outpre + "_bicov.txt", ios::out | ios::trunc);
        tricov.open("PloidyFrost_output/" + outpre + "_tricov.txt", ios::out | ios::trunc);
        tetracov.open("PloidyFrost_output/" + outpre + "_tetracov.txt", ios::out | ios::trunc);
        pentacov.open("PloidyFrost_output/" + outpre + "_pentacov.txt", ios::out | ios::trunc);
        ofstream s_var;
        s_var.open("PloidyFrost_output/" + outpre + "_alignseq.txt", ios::out | ios::trunc);
        if (!allfre.is_open() || !bifre.is_open() || !trifre.is_open() ||
            !tetrafre.is_open() || !bicov.is_open() || !tricov.is_open() || !tetracov.is_open() || !s_var.is_open() ||
            !pentafre.is_open() || !pentacov.is_open())
        {
            cout << "CDBG:: PloidyEstimation():Open file error" << endl;
            exit(EXIT_FAILURE);
        }
        vector<thread> workers;
        unitigIterator<MyUnitig> all_iter = cdbg.begin();
        mutex mtx_uni;
        mutex mtx_file;
        mutex mtx_bicov;
        mutex mtx_tricov;
        mutex mtx_tetracov;
        mutex mtx_allfre;
        mutex mtx_pentacov;
        mutex mtx_var;
        size_t nb_unitig_processed = 0;

        for (size_t t = 0; t < thr; ++t)
        {
            workers.emplace_back([&]
                                 {
                unitigIterator<MyUnitig> iter;
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
                        cout << "CDBG::PloidyEstimation(): Processed " << nb_unitig_processed << " unitigs " << endl;
                    }
                    UnitigMap<MyUnitig> u(*iter);
                    MyUnitig *u_data = iter->getData();
                    if (u_data->is_both_visited())
                    {
                        lock_guard<std::mutex> lock(mtx_uni);
                        if (all_iter != cdbg.end())
                        {
                            iter = all_iter;
                            all_iter++;
                            continue;
                        }
                        else
                        {
                            return;
                        }
                    }
                    while (!u_data->is_both_visited())
                    {
                        if (u_data->is_plus_visited() != true)
                        {                
                            u.strand = true;
                            if(u_data->isComplex(u.strand))
                            {
                                u_data->set_plus_visited();
                                continue;
                            }
                        }
                        else if (u_data->is_minus_visited() != true)
                        {
                            u.strand = false;
                            if(u_data->isComplex(u.strand))
                            {
                                u_data->set_minus_visited();
                                break;
                            }
                        }
                        else
                        {
                            break;
                        }
                        double sum = 0;
                        bool isStrict = u_data->isStrict(u.strand);
                        UnitigMap<MyUnitig> exit_uni=NULL;
                        
                        pair<double, uint> core = readCov(u);
                        if(isStrict)
                        {                            
                            exit_uni = (*u.getSuccessors().begin()->getSuccessors().begin());
                            if(u.referenceUnitigToString().compare(exit_uni.referenceUnitigToString())<0){
                                if (u.strand)
                                {
                                    u.getData()->set_plus_visited();
                                }
                                else
                                {
                                    u.getData()->set_minus_visited();
                                }
                                continue;
                            }
                            vector<UnitigMap<MyUnitig>> unitig_vec;
                       
                            bool flag=true;
                            vector<double> cov;
                            uint8_t allele_count = 0;
                            for (const auto &uu : u.getSuccessors())
                            {
                                unitig_vec.push_back(uu);
                                pair<double, uint> inside = readCov(uu); 
                                if (inside.second > lower && inside.second < upper)
                                {
                                    cov.push_back(inside.first);
                                    sum += inside.first;
                                    allele_count++;
                                }
                                else
                                {
                                    flag = false;
                                    break;
                                }
                            }
                        if (flag == true)
                        {
                        sortSeq_simple(cov,unitig_vec,0,cov.size()-1);
                        SeqAlign seqalign(match,mismatch,gap);
                        uint variantNum = 0;
                        vector<uint> var_site;
                        vector<uint> snp_pos;
                        vector<uint> indel_pos;
                        vector<uint> indel_len;
                        vector<vector<unsigned short>> partition;
                        vector<string> str_vec;
                
                        for (auto uni : unitig_vec)
                        {
                            str_vec.push_back(uni.mappedSequenceToString());
                        }
                               
                                seqalign.SequenceAlignment(str_vec, snp_pos, indel_pos, partition, indel_len);
                                if (str_vec.size() != 0)
                                {
                                    coreNum.fetch_add(1);
                                    coreCov.fetch_add((size_t)core.first);
                                    stringstream var_info_stream; 
                                    size_t var_count=var_count_all.fetch_add(1);
                                    for (auto s : str_vec)
                                    {
                                        var_info_stream << var_count << "\t" <<1<<"\t"<<u.getData()->get_id()<<"\t"<<exit_uni.getData()->get_id()<< "\t" << s << endl;
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
                                        vector<unsigned short> part = partition[var_site[i]];
                                        unsigned short maxnum = *max_element(part.begin(),part.end());
                                        vector<double> temp_cov(maxnum, 0);
                                        uint var_distance = 0;
                                        if (i == 0)
                                        {
                                            if (i != var_site.size() - 1)
                                            {
                                                var_distance = min((size_t)(var_site[i + 1] - var_site[i] - 1) , u.size );
                                            }
                                            else
                                            {
                                                var_distance = min(u.size , exit_uni.size );
                                            }
                                        }
                                        else if (i == var_site.size() - 1)
                                        {
                                            var_distance = min((size_t)(var_site[i] - var_site[i - 1] - 1) , exit_uni.size );
                                        }
                                        else
                                        {
                                            var_distance = min(var_site[i] - var_site[i - 1] - 1,var_site[i + 1] - var_site[i] - 1);
                                        }
                                        for (int j = 0; j < part.size(); j++)
                                        {
                                            temp_cov[part[j] - 1] += cov[j];
                                        }
                                        for (auto c : temp_cov)
                                        {
                                            cov_info << c << "\t";
                                            fre_info << (double)c / sum << "\n";
                                        }
                                        if (find(indel_pos.begin(), indel_pos.end(), var_site[i]) != indel_pos.end())
                                        {
                                            ++indel;
                                            cov_info << 1 << "\t" << indel_len[indel - 1] << "\t" << var_count << "\t" << var_site.size() << "\t" << var_distance << "\t" << endl;

                                        }
                                        else
                                        {
                                            cov_info << 1 << "\t" << "0\t" << var_count << "\t" << var_site.size() << "\t" << var_distance << "\t"<<endl;
                                        
                                        
                                        }
                                            switch (maxnum)
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
                                        allfre<<bifre_info.str()<<trifre_info.str()<<tetrafre_info.str()<<pentafre_info.str();              
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
                            while(exit_uni.getData()->get_id()!=u.getData()->get_bubble_id(u.strand))
                            {
                                exit_uni=*exit_uni.getSuccessors().begin();                    
                            }                           
                            if(u.referenceUnitigToString().compare(exit_uni.referenceUnitigToString())<0){
                                if (u.strand)
                                {
                                    u.getData()->set_plus_visited();
                                }
                                else
                                {
                                    u.getData()->set_minus_visited();
                                }
                                continue;
                            }
                        
                            
                            vector<string> str_vec;
                            stack<UnitigMap<MyUnitig>> major;
                            stack<UnitigMap<MyUnitig>> minor;
                            string bubbleStr = "";
                            int k = cdbg.getK();
                            minor.push(u);
                            while (!minor.empty())
                            {
                                UnitigMap<MyUnitig> umi = minor.top();
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
                            uint variantNum = 0;
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
                                    var_info_stream << var_count << "\t" <<0<<"\t"<<u.getData()->get_id()<<"\t"<<exit_uni.getData()->get_id()<< "\t" << s << endl;
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
                                uint pre_var_pos = 0;
                                uint indel = 0;
                                for (uint i = 0; i < var_site.size(); i++)
                                {
                                    stringstream cov_info; 
                                    stringstream fre_info; 
                                    uint var_distance = 0;
                                    unsigned short maxnum = *max_element(partition[var_site[i]].begin(),partition[var_site[i]].end());
                                    vector<string> k_length_str(str_vec.size());
                                    vector<set<string>> string_set(maxnum); 
                                    if (i == 0)
                                    {
                                        if (i != var_site.size() - 1)
                                        {
                                            var_distance = min((size_t)(var_site[i + 1] - var_site[i] - 1) , u.size );
                                        }
                                        else
                                        {
                                            var_distance = min(u.size , exit_uni.size );
                                        }
                                    }
                                    else if (i == var_site.size() - 1)
                                    {
                                        var_distance = min((size_t)(var_site[i] - var_site[i - 1] - 1) , exit_uni.size );
                                    }
                                    else
                                    {
                                        var_distance = min(var_site[i] - var_site[i - 1] - 1,var_site[i + 1] - var_site[i] - 1);
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
                                        vector<double> temp_cov(maxnum, 0);
                                        double sum = 0;
                                        bool flag_site_cov = true;
                                        for (int string_i = 0; string_i < string_set.size() && flag_site_cov; string_i++)
                                        {
                                            for (auto s : string_set[string_i])
                                            {
                                                pair<double, bool> cov_res = readCov(s, lower, upper);
                                                if (cov_res.second == false)
                                                {
                                                    flag_site_cov = false;
                                                    break;
                                                }
                    
                                                temp_cov[string_i] += cov_res.first;
                                            }
                                            sum += temp_cov[string_i];
                                        }
                                        if (flag_site_cov == false)
                                        {
                                            continue;
                                        }
                              
                                        for (auto tc : temp_cov)
                                        {
                                            cov_info << tc << "\t";
                                            fre_info << tc / sum << endl;
                                        }
                                        cov_info << 0 << "\t" << indel_len[indel - 1] << "\t" << var_count << "\t" << var_site.size() << "\t" << var_distance << "\t" << endl;
                                        
                                        
                                    }
                                    
                                    else
                                    {
                                        
                                        if (indel > 0)
                                {
                                    
                                    uint indel_len_sum = 0; 
                                    for (int len_i = 0; len_i < indel; len_i++)
                                    {
                                        indel_len_sum += indel_len[len_i];
                                    }
                                    int start=var_site[i]- cdbg.getK() - indel_len_sum + 1;
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
                                        vector<double> temp_cov(maxnum, 0);
                                        double sum = 0;
                                        bool flag_site_cov = true;
                                        for (int string_i = 0; string_i < string_set.size() && flag_site_cov; string_i++)
                                        {
                                            for (auto s : string_set[string_i])
                                            {
                                                pair<double, bool> cov_res = readCov(s, lower, upper);
                                                if (cov_res.second == false)
                                                {
                                                    flag_site_cov = false;
                                                    break;
                                                }
                                                temp_cov[string_i] += cov_res.first;
                                            }
                                            sum += temp_cov[string_i];
                                        }
                                        if (flag_site_cov == false)
                                        {
                                            continue;
                                        }
                                     
                                        for (auto tc : temp_cov)
                                        {
                                            cov_info << tc << "\t";
                                            fre_info << tc / sum << endl;
                                        }
                                        cov_info << 0 << "\t" << "0\t" << var_count << "\t" << var_site.size() << "\t" << var_distance <<"\t"<< endl;
                                        
                                        
                                    }
                                            switch (maxnum)
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
                        if (u.strand==true)
                        {
                            u_data->set_plus_visited();
                            if(exit_uni.strand==true)
                            {
                                exit_uni.getData()->set_minus_visited();
                            }
                            else
                            {
                                exit_uni.getData()->set_plus_visited();
                            }
                        }
                        else
                        {
                            u_data->set_minus_visited();
                            if(exit_uni.strand==true)
                            {
                                exit_uni.getData()->set_minus_visited();
                            }
                            else
                            {
                                exit_uni.getData()->set_plus_visited();
                            }
                        }
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
        for (auto &w : workers)
            w.join();
        time_t end_time = time(NULL);
        cout << "CDBG::PloidyEstimation(): Cpu time : "
             << (double)(clock() - start_clock) / CLOCKS_PER_SEC << "s" << endl;
        cout << "CDBG::PloidyEstimation(): Real time : "
             << (double)difftime(end_time, start_time) << "s" << endl;
        cout << "CDBG::PloidyEstimation(): Alleles in SuperBubbles  :\t"
             << "2 :" << allele2 << "\t"
             << "3 :" << allele3 << "\t"
             << "4 :" << allele4 << "\t"
             << "5 :" << allele5 << endl;
        s_var.close();
        bifre.close();
        trifre.close();
        tetrafre.close();
        bicov.close();
        tricov.close();
        tetracov.close();
        pentacov.close();
        pentafre.close();
        allfre.close();
        const int avg = coreCov / coreNum;
        cout << "CDBG::PloidyEstimation(): Sites' Average Coverage:" << avg << endl;
    }
    else
    {
        ploidyEstimation_ptr(outpre, lower, upper);
    }
}
void CDBG::extractSuperBubble_multithread_ptr(const UnitigMap<MyUnitig> &s, mutex &m)
{
    bool flag_cycle = false;
    bool flag_tip = false;
    bool flag_large_cycle = false;
    vector<UnitigMap<MyUnitig>> vertices_visit;
    pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> p;
    vector<UnitigMap<MyUnitig>> vec_km_seen;
    map<size_t, uint8_t> state_map;
    map<size_t, bool> strand_map;
    unordered_set<UnitigMap<MyUnitig>, UnitigMapHash<MyUnitig>> cycle_unitig_set;
    UnitigMap<MyUnitig> v(s);
    vertices_visit.push_back(v);
    vec_km_seen.push_back(v);
    while (!vertices_visit.empty())
    {
        v = vertices_visit.back();
        vertices_visit.pop_back();
        state_map[v.getData()->get_id()] = 0x01;
        strand_map[v.getData()->get_id()] = v.strand;
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
                if (state_map.find(u.getData()->get_id()) == state_map.end())
                {
                    vec_km_seen.push_back(u);
                    strand_map[u.getData()->get_id()] = u.strand;
                    state_map[u.getData()->get_id()] = 0x02;
                    bool all_predecessor_visited = true;
                    for (const auto &predecessor : u.getPredecessors())
                    {
                        if (state_map.find(predecessor.getData()->get_id()) != state_map.end())
                        {
                            if (state_map[predecessor.getData()->get_id()] != 0x01)
                            {
                                all_predecessor_visited = false;
                            }
                            if (strand_map[predecessor.getData()->get_id()] != predecessor.strand)
                            {
                                flag_cycle = true;
                                cycle_unitig_set.insert(u);
                                cycle_unitig_set.insert(predecessor);
                            }
                        }
                        else
                            all_predecessor_visited = false;
                    }
                    if (all_predecessor_visited)
                        vertices_visit.push_back(u);
                }
                else if (state_map[u.getData()->get_id()] != 0x01)
                {
                    if (strand_map[u.getData()->get_id()] != u.strand)
                    {
                        flag_cycle = true;
                        cycle_unitig_set.insert(u);
                        cycle_unitig_set.insert(v);
                    }
                    state_map[u.getData()->get_id()] = 0x02;
                    bool all_predecessor_visited = true;
                    for (const auto &predecessor : u.getPredecessors())
                    {
                        if (state_map.find(predecessor.getData()->get_id()) != state_map.end())
                        {
                            if (state_map[predecessor.getData()->get_id()] != 0x01)
                            {
                                all_predecessor_visited = false;
                            }
                            if (strand_map[predecessor.getData()->get_id()] != predecessor.strand)
                            {
                                flag_cycle = true;
                                cycle_unitig_set.insert(u);
                                cycle_unitig_set.insert(predecessor);
                            }
                        }
                        else
                            all_predecessor_visited = false;
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
                    if (state_map[cucm.getData()->get_id()] == 0x02)
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
            if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
            {
                MyUnitig *ex = ucm.getData()->get_plus();
                if (ex->get_plus() == ucm.getData())
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->set_plus_self();
            if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
            {
                MyUnitig *ex = ucm.getData()->get_minus();
                if (ex->get_plus() == ucm.getData())
                {
                    ex->set_plus_self();
                }
                else
                {
                    ex->set_minus_self();
                }
            }
            ucm.getData()->set_minus_self();
            ucm.getData()->set_non_super();
        }
        if (s.strand == true)
        {
            s.getData()->set_plus_self();
        }
        else
        {
            s.getData()->set_minus_self();
        }
    }
    return;
}
void CDBG::setNoBubble_multithread_ptr_cycle(const vector<UnitigMap<MyUnitig>> &vec_km_seen, mutex &x, pair<UnitigMap<MyUnitig>, UnitigMap<MyUnitig>> &p)
{
    lock_guard<std::mutex> lock(x);
    for (const auto &ucm : vec_km_seen)
    {

        if (ucm.getData()->get_plus() != NULL && ucm.getData()->get_plus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_plus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_plus_self();

        if (ucm.getData()->get_minus() != NULL && ucm.getData()->get_minus() != ucm.getData())
        {
            MyUnitig *ex = ucm.getData()->get_minus();
            if (ex->get_plus() == ucm.getData())
            {
                ex->set_plus_self();
            }
            else
            {
                ex->set_minus_self();
            }
        }
        ucm.getData()->set_minus_self();

        ucm.getData()->set_non_super();
    }
    if (p.first.strand == true)
    {
        p.first.getData()->set_plus_self();
    }
    else
    {
        p.first.getData()->set_minus_self();
    }
    if (p.second.strand == false)
    {
        p.second.getData()->set_plus_self();
    }
    else
    {
        p.second.getData()->set_minus_self();
    }
}
