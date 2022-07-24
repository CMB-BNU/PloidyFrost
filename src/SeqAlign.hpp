#include <iostream>
#include <vector>
#include <set>
using namespace std;
struct MatrixUnit;
struct AlignUnit;
class SeqAlign
{
public:
    SeqAlign(double &m, double &d, double &g) : MATCH(m), DIS_MATCH(d), GAP(g) {}
    SeqAlign() : MATCH(2), DIS_MATCH(-1), GAP(-3) {}
    vector<string> compareStrPair(const vector<vector<string>> &str_pair, vector<uint> &max_snp_pos, vector<uint> &max_indel_pos, vector<vector<unsigned short>> &max_num_all, vector<uint> &indel_len_vec);
    AlignUnit variantAnalyze(const string &A, const string &B);
    vector<AlignUnit> traceback(vector<vector<MatrixUnit>> matrix, string str1, string str2);
    vector<AlignUnit> needlemanWunch(const string &A, const string &B);
    void SequenceAlignment(vector<string> &str, vector<uint> &max_snp_pos, vector<uint> &max_indel_pos, vector<vector<unsigned short>> &max_num_all, vector<uint> &indel_len_vec);

private:
    double MATCH;
    double DIS_MATCH;
    double GAP;
};
struct MatrixUnit
{
    uint8_t Up = 0;
    uint8_t LeftUp = 0;
    uint8_t Left = 0;
    long score = 0;
};
struct AlignUnit
{
    string str1 = "";
    string str2 = "";
    vector<uint> gap_pos;
    vector<pair<uint8_t, pair<string, string>>> str_pairs;
    vector<pair<string, string>> str_vec;
    long score = 0;
    vector<uint> pos;
    uint snp = 0;
    uint indel = 0;
    int tag;
    size_t min_distance = 0;
    long operator-(const AlignUnit &x)
    {
        if (score == x.score)
        {
            if (pos.size() == x.pos.size())
            {
                if (indel == x.indel)
                {
                    return 0;
                }
                else
                {
                    return long(x.indel) - indel;
                }
            }
            else
            {
                return long(x.pos.size()) - pos.size();
            }
        }
        else
        {
            return score > x.score ? 1 : -1;
        }
    }
};
