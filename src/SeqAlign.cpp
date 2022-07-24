#include "SeqAlign.hpp"
#include <climits>
#include <string.h>
#include <algorithm>
typedef unsigned int uint;
using namespace std;

vector<string> SeqAlign::compareStrPair(const vector<vector<string>> &str_pair, vector<uint> &max_snp_pos, vector<uint> &max_indel_pos, vector<vector<unsigned short>> &max_num_all, vector<uint> &indel_len_vec)
{
    auto compute_dis = [=](vector<uint> v) -> size_t
    {
        size_t count = 0;
        if (v.size() > 0)
        {
            if (v.size() == 1)
            {
                int left = int(v[0]);
                int right = int(str_pair.back().back().length() - v[0]) - 1;
                if (left > right)
                {
                    left += 1;
                    count = left;
                }
                else
                {
                    count = right;
                }
            }
            else
            {
                count = v[0];
                for (int i = 1; i < v.size(); i++)
                    count = min(int(v[i] - v[i - 1] - 1), int(count));
                count = min(count, str_pair.back().back().length() - v[v.size() - 1] - 1);
            }
        }
        return count;
    };
    vector<string> max_pair;
    int snp_dis = INT_MAX;
    int indel_dis = INT_MAX;
    int snp_count = INT_MAX / 2;
    int indel_count = INT_MAX / 2;
    int all_dis = INT_MAX;
    int site_l = -1;
    int site_r = -1;
    for (int i = 0; i < str_pair.size(); i++)
    {
        vector<uint> snp_pos;
        vector<uint> indel_pos;
        vector<uint> indel_len;
        vector<vector<unsigned short>> num_all;
        bool INDEL = false;
        uint8_t indel = 0;
        uint8_t snp = 0;
        for (int j = 0; j < str_pair[i].back().length(); j++)
        {
            set<char> char_set;
            vector<unsigned short> num(str_pair[i].size(), 0);
            for (int k = 0; k < str_pair[i].size(); k++)
            {
                char_set.insert(str_pair[i][k][j]);
            }
            if (char_set.size() > 1)
            {
                if (char_set.find('-') == char_set.end())
                {
                    if (INDEL == true)
                    {
                        indel_len.push_back(j - indel_pos[indel - 1]);
                        INDEL = false;
                    }
                    snp_pos.push_back(j);
                    snp++;
                    int count_snp = 0;
                    for (int ki = 0; ki < str_pair[i].size(); ki++)
                    {
                        bool is_same = false;
                        for (int kj = 0; kj < ki; kj++)
                        {
                            if (str_pair[i][kj][j] == str_pair[i][ki][j])
                            {
                                is_same = true;
                                num[ki] = num[kj];
                                break;
                            }
                        }
                        if (is_same == false)
                        {
                            num[ki] = ++count_snp;
                        }
                    }
                }
                else
                {
                    bool old_indel = true;
                    if (INDEL == true)
                    {
                        for (int kj = 0; kj < str_pair[i].size(); kj++)
                        {
                            if ((str_pair[i][kj][j] == '-' && str_pair[i][kj][j - 1] != '-') || (str_pair[i][kj][j] != '-' && str_pair[i][kj][j - 1] == '-'))
                            {
                                old_indel = false;
                                break;
                            }
                        }
                        if (old_indel == false)
                        {
                            indel_len.push_back(j - indel_pos[indel - 1]);
                            ++indel;
                            indel_pos.push_back(j);
                        }
                    }
                    else
                    {
                        old_indel = false;
                        ++indel;
                        indel_pos.push_back(j);
                        INDEL = true;
                    }
                    if (old_indel == false || char_set.size() > 2)
                    {
                        int count_char = 0;
                        for (int ki = 0; ki < str_pair[i].size(); ki++)
                        {
                            bool is_same = false;
                            for (int kii = 0; kii < ki; kii++)
                            {
                                if (str_pair[i][kii][j] == str_pair[i][ki][j])
                                {
                                    is_same = true;
                                    num[ki] = num[kii];
                                    break;
                                }
                                else
                                {
                                    is_same = false;
                                }
                            }
                            if (is_same == false)
                            {
                                num[ki] = ++count_char;
                            }
                        }
                    }
                }
            }
            else
            {
                if (INDEL == true)
                {
                    indel_len.push_back(j - indel_pos[indel - 1]);
                    INDEL = false;
                }
            }
            num_all.push_back(num);
        }
        bool flag = false;
        if (snp + indel < snp_count + indel_count)
            flag = true;
        else if (snp + indel == snp_count + indel_count)
            if (indel < indel_count)
                flag = true;
            else if (indel == indel_count)
            {
                size_t now_indel_dis = compute_dis(indel_pos);
                if (now_indel_dis > indel_dis)
                    flag = true;
                else if (now_indel_dis == indel_dis)
                {
                    size_t now_snp_dis = compute_dis(snp_pos);
                    if (now_snp_dis > snp_dis)
                        flag = true;
                    else if (now_snp_dis == snp_dis)
                    {
                        vector<uint> temp_vec;
                        temp_vec.resize(snp_pos.size() + indel_pos.size());
                        merge(snp_pos.begin(), snp_pos.end(), indel_pos.begin(), indel_pos.end(), temp_vec.begin());
                        size_t now_all_dis = compute_dis(temp_vec);
                        if (now_all_dis > all_dis)
                            flag = true;
                        else if (now_all_dis == all_dis)
                        {
                            int now_site_l = int(temp_vec[0]);
                            int now_site_r = int(temp_vec.back());
                            if (now_site_l > site_l || now_site_r > site_r)
                            {
                                flag = true;
                            }
                            else if (now_site_l == site_l && now_site_r == site_r)
                            {
                                for (int m = 0; m < str_pair[i].size(); m++)
                                {
                                    if (strcmp(str_pair[i][m].c_str(), max_pair[m].c_str()) > 0)
                                    {
                                        all_dis = now_all_dis;
                                        site_l = now_site_l;
                                        site_r = now_site_r;
                                        snp_count = snp;
                                        indel_count = indel;
                                        snp_dis = now_snp_dis;
                                        indel_dis = now_indel_dis;
                                        max_pair = str_pair[i];
                                        max_snp_pos = snp_pos;
                                        max_indel_pos = indel_pos;
                                        max_num_all = num_all;
                                        indel_len_vec = indel_len;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        if (flag == true)
        {
            vector<uint> temp_vec;
            temp_vec.resize(snp_pos.size() + indel_pos.size());
            merge(snp_pos.begin(), snp_pos.end(), indel_pos.begin(), indel_pos.end(), temp_vec.begin());
            all_dis = compute_dis(temp_vec);
            site_l = max(int(site_l), int(temp_vec.size() > 0 ? temp_vec[0] : -1));
            site_r = max(int(site_r), int(temp_vec.size() > 0 ? temp_vec.back() : -1));
            snp_count = snp;
            indel_count = indel;
            snp_dis = compute_dis(snp_pos);
            indel_dis = compute_dis(indel_pos);
            max_pair = str_pair[i];
            max_snp_pos = snp_pos;
            max_indel_pos = indel_pos;
            max_num_all = num_all;
            indel_len_vec = indel_len;
        }
    }
    return max_pair;
}
AlignUnit SeqAlign::variantAnalyze(const string &A, const string &B)
{
    auto score_func = [&](char a, char b) -> double
    {
        if (a == '-' || b == '-')
            return GAP;
        else if (a == b)
            return MATCH;
        else
            return DIS_MATCH;
    };
    AlignUnit au;
    au.score = 0;
    au.str1 = A;
    au.str2 = B;
    uint8_t flag = 0;
    for (size_t i = 0; i < A.length(); i++)
    {
        au.score += score_func(A[i], B[i]);
        if (A[i] != B[i])
        {
            if (A[i] == '-')
            {
                if (flag != 1)
                {
                    flag = 1;
                    au.indel++;
                    au.pos.push_back(i);
                }
            }
            else if (B[i] == '-')
            {
                if (flag != 2)
                {
                    flag = 2;
                    au.indel++;
                    au.pos.push_back(i);
                }
            }
            else
            {
                au.snp++;
                flag = 0;
                au.pos.push_back(i);
            }
        }
        else
        {
            flag = 0;
        }
    }
    if (au.pos.size() > 0)
    {
        if (au.pos.size() == 1)
        {
            au.min_distance = min(int(au.pos[0]), int(A.length() - au.pos[0] - 1));
        }
        else
        {
            au.min_distance = int(au.pos[0]);
            for (int i = au.pos.size() - 1; i >= 1; i--)
            {
                au.min_distance = min(int(au.pos[i] - au.pos[i - 1] - 1), int(au.min_distance));
            }
            au.min_distance = min(int(A.length() - au.pos[0] - 1), int(au.min_distance));
        }
    }
    return au;
}
vector<AlignUnit> SeqAlign::traceback(vector<vector<MatrixUnit>> matrix, string str1, string str2)
{
    vector<AlignUnit> au_vec;
    vector<pair<size_t, size_t>> stack;
    string resA = "";
    string resB = "";
    size_t indel1 = 0;
    size_t indel2 = 0;
    size_t indel1_max = 5;
    size_t indel2_max = 5;
    stack.push_back(pair<size_t, size_t>(str1.length(), str2.length()));
    vector<vector<MatrixUnit>> matrix_temp = matrix;
    vector<unsigned int> gap_pos;
    while (stack.size())
    {
        pair<size_t, size_t> p = stack.back();
        if (p.first == 0 && p.second == 0 && indel1 <= indel1_max && indel2 <= indel2_max)
        {
            string res_temp = resA;
            for (int j = 0; j < gap_pos.size(); j++)
            {
                res_temp[gap_pos[j] + gap_pos.size() - j - 1] = '-';
            }
            AlignUnit au = variantAnalyze(res_temp, resB);
            au.gap_pos = gap_pos;
            if (au_vec.size() > 0)
            {
                AlignUnit au_back = au_vec.back();
                int diff = au_back - au;
                if (diff == 0)
                {
                    au_vec.push_back(au);
                    indel1_max = indel1;
                    indel2_max = indel2;
                }
                else if (diff < 0)
                {
                    au_vec.clear();
                    au_vec.push_back(au);
                    indel1_max = indel1;
                    indel2_max = indel2;
                }
            }
            else
            {
                au_vec.push_back(au);
                indel1_max = indel1;
                indel2_max = indel2;
            }
        }
        if (matrix_temp[p.first][p.second].Left)
        {
            if (indel1 < indel1_max)
            {
                if (resA == "" || resA[0] != '+')
                {
                    ++indel1;
                }
                stack.push_back(pair<size_t, size_t>(p.first, p.second - 1));
                resA = "+" + resA;
                gap_pos.push_back(p.first);
                resB = str2[p.second - 1] + resB;
            }
            else if (indel1 == indel1_max)
            {
                if (resA[0] != '+')
                {
                    matrix[p.first][p.second].Left = 0;
                    matrix_temp[p.first][p.second].Left = 0;
                    continue;
                }
                else
                {
                    stack.push_back(pair<size_t, size_t>(p.first, p.second - 1));
                    resA = "+" + resA;
                    gap_pos.push_back(p.first);
                    resB = str2[p.second - 1] + resB;
                }
            }
            else
            {
                matrix[p.first][p.second].Left = 0;
                matrix_temp[p.first][p.second].Left = 0;
                continue;
            }
            matrix_temp[p.first][p.second].Left = 0;
        }
        else if (matrix_temp[p.first][p.second].Up)
        {
            if (indel2 < indel2_max)
            {
                if (resB == "" || resB[0] == '-')
                {
                    ++indel2;
                }
                stack.push_back(pair<size_t, size_t>(p.first - 1, p.second));
                resA = str1[p.first - 1] + resA;
                resB = "-" + resB;
            }
            else if (indel2 == indel2_max)
            {
                if (resB[0] != '-')
                {
                    matrix_temp[p.first][p.second].Up = 0;
                    matrix[p.first][p.second].Up = 0;
                    continue;
                }
                stack.push_back(pair<size_t, size_t>(p.first - 1, p.second));
                resA = str1[p.first - 1] + resA;
                resB = "-" + resB;
            }
            else
            {
                matrix_temp[p.first][p.second].Up = 0;
                matrix[p.first][p.second].Up = 0;
                continue;
            }
            matrix_temp[p.first][p.second].Up = 0;
        }
        else if (matrix_temp[p.first][p.second].LeftUp)
        {
            stack.push_back(pair<size_t, size_t>(p.first - 1, p.second - 1));
            resA = str1[p.first - 1] + resA;
            resB = str2[p.second - 1] + resB;
            matrix_temp[p.first][p.second].LeftUp = 0;
        }
        else
        {
            if (resA.length() == 0)
            {
                break;
            }
            stack.pop_back();
            matrix_temp[p.first][p.second] = matrix[p.first][p.second];
            if (resA[0] == '+')
            {
                if (resA.length() >= 2)
                {
                    if (resA[1] != '+')
                    {
                        --indel1;
                    }
                }
                else
                {
                    --indel1;
                }
            }
            if (resB[0] == '-')
            {
                if (resB.length() >= 2)
                {
                    if (resB[1] != '-')
                    {
                        --indel2;
                    }
                }
                else
                {
                    --indel2;
                }
            }
            if (resA[0] == '+')
            {
                gap_pos.pop_back();
            }
            resA = resA.substr(1, resA.length() - 1);
            resB = resB.substr(1, resB.length() - 1);
        }
    }

    return au_vec;
}

vector<AlignUnit> SeqAlign::needlemanWunch(const string &A, const string &B)
{
    size_t m = A.length();
    size_t n = B.length();
    vector<vector<MatrixUnit>> matrix(m + 1);
    matrix[0].resize(n + 1);
    for (size_t i = 1; i <= m; i++)
    {
        matrix[i].resize(n + 1);
        matrix[i][0].score = GAP * i;
        matrix[i][0].Up = 1;
    }
    for (size_t j = 1; j <= n; j++)
    {
        matrix[0][j].score = GAP * j;
        matrix[0][j].Left = 1;
    }
    int up_score, leftup_score, left_score, max_score;
    auto score_func = [&](char a, char b) -> double
    {
        if (a == b)
            return MATCH;
        else if (a == '-' || b == '-')
            return GAP;
        else
            return DIS_MATCH;
    };

    for (size_t i = 1; i <= m; i++)
    {
        for (size_t j = 1; j <= n; j++)
        {
            up_score = matrix[i - 1][j].score + GAP;
            if (matrix[i - 1][j].Up == 1)
            {
                up_score += 1;
            }
            leftup_score = matrix[i - 1][j - 1].score + score_func(A[i - 1], B[j - 1]);
            if (matrix[i - 1][j - 1].LeftUp == 1)
            {
                leftup_score += 1;
            }
            left_score = matrix[i][j - 1].score + GAP;
            if (matrix[i][j - 1].Left == 1)
            {
                left_score += 1;
            }
            max_score = max(max(up_score, leftup_score), left_score);
            if (max_score == left_score && i != m && A[i] == '-')
            {
                left_score = INT_MIN;
                max_score = up_score > leftup_score ? up_score : leftup_score;
            }
            matrix[i][j].score = max_score;
            if (up_score == max_score)
            {
                matrix[i][j].Up = 1;
            }
            if (leftup_score == max_score)
            {
                matrix[i][j].LeftUp = 1;
            }
            if (left_score == max_score)
            {
                matrix[i][j].Left = 1;
            }
        }
    }
    return traceback(matrix, A, B);
}
void SeqAlign::SequenceAlignment(vector<string> &str, vector<uint> &max_snp_pos, vector<uint> &max_indel_pos, vector<vector<unsigned short>> &max_num_all, vector<uint> &indel_len_vec)
{
    vector<AlignUnit> align_vec = needlemanWunch(str[0], str[1]);
    vector<vector<string>> str_pairs(align_vec.size());
    for (int i = 0; i < str_pairs.size(); i++)
    {
        str_pairs[i].push_back(align_vec[i].str1);
        str_pairs[i].push_back(align_vec[i].str2);
    }
    for (int i = 2; i < str.size(); i++)
    {
        vector<vector<string>> temp_pairs = str_pairs;
        str_pairs.clear();
        int max_score = INT_MIN;
        for (int k = 0; k < temp_pairs.size(); k++)
        {
            int max_score_k = 0;
            vector<AlignUnit> align_temp = needlemanWunch(temp_pairs[k][0], str[i]);
            vector<vector<string>> str_pair_vec_all(align_temp.size());
            vector<int> valid_au_pos;
            for (int vp = 0; vp < align_temp.size(); vp++)
            {
                valid_au_pos.push_back(vp);
                str_pair_vec_all[vp].push_back(align_temp[vp].str1);
            }
            for (int j = 1; j < i; j++)
            {
                int max_score_j = INT_MIN;
                AlignUnit au_max;
                au_max.score = INT_MIN;
                vector<int> valid_au_pos_j;
                for (int c = 0; c < valid_au_pos.size(); c++)
                {
                    uint pre = 0;
                    string temp_str = "";
                    if (align_temp[valid_au_pos[c]].gap_pos.size())
                    {
                        for (int s = align_temp[valid_au_pos[c]].gap_pos.size() - 1; s >= 0; s--)
                        {
                            temp_str = temp_str + temp_pairs[k][j].substr(pre, align_temp[valid_au_pos[c]].gap_pos[s] - pre) + "-";
                            pre = align_temp[valid_au_pos[c]].gap_pos[s];
                        }
                        temp_str = temp_str + temp_pairs[k][j].substr(pre);
                    }
                    else
                    {
                        temp_str = temp_pairs[k][j];
                    }
                    AlignUnit au = variantAnalyze(temp_str, align_temp[valid_au_pos[c]].str2);
                    int diff = au - au_max;
                    if (diff > 0)
                    {
                        au_max = au;
                        max_score_j = au_max.score;
                        valid_au_pos_j.clear();
                        valid_au_pos_j.push_back(valid_au_pos[c]);
                        str_pair_vec_all[valid_au_pos[c]].push_back(temp_str);
                    }
                    else if (diff == 0)
                    {
                        max_score_j = au_max.score;
                        valid_au_pos_j.push_back(valid_au_pos[c]);
                        str_pair_vec_all[valid_au_pos[c]].push_back(temp_str);
                    }
                }
                valid_au_pos.clear();
                valid_au_pos = valid_au_pos_j;
                max_score_k += max_score_j;
            }
            if (max_score_k > max_score)
            {
                max_score = max_score_k;
                str_pairs.clear();
                for (int vp = 0; vp < valid_au_pos.size(); vp++)
                {
                    str_pair_vec_all[valid_au_pos[vp]].push_back(align_temp[valid_au_pos[vp]].str2);
                    str_pairs.push_back(str_pair_vec_all[valid_au_pos[vp]]);
                }
            }
            else if (max_score_k == max_score)
            {
                for (int vp = 0; vp < valid_au_pos.size(); vp++)
                {
                    str_pair_vec_all[valid_au_pos[vp]].push_back(align_temp[valid_au_pos[vp]].str2);
                    str_pairs.push_back(str_pair_vec_all[valid_au_pos[vp]]);
                }
            }
        }
    }
    str = compareStrPair(str_pairs, max_snp_pos, max_indel_pos, max_num_all, indel_len_vec);
}
