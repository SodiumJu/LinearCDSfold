#include <algorithm>
#include <cstring>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <sys/time.h>
#include <stack>
#include <set>
#include <string>
#include <tuple>
#include "utility.h"
#include <unordered_map>
#include <vector>

using namespace std;

#define VALUE_MIN std::numeric_limits<float>::lowest()

const float skip_sc = 0; 
const int gap = 3;
float beta = 1;
float lambda = 1000;

#define NOTON 7 // NUM_OF_TYPE_OF_NUCS
#define gap 3 // NUM_OF_SMALLEST_DIS_CAN_PAIR
// #define timer

bool _allowed_pairs[NOTON][NOTON];
bool xcon_table[3][2][NOTON][NOTON];

// #define BASE(x) ((x=='A'? 0 : (x=='C'? 1 : (x=='G'? 2 : (x=='U'?3: (x=='u'?4:(x=='g'?5:6)))))))
// #define reBASE(x) ((x==0? 'A' : (x==1? 'C' : (x==2? 'G' : (x==3?'U': (x==4?'u':(x==5?'g':'*')))))))

enum Manner {
    NONE=0,
    MANNER_NtoS,
    MANNER_PtoPS,
    MANNER_M1_M1toM2,
    MANNER_P_StoM1,
    MANNER_PtoM1,
    MANNER_M2toP,
    MANNER_S_M2toP,
    MANNER_S_PStoP,
    MANNER_PStoP,
    MANNER_HtoP,
    MANNER_P_StoPS,
    MANNER_StoH,
    MANNER_StoS,
    MANNER_PStoC,
    MANNER_C_PStoC,
    MANNER_C_StoC,
    MANNER_M2_M1toM2,
    MANNER_StoC,
    MANNER_PtoP
};

std::unordered_map<int, std::string> reManner = {
        {0,"NONE"},
        {1,"MANNER_NtoS"},
        {2,"MANNER_PtoPS"},
        {3,"MANNER_M1_M1toM2"},
        {4,"MANNER_P_StoM1"},
        {5,"MANNER_PtoM1"},
        {6,"MANNER_M2toP"},
        {7,"MANNER_S_M2toP"},
        {8,"MANNER_S_PStoP"},
        {9,"MANNER_PStoP"},
        {10,"MANNER_HtoP"},
        {11,"MANNER_P_StoPS"},
        {12,"MANNER_StoH"},
        {13,"MANNER_StoS"},
        {14,"MANNER_PStoC"},
        {15,"MANNER_C_PStoC"},
        {16,"MANNER_C_StoC"},
        {17,"MANNER_M2_M1toM2"},
        {18,"MANNER_StoC"},
        {19,"MANNER_PtoP"}

    };

struct State {
    float score=VALUE_MIN;
    int index_1=-1;
    int index_2=-1;
    Manner MANNER=NONE;
    vector<int> start_codon = {-1,-1};
    vector<int> end_codon = {-1,-1};
    int length = 0;
    float CAI = 0;
};

struct State_pkg {

    std::vector<std::unordered_map<int, State*>> bestS; // for skipping
    std::vector<std::unordered_map<int, State*>> bestC;
    std::vector<std::unordered_map<int, State*>> bestH;
    std::vector<std::unordered_map<int, State*>> bestP;
    std::vector<std::unordered_map<int, State*>> bestPS;
    std::vector<std::unordered_map<int, State*>> bestM1;
    std::vector<std::unordered_map<int, State*>> bestM2;

};



#define BASE(x) ((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: (x=='u'?5:(x=='g'?6:0)))))))
#define reBASE(x) ((x==1? 'A' : (x==2? 'C' : (x==3? 'G' : (x==4?'U': (x==5?'u':(x==6?'g':'*')))))))

class InputParser{
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

void initialize(){
    _allowed_pairs[BASE('A')][BASE('U')] = true;
    _allowed_pairs[BASE('A')][BASE('u')] = true;
    _allowed_pairs[BASE('U')][BASE('A')] = true;
    _allowed_pairs[BASE('u')][BASE('A')] = true;
    _allowed_pairs[BASE('C')][BASE('G')] = true;
    _allowed_pairs[BASE('G')][BASE('C')] = true;
    _allowed_pairs[BASE('C')][BASE('g')] = true;
    _allowed_pairs[BASE('g')][BASE('C')] = true;
    _allowed_pairs[BASE('G')][BASE('U')] = true;
    _allowed_pairs[BASE('g')][BASE('U')] = true;
    _allowed_pairs[BASE('G')][BASE('u')] = true;
    _allowed_pairs[BASE('g')][BASE('u')] = true;
    // _allowed_pairs[BASE('G')][BASE('U')] = true;
    // _allowed_pairs[BASE('U')][BASE('G')] = true;

    xcon_table[0][0][BASE('C')][BASE('G')] = true; //R
    xcon_table[0][0][BASE('C')][BASE('g')] = true; //R
    xcon_table[0][1][BASE('G')][BASE('C')] = true; //R
    xcon_table[0][1][BASE('G')][BASE('U')] = true; //R
    xcon_table[0][0][BASE('A')][BASE('g')] = true; //R
    xcon_table[0][1][BASE('g')][BASE('A')] = true; //R
    xcon_table[0][1][BASE('g')][BASE('G')] = true; //R

    xcon_table[1][0][BASE('C')][BASE('U')] = true; //L
    xcon_table[1][0][BASE('C')][BASE('u')] = true; //L
    xcon_table[1][1][BASE('U')][BASE('C')] = true; //L
    xcon_table[1][1][BASE('U')][BASE('U')] = true; //L
    xcon_table[1][0][BASE('U')][BASE('u')] = true; //L
    xcon_table[1][1][BASE('u')][BASE('A')] = true; //L
    xcon_table[1][1][BASE('u')][BASE('G')] = true; //L

    xcon_table[2][0][BASE('U')][BASE('C')] = true; //S
    xcon_table[2][0][BASE('A')][BASE('G')] = true; //S
    xcon_table[2][1][BASE('C')][BASE('U')] = true; //S
    xcon_table[2][1][BASE('C')][BASE('C')] = true; //S
    xcon_table[2][1][BASE('C')][BASE('A')] = true; //S
    xcon_table[2][1][BASE('C')][BASE('G')] = true; //S
    xcon_table[2][1][BASE('G')][BASE('U')] = true; //S
    xcon_table[2][1][BASE('G')][BASE('C')] = true; //S
    //1:R
    //xcon_table[1]
}

vector<string> codonset={"AUU","AUA","AUC","CUA","CUC","CUG","CUU","UUA","UUG","GUU","GUA","GUC","GUG","UUU","UUC","AUG","UGU","UGC","GCA","GCC","GCG","GCU","GGU","GGC","GGA","GGG","CCU","CCC","CCA","CCG","ACU","ACC","ACA","ACG","UCU","UCC","UCA","UCG","AGU","AGC","UAU","UAC","UGG","CAA","CAG","AAU","AAC","CAU","CAC","GAA","GAG","GAU","GAC","AAA","AAG","CGU","CGC","CGA","CGG","AGA","AGG","UAA","UAG","UGA"};

std::unordered_map<int, float> codonset_CAI;

inline int cal_CAI_index(int base_0 , int base_1, int base_2){
    return base_0*100 + base_1*10 + base_2;
}

inline void initialize_CAI_table(vector<float> cai_vector){
    for (int i = 0; i<codonset.size(); i++){
        // int tmp = BASE(codonset[i][0])*100+BASE(codonset[i][1])*10+BASE(codonset[i][2]);
        if(codonset[i] == "CUA" || codonset[i] == "CUG" || codonset[i] == "UUA" || codonset[i] == "UUG"){
            int tmp = cal_CAI_index(BASE(codonset[i][0]), BASE('u'), BASE(codonset[i][2]));
            codonset_CAI[tmp]=cai_vector[i];
        }else if(codonset[i] == "CGA" || codonset[i] == "CGG" || codonset[i] == "AGA" || codonset[i] == "AGG"){
            int tmp = cal_CAI_index(BASE(codonset[i][0]), BASE('g'), BASE(codonset[i][2]));
            codonset_CAI[tmp]=cai_vector[i];
        }

        int tmp = cal_CAI_index(BASE(codonset[i][0]), BASE(codonset[i][1]), BASE(codonset[i][2]));

        codonset_CAI[tmp]=cai_vector[i];

    }    
}

std::unordered_map<std::string, std::string>  Amino_label = {
        {"A","GCN"},{"C","UGY"},{"D","GAY"},{"E","GAR"},{"F","UUY"},
        {"G","GGN"},{"H","CAY"},{"I","AUH"},{"K","AAR"},{"L","YVN"},
        {"M","AUG"},{"N","AAY"},{"P","CCN"},{"Q","CAR"},{"R","MON"},
        {"S","WSN"},{"T","ACN"},{"V","GUN"},{"W","UGG"},{"Y","UAY"},
        {"*","NNN"}
    };

std::unordered_map<std::string, std::string>  re_Amino_label = {
        {"GCU","A"}, {"GCC","A"}, {"GCA","A"}, {"GCG","A"},
        {"CGU","R"}, {"CGC","R"}, {"CGA","R"}, {"CGG","R"}, {"AGA","R"}, {"AGG","R"},
        {"AAU","N"}, {"AAC","N"},
        {"GAU","D"}, {"GAC","D"},
        {"UGU","C"}, {"UGC","C"},
        {"CAA","Q"}, {"CAG","Q"},
        {"GAA","E"}, {"GAG","E"},
        {"GGU","G"}, {"GGC","G"}, {"GGA","G"}, {"GGG","G"},
        {"CAU","H"}, {"CAC","H"},
        {"AUU","I"}, {"AUC","I"}, {"AUA","I"},
        {"CUU","L"}, {"CUC","L"}, {"CUA","L"}, {"CUG","L"}, {"UUA","L"}, {"UUG","L"},
        {"AAA","K"}, {"AAG","K"},
        {"AUG","M"},
        {"UUU","F"}, {"UUC","F"},
        {"CCU","P"}, {"CCC","P"}, {"CCA","P"}, {"CCG","P"},
        {"UCU","S"}, {"UCC","S"}, {"UCA","S"}, {"UCG","S"}, {"AGU","S"}, {"AGC","S"},
        {"ACU","T"}, {"ACC","T"}, {"ACA","T"}, {"ACG","T"},
        {"UGG","W"},
        {"UAU","Y"}, {"UAC","Y"},
        {"GUU","V"}, {"GUC","V"}, {"GUA","V"}, {"GUG","V"}

    };


std::unordered_map<char, vector<int>> Base_table={
        {'N',{BASE('A'),BASE('C'),BASE('G'),BASE('U')}},
        {'A',{BASE('A')}},
        {'U',{BASE('U')}},
        {'G',{BASE('G')}},
        {'C',{BASE('C')}},
        {'W',{BASE('A'),BASE('U')}},
        {'S',{BASE('G'),BASE('C')}},
        {'M',{BASE('A'),BASE('C')}},
        // {'K',{BASE('G'),BASE('U')}},
        {'R',{BASE('A'),BASE('G')}},
        {'Y',{BASE('C'),BASE('U')}},
        {'H',{BASE('A'),BASE('C'),BASE('U')}},
        {'V',{BASE('U'),BASE('u')}},
        {'O',{BASE('G'),BASE('g')}}
    };



std::vector<std::pair<float, int>> scores;


inline int cal_base_index(int i, int base_i, int base_j){
    return i*100+base_i*10+base_j;
}



//...[i][j]...
bool val_con(vector<int>& con_seq, int base_i, int base_j, int pos_i){
    //j must be i+1
    int pos_j=pos_i+1;
    if((pos_i/3)==(pos_j/3)){
        if((con_seq[pos_i] == con_seq[pos_j]) & (con_seq[pos_i]!=3)){
            int ami_index = pos_i%3;
            if(xcon_table[con_seq[pos_i]][ami_index][base_i][base_j]){
                return true;
            }else{
                return false;
            }
        }else{
            return true;
        }   
    }else{
        return true;
    }
    
}

unsigned long quickselect_partition(vector<pair<float, int>>& scores, unsigned long lower, unsigned long upper) {
    float pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);

    }
    return upper;
}


// in-place quick-select
float quickselect(vector<pair<float, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}

float beam_prune(int beam, vector<int>& con_seq, string& rna_seq, std::unordered_map<int, State*> &beamstep, std::vector<std::unordered_map<int, State*>>& bestH, bool skip) {
    scores.clear();

    if (beamstep.size() <= beam) return VALUE_MIN;

    for (auto &item : beamstep) {
        int ind = item.first;
        int i = ind/100;
        int base_i = (ind%100)/10;
        int base_j = ind%10;

        State* cand = item.second;
        int k = i - 1;

        float newscore;

        if(skip == false){

            if(k>=0){
                vector<int> k_base_list=Base_table.at(rna_seq[k]);
                vector<int> f_base_list=Base_table.at(rna_seq[0]);

                int base_f;
                int base_k;

                int base_index;

                
                float tmp_newscore = VALUE_MIN;
                float best_newscore = VALUE_MIN;
                bool have_previous=false;
                // if(skip==false){
                for(int fi=0; fi<f_base_list.size(); fi++){
                    for(int ck=0; ck<k_base_list.size(); ck++){
                        base_f=f_base_list[fi];
                        base_k=k_base_list[ck];
                        base_index = cal_base_index(0,base_f,base_k);
                        if(val_con(con_seq, base_k, base_i, k) & bestH[k].count(base_index)){
                            //tmp_new_score = bestH[k][base_index]->score + cand->score;

                            tmp_newscore = bestH[k][base_index]->score;
                            have_previous = true;
                            if(tmp_newscore > best_newscore){
                                best_newscore = tmp_newscore;
                            }
                        }
                    }
                }


                // }else{
                //     have_previous = false;
                // }
                
                if(have_previous){
                    newscore = best_newscore + cand->score;
                }else{
                    newscore = cand->score;
                }
            }else{
                newscore = cand->score;
            }
        }else{
            newscore = i;
        }
        // cout << newscore << " (" << cand->score << ")" << " ";
        scores.push_back(make_pair(newscore, ind));

        // lisiz: for _V, avoid -inf-int=+inf
        // if ((k >= 0) && (bestC[k].score == VALUE_MIN)) newscore = VALUE_MIN;
        // else newscore = (k >= 0 ? bestC[k].score : 0) + cand.score;
        
    }
    // cout << endl;
    
    float threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void ami_to_rna(vector<string>& rna_seq_list,string& ami_seq){
    string rna_seq;
    string Amino1;
    string codon;
    for(int i=0;i<ami_seq.size();i++){
        Amino1 = ami_seq[i];
        codon=Amino_label.at(Amino1);
        for (int j=0; j<codon.size();j++){
            rna_seq.push_back(codon[j]);
        }       
    }
    rna_seq_list.push_back(rna_seq);
}

void add_con_seq(vector<vector<int>>& con_seq_list,string& ami_seq){
    vector<int> con_seq;
    char Amino1; 
    for(int i=0; i<ami_seq.size(); i++){
        Amino1 = ami_seq[i];
        if(Amino1=='R'){
            for(int j=0; j<3; j++)
                con_seq.push_back(0); //R==0
        }else if(Amino1=='L'){
            for(int j=0; j<3; j++)
                con_seq.push_back(1); //L==1
        }else if(Amino1=='S'){
            for(int j=0; j<3; j++)
                con_seq.push_back(2); //S==2
        }else{
            for(int j=0; j<3; j++)
                con_seq.push_back(3); //3
        }
    }
    con_seq_list.push_back(con_seq);
}

void check_valid_ami(int seq_size, string result_seq,string& ami_seq){
    string sub_str;
    string index;
    for(int i=0; i<seq_size; i++){
        sub_str =  result_seq.substr(3*i,3);
        // cout << sub_str << endl;
        if(re_Amino_label.count(sub_str)){
            index = re_Amino_label.at(sub_str);    
        }else{
            // cout << sub_str << endl;
            break;
        }
        
        // cout << index << endl;
        
        ami_seq.push_back(index[0]);
        
    }
}

int get_parentheses(char* result_seq, char* result, string& rna_seq, State_pkg& state_pkg){

    State* state = new State;
    std::vector<std::unordered_map<int, State*>>& bestS = state_pkg.bestS; // for skipping
    std::vector<std::unordered_map<int, State*>>& bestC = state_pkg.bestC;
    std::vector<std::unordered_map<int, State*>>& bestH = state_pkg.bestH;
    std::vector<std::unordered_map<int, State*>>& bestP = state_pkg.bestP;
    std::vector<std::unordered_map<int, State*>>& bestPS = state_pkg.bestPS;
    std::vector<std::unordered_map<int, State*>>& bestM1 = state_pkg.bestM1;
    std::vector<std::unordered_map<int, State*>>& bestM2 = state_pkg.bestM2;

    int seq_length = rna_seq.size();
    memset(result, '*', seq_length);
    result[seq_length] = 0;
    result_seq[seq_length] = 0;
    stack<tuple<int, int, int, State*>> stk;
    int best_sc=VALUE_MIN;
    int best_index=-1;
    float tmp_sc=0;
    int tmp_index=-1;
    int base_i;
    int base_j;
    vector<int> f_base_list=Base_table.at(rna_seq[0]);//i list
    vector<int> l_base_list=Base_table.at(rna_seq[seq_length-1]);//j list

    for(int ci=0; ci<f_base_list.size(); ci++){
        for(int cj=0; cj<l_base_list.size(); cj++){
            base_i=f_base_list[ci];
            base_j=l_base_list[cj];
            tmp_index=cal_base_index(0, base_i, base_j);
            // cout << "tmp_index:" << tmp_index << endl;
            // cout << ci << endl;
            if(bestC[seq_length-1].count(tmp_index)){
                tmp_sc=bestC[seq_length-1][tmp_index]->score;
                // cout << "tmp_sc:" << tmp_sc << endl;
                if(tmp_sc>best_sc){
                    best_index=tmp_index;
                    best_sc=tmp_sc;
                }
            }
            
        }
    }

    if(best_index!=-1){
        int nuc_index = best_index%100;
        int f_base = nuc_index/10;
        int l_base = nuc_index%10;
        result_seq[0] = reBASE(f_base);
        result_seq[seq_length-1] = reBASE(l_base);

        // cout << best_sc << endl;
        stk.push(make_tuple(0, seq_length-1, best_index, bestC[seq_length-1][best_index]));
        // cout << bestH[seq_length-1].count(best_index) << endl;
        // cout << bestH[seq_length-1][best_index]->case_ << endl;
        // cout <<  "MANNER_StoS:" << MANNER_StoS <<endl;
        // cout <<  "MANNER_NtoS:" << MANNER_NtoS <<endl;
        while ( !stk.empty() ) {
            tuple<int, int, int, State*> top = stk.top();
            int i = get<0>(top), j = get<1>(top);

            // cout << "i: " << i << " j: " << j <<endl;


            int c_Index = get<2>(top);

            // cout << "c_Index: " << c_Index <<endl;

            state = get<3>(top);

            nuc_index = c_Index%100;
            f_base = nuc_index/10;
            l_base = nuc_index%10;

            Manner MANNER = state->MANNER;
            // cout << "MANNER: " << reManner[MANNER] <<endl;

            // cout << "---" << endl;
            // cout << i << " " << j << endl;
            // cout << reManner[MANNER] << " " << reBASE(f_base) << " " << state->score << " " << reBASE(l_base) << endl;

            int Index_1 = state->index_1;
            int Index_2 = state->index_2;

            // cout << "Index_1: " << Index_1 << " Index_2: " << Index_2 <<endl;

            int nuc_index_1 = Index_1%100;
            int s_1 = Index_1/100;
            int f_base_1 = nuc_index_1/10;
            int l_base_1 = nuc_index_1%10;

            int nuc_index_2 = Index_2%100;
            int s_2 = Index_2/100;
            int f_base_2 = nuc_index_2/10;
            int l_base_2 = nuc_index_2%10;

            // cout <<"--------------------------" <<endl;

            stk.pop();

            switch (MANNER){
                case MANNER_C_StoC:

                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestC[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestS[j][Index_2]));
                    break;

                case MANNER_C_PStoC:

                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestC[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestPS[j][Index_2]));
                    break;

                case MANNER_PStoC:

                    stk.push(make_tuple(s_1, j, Index_1, bestPS[j][Index_1]));
                    break;

                case MANNER_StoS:

                    result[j]='.';
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestS[j-1][Index_1]));
                    break;

                case MANNER_StoH:

                    result[i]='(';
                    result[j]=')';
                    // cout << i << "(" << state->score << ")" << j << endl;
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestS[j-1][Index_1]));
                    break;

                case MANNER_P_StoPS:

                    // cout << Index_1 << " " << Index_2 << endl;
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestP[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestS[j][Index_2]));
                    break;

                case MANNER_HtoP:

                    stk.push(make_tuple(s_1, j, Index_1, bestH[j][Index_1]));
                    break;

                case MANNER_PStoP:

                    result[i]='(';
                    result[j]=')';
                    // cout << i << "(" << state->score << ")" << j << endl;
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestPS[j-1][Index_1]));
                    break;

                case MANNER_PtoP:

                    result[i]='(';
                    result[j]=')';
                    // cout << i << "(" << state->score << ")" << j << endl;
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestP[j-1][Index_1]));
                    break;

                case MANNER_S_PStoP:

                    result[i]='(';
                    result[j]=')';
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestS[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j-1, Index_2, bestPS[j-1][Index_2]));
                    break;

                case MANNER_S_M2toP:

                    result[i]='(';
                    result[j]=')';
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestS[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j-1, Index_2, bestM2[j-1][Index_2]));
                    break;

                case MANNER_M2toP:

                    result[i]='(';
                    result[j]=')';
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestM2[j-1][Index_1]));
                    break;

                case MANNER_PtoM1:

                    stk.push(make_tuple(s_1, j, Index_1, bestP[j][Index_1]));
                    break;

                case MANNER_P_StoM1:

                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestP[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestS[j][Index_2]));
                    break;
                
                case MANNER_M1_M1toM2:

                    // cout << Index_1 << " " << Index_2 << endl;
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestM1[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestM1[j][Index_2]));
                    // cout << "next" << endl;
                    break;

                case MANNER_M2_M1toM2:

                    // cout << Index_1 << " " << Index_2 << endl;
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestM2[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestM1[j][Index_2]));
                    // cout << "next" << endl;
                    break;

                case MANNER_PtoPS:

                    stk.push(make_tuple(s_1, j, Index_1, bestP[j][Index_1]));
                    break;

                case MANNER_NtoS:

                    result[j]='.';
                    result_seq[j]=reBASE(l_base);
                    break;

                case MANNER_StoC:

                    stk.push(make_tuple(s_1, j, Index_1, bestS[j][Index_1]));
                    break;

                default:
                    cout << "Something error" << endl;
            }

        }
    }
    return best_sc;
}

void check_parentheses(int i, int j, int current_index, State* current_state, string& rna_seq, State_pkg& state_pkg){

    State* state = new State;
    std::vector<std::unordered_map<int, State*>>& bestS = state_pkg.bestS; // for skipping
    std::vector<std::unordered_map<int, State*>>& bestC = state_pkg.bestC;
    std::vector<std::unordered_map<int, State*>>& bestH = state_pkg.bestH;
    std::vector<std::unordered_map<int, State*>>& bestP = state_pkg.bestP;
    std::vector<std::unordered_map<int, State*>>& bestPS = state_pkg.bestPS;
    std::vector<std::unordered_map<int, State*>>& bestM1 = state_pkg.bestM1;
    std::vector<std::unordered_map<int, State*>>& bestM2 = state_pkg.bestM2;

    int seq_length = rna_seq.size();
    char result[seq_length+1];
    char result_seq[seq_length+1];

    memset(result, '*', seq_length);
    memset(result_seq, '*', seq_length);
    result[seq_length] = 0;
    result_seq[seq_length] = 0;
    stack<tuple<int, int, int, State*>> stk;

    if(current_state->score!=-1){
        int nuc_index = current_index%100;
        int f_base = nuc_index/10;
        int l_base = nuc_index%10;
        result_seq[i] = reBASE(f_base);
        result_seq[j] = reBASE(l_base);


        // cout << best_sc << endl;
        stk.push(make_tuple(i, j, current_index, current_state));
        // cout << bestH[seq_length-1].count(best_index) << endl;
        // cout << bestH[seq_length-1][best_index]->case_ << endl;
        // cout <<  "MANNER_StoS:" << MANNER_StoS <<endl;
        // cout <<  "MANNER_NtoS:" << MANNER_NtoS <<endl;

        while ( !stk.empty() ) {
            tuple<int, int, int, State*> top = stk.top();
            int i = get<0>(top), j = get<1>(top);

            // cout << "i: " << i << " j: " << j <<endl;


            int c_Index = get<2>(top);

            // cout << "c_Index: " << c_Index <<endl;

            state = get<3>(top);

            nuc_index = c_Index%100;
            f_base = nuc_index/10;
            l_base = nuc_index%10;

            Manner MANNER = state->MANNER;


            cout << reManner[MANNER] << " " << reBASE(f_base) << " " << state->score << " " << reBASE(l_base) << endl;
            // cout << "MANNER: " << reManner[MANNER] <<endl;

            int Index_1 = state->index_1;
            int Index_2 = state->index_2;

            // cout << "Index_1: " << Index_1 << " Index_2: " << Index_2 <<endl;

            int nuc_index_1 = Index_1%100;
            int s_1 = Index_1/100;
            int f_base_1 = nuc_index_1/10;
            int l_base_1 = nuc_index_1%10;

            int nuc_index_2 = Index_2%100;
            int s_2 = Index_2/100;
            int f_base_2 = nuc_index_2/10;
            int l_base_2 = nuc_index_2%10;

            // cout <<"--------------------------" <<endl;

            stk.pop();

            switch (MANNER){
                case MANNER_C_StoC:

                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestC[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestS[j][Index_2]));
                    break;

                case MANNER_C_PStoC:

                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestC[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestPS[j][Index_2]));
                    break;

                case MANNER_PStoC:

                    stk.push(make_tuple(s_1, j, Index_1, bestPS[j][Index_1]));
                    break;

                case MANNER_StoS:

                    result[j]='.';
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestS[j-1][Index_1]));
                    break;

                case MANNER_StoH:

                    result[i]='(';
                    result[j]=')';
                    // cout << i << "(" << state->score << ")" << j << endl;
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestS[j-1][Index_1]));
                    break;

                case MANNER_P_StoPS:

                    // cout << Index_1 << " " << Index_2 << endl;
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestP[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestS[j][Index_2]));
                    break;

                case MANNER_HtoP:

                    stk.push(make_tuple(s_1, j, Index_1, bestH[j][Index_1]));
                    break;

                case MANNER_PStoP:

                    result[i]='(';
                    result[j]=')';
                    // cout << i << "(" << state->score << ")" << j << endl;
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestPS[j-1][Index_1]));
                    break;

                case MANNER_PtoP:

                    result[i]='(';
                    result[j]=')';
                    // cout << i << "(" << state->score << ")" << j << endl;
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestP[j-1][Index_1]));
                    break;

                case MANNER_S_PStoP:

                    result[i]='(';
                    result[j]=')';
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestS[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j-1, Index_2, bestPS[j-1][Index_2]));
                    break;

                case MANNER_S_M2toP:

                    result[i]='(';
                    result[j]=')';
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestS[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j-1, Index_2, bestM2[j-1][Index_2]));
                    break;

                case MANNER_M2toP:

                    result[i]='(';
                    result[j]=')';
                    result_seq[i]=reBASE(f_base);
                    result_seq[j]=reBASE(l_base);
                    stk.push(make_tuple(s_1, j-1, Index_1, bestM2[j-1][Index_1]));
                    break;

                case MANNER_PtoM1:

                    stk.push(make_tuple(s_1, j, Index_1, bestP[j][Index_1]));
                    break;

                case MANNER_P_StoM1:

                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestP[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestS[j][Index_2]));
                    break;
                
                case MANNER_M1_M1toM2:

                    // cout << Index_1 << " " << Index_2 << endl;
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestM1[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestM1[j][Index_2]));
                    // cout << "next" << endl;
                    break;

                case MANNER_M2_M1toM2:

                    // cout << Index_1 << " " << Index_2 << endl;
                    stk.push(make_tuple(s_1, s_2-1, Index_1, bestM2[s_2-1][Index_1]));
                    stk.push(make_tuple(s_2, j, Index_2, bestM1[j][Index_2]));
                    // cout << "next" << endl;
                    break;

                case MANNER_PtoPS:

                    stk.push(make_tuple(s_1, j, Index_1, bestP[j][Index_1]));
                    break;

                case MANNER_NtoS:

                    result[j]='.';
                    result_seq[j]=reBASE(l_base);
                    break;

                case MANNER_StoC:

                    stk.push(make_tuple(s_1, j, Index_1, bestS[j][Index_1]));
                    break;

                default:
                    cout << "Something error" << endl;
            }

        }
    }
    cout << result << endl;
    cout << result_seq << endl;
}






void update_if_better(unordered_map<int, State*>& beamstepU, int base_index, float newscore, int index_1, int index_2, Manner MANNER,
    vector<int>& new_start_codon, vector<int>& new_end_codon){
    
    if(!beamstepU.count(base_index)){ //if the hash does not exists

        State* state = new State;
        state->score = newscore;
        state->index_1 = index_1;
        state->index_2 = index_2;
        state->MANNER = MANNER;

        //CAI
        state->start_codon = new_start_codon;
        state->end_codon = new_end_codon;
        //CAI

        beamstepU[base_index]=state;
        state = nullptr;

    }else{
        
        if(newscore >= beamstepU[base_index]->score){

            State* state = new State;
            state->score = newscore;
            state->index_1 = index_1;
            state->index_2 = index_2;
            state->MANNER = MANNER;

            //CAI
            state->start_codon = new_start_codon;
            state->end_codon = new_end_codon;
            //CAI

            beamstepU[base_index]=state;
            state = nullptr;

        }
    }
}

void update_if_better(unordered_map<int, State*>& beamstepU, int base_index, float newscore, int index_1, Manner MANNER,
    vector<int>& new_start_codon, vector<int>& new_end_codon){
    
    if(!beamstepU.count(base_index)){ //if the hash does not exists

        State* state = new State;
        state->score = newscore;
        state->index_1 = index_1;
        state->MANNER = MANNER;

        //CAI
        state->start_codon = new_start_codon;
        state->end_codon = new_end_codon;
        //CAI

        beamstepU[base_index]=state;
        state = nullptr;

    }else{
        
        if(newscore >= beamstepU[base_index]->score){

            State* state = new State;
            state->score = newscore;
            state->index_1 = index_1;
            state->MANNER = MANNER;

            //CAI
            state->start_codon = new_start_codon;
            state->end_codon = new_end_codon;
            //CAI
            
            beamstepU[base_index]=state;
            state = nullptr;

        }
    }
}

void two_part_CAI(float& new_CAI, vector<int>& new_start_codon, vector<int>& new_end_codon,
    int i ,int j ,int k ,State* state_k_pi , State* state_i_j){
    //[k       pi] [i            j]
    int codon_index;
    new_CAI = 0;

    new_start_codon = {-1,-1};
    new_end_codon = {-1,-1};

    if((j-i == 0) & (i-1-k != 0)){ //length i-j == 1
        //[0][1] | [2]

        if(i%3 == 2){
            codon_index = cal_CAI_index(state_k_pi->end_codon[0], state_k_pi->end_codon[1], state_i_j->start_codon[0]);
            new_CAI = codonset_CAI[codon_index];
        }

        new_end_codon[1] = state_i_j->end_codon[1];
        new_end_codon[0] = state_k_pi->end_codon[1];
        new_start_codon = state_k_pi->start_codon;

    }else if((j-i != 0) & (i-1-k == 0)){ //length i-1-k == 1
        //[0] | [1][2]
        if((i-1)%3 == 0){
            codon_index = cal_CAI_index(state_k_pi->end_codon[1], state_i_j->start_codon[0], state_i_j->start_codon[1]);
            new_CAI = codonset_CAI[codon_index];
        }

        new_end_codon = state_i_j->end_codon;
        new_start_codon[0] = state_k_pi->start_codon[0];
        new_start_codon[1] = state_i_j->start_codon[0];

    }else{
        if(i%3 == 0){
            new_CAI = 0;
            
        }else if(i%3 == 1){ //i%3!=0

            //[0][1][2] [3] | [4][5] [6][7][8]
            codon_index = cal_CAI_index(state_k_pi->end_codon[1], state_i_j->start_codon[0], state_i_j->start_codon[1]);
            new_CAI = codonset_CAI[codon_index];
        }else if(i%3 == 2){ //i%3!=0
            //[0][1][2] [3][4] | [5] [6][7][8]
            // if( i == 17 && j == 44 && k == 3)
            //     cout << reBASE(state_k_pi->end_codon[0]) << reBASE(state_k_pi->end_codon[1]) << reBASE(state_i_j->start_codon[0]) << endl;
            codon_index = cal_CAI_index(state_k_pi->end_codon[0], state_k_pi->end_codon[1], state_i_j->start_codon[0]);
            new_CAI = codonset_CAI[codon_index];   
        }
        new_end_codon = state_i_j->end_codon;
        new_start_codon = state_k_pi->start_codon;
    }
}

void bracket_CAI(float& new_CAI, vector<int>& new_start_codon, vector<int>& new_end_codon,
    int i ,int j ,State* state_ni_pj, int base_i, int base_j){

    new_CAI = 0;
    int codon_index;
   
    //CAI
    //[0 1 2]
    //len of state_ni_pj >=3
    if(j%3 == 2){
        codon_index = cal_CAI_index(state_ni_pj->end_codon[0], state_ni_pj->end_codon[1], base_j); // ...[0][1][j]
        new_CAI = new_CAI + codonset_CAI[codon_index];   
    }

    if(i%3 == 0){
        codon_index = cal_CAI_index(base_i, state_ni_pj->start_codon[0], state_ni_pj->start_codon[1]); // [i][0][1]...
        new_CAI = new_CAI + codonset_CAI[codon_index];
    }

    { //update start/end codon 
        new_end_codon[0]=state_ni_pj->end_codon[1]; //0 is outside pos.
        new_end_codon[1]=base_j;

        new_start_codon[1]=state_ni_pj->start_codon[0]; //0 is outside pos.
        new_start_codon[0]=base_i;
    }
    // update_if_better(beamstepH, base_index, newscore, ni_pj_index,
    //     MANNER_StoH);
}

void bracket_two_CAI(float& new_CAI, vector<int>& new_start_codon, vector<int>& new_end_codon,
    int i ,int j , int k, State* state_nk_i, State* state_ni_pj, int base_i, int base_j, int base_k){
    new_CAI = 0;
    int codon_index;

    // [k]            [j]
    //  ( [0][1] | [2] )
    if(j%3 == 2){

        codon_index = cal_CAI_index(state_ni_pj->end_codon[0], state_ni_pj->end_codon[1], base_j);       
        new_CAI = new_CAI + codonset_CAI[codon_index];     
    }

    if(i-(k+1) == 0){
        if(i%3 == 0){
            //    i
            // ( [0] | [1][2] )
            codon_index = cal_CAI_index(state_nk_i->end_codon[1], state_ni_pj->start_codon[0], state_ni_pj->start_codon[1]);
            new_CAI = new_CAI + codonset_CAI[codon_index];
        }

        if(i%3 == 1){
            // k  i
            // ( [1] | [2][3] )
            codon_index = cal_CAI_index(base_k, state_nk_i->end_codon[1], state_ni_pj->start_codon[0]);
            new_CAI = new_CAI + codonset_CAI[codon_index];
        }
         // k  i
         // ( [2] | [3][4] )
    }else{
        //[i]   [ni]
        //[0] | [1][2] [3][4][5] [6][7][8]
        if(i%3 == 0){
            codon_index = cal_CAI_index(state_nk_i->end_codon[1], state_ni_pj->start_codon[0], state_ni_pj->start_codon[1]);
            new_CAI = new_CAI + codonset_CAI[codon_index];
        //   [i]   [ni]
        //[0][1] | [2] [3][4][5] [6][7][8]
        }else if(i%3 == 1){
            codon_index = cal_CAI_index(state_nk_i->end_codon[0], state_nk_i->end_codon[1], state_ni_pj->start_codon[0]);
            new_CAI = new_CAI + codonset_CAI[codon_index];       
        }

        //[k]   [nk][i]
        //[0] | [1] [2] | 
        if(k%3 == 0){
            codon_index = cal_CAI_index(base_k, state_nk_i->start_codon[0], state_nk_i->start_codon[1]);
            new_CAI = new_CAI + codonset_CAI[codon_index];  
        }

    }
    
    { //update start/end codon 
        new_end_codon[1]=base_j;
        new_end_codon[0]=state_ni_pj->end_codon[1];

        new_start_codon[1]=state_nk_i->start_codon[0]; //0 is outside pos.
        new_start_codon[0]=base_k;
    }
}

void skip_CAI(float& new_CAI, vector<int>& new_start_codon, vector<int>& new_end_codon,
    int ni ,int j ,State* state_ni_pj, int base_j){
    new_CAI = 0;
    
    int codon_index;

    //CAI
    if((j%3 == 2) & (j-1-ni != 0)){
        codon_index = cal_CAI_index(state_ni_pj->end_codon[0], state_ni_pj->end_codon[1], base_j);
        //   [i][0][1]...[0][1][j]
        new_CAI = codonset_CAI[codon_index];

    }
    { //update start/end codon 
        new_end_codon[0]=state_ni_pj->end_codon[1]; //   [i][0][1]...[0][1][j]
        new_end_codon[1]=base_j;

        if(j-1-ni == 0){
            new_start_codon[0]=state_ni_pj->start_codon[0];
            new_start_codon[1]=base_j;
        }else{
            new_start_codon = state_ni_pj->start_codon ;
        }
    }
    //CAI    
}

inline float cal_newscore(float score, float new_CAI){
    return beta*score + lambda*new_CAI;
}

void mv_CDSfold(int beamsize, string& rna_seq, vector<int>& con_seq, State_pkg& state_pkg){

#ifdef timer
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::time_point<std::chrono::system_clock> start1, end1;
#endif    

    std::vector<std::unordered_map<int, State*>>& bestS = state_pkg.bestS; // for skipping
    std::vector<std::unordered_map<int, State*>>& bestC = state_pkg.bestC;
    std::vector<std::unordered_map<int, State*>>& bestH = state_pkg.bestH;
    std::vector<std::unordered_map<int, State*>>& bestP = state_pkg.bestP;
    std::vector<std::unordered_map<int, State*>>& bestPS = state_pkg.bestPS;
    std::vector<std::unordered_map<int, State*>>& bestM1 = state_pkg.bestM1;
    std::vector<std::unordered_map<int, State*>>& bestM2 = state_pkg.bestM2;

    //variables 

    float threshold;

    bestH.clear();
    int seq_length=rna_seq.size();

    bestS.resize(seq_length);
    bestC.resize(seq_length);
    bestH.resize(seq_length);
    bestP.resize(seq_length);
    bestPS.resize(seq_length);
    bestM1.resize(seq_length);
    bestM2.resize(seq_length);

    int base;
    int base_i;
    int base_j;
    int base_pj;
    int base_pi;
    int base_ni;
    int base_k;
    int base_nk;
    int base_q;
    int base_nq;

    int base_jnext;

    int pi_index;
    int nuc_index;
    int base_index;
    int pre_index;
    int ppre_index;
    int pj_index;

    int ni;
    int i;
    int k;
    int nk;
    int pi;
    int p;
    int q;
    int nq;
    int pj;

    int inner_pairs;

    vector<int> j_base_list;
    vector<int> pj_base_list;
    vector<int> i_base_list;
    vector<int> jnext_base_list;
    vector<int> k_base_list;

    vector<int> next_pair[NOTON];
    //variables

    float newscore;

    float new_CAI;
    vector<int> new_start_codon;
    vector<int> new_end_codon;
    //1st initialized j loop
    for(int j=0; j<seq_length; j++){
        vector<int> j_base_list=Base_table.at(rna_seq[j]);
        for(int cj=0; cj<j_base_list.size(); cj++){
            base_j = j_base_list[cj];
            base_index = cal_base_index(j, base_j, base_j); //int i, int base_i, int base_j
            
            State* state = new State;
            state->score = skip_sc;
            state->index_1 = base_index;
            state->MANNER = MANNER_NtoS;
            state->length = 1;
            state->end_codon[1] =base_j; //[][1]
            state->start_codon[0] =base_j; //[0][]

            bestS[j][base_index] = state;
            state = nullptr;

        }
    }

    //main j loop
    for(int j=0; j<seq_length; j++){

#ifdef timer
        start = std::chrono::system_clock::now();
#endif

        unordered_map<int, State*>& beamstepS = bestS[j];
        unordered_map<int, State*>& beamstepH = bestH[j];
        unordered_map<int, State*>& beamstepP = bestP[j];
        unordered_map<int, State*>& beamstepPS = bestPS[j];
        unordered_map<int, State*>& beamstepM1 = bestM1[j];
        unordered_map<int, State*>& beamstepM2 = bestM2[j];
        unordered_map<int, State*>& beamstepC = bestC[j];

        
        j_base_list=Base_table.at(rna_seq[j]);//j list

        //beam of S
        //S=S=.
        //H= (+ S[...] +)

#ifdef timer
        start1 = std::chrono::system_clock::now();
#endif
        pj=j-1;
        if(pj>=0){

            unordered_map<int, State*>& pj_stepS = bestS[j-1];

            for (auto &pj_item : pj_stepS) {
                int ni_pj_index = pj_item.first;
                State* state_ni_pj = pj_item.second;
                ni = ni_pj_index/100;
                nuc_index = ni_pj_index%100;
                base_ni = nuc_index/10;
                base_pj = nuc_index%10;
                i=ni-1;


                for(int cj=0; cj<j_base_list.size(); cj++){
                    base_j = j_base_list[cj];

                    //S[j]=S[ni...j-1]+j[.]
                    if(val_con(con_seq, base_pj, base_j, j-1)){ // validate j and j-1
                                       
                        base_index = cal_base_index(ni, base_ni, base_j);

                        skip_CAI(new_CAI, new_start_codon, new_end_codon,
                            ni ,j ,state_ni_pj, base_j);
                        newscore = state_ni_pj->score;
                        newscore = newscore + cal_newscore(0, new_CAI);
                        // newscore = state_ni_pj->score + new_CAI;

                        update_if_better(beamstepS, base_index, newscore, ni_pj_index, MANNER_StoS, new_start_codon, new_end_codon);
                    } //end
                    
                    //H[j]=i[(]+S[ni...j-1]+j[)]
                    if((i>=0) & (j-i>gap)){
                        i_base_list = Base_table.at(rna_seq[i]);
                        for(int ci = 0; ci<i_base_list.size(); ci++){  
                            base_i = i_base_list[ci];

                            if (_allowed_pairs[base_i][base_j]){ //i can pair with j
                                /*cout << "i: " << i << " j: " << j << endl;
                                cout << "base_i:" << reBASE(base_i) << " base_j:" << reBASE(base_j) << endl;
                                cout << "-------------------------------------------------" << endl;*/
                                if(val_con(con_seq, base_i, base_ni, i) & val_con(con_seq, base_pj, base_j, j-1)){
                                    //now

                                    bracket_CAI(new_CAI, new_start_codon, new_end_codon,
                                        i ,j ,state_ni_pj, base_i, base_j);
                                    base_index = cal_base_index(i, base_i, base_j); //int i, int base_i, int base_j
                                    newscore = state_ni_pj->score;
                                    newscore = newscore + cal_newscore((- v_score_hairpin(i, j, base_i, base_ni, base_pj, base_j)), new_CAI);
                                    // newscore = - v_score_hairpin(i, j, base_i, base_ni, base_pj, base_j) + state_ni_pj->score + new_CAI;

                                    update_if_better(beamstepH, base_index, newscore, ni_pj_index,
                                     MANNER_StoH, new_start_codon, new_end_codon);

                                }    
                            }
                        }  
                    } //end
                }
            }//end of pj_stepS
        }
        if(beamsize != 0){
            threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepS, bestC, true);
            threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepH, bestC, false);
        }
        

#ifdef timer
        end1 = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds1 = end1 - start1;

        cout << "S elapsed time: " << elapsed_seconds1.count() << "s" << endl;
#endif

        //PS=P[(a)]+S[...]
        //PS[k...j]=P[k...pi]+S[i...j]
        //[0][1][2] [3][4][5] [6][7][8]
#ifdef timer
        start1 = std::chrono::system_clock::now();
#endif

        for (auto &j_item : beamstepS) {

            int i_j_index = j_item.first;
            State* state_i_j = j_item.second;
            i = i_j_index/100;
            nuc_index = i_j_index%100;
            base_i = nuc_index/10;
            base_j = nuc_index%10;

            pi=i-1;
            if(pi>=0){

                // unordered_map<int, State*>& pi_stepP = bestP[pi];
                unordered_map<int, State*>& pi_stepP = bestP[pi];

                for (auto &pi_item : pi_stepP) {

                    int k_pi_index = pi_item.first;
                    State* state_k_pi = pi_item.second;
                    k = k_pi_index/100;
                    nuc_index = k_pi_index%100;
                    base_k = nuc_index/10;
                    base_pi = nuc_index%10;

                    if(val_con(con_seq, base_pi, base_i, pi)){
                        two_part_CAI(new_CAI, new_start_codon, new_end_codon, i, j, k ,state_k_pi ,state_i_j);
                        newscore = state_k_pi->score + state_i_j->score;
                        newscore = newscore + cal_newscore(0, new_CAI);
                        // newscore = state_k_pi->score + state_i_j->score + new_CAI;
                        base_index = cal_base_index(k,base_k,base_j);

                        update_if_better(beamstepPS, base_index, newscore, k_pi_index, i_j_index, 
                            MANNER_P_StoPS, new_start_codon, new_end_codon);    
                    }
                    
                     
                    
                }  
            }   
        } //end
        // end of beam S

#ifdef timer
        end1 = std::chrono::system_clock::now();

        elapsed_seconds1 = end1 - start1;

        cout << "PS elapsed time: " << elapsed_seconds1.count() << "s" << endl;
#endif

        // beam of H
        //P=H

#ifdef timer
        start1 = std::chrono::system_clock::now();
#endif

        for (auto &j_item : beamstepH) {

            int i_j_index = j_item.first;
            State* state_i_j = j_item.second;
            i = i_j_index/100;
            nuc_index = i_j_index%100;
            base_i = nuc_index/10;
            base_j = nuc_index%10;
            base_index = cal_base_index(i,base_i,base_j);
            newscore = state_i_j->score;

            new_start_codon = state_i_j->start_codon;
            new_end_codon = state_i_j->end_codon;

            update_if_better(beamstepP, base_index, newscore, i_j_index,
                MANNER_HtoP, new_start_codon, new_end_codon);
        }//end P=H
        //end of beam H

#ifdef timer
        end1 = std::chrono::system_clock::now();

        elapsed_seconds1 = end1 - start1;

        cout << "H elapsed time: " << elapsed_seconds1.count() << "s" << endl;
#endif

#ifdef timer
        start1 = std::chrono::system_clock::now();
#endif

        //P=(+PS[(...)...]+)
        //P=(+S[...]+PS+)
        //beam of PS
        if(j-1>=0){

            unordered_map<int, State*>& pj_stepPS = bestPS[j-1];

            for (auto &pj_item : pj_stepPS) {

                int ni_pj_index = pj_item.first;
                // printf("%d\n", i);
                State* state_ni_pj = pj_item.second;
                ni = ni_pj_index/100;
                nuc_index = ni_pj_index%100;
                base_ni = nuc_index/10;
                base_pj = nuc_index%10;
                i = ni - 1;

                if(i>=0){

                    int pre_k_pi_index = state_ni_pj->index_1; //P k...pi
                    int pre_i_j_index = state_ni_pj->index_2; //S i...j

                    // cout << "pre index:" << endl;
                    if(pre_i_j_index==-1){ //PS[p(...q)]
                        nq=j;
                        q=j-1;
                        base_q = base_pj;
                        base_nq = base_j;

                    }else{
                        nq=pre_i_j_index/100; //i in S
                        q=nq-1;

                        nuc_index = pre_k_pi_index%100;
                        base_q = nuc_index%10;
                        nuc_index = pre_i_j_index%100;
                        base_nq = nuc_index/10;

                    }

                    for(int cj=0; cj<j_base_list.size(); cj++){
                        // if(i>=0){
                            i_base_list=Base_table.at(rna_seq[i]);//j list
                            for(int ci=0; ci<i_base_list.size(); ci++){
                                base_j = j_base_list[cj];
                                base_i = i_base_list[ci];
                                if (_allowed_pairs[base_i][base_j]){
                                    if( val_con(con_seq, base_pj, base_j, j-1) & val_con(con_seq, base_i, base_ni, i) ){
                                        //P=(+PS[(...)...]+)

                                        bracket_CAI(new_CAI, new_start_codon, new_end_codon,
                                            i ,j ,state_ni_pj, base_i, base_j);
                                        newscore = state_ni_pj->score;
                                        newscore = newscore + cal_newscore((- v_score_single(i, j, ni, q, base_i, base_ni, base_pj, base_j, base_i, base_ni, base_q, base_nq)), new_CAI);
                                //         newscore = - v_score_single(i, j, ni, q,
                                // base_i, base_ni, base_pj, base_j,
                                // base_i, base_ni, base_q, base_nq)+state_ni_pj->score + new_CAI;
                                        
                                        base_index = cal_base_index(i,base_i,base_j);
                                        update_if_better(beamstepP, base_index, newscore, ni_pj_index,
                                            MANNER_PStoP, new_start_codon, new_end_codon);
                                        

                                    }  
                                }
                                   
                            }    
                        // }
                    }
                    //end P[i...j]=i[(]+PS[ni(A)...pj]+j[)]

                    // P[i...j]=k[(]+S[nk...i]+PS[ni(A)...pj]+j[)]

                    unordered_map<int, State*>& i_stepS = bestS[i];

                    for (auto &i_item : i_stepS) {
                        int nk_i_index = i_item.first;
                        State* state_nk_i = i_item.second;
                        nk = nk_i_index/100;
                        nuc_index = nk_i_index%100;
                        base_nk = nuc_index/10;
                        base_i = nuc_index%10;

                        k = nk - 1;

                        if(k >= 0){

                            k_base_list=Base_table.at(rna_seq[k]);//j list

                            for(int ck=0; ck<k_base_list.size(); ck++){
                                for(int cj=0; cj<j_base_list.size(); cj++){

                                    base_k = k_base_list[ck];
                                    base_j = j_base_list[cj];
                                    if (_allowed_pairs[base_k][base_j]){
                                        if( val_con(con_seq, base_pj, base_j, j-1) & val_con(con_seq, base_i, base_ni, i) & val_con(con_seq, base_k, base_nk, k) ){
                                            //P=(+S[...]+PS+)

                                            bracket_two_CAI(new_CAI, new_start_codon, new_end_codon,
                                                i ,j ,k, state_nk_i, state_ni_pj, base_i, base_j, base_k);
                                            newscore = state_nk_i->score + state_ni_pj->score;
                                            newscore = newscore + cal_newscore((- v_score_single(k, j, ni, q, base_k, base_nk, base_pj, base_j, base_i, base_ni, base_q, base_nq)), new_CAI);
                                            // newscore = - v_score_single(k, j, ni, q,
                                            //     base_k, base_nk, base_pj, base_j,
                                            //     base_i, base_ni, base_q, base_nq)+state_ni_pj->score + new_CAI;
                                            base_index = cal_base_index(k,base_k,base_j);

                                            update_if_better(beamstepP, base_index, newscore, nk_i_index, ni_pj_index,
                                                MANNER_S_PStoP, new_start_codon, new_end_codon);

                                        }
                                    }
                                    


                                }


                            }    
                        }
                        
                    }
                    //end P[i...j]=k[(]+S[nk...i]+PS[ni(A)...pj]+j[)]
                }
            }

            unordered_map<int, State*>& pj_stepP = bestP[j-1];
            
            for (auto &pj_item : pj_stepP) {

                int ni_pj_index = pj_item.first;
                // printf("%d\n", i);
                State* state_ni_pj = pj_item.second;
                ni = ni_pj_index/100;
                nuc_index = ni_pj_index%100;
                base_ni = nuc_index/10;
                base_pj = nuc_index%10;
                i = ni - 1;

                if(i>=0){

                    nq=j;
                    q=j-1;
                    base_q = base_pj;
                    base_nq = base_j;

                    for(int cj=0; cj<j_base_list.size(); cj++){
                        i_base_list=Base_table.at(rna_seq[i]);//j list
                        for(int ci=0; ci<i_base_list.size(); ci++){
                            base_j = j_base_list[cj];
                            base_i = i_base_list[ci];
                            if (_allowed_pairs[base_i][base_j]){
                                if( val_con(con_seq, base_pj, base_j, j-1) & val_con(con_seq, base_i, base_ni, i) ){
                                    //P[i...j]=i[(]+P[ni(A)...pj]+j[)]

                                    bracket_CAI(new_CAI, new_start_codon, new_end_codon,
                                        i ,j ,state_ni_pj, base_i, base_j);
                                    newscore = state_ni_pj->score;
                                    newscore = newscore + cal_newscore((- v_score_single(i, j, ni, q, base_i, base_ni, base_pj, base_j, base_i, base_ni, base_q, base_nq)), new_CAI);
                            //         newscore = - v_score_single(i, j, ni, q,
                            // base_i, base_ni, base_pj, base_j,
                            // base_i, base_ni, base_q, base_nq)+state_ni_pj->score + new_CAI;

                                    base_index = cal_base_index(i,base_i,base_j);

                                    update_if_better(beamstepP, base_index, newscore, ni_pj_index,
                                        MANNER_PtoP, new_start_codon, new_end_codon);

                                }  
                            }
                                   
                        }    
                    }
                }
            }

        }
        // end of beam of PS

#ifdef timer
        end1 = std::chrono::system_clock::now();

        elapsed_seconds1 = end1 - start1;

        cout << "P PS elapsed time: " << elapsed_seconds1.count() << "s" << endl;
#endif        

#ifdef timer
        start1 = std::chrono::system_clock::now();
#endif

        //beam of Multi
        //P=(+S[...]+M2[(a)...(b)...]+)
        if(j-1>0){

            unordered_map<int, State*>& pj_stepM2 = bestM2[j-1];
            
            for (auto &pj_item : pj_stepM2) {

                int ni_pj_index = pj_item.first;
                State* state_ni_pj = pj_item.second;
                ni = ni_pj_index/100;
                nuc_index = ni_pj_index%100;
                base_ni = nuc_index/10;
                base_pj = nuc_index%10;

                i = ni - 1;

                if(i >= 0){
                    //P=[(]+S+M2+[)]
                    unordered_map<int, State*>& i_stepS = bestS[i];
                    for (auto &i_item : i_stepS) {
                        int nk_i_index = i_item.first;
                        State* state_nk_i = i_item.second;
                        nk = nk_i_index/100;
                        nuc_index = nk_i_index%100;
                        base_nk = nuc_index/10;
                        base_i = nuc_index%10;
                        k = nk - 1;
                        if(k>=0){
                            k_base_list=Base_table.at(rna_seq[k]);
                            for(int cj=0; cj<j_base_list.size(); cj++){
                                for(int ck=0; ck<k_base_list.size(); ck++){
                                    base_j=j_base_list[cj];
                                    base_k=k_base_list[ck];
                                    if(_allowed_pairs[base_k][base_j]){
                                        if( val_con(con_seq, base_k, base_nk, k) & val_con(con_seq, base_i, base_ni, i) & val_con(con_seq, base_pj, base_j, j-1) ){

                                            bracket_two_CAI(new_CAI, new_start_codon, new_end_codon,
                                                i ,j ,k, state_nk_i, state_ni_pj, base_i, base_j, base_k);
                                            newscore = state_nk_i->score + state_ni_pj->score;
                                            newscore = newscore + cal_newscore((- v_score_multi(base_k, base_j)),new_CAI);
                                            // newscore = - v_score_multi(base_k, base_j) + state_ni_pj->score + new_CAI;

                                            base_index = cal_base_index(k,base_k,base_j);

                                            update_if_better(beamstepP, base_index, newscore, nk_i_index, ni_pj_index,
                                            MANNER_S_M2toP, new_start_codon, new_end_codon);

                                        }
                                    }
                                }
                            }   
                        }
                        

                    } //end P=[(]+S+M2+[)]

                    ////P=(+M2[(a)...(b)...]+)
                    i_base_list=Base_table.at(rna_seq[i]);
                    for(int cj=0; cj<j_base_list.size(); cj++){
                        for(int ci=0; ci<i_base_list.size(); ci++){
                            base_j = j_base_list[cj];
                            base_i = i_base_list[ci];
                            if(_allowed_pairs[base_i][base_j]){
                                if( val_con(con_seq, base_i, base_ni, i) & val_con(con_seq, base_pj, base_j, j-1) ){

                                    bracket_CAI(new_CAI, new_start_codon, new_end_codon,
                                        i ,j ,state_ni_pj, base_i, base_j);
                                    newscore = state_ni_pj->score;
                                    newscore = newscore + cal_newscore((- v_score_multi(base_i, base_j)), new_CAI);
                                    // newscore = - v_score_multi(base_i, base_j) + state_ni_pj->score + new_CAI;

                                    base_index = cal_base_index(i,base_i,base_j);

                                    update_if_better(beamstepP, base_index, newscore, ni_pj_index,
                                    MANNER_M2toP, new_start_codon, new_end_codon);
                                }    
                            }
                        }
                    }
                }   
            }
        }
        if(beamsize !=  0){
            threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepP, bestC, false);    
        }
        

#ifdef timer
        end1 = std::chrono::system_clock::now();

        elapsed_seconds1 = end1 - start1;

        cout << "Mt elapsed time: " << elapsed_seconds1.count() << "s" << endl;
#endif

        //end of beam Multi

#ifdef timer
        start1 = std::chrono::system_clock::now();
#endif

        //beam of M1
        //M1=P[(a)]
        //M1=P
        for (auto &j_item : beamstepP) {
            int i_j_index = j_item.first;
            State* state_i_j = j_item.second;
            i = i_j_index/100;
            nuc_index = i_j_index%100;
            base_i = nuc_index/10;
            base_j = nuc_index%10;
            base_index = cal_base_index(i,base_i,base_j);
            
            new_start_codon = state_i_j->start_codon;
            new_end_codon = state_i_j->end_codon;

            newscore = state_i_j->score;
            newscore = newscore + cal_newscore((- v_score_M1(base_i, base_j)), 0);
            // newscore = - v_score_M1(base_i, base_j)+state_i_j->score;

            update_if_better(beamstepM1, base_index, newscore, i_j_index, MANNER_PtoM1,
                new_start_codon, new_end_codon);
        }//end M1=P

        //M1=P[(a)]+S[...]
        //M1=P+S
        for (auto &j_item : beamstepS) {
            int i_j_index = j_item.first;
            State* state_i_j = j_item.second;
            i = i_j_index/100;
            nuc_index = i_j_index%100;
            base_i = nuc_index/10;
            base_j = nuc_index%10;
            pi=i-1;
            if(pi>=0){
                unordered_map<int, State*>& pi_stepP = bestP[pi];
                for (auto &pi_item : pi_stepP) {
                    int k_pi_index = pi_item.first;
                    State* state_k_pi = pi_item.second;
                    k = k_pi_index/100;
                    nuc_index = k_pi_index%100;
                    base_k = nuc_index/10;
                    base_pi = nuc_index%10;

                    if( val_con(con_seq, base_pi, base_i, pi) ){

                        two_part_CAI(new_CAI, new_start_codon, new_end_codon,
                            i ,j ,k ,state_k_pi ,state_i_j);

                        newscore = state_k_pi->score + state_i_j->score;
                        newscore = newscore + cal_newscore((- v_score_M1(base_k, base_pi)), new_CAI);
                        // newscore = - v_score_M1(base_k, base_pi)+state_k_pi->score + new_CAI;

                        base_index = cal_base_index(k,base_k,base_j);

                        update_if_better(beamstepM1, base_index, newscore, k_pi_index, i_j_index,
                                        MANNER_P_StoM1, new_start_codon, new_end_codon);
                    }

                }    
            }    
        }//end M1=P+S
        if(beamsize != 0){
            threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepM1, bestC, false);    
        }
        
        //end of beam M1

#ifdef timer
        end1 = std::chrono::system_clock::now();

        elapsed_seconds1 = end1 - start1;

        cout << "M1 elapsed time: " << elapsed_seconds1.count() << "s" << endl;
#endif





#ifdef timer
        start1 = std::chrono::system_clock::now();
#endif

        //beam of M2

        //M2=M1[(a)...]+M1[(b)...]
        //M2=M1+M1
        for (auto &j_item : beamstepM1) {
            int i_j_index = j_item.first;
            State* state_i_j = j_item.second;
            i = i_j_index/100;
            nuc_index = i_j_index%100;
            base_i = nuc_index/10;
            base_j = nuc_index%10;

            pi=i-1;
            if(pi>=0){
               unordered_map<int, State*>& pi_stepM1 = bestM1[pi];
               for (auto &pi_item : pi_stepM1) {
                    int k_pi_index = pi_item.first;
                    State* state_k_pi = pi_item.second;
                    k = k_pi_index/100;
                    nuc_index = k_pi_index%100;
                    base_k = nuc_index/10;
                    base_pi = nuc_index%10;

                    if( val_con(con_seq, base_pi, base_i, pi) ){

                        two_part_CAI(new_CAI, new_start_codon, new_end_codon,
                            i ,j ,k ,state_k_pi ,state_i_j);

                        newscore = state_i_j->score + state_k_pi->score;
                        newscore = newscore + cal_newscore(0, new_CAI);
                        // newscore = state_i_j->score + state_k_pi->score + new_CAI;

                        base_index = cal_base_index(k,base_k,base_j);

                        update_if_better(beamstepM2, base_index, newscore, k_pi_index, i_j_index,
                                        MANNER_M1_M1toM2, new_start_codon, new_end_codon);
                    }
                }

                //M2=M2+M1
                unordered_map<int, State*>& pi_stepM2 = bestM2[pi];
                for (auto &pi_item : pi_stepM2) {
                    int k_pi_index = pi_item.first;
                    State* state_k_pi = pi_item.second;
                    k = k_pi_index/100;
                    nuc_index = k_pi_index%100;
                    base_k = nuc_index/10;
                    base_pi = nuc_index%10;

                    if( val_con(con_seq, base_pi, base_i, pi) ){
                        
                        two_part_CAI(new_CAI, new_start_codon, new_end_codon,
                            i , j, k, state_k_pi ,state_i_j);

                        newscore = state_i_j->score + state_k_pi->score;
                        newscore = newscore + cal_newscore(0, new_CAI);
                        // newscore = state_i_j->score + state_k_pi->score + new_CAI;

                        base_index = cal_base_index(k,base_k,base_j);
                        update_if_better(beamstepM2, base_index, newscore, k_pi_index, i_j_index,
                                        MANNER_M2_M1toM2, new_start_codon, new_end_codon);
                    }
                } 
            }    
        }//end M2=M1+M1
        if(beamsize != 0){
            threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepM2, bestC, false);  
        }
        
        //end beam M2

#ifdef timer
        end1 = std::chrono::system_clock::now();

        elapsed_seconds1 = end1 - start1;

        cout << "M2 elapsed time: " << elapsed_seconds1.count() << "s" << endl;
#endif

        //beam of P

        //PS=P
        for (auto &j_item : beamstepP) {
            int i_j_index = j_item.first;
            State* state_i_j = j_item.second;
            i = i_j_index/100;
            nuc_index = i_j_index%100;
            base_i = nuc_index/10;
            base_j = nuc_index%10;
            base_index = cal_base_index(i,base_i,base_j);

            newscore = state_i_j->score;

            new_start_codon = state_i_j->start_codon;
            new_end_codon = state_i_j->end_codon;

            update_if_better(beamstepPS, base_index, newscore, i_j_index,
                MANNER_PtoPS, new_start_codon, new_end_codon);
        }//end PS=P
        if(beamsize != 0){
            threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepPS, bestC, false);    
        }
        
        //end beam P



#ifdef timer
        start1 = std::chrono::system_clock::now();
#endif

        //beam of C
        for (auto &j_item : beamstepPS) {
            int i_j_index = j_item.first;
            State* state_i_j = j_item.second;
            i = i_j_index/100;
            nuc_index = i_j_index%100;
            base_i = nuc_index/10;
            base_j = nuc_index%10;

            pi = i-1;
            //C=S[...]+PS[(a)...]
            if(pi>=0){
                unordered_map<int, State*>& pi_stepC = bestC[pi];
                for (auto &pi_item : pi_stepC) {
                    int k_pi_index = pi_item.first;
                    State* state_k_pi = pi_item.second;
                    k = k_pi_index/100;
                    nuc_index = k_pi_index%100;
                    base_k = nuc_index/10;
                    base_pi = nuc_index%10;
                    
                    if( val_con(con_seq, base_pi, base_i, pi) & (k==0) ){

                        two_part_CAI(new_CAI, new_start_codon, new_end_codon,
                            i , j, k, state_k_pi ,state_i_j);

                        base_index = cal_base_index(0,base_k,base_j);

                        newscore = state_k_pi->score + state_i_j->score;
                        newscore = newscore + cal_newscore(0, new_CAI);
                        // newscore = state_i_j->score + new_CAI;

                        update_if_better(beamstepC, base_index, newscore, k_pi_index, i_j_index,
                            MANNER_C_PStoC, new_start_codon, new_end_codon);

                    }

                }
            }else if(i==0){
                base_index = cal_base_index(i,base_i,base_j);
                newscore = state_i_j->score;

                new_start_codon = state_i_j->start_codon;
                new_end_codon = state_i_j->end_codon;

                update_if_better(beamstepC, base_index, newscore, i_j_index,
                            MANNER_PStoC, new_start_codon, new_end_codon);    
            }

            
        }

        if(j>0){
            unordered_map<int, State*>& pj_stepC = bestC[j-1];
            for (auto &pj_item : pj_stepC) {
                int i_pj_index = pj_item.first;
                State* state_i_pj = pj_item.second;
                i = i_pj_index/100;
                nuc_index = i_pj_index%100;
                base_i = nuc_index/10;
                base_pj = nuc_index%10;

                for(int cj=0; cj<j_base_list.size(); cj++){
                    base_j = j_base_list[cj];

                    if(val_con(con_seq, base_pj, base_j, j-1)){
                        //C = C[]+[.]
                        int base_index = cal_base_index(i,base_i,base_j);

                        skip_CAI(new_CAI, new_start_codon, new_end_codon,
                            i ,j ,state_i_pj, base_j);
                        

                        newscore = state_i_pj->score;
                        newscore = newscore + cal_newscore(0, new_CAI);
                        // newscore = state_i_pj->score + new_CAI;
                        
                        int s_index = cal_base_index(j,base_j,base_j);
                        update_if_better(beamstepC, base_index, newscore, i_pj_index, s_index,
                        MANNER_C_StoC, new_start_codon, new_end_codon);
                    }
                }
            }
        }else if(j==0){
            unordered_map<int, State*>& j_stepS = beamstepS;
            for (auto &j_item : j_stepS){
                int i_j_index = j_item.first;
                State* state_i_j = j_item.second;
                newscore=0;

                new_start_codon = state_i_j->start_codon;
                new_end_codon = state_i_j->end_codon;

                update_if_better(beamstepC, i_j_index, newscore, i_j_index,
                MANNER_StoC, new_start_codon, new_end_codon);
            }
        }




#ifdef timer
        end1 = std::chrono::system_clock::now();

        elapsed_seconds1 = end1 - start1;

        cout << "C elapsed time: " << elapsed_seconds1.count() << "s" << endl;
#endif





        //end of beam C




// #ifdef timer
//         start1 = std::chrono::system_clock::now();
// #endif

//         //beam pruning
//         // cout << "--------------------------------------------------------" << endl;
        // cout << "Number of states:" << beamstepPS.size() << endl;
        
//         threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepS, bestC);
//         threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepH, bestC);
//         threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepP, bestC);
//         threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepPS, bestC);
//         threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepM1, bestC);
//         threshold = beam_prune(beamsize, con_seq, rna_seq, beamstepM2, bestC);

// #ifdef timer
//         end1 = std::chrono::system_clock::now();

//         elapsed_seconds1 = end1 - start1;

//         cout << "Beam pruning elapsed time: " << elapsed_seconds1.count() << "s" << endl;
// #endif



#ifdef verbose
        cout << j << endl;  
        cout << "S Number of states (pruning):" << beamstepS.size() << endl;
        cout << "H Number of states (pruning):" << beamstepH.size() << endl;
        cout << "P Number of states (pruning):" << beamstepP.size() << endl;
        cout << "PS Number of states (pruning):" << beamstepPS.size() << endl;
        cout << "M1 Number of states (pruning):" << beamstepM1.size() << endl;
        cout << "M2 Number of states (pruning):" << beamstepM2.size() << endl;
        cout << "C Number of states (pruning):" << beamstepC.size() << endl;
        cout << "--------------------------------------------------------" << endl;
       
        for (auto x : beamstepS){
            cout << x.first << " " << x.second->score << endl;
        }
#endif

        // cout << "---------------------------" << endl;
        // cout << "j " << j << endl;
/*        if(j==65){
            for (auto x : beamstepP){
                cout << x.first << " " << x.second->score << endl;
                i =  x.first/100;
                if(i==8){
                    check_parentheses(i, j, x.first, beamstepP[x.first], rna_seq, state_pkg);
                }
            }    
        }*/
        

#ifdef timer
        end = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds = end - start;
        cout << j << " step "<<  "elapsed time: " << elapsed_seconds.count() << "s" << endl;
#endif




    } //end of j

    // for (auto x : bestC[seq_length-1]){
    //     cout << x.first << " " << x.second->score << endl;
    // }
}

inline void log_CAI(vector<float>& cai_vector, vector<float>& log_cai_vector){

    for(int i = 0; i < cai_vector.size(); i++){
        log_cai_vector.push_back(log10f(cai_vector[i]));
    }
}

inline float cal_CAI(char* result_seq, int size){
    float CAI=0;
    int cai_index;
    for(int i = 0; i < size; i++){
        cai_index = cal_base_index(BASE(result_seq[3*i]),BASE(result_seq[3*i+1]),BASE(result_seq[3*i+2]));
        CAI = CAI + codonset_CAI[cai_index];
    }
    return pow(10.0, CAI/size);
}

int main(int argc, char** argv){
    std::chrono::time_point<std::chrono::system_clock> start, end;
    bool fasta = false;
    int beamsize=100;
    string seq="";
    ifstream fasta_file;
    ifstream cai_file;

    string rna_seq,ami_seq;
    vector<int> con_seq;
    vector<string> rna_seq_list, rna_name_list, ami_seq_list, ami_name_list;
    vector<vector<int> > con_seq_list;
    std::vector<float> cai_vector;
    std::vector<float> log_cai_vector;

    initialize();

    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        cout << "./LinearCDSfold [-f] [fasta_file_path] [-b] [beam_size] [-cai] [cai_file_path]" << endl;
    }
    const std::string &beamsize_ = input.getCmdOption("-b");
    if (!beamsize_.empty()){
        // cout << "old beam size: " << beamsize << endl;
        beamsize = stoi(beamsize_);
        cout << "beam size: " << beamsize << endl;
    }

    const std::string &cai_file_str = input.getCmdOption("-cai");
    if (!cai_file_str.empty()){

        cai_file.open(cai_file_str);
        std::string line;
        while (std::getline(cai_file, line)) {
            std::string token;
            int pos = 0;
            
            while ((pos = line.find(",")) != std::string::npos) {
                token = line.substr(0, pos);
                cai_vector.push_back(std::stof(token));
                line.erase(0, pos + 1);
            }
            
            cai_vector.push_back(std::stof(line));
        }
        log_CAI(cai_vector, log_cai_vector);
        initialize_CAI_table(log_cai_vector);
        // initialize_CAI_table(cai_vector);
    }

    const std::string &lambda_ = input.getCmdOption("-lambda");
    if (!lambda_.empty()){
        // cout << "old beam size: " << beamsize << endl;
        lambda = stof(lambda_);
        lambda = lambda*100;
        cout << "lambda: " << lambda << endl;
    }

    const std::string &beta_ = input.getCmdOption("-beta");
    if (!beta_.empty()){
        // cout << "old beam size: " << beamsize << endl;
        beta = stof(beta_);
        cout << "beta: " << beta << endl;
    }

    const std::string &file = input.getCmdOption("-f");
    if (!file.empty()){
        fasta_file.open(file);

        while (fasta_file >> seq){
            // cout << seq << endl;
            if (seq.empty()) continue;
            else if (seq[0] == '>' or seq[0] == ';'){
                if (!ami_seq.empty())
                    ami_seq_list.push_back(ami_seq);
                ami_seq.clear();
                continue;
            }else{
                rtrim(seq);
                ami_seq += seq;
            }
        }
        if (!ami_seq.empty())
            ami_seq_list.push_back(ami_seq);
        cout << "input file: " << file << endl;

    }else{
        for (string seq; getline(cin, seq);){
            if (seq.empty()) continue;
            if (!isalpha(seq[0])){
                printf("Unrecognized sequence: %s\n", seq.c_str());
                continue;
            }
            ami_seq_list.push_back(seq);
        }
    }

    for(int i = 0; i < ami_seq_list.size(); i++){

#ifdef timer
        start = std::chrono::system_clock::now();
#endif

        ami_seq = ami_seq_list[i];
        transform(ami_seq.begin(), ami_seq.end(), ami_seq.begin(), ::toupper);
        ami_to_rna(rna_seq_list,ami_seq);
        add_con_seq(con_seq_list,ami_seq);
        
        rna_seq = rna_seq_list[i];
        con_seq = con_seq_list[i];

        cout << "input: ";
        
        printf("%s\n", ami_seq.c_str());
        // printf("%s\n", rna_seq.c_str());

        char result[rna_seq.size()+1];
        char result_seq[rna_seq.size()+1];

        string check_ami_seq;
        std::vector<std::unordered_map<int, State*>> bestH;
        State_pkg state_pkg;
        mv_CDSfold(beamsize, rna_seq, con_seq, state_pkg);

#ifdef timer
        end = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds = end - start;

        cout << "LCDS elapsed time: " << elapsed_seconds.count() << "s" << endl;
#endif

        /*get_parentheses(result_seq, result, rna_seq, state_pkg);
        cout << result << endl;
        cout << result_seq << endl;*/


#ifdef timer
        start = std::chrono::system_clock::now();
#endif

        int best_score = get_parentheses(result_seq, result, rna_seq, state_pkg);

#ifdef timer
        end = std::chrono::system_clock::now();

        elapsed_seconds = end - start;
        
        cout << "Back parentheses elapsed time: " << elapsed_seconds.count() << "s" << endl;
#endif

        // transform(result_seq.begin(), result_seq.end(), result_seq.begin(), ::toupper);
        for(int i=0; i<rna_seq.size(); i++){
            result_seq[i]=toupper(result_seq[i]);
        }
        check_valid_ami(ami_seq.size(), result_seq,check_ami_seq);
        // cout << check_ami_seq << endl;

        //check correct
        bool codon_coorect = true;

        for(int i = 0; i<ami_seq.size(); i++ ){
            if(ami_seq[i] != check_ami_seq[i]){
                cout << "codon " << i << " ";
                cout << result_seq[i*3] << result_seq[i*3+1] << result_seq[i*3+2] << " error" << endl;
                codon_coorect = false;
            }
        }

        if(!codon_coorect)
            cout << "ami codon error" << endl;
        else
            cout << "ami codon correct" << endl;

        //print results
        cout << "output:" << check_ami_seq << endl;
        cout << result << endl;
        cout << result_seq << endl;
        cout << "score: " << float(best_score)/-100 << endl;

        cout << "log(CAI): " << cal_CAI(result_seq, ami_seq.size()) << endl;


    }

    
        
    return 0;
}
