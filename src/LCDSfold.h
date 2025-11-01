#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <limits>
#include <cctype>
#include <chrono>
#include <ctime>
#include <math.h>
#include <queue>

#include "functions.h"
#include "utility.h"
#include "beamprune.h"
#include "ac_matcher.h"

using namespace std;
using namespace oligodesign;
#define GET_ACGU_NUM_V(x) ((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: 0)))))
bool debug = false;

// #define debug

template <typename T>
T BackTrack(string& rna_solution, string& structure_solution, string& rna_seq, AllTables<T>& alltables, vector<int>& con_seq){
    // cout<<"start back tracking"<<endl;
    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    
    int seq_length = rna_seq.size();
    rna_solution.assign(seq_length, '.');
    structure_solution.assign(seq_length, '.');

    vector<int> first_nuclist = Base_table.at(rna_seq[0]);//i list
    vector<int> last_nuclist = Base_table.at(rna_seq[seq_length-1]);//j list

    State<T> state;
    T maxscore = VALUE_MIN<T>();
    int maxindex;
    for(int ci = 0; ci < first_nuclist.size(); ci++){
        for(int cj = 0; cj < last_nuclist.size(); cj++){
            int nuci = first_nuclist[ci];
            int nucj = last_nuclist[cj];
            int index = GetIndex(0, nuci, nucj);
            if(bestF[seq_length-1].count(index)){
                T newscore = bestF[seq_length-1][index].score;
                if (maxscore < newscore){
                    maxscore = newscore;
                    maxindex = index;
                    state = bestF[seq_length-1][index];}}
    }}

    stack<tuple<int, int, State<T>>> stk;
    stk.push(make_tuple(seq_length-1, maxindex, state));
    while ( !stk.empty() ) {

        tuple<int, int, State<T>> top = stk.top();
        int j = get<0>(top);
        int cindex = get<1>(top);
        state = get<2>(top);
        int index1 = state.index_1;
        int index2 = state.index_2;
        int index3 = state.index_3;
        int i, nuci, nuci_pair,nucj, nucx,len;
        Manner MANNER = state.MANNER;
        State<T> next_state, next_state1, next_state2, next_state3;
        // cout<<mannerToString(MANNER)<<endl;
        if(MANNER == MANNER_C_StoCS)//make_tuple(i, nuci, nuci_pair,nucj, nucx, len);
            std::tie(i, nuci, nuci_pair,nucj, nucx, len) = GetIndexTupleCS(cindex);
        else{
            std::tie(i, nuci, nucj) = GetIndexTuple(cindex);
        }

        
        stk.pop();
        int k, nuck, nuck_pair,nucl, posi, s, nucs, nucs_end;
        int struct_score,p,q,q1,p_1,l1,l2;
        switch (MANNER){
            
            case NONE:
                break;
            case MANNER_NONEtoN:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_NONEtoF:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_NONEtoS:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_N_EtoN:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestN[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_NtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestN[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_S_HP:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                for(int x = i + 1; x < j ; x ++ ){
                    rna_solution[x]=sp_backtrack[cindex][x - i];
                    structure_solution[x]='.';
                }
                break;
            case MANNER_S_EtoS:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                std::tie(k, nuck, nucl) = GetIndexTuple(index1);
                next_state = bestS[j-1][(j-1)-k+1][index1];//consider len of single
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_CtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestC[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_S_CtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestS[k-1][(k-1)-(i+1)+1][index1];//consider len of single
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j-1][index2];
                stk.push(make_tuple(j-1 ,index2, next_state2));
                break;
            case MANNER_CStoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state1 = bestCS[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state1));
                break;
            case MANNER_C_StoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                std::tie(q1, nucs, nuck) = GetIndexTuple(index2);
                q = q1-1;
                next_state1 = bestC[q][index1];
                stk.push(make_tuple(q ,index1, next_state1));
                next_state2 = bestS[j-1][(j-1)-(q1)+1][index2];//consider len of single
                stk.push(make_tuple(j-1 ,index2, next_state2));
                break;
            case MANNER_S_CStoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state2 = bestCS[j-1][index2];
                stk.push(make_tuple(j-1 ,index2, next_state2));
                std::tie(k, nuck, nuck_pair,nucl, nucx,posi) = GetIndexTupleCS(index2);//from CS to get S
                next_state1 = bestS[k-1][(k-1)-(i+1)+1][index1];//consider len of single
                stk.push(make_tuple(k-1 ,index1, next_state1));
                break;
            case MANNER_S_C_StoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                std::tie(q1,nuck,nucs) = GetIndexTuple(index3);//nuc is not important
                l1 = (j-1) - q1 + 1;
                next_state3 = bestS[j-1][l1][index3];//consider len of single
                stk.push(make_tuple(j-1 ,index3, next_state3));
                q = q1-1;

                next_state2 = bestC[q][index2];
                stk.push(make_tuple(q ,index2, next_state2));
                std::tie(p,nuck,nucs) = GetIndexTuple(index2);//nuc is not important
                p_1 = p-1;

                l2 = p_1 - (i+1)+ 1;
                next_state1 = bestS[p_1][l2][index1];//consider len of single
                stk.push(make_tuple(p_1 ,index1, next_state1));
                break;
            case MANNER_C_StoCS:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestC[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestS[j][len][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_CtoM1:
                next_state = bestC[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_M1_CtoM2:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestM1[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_M1_EtoM1:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestM1[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_M2toM1:
                next_state = bestM2[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_M2toMulti:
                next_state = bestM2[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_S_M2toMulti:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                std::tie(s, nucs, nucs_end) = GetIndexTuple(index1);
                next_state1 = bestS[k-1][(k-1)-s+1][index1];//consider len of single
                // next_state1 = bestN[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestM2[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_Multi_EtoMulti:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestMulti[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_MultitoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestMulti[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));

                // if(debug){

                //     struct_score += - v_score_multi(GET_ACGU_NUM_V(rna_solution[i]),GET_ACGU_NUM_V(rna_solution[j]));
                //     printf("Multi loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, rna_solution[i], rna_solution[j], struct_score / -100.0);

                // }
                break;
            case MANNER_F_EtoF:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestF[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_CtoF:
                next_state = bestC[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_F_CtoF:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestF[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            default:
                cout << "Something error" << endl;
            
        }
        //cout<<mannerToString(MANNER)<<", "<<state.score/lambda<<", "<<i+1<<", "<<j+1<<endl;
        
    }


    return maxscore;
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, Manner MANNER){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, -1, -1,MANNER);
        stateMap[index] = newState;
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, int preindex_2, Manner MANNER){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, preindex_2, -1,MANNER);
        stateMap[index] = newState;
    }
}

//new update fun()
template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, Manner MANNER, int last_pair_pos){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, -1, -1,MANNER, last_pair_pos);
        stateMap[index] = newState;
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, int preindex_2, Manner MANNER, int last_pair_pos){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, preindex_2, -1, MANNER, last_pair_pos);
        stateMap[index] = newState;
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, int preindex_2, int preindex_3, Manner MANNER){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, preindex_2, preindex_3, MANNER);
        stateMap[index] = newState;
    }
}
//end new update fun()
//end new update fun()

inline bool isLastNuc(int pos){
    if(pos%3 == 2)
        return true;
    return false;
}

inline double getCAI(char Amino_label, int nuc){
    return CodonSetCAIMap[Amino_label][nuc];
}

inline double getCAI_DN(char Amino_label, int nuc){
    return CodonSetCAIMap_DN[Amino_label][nuc];
}

inline bool compare(double x, double y, double epsilon = 1.0E-5) {
    if(fabs(x - y) < epsilon) return true;
    return false;
}

template <typename T>
void initialize_Special_HP_DN(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq){
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    ACMatcher<char,256> matcher;
            
    for(auto sp_loop : sp_loops_vec){
        string s = sp_loop;
        matcher.AddString(s.c_str(),s.length());
    }
            
    if (!matcher.MakeTree()) {
        std::cout << "Failed to build the Aho-Corasick tree" << std::endl;
        return ;
    }

    
    unordered_map<int,unordered_set<string>> global_sp;//<index, all special hairpins could appear>
    for(int i = 0;i < ami_seq.length()-3;i++){
        vector<string> &Amino_1 = Amino_to_nucs[ami_seq[i]];
        vector<string> &Amino_2 = Amino_to_nucs[ami_seq[i+1]];
        vector<string> &Amino_3 = Amino_to_nucs[ami_seq[i+2]];
        vector<string> &Amino_4 = Amino_to_nucs[ami_seq[i+3]];
        for(auto &A1 : Amino_1)
            for(auto &A2 : Amino_2)
                for(auto &A3 : Amino_3){
                    for(auto &A4 : Amino_4){
                        std::string sub_seq = A1 + A2 + A3 + A4;
                        //cout<<sub_seq<<endl;
                        unordered_map<int,int> sp;
                        matcher.SearchPositions(sub_seq.c_str(),i*3,0,sub_seq.length()-1,sp);
                        for(auto ss:sp){
                            string special_hp = sub_seq.substr(ss.first - i*3 - ss.second +1,ss.second);
                            if(global_sp[ss.first].count(special_hp) == 0)
                                global_sp[ss.first].insert(special_hp);
                            //cout<<"position "<<ss.first<<" have sp loop :"<<sub_seq.substr(ss.first - i*3 - ss.second +1,ss.second)<<endl;
                        }
                    }
                }
        
    }
    
    for(auto g:global_sp){
        // cout<<"position "<<g.first<<" have sp loop(s) : ";
        for(auto h : g.second){
            vector<int> last_nuc ,first_nuc;
            int last_idx = g.first,first_idx = last_idx - h.length() +1;
            unordered_map<int,T> first_CAI;
            // cout<<"{"<<h <<", "<<sp_loops[h]<< ", ";
            if(isLastNuc(last_idx)){//third place in a codon
                string codon_ = h.substr(h.length()-3,3);
                last_nuc.push_back(getLastExtendedNuc(re_Amino_label[codon_][0],codon_));
                //cout<<CheckReBASE(getLastExtendedNuc(re_Amino_label[codon_][0],codon_))<<endl;//Not finish
            }
            else{
                vector<int> &nuc_list = Base_table.at(rna_seq[last_idx]);
                vector<int> &nuc_list2 = Base_table.at(rna_seq[last_idx-1]);
                for(auto &nuc : nuc_list){
                    if(reBASE(nuc) == h[h.length()-1]){
                        for(auto &nuc2 : nuc_list2){
                            if(reBASE(nuc2) == h[h.length()-2] && IsLegal(con_seq,nuc2,nuc,last_idx -1))
                                last_nuc.push_back(nuc);
                        }
                    }
                    
                }
            }

            if(first_idx % 3 == 1){//second place in a codon
                vector<int> &nuc_list = Base_table.at(rna_seq[first_idx]);
                vector<int> &nuc_list2 = Base_table.at(rna_seq[first_idx+1]);
                for(auto &nuc : nuc_list){
                    if(reBASE(nuc) == h[0]){
                        for(auto &nuc2 : nuc_list2){
                            if(reBASE(nuc2) == h[1] && IsLegal(con_seq,nuc,nuc2,first_idx)){
                                first_nuc.push_back(nuc);
                                first_CAI[nuc] = getCAI_DN(ami_seq[(first_idx+1)/3], nuc2);
                            }
                                
                        }
                    }
                }
            }
            else{
                vector<int> &nuc_list = Base_table.at(rna_seq[first_idx]);
                for(auto &nuc : nuc_list){
                    if(reBASE(nuc) == h[0]){
                        first_nuc.push_back(nuc);
                    }
                }
            }
            
            std::unordered_map<int, State<T>>& bestC_j = bestC[last_idx];
            for(auto &nuci : first_nuc){
                for(auto &nucj : last_nuc){
                    T newscore = -lambda * sp_loops[h];
                    for(int x = first_idx ; x <= last_idx ; x++){
                        if(isLastNuc(x) && (1-lambda)){
                            if(x - first_idx == 0){//special case1
                                newscore += getCAI_DN(ami_seq[x/3], nuci);
                            }
                            else if(x - first_idx == 1){//special case2
                                newscore += first_CAI[nuci];
                            }
                            else{
                                string codon_ = h.substr(x-first_idx-2,3);
                                int nucx = getLastExtendedNuc(re_Amino_label[codon_][0],codon_);
                                newscore += getCAI_DN(ami_seq[x/3], nucx);
                            } 
                        }
                            
                    }
                    int index = GetIndex(first_idx, nuci, nucj);
                    update(bestC_j, index, newscore, -1, MANNER_S_HP);
                    if(bestC_j[index].score == newscore)
                        sp_backtrack[index] = h;
                }
            }
        }        
    }
}

template<typename T>
void initialize_Special_HP_LD(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq){
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    ACMatcher<char,256> matcher;
            
    for(auto sp_loop : sp_loops_vec){
        string s = sp_loop;
        matcher.AddString(s.c_str(),s.length());
    }
            
    if (!matcher.MakeTree()) {
        std::cout << "Failed to build the Aho-Corasick tree" << std::endl;
        return ;
    }

    
    unordered_map<int,unordered_set<string>> global_sp;//<index, all special hairpins could appear>
    for(int i = 0;i < ami_seq.length()-3;i++){
        vector<string> &Amino_1 = Amino_to_nucs[ami_seq[i]];
        vector<string> &Amino_2 = Amino_to_nucs[ami_seq[i+1]];
        vector<string> &Amino_3 = Amino_to_nucs[ami_seq[i+2]];
        vector<string> &Amino_4 = Amino_to_nucs[ami_seq[i+3]];
        for(auto &A1 : Amino_1)
            for(auto &A2 : Amino_2)
                for(auto &A3 : Amino_3){
                    for(auto &A4 : Amino_4){
                        std::string sub_seq = A1 + A2 + A3 + A4;
                        //cout<<sub_seq<<endl;
                        unordered_map<int,int> sp;
                        matcher.SearchPositions(sub_seq.c_str(),i*3,0,sub_seq.length()-1,sp);
                        for(auto ss:sp){
                            string special_hp = sub_seq.substr(ss.first - i*3 - ss.second +1,ss.second);
                            if(global_sp[ss.first].count(special_hp) == 0)
                                global_sp[ss.first].insert(special_hp);
                            //cout<<"position "<<ss.first<<" have sp loop :"<<sub_seq.substr(ss.first - i*3 - ss.second +1,ss.second)<<endl;
                        }
                    }
                }
        
    }
    
    for(auto g:global_sp){
        // cout<<"position "<<g.first<<" have sp loop(s) : ";
        for(auto h : g.second){
            vector<int> last_nuc ,first_nuc;
            int last_idx = g.first,first_idx = last_idx - h.length() +1;
            unordered_map<int,T> first_CAI;
            // cout<<"{"<<h <<", "<<sp_loops[h]<< ", ";
            if(isLastNuc(last_idx)){//third place in a codon
                string codon_ = h.substr(h.length()-3,3);
                last_nuc.push_back(getLastExtendedNuc(re_Amino_label[codon_][0],codon_));
                //cout<<CheckReBASE(getLastExtendedNuc(re_Amino_label[codon_][0],codon_))<<endl;//Not finish
            }
            else{
                vector<int> &nuc_list = Base_table.at(rna_seq[last_idx]);
                vector<int> &nuc_list2 = Base_table.at(rna_seq[last_idx-1]);
                for(auto &nuc : nuc_list){
                    if(reBASE(nuc) == h[h.length()-1]){
                        for(auto &nuc2 : nuc_list2){
                            if(reBASE(nuc2) == h[h.length()-2] && IsLegal(con_seq,nuc2,nuc,last_idx -1))
                                last_nuc.push_back(nuc);
                        }
                    }
                    
                }
            }

            if(first_idx % 3 == 1){//second place in a codon
                vector<int> &nuc_list = Base_table.at(rna_seq[first_idx]);
                vector<int> &nuc_list2 = Base_table.at(rna_seq[first_idx+1]);
                for(auto &nuc : nuc_list){
                    if(reBASE(nuc) == h[0]){
                        for(auto &nuc2 : nuc_list2){
                            if(reBASE(nuc2) == h[1] && IsLegal(con_seq,nuc,nuc2,first_idx)){
                                first_nuc.push_back(nuc);
                                first_CAI[nuc] = getCAI(ami_seq[(first_idx+1)/3], nuc2);
                            }
                                
                        }
                    }
                }
            }
            else{
                vector<int> &nuc_list = Base_table.at(rna_seq[first_idx]);
                for(auto &nuc : nuc_list){
                    if(reBASE(nuc) == h[0]){
                        first_nuc.push_back(nuc);
                    }
                }
            }
            std::unordered_map<int, State<T>>& bestC_j = bestC[last_idx];
            for(auto &nuci : first_nuc){
                for(auto &nucj : last_nuc){
                    T newscore = -1 * sp_loops[h];
                    for(int x = first_idx ; x <= last_idx ; x++){
                        if(isLastNuc(x) && (lambda)){
                            if(x - first_idx == 0){//special case1
                                newscore += getCAI(ami_seq[x/3], nuci);
                            }
                            else if(x - first_idx == 1){//special case2
                                newscore += first_CAI[nuci];
                            }
                            else{
                                string codon_ = h.substr(x-first_idx-2,3);
                                int nucx = getLastExtendedNuc(re_Amino_label[codon_][0],codon_);
                                newscore += getCAI(ami_seq[x/3], nucx);
                            } 
                        }
                            
                    }
                    int index = GetIndex(first_idx, nuci, nucj);
                    update(bestC_j, index, newscore, -1, MANNER_S_HP);
                    if(bestC_j[index].score == newscore)
                        sp_backtrack[index] = h;
                }
            }
        }        
    }
}

template<typename T>
void LCDSfoldCAI_LD_exact(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq){
    //LinearDesign
    int seq_length = rna_seq.size();
    T newscore;
    int index;
    int preindex;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    vector<int> nuc0_list = Base_table.at(rna_seq[0]); //j list
    for(int c0 = 0; c0 < nuc0_list.size(); ++c0){
        int nuc0 = nuc0_list[c0];
        index = GetIndex(0, nuc0, nuc0);
        update(bestN[0], index, 0, -1, MANNER_NONEtoN);
        update(bestF[0], index, 0, -1, MANNER_NONEtoF);
        update(bestS[0][1], index, 0, -1, MANNER_NONEtoS);
    }
    
    for(int j = 1; j < seq_length; ++j){// 1 ~ seq_length
        Processing(seq_length);
        std::unordered_map<int, State<T>>& bestN_j = bestN[j];
        std::vector<std::unordered_map<int, State<T>>>& bestS_j = bestS[j];
        std::unordered_map<int, State<T>>& bestF_j = bestF[j];
        std::unordered_map<int, State<T>>& bestC_j = bestC[j];
        std::unordered_map<int, State<T>>& bestCS_j = bestCS[j];
        std::unordered_map<int, State<T>>& bestM1_j = bestM1[j];
        std::unordered_map<int, State<T>>& bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>>& bestMulti_j = bestMulti[j];
        
        vector<int> nucj_list = Base_table.at(rna_seq[j]); //j list
        int j_1 = j - 1;
        for (auto &nucj : nucj_list){// all possible nuc at position j 
            
            // NONE->S
            // NONE->N
            index = GetIndex(j, nucj, nucj);
            newscore = 0;

            if(isLastNuc(j) && (lambda))
                newscore += getCAI(ami_seq[j/3], nucj);
            
            update(bestS_j[1], index, newscore, -1, MANNER_NONEtoS);// initialize bestS[j]
            update(bestN_j, index, newscore, -1, MANNER_NONEtoN);// initialize bestN[j] 

            /*all N -> another state*/
            std::unordered_map<int, State<T>>& bestN_j_1 = bestN[j_1];
            for (auto &itemN_j_1 : bestN_j_1) {// all different states -> N at position j-1
                int index_j_1 = itemN_j_1.first;// key
                State<T> stateN_j_1 = itemN_j_1.second;// state
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);// from which position i and nuc of i and j-1
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateN_j_1.score;
                    
                    if(isLastNuc(j) && (lambda))
                        newscore += getCAI(ami_seq[j/3], nucj);
                    // N -> N
                    index = GetIndex(i, nuci, nucj);
                    update(bestN_j, index, newscore, index_j_1, MANNER_N_EtoN);

                    int i_1 = i - 1;//the continuous Non-pairing is from i to j-1, so if we want to transfer N -> C ,the pairing must be i-1 and j (loop : i ~ j-1)
                    if((i_1 >= 0) && (j - i_1 > HAIRPIN_GAP)){//consider hairpin(gap > 3)
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                            if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                                if(_allowed_pairs[nuci_1][nucj]){
                                    //because N -> C means Non-paired to Closed, need to consider hairpin structure
                                    newscore = -v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)) + stateN_j_1.score;
                                    if(isLastNuc(j) && (lambda))
                                        newscore += getCAI(ami_seq[j/3], nucj);
                                    if(isLastNuc(i_1) && (lambda))
                                        newscore += getCAI(ami_seq[i_1/3], nuci_1);
                                    //N->C
                                    index = GetIndex(i_1, nuci_1, nucj);
                                    update(bestC_j, index, newscore, index_j_1, MANNER_NtoC);
                            }}
                }}}
            }
        }

        /*all S -> another state*/
        for (auto &nucj : nucj_list){
            int counter = 0;
            for(int l = 1; l <= min(SINGLE_MAX_LEN-1,j) ; l++){
                std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1][l];
                for (auto &itemS_j_1 : bestS_j_1) {
                    int index_j_1 = itemS_j_1.first;
                    State<T> stateS_j_1 = itemS_j_1.second;
                    int i, nucj_1, nuci;
                    std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                    if(IsLegal(con_seq, nucj_1, nucj, j_1) && (j - i) < SINGLE_MAX_LEN){
                        //S means single length extend
                        newscore = stateS_j_1.score;
                        if(isLastNuc(j) && (lambda))
                            newscore += getCAI(ami_seq[j/3], nucj);
                        // S -> S
                        index = GetIndex(i, nuci, nucj);
                        update(bestS_j[l+1], index, newscore, index_j_1, MANNER_S_EtoS);
                    }
                }
            }
        }

        for (auto &nucj : nucj_list){
            for(int l1 = min(SINGLE_MAX_LEN,j-4) ; l1 >= 1; --l1){
                std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1][l1];
                for (auto &itemS_j_1 : bestS_j_1) {
                    int index_j_1 = itemS_j_1.first;
                    State<T> stateS_j_1 = itemS_j_1.second;
                    int q1, nucq1, nucj_1;
                    std::tie(q1, nucq1, nucj_1) = GetIndexTuple(index_j_1);
                    int q = q1 - 1;
                    
                    if(q >= 0 &&  IsLegal(con_seq, nucj_1, nucj, j_1)){
                        std::unordered_map<int, State<T>>& bestC_q = bestC[q];    
                        for (auto &itemC_q : bestC_q) {
                            int index_q = itemC_q.first;
                            State<T> stateC_q = itemC_q.second;
                            int p, nucp, nucq;    
                            std::tie(p, nucp, nucq) = GetIndexTuple(index_q);
                            int p_1 = p - 1;

                            if(p_1 >= 0 && IsLegal(con_seq, nucq, nucq1, q)){
                                //C_StoS (right bulge)
                                int i = p_1;
                                vector<int> nuci_list = Base_table.at(rna_seq[i]); //i list
                                for(auto &nuci : nuci_list){
                                    if(IsLegal(con_seq, nuci, nucp, i) && _allowed_pairs[nuci][nucj]){
                                        //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                                        newscore = -v_score_single(i,j,p,q,EnrBASE(nuci),EnrBASE(nucp),EnrBASE(nucj_1),EnrBASE(nucj),EnrBASE(nuci),EnrBASE(nucp),EnrBASE(nucq),EnrBASE(nucq1))
                                        + stateS_j_1.score + stateC_q.score;  

                                        if(isLastNuc(j) && (lambda))
                                            newscore += getCAI(ami_seq[j/3], nucj);
                                        if(isLastNuc(i) && (lambda))
                                            newscore += getCAI(ami_seq[i/3], nuci);

                                        index = GetIndex(i, nuci, nucj);
                                        update(bestC_j, index, newscore, index_q, index_j_1, MANNER_C_StoC);
                                    }
                                }//end nuci list


                                //S_C_StoC (internal loop)
                                for(int l2 = 1; l2 <= min(SINGLE_MAX_LEN-l1 , p_1) ; ++l2){
                                    std::unordered_map<int, State<T>>& bestS_p_1 = bestS[p_1][l2];
                                    for (auto &itemS_p_1 : bestS_p_1) {
                                        int index_p_1 = itemS_p_1.first;
                                        State<T> stateS_p_1 = itemS_p_1.second;
                                        int i1, nuci1, nucp_1;
                                        std::tie(i1, nuci1, nucp_1) = GetIndexTuple(index_p_1);
                                        int i = i1 - 1;
                                        if(i >= 0 && IsLegal(con_seq, nucp_1, nucp, p_1)){
                                            vector<int> nuci_list = Base_table.at(rna_seq[i]); //i list
                                            for(auto &nuci : nuci_list){
                                                if(IsLegal(con_seq, nuci, nuci1, i) && _allowed_pairs[nuci][nucj]){
                                                    //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                                                    newscore = -v_score_single(i,j,p,q,EnrBASE(nuci),EnrBASE(nuci1),EnrBASE(nucj_1),EnrBASE(nucj),EnrBASE(nucp_1),EnrBASE(nucp),EnrBASE(nucq),EnrBASE(nucq1))
                                                    + stateS_j_1.score + stateC_q.score + stateS_p_1.score;  

                                                    if(isLastNuc(j) && (lambda))
                                                        newscore += getCAI(ami_seq[j/3], nucj);
                                                    if(isLastNuc(i) && (lambda))
                                                        newscore += getCAI(ami_seq[i/3], nuci);

                                                    index = GetIndex(i, nuci, nucj);
                                                    update(bestC_j, index, newscore, index_p_1, index_q, index_j_1, MANNER_S_C_StoC);
                                                }
                                            }//end nuci list
                                        
                                        }
                                        
                                    }//end bestS[p_1][l2]
                                }
                            }
                        }//end bestC_q
                    }
            
                }//end bestS[j_1]
            }
        }//end of nucj
        
        for (auto &nucj : nucj_list){

            std::unordered_map<int, State<T>>& bestC_j_1 = bestC[j_1];
            for (auto &itemC_j_1 : bestC_j_1) {
                int index_j_1 = itemC_j_1.first;
                State<T> stateC_j_1 = itemC_j_1.second;
                int i1, nucj_1, nuci1;
                std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);
                int i = i1 - 1;

                if(IsLegal(con_seq, nucj_1, nucj, j_1) && (i >= 0) && (j - i > (HAIRPIN_GAP + 2))){
                    
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i-1 list
                    for(auto &nuci : nuci_list){
                        
                        if(_allowed_pairs[nuci][nucj] && IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j_1)){
                            
                            newscore = - v_score_single(i, j, i1, j_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score;

                            if(isLastNuc(j) && (lambda))
                                newscore += getCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && (lambda))
                                newscore += getCAI(ami_seq[i/3], nuci);
                            
                            //C->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CtoC);
                            
                            
                    }}

                    int i_1 = i - 1;
                    
                    if(i_1 >= 0){
                        for(int l = 1;l<=min(SINGLE_MAX_LEN,i);l++){
                            std::unordered_map<int, State<T>>& bestS_i = bestS[i][l];
                            for (auto &itemS_i : bestS_i) {
                                int index_i = itemS_i.first;
                                State<T> stateS_i = itemS_i.second;
                                int k1, nuck1, nuci;
                                std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                                int k = k1 - 1;
                                if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i)){
                                    vector<int> nuck_list = Base_table.at(rna_seq[k]); //k list
                                    for(auto &nuck : nuck_list){
                                        
                                        if(_allowed_pairs[nuck][nucj] && IsLegal(con_seq, nuck, nuck1, k)){
                                            
                                            newscore = - v_score_single(k, j, i1, j_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score + stateS_i.score;

                                            if(isLastNuc(j) && (lambda))
                                                newscore += getCAI(ami_seq[j/3], nucj);
                                            if(isLastNuc(k) && (lambda))
                                                newscore += getCAI(ami_seq[k/3], nuck);
                                            
                                            //S+C->C
                                            index = GetIndex(k, nuck, nucj);
                                            update(bestC_j, index, newscore, index_i, index_j_1, MANNER_S_CtoC);

                            }}}}
                        }
                        
                    }
                    
            }} //Cj_1
            
        } //nucj
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestMulti_j_1 = bestMulti[j_1];
            for (auto &itemMulti_j_1 : bestMulti_j_1) { 
                int index_j_1 = itemMulti_j_1.first;
                State<T> stateMulti_j_1 = itemMulti_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    
                    int last_pair_pos_Multi = stateMulti_j_1.last_pair_pos;
                    if(j - last_pair_pos_Multi <= SINGLE_MAX_LEN){
                        
                        newscore = stateMulti_j_1.score;

                        if(isLastNuc(j) && lambda)
                            newscore += getCAI(ami_seq[j/3], nucj);

                        //Multi->Multi
                        index = GetIndex(i, nuci, nucj);
                        update(bestMulti_j, index, newscore, index_j_1, MANNER_Multi_EtoMulti, last_pair_pos_Multi);
                    }
                    
                    int i_1 = i - 1;
                    if(i_1 >= 0){
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                        if(IsLegal(con_seq, nuci_1, nuci, i_1) && _allowed_pairs[nuci_1][nucj]){
                            newscore = - v_score_multi(EnrBASE(nuci_1), EnrBASE(nucj)) + stateMulti_j_1.score;
                            
                            if(isLastNuc(j) && (lambda))
                                newscore += getCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i_1) && (lambda))
                                newscore += getCAI(ami_seq[i_1/3], nuci_1);

                            //Multi->C
                            index = GetIndex(i_1, nuci_1, nucj);
                            //newscore = INT32_MIN;//test
                            update(bestC_j, index, newscore, index_j_1, MANNER_MultitoC);
                        }}
                    }
                }
            }
        } //nucj
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);

        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = - v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score;

            //C->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_CtoM1, j);

            int i_1 = i - 1;
            if(i_1 >= HAIRPIN_GAP + 1){
                std::unordered_map<int, State<T>>& bestM1_i_1 = bestM1[i_1];
                for (auto &itemM1_i_1 : bestM1_i_1) {
                    int index_i_1 = itemM1_i_1.first;
                    State<T> stateM1_i_1 = itemM1_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = - v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score + stateM1_i_1.score;
                        // M1+C->M2
                        index = GetIndex(k, nuck, nucj);
                        update(bestM2_j, index, newscore, index_i_1 ,index_j, MANNER_M1_CtoM2, j);
                    }
                }
            }
        }
        
        BeamPrune(con_seq, rna_seq, bestM2_j, bestF, false);

        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestM1_j_1 = bestM1[j_1];
            for (auto &itemM1_j_1 : bestM1_j_1) {
                int index_j_1 = itemM1_j_1.first;
                State<T> stateM1_j_1 = itemM1_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);

                int last_pair_pos_M1 = stateM1_j_1.last_pair_pos;
                assert(last_pair_pos_M1 != -1);
                if(j - last_pair_pos_M1 <= SINGLE_MAX_LEN && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateM1_j_1.score;

                    if(isLastNuc(j) && (lambda))
                        newscore += getCAI(ami_seq[j/3], nucj);

                    //M1->M1
                    index = GetIndex(i, nuci, nucj);
                    update(bestM1_j, index, newscore, index_j_1, MANNER_M1_EtoM1, last_pair_pos_M1);
                }
            }
        } //nucj

        for (auto &itemM2_j : bestM2_j) {
            int index_j = itemM2_j.first;
            State<T> stateM2_j = itemM2_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = stateM2_j.score;
            int last_pair_pos_M2 = stateM2_j.last_pair_pos;

            //M2->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_M2toM1, last_pair_pos_M2);
            //M2->Multi
            update(bestMulti_j, index_j, newscore, index_j, MANNER_M2toMulti, last_pair_pos_M2);
            int i_1 = i - 1;
            if(i_1 >= 0){
                for(int l = 1; l <= min(SINGLE_MAX_LEN,i_1); l++){
                    std::unordered_map<int, State<T>>& bestS_i_1 = bestS[i_1][l];
                    for (auto &itemS_i_1 : bestS_i_1) {
                        int index_i_1 = itemS_i_1.first;
                        State<T> stateS_i_1 = itemS_i_1.second;
                        int k, nuck, nuci_1;
                        std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                        if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                            newscore = stateM2_j.score + stateS_i_1.score;
                            //S+M2->Multi
                            index = GetIndex(k, nuck, nucj);
                            update(bestMulti_j, index, newscore, index_i_1, index_j, MANNER_S_M2toMulti, last_pair_pos_M2);
                        }
                    }
                }
                
                
            }
            
        }


        {//conclusion!
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestF_j_1 = bestF[j_1];
            for (auto &itemF_j_1 : bestF_j_1) {
                int index_j_1 = itemF_j_1.first;
                State<T> stateF_j_1 = itemF_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateF_j_1.score;
                    if(isLastNuc(j) && (lambda))
                        newscore += getCAI(ami_seq[j/3], nucj);
                    //F->F
                    index = GetIndex(i, nuci, nucj);
                    update(bestF_j, index, newscore, index_j_1, MANNER_F_EtoF);
                }
            }
        }
        
        
        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);

            if(i == 0 ){
                newscore = stateC_j.score - v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                //C->F
                update(bestF_j, index_j, newscore, index_j, MANNER_CtoF);
            }

            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestF_i_1 = bestF[i_1];
                for (auto &itemF_i_1 : bestF_i_1) {
                    int index_i_1 = itemF_i_1.first;
                    State<T> stateF_i_1 = itemF_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if((k == 0 || j == seq_length-1) && IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateF_i_1.score + stateC_j.score - v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                        //F+C->F
                        index = GetIndex(k, nuck, nucj);
                        update(bestF_j, index, newscore, index_i_1, index_j, MANNER_F_CtoF);
            }}}
        }
        
        //BeamPrune(con_seq, rna_seq, bestF_j, bestF, false);
        }

    } //jend
}

template<typename T>
void LCDSfoldCAI_LD_beam(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq){
    //LinearDesign
    int seq_length = rna_seq.size();
    T newscore;
    int index;
    int preindex;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    vector<int> nuc0_list = Base_table.at(rna_seq[0]); //j list
    for(int c0 = 0; c0 < nuc0_list.size(); ++c0){
        int nuc0 = nuc0_list[c0];
        index = GetIndex(0, nuc0, nuc0);
        update(bestN[0], index, 0, -1, MANNER_NONEtoN);
        update(bestF[0], index, 0, -1, MANNER_NONEtoF);
        update(bestS[0][1], index, 0, -1, MANNER_NONEtoS);
    }
    
    for(int j = 1; j < seq_length; ++j){// 1 ~ seq_length
        Processing(seq_length);
        std::unordered_map<int, State<T>>& bestN_j = bestN[j];
        std::vector<std::unordered_map<int, State<T>>>& bestS_j = bestS[j];
        std::unordered_map<int, State<T>>& bestF_j = bestF[j];
        std::unordered_map<int, State<T>>& bestC_j = bestC[j];
        std::unordered_map<int, State<T>>& bestCS_j = bestCS[j];
        std::unordered_map<int, State<T>>& bestM1_j = bestM1[j];
        std::unordered_map<int, State<T>>& bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>>& bestMulti_j = bestMulti[j];
        
        vector<int> nucj_list = Base_table.at(rna_seq[j]); //j list
        int j_1 = j - 1;
        for (auto &nucj : nucj_list){// all possible nuc at position j 
            
            // NONE->S
            // NONE->N
            index = GetIndex(j, nucj, nucj);
            newscore = 0;

            if(isLastNuc(j) && (lambda))
                newscore += getCAI(ami_seq[j/3], nucj);
            
            update(bestS_j[1], index, newscore, -1, MANNER_NONEtoS);// initialize bestS[j]
            update(bestN_j, index, newscore, -1, MANNER_NONEtoN);// initialize bestN[j] 

            /*all N -> another state*/
            std::unordered_map<int, State<T>>& bestN_j_1 = bestN[j_1];
            for (auto &itemN_j_1 : bestN_j_1) {// all different states -> N at position j-1
                int index_j_1 = itemN_j_1.first;// key
                State<T> stateN_j_1 = itemN_j_1.second;// state
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);// from which position i and nuc of i and j-1
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    //because N -> N means Non-pairing to Non-pairing, only need to consider CAI(at position j)
                    newscore = stateN_j_1.score;
                    
                    if(isLastNuc(j) && (lambda))
                        newscore += getCAI(ami_seq[j/3], nucj);
                    // N -> N
                    index = GetIndex(i, nuci, nucj);
                    update(bestN_j, index, newscore, index_j_1, MANNER_N_EtoN);

                    int i_1 = i - 1;//the continuous Non-pairing is from i to j-1, so if we want to transfer N -> C ,the pairing must be i-1 and j (loop : i ~ j-1)
                    if((i_1 >= 0) && (j - i_1 > HAIRPIN_GAP)){//consider hairpin(gap > 3)
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                            if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                                if(_allowed_pairs[nuci_1][nucj]){
                                    //because N -> C means Non-pairing to Closing, need to consider hairpin structure
                                    newscore = -v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)) + stateN_j_1.score;
                                    if(isLastNuc(j) && (lambda))
                                        newscore += getCAI(ami_seq[j/3], nucj);
                                    if(isLastNuc(i_1) && (lambda))
                                        newscore += getCAI(ami_seq[i_1/3], nuci_1);
                                    //N->C
                                    index = GetIndex(i_1, nuci_1, nucj);
                                    update(bestC_j, index, newscore, index_j_1, MANNER_NtoC);
                            }}
                }}}
            }
        }
        BeamPrune(con_seq, rna_seq, bestN_j, bestF, !(lambda));
         /*all S -> another state*/
        
        for (auto &nucj : nucj_list){
            int counter = 0;
            for(int l = 1; l <= min(SINGLE_MAX_LEN-1,j) ; l++){
                std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1][l];
                for (auto &itemS_j_1 : bestS_j_1) {
                    int index_j_1 = itemS_j_1.first;
                    State<T> stateS_j_1 = itemS_j_1.second;
                    int i, nucj_1, nuci;
                    std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                    if(IsLegal(con_seq, nucj_1, nucj, j_1) && (j - i) < SINGLE_MAX_LEN){
                        //S means single length extend
                        newscore = stateS_j_1.score;
                        if(isLastNuc(j) && (lambda))
                            newscore += getCAI(ami_seq[j/3], nucj);
                        // S -> S
                        index = GetIndex(i, nuci, nucj);
                        update(bestS_j[l+1], index, newscore, index_j_1, MANNER_S_EtoS);
                    }
                }
            }
        }
        
        // for(int l = 1; l <= (SINGLE_MAX_LEN,j) ; l++)
        //     BeamPrune(con_seq, rna_seq, bestS_j[l], bestF, !(1-lambda));
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestCS_j_1 = bestCS[j_1];
            for (auto &itemCS_j_1 : bestCS_j_1) {
                int index_j_1 = itemCS_j_1.first;
                State<T> stateCS_j_1 = itemCS_j_1.second;
                int i1, nucj_1, nuci1, nuci1_pair, nucx,len;
                std::tie(i1, nuci1, nuci1_pair, nucj_1, nucx, len) = GetIndexTupleCS(index_j_1);//make_tuple(i, nuci, nucj, nucx, len)
                int i = i1 - 1;
                if(i >= 0 && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i list
                    for(auto &nuci : nuci_list){
                        if(IsLegal(con_seq, nuci, nuci1, i) && _allowed_pairs[nuci][nucj]){
                            //CS means the closing pair (nucleotide p,nucleotide q) plus all POSSIBLE single extend.
                            int pre_index_2 = stateCS_j_1.index_2;
                            int pre_index_1 = stateCS_j_1.index_1;
                            int s, nucs, nucl;
                            std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                            int s_1 = s - 1;
                            int a, nuca, nucs_1;
                            std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);
                            //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                            newscore = -v_score_single(i, j, i1, s_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                                EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score;
                            
                            if(isLastNuc(j) && (lambda))
                                newscore += getCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && (lambda))
                                newscore += getCAI(ami_seq[i/3], nuci);
                            
                            //CS->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CStoC);//left bulge
                        }
                    }
                    //cout<<"solve state S + CS to C"<<endl;
                    
                    for(int l = 1; l <= min(SINGLE_MAX_LEN-len,i) ; l++){
                        std::unordered_map<int, State<T>>& bestS_i = bestS[i][l];
                        for (auto &itemS_i : bestS_i) {
                            int index_i = itemS_i.first;
                            State<T> stateS_i = itemS_i.second;
                            int k1, nuck1, nuci;
                            std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                            int k = k1 - 1;
                            
                            if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i) && i - k + len <= SINGLE_MAX_LEN ){
                                vector<int> nuck_list = Base_table.at(rna_seq[k]); //i list
                                for(auto &nuck : nuck_list){
                                    if(IsLegal(con_seq, nuck, nuck1, k) && _allowed_pairs[nuck][nucj]){
                                        int pre_index_2 = stateCS_j_1.index_2;
                                        int pre_index_1 = stateCS_j_1.index_1;
                                        int s, nucs, nucl;
                                        std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                                        int s_1 = s - 1;
                                        int a, nuca, nucs_1;
                                        std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);
                                        //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                                        newscore = - v_score_single(k, j, i1, s_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score + stateS_i.score;
                                        
                                        if(isLastNuc(j) && (lambda))
                                            newscore += getCAI(ami_seq[j/3], nucj);
                                        if(isLastNuc(k) && (lambda))
                                            newscore += getCAI(ami_seq[k/3], nuck);
                                        
                                        //S+CS->C
                                        index = GetIndex(k, nuck, nucj);
                                        update(bestC_j, index, newscore, index_i ,index_j_1, MANNER_S_CStoC);
                                    }
                                }
                            }
                            
                            
                        }
                    }
                    
                }
            }
        }
        
        // BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        for(int l = min(SINGLE_MAX_LEN,j-4) ; l >= 1; --l){
            for (auto &itemS_j : bestS_j[l]) {
                int index_j = itemS_j.first;
                State<T> stateS_j = itemS_j.second;
                int i, nucj, nuci;
                std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
                int i_1 = i - 1;
                
                if(i_1 >= 0){
                    std::unordered_map<int, State<T>>& bestC_i_1 = bestC[i_1];
                    
                    for (auto &itemC_i_1 : bestC_i_1) {
                        int index_i_1 = itemC_i_1.first;
                        State<T> stateC_i_1 = itemC_i_1.second;
                        int k, nuck, nuci_1;    
                        std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);

                        if(k >= 0 && IsLegal(con_seq, nuci_1, nuci, i_1)){
                            newscore = stateS_j.score + stateC_i_1.score;
                            // C+S->CS
                            // GetIndexCS(int i, int nuci, int nuci_pair,int nucj, int nucx ,int len)
                            index = GetIndexCS(k, nuck, nuci_1,nucj, nuci,l);
                            update(bestCS_j, index, newscore, index_i_1, index_j, MANNER_C_StoCS);
                        }
                    }
                }
        
            }
        }

        BeamPrune(con_seq, rna_seq, bestCS_j, bestF, false);
        
        for (auto &nucj : nucj_list){

            std::unordered_map<int, State<T>>& bestC_j_1 = bestC[j_1];
            for (auto &itemC_j_1 : bestC_j_1) {
                int index_j_1 = itemC_j_1.first;
                State<T> stateC_j_1 = itemC_j_1.second;
                int i1, nucj_1, nuci1;
                std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);
                int i = i1 - 1;

                
                
                if(IsLegal(con_seq, nucj_1, nucj, j_1) && (i >= 0) && (j - i > (HAIRPIN_GAP + 2))){
                    
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i-1 list
                    for(auto &nuci : nuci_list){
                        
                        if(_allowed_pairs[nuci][nucj] && IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j_1)){
                            
                            newscore = - v_score_single(i, j, i1, j_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score;

                            if(isLastNuc(j) && (lambda))
                                newscore += getCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && (lambda))
                                newscore += getCAI(ami_seq[i/3], nuci);
                            
                            //C->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CtoC);
                            
                            
                    }}

                    int i_1 = i - 1;
                    
                    if(i_1 >= 0){
                        for(int l = 1;l<=min(SINGLE_MAX_LEN,i);l++){
                            std::unordered_map<int, State<T>>& bestS_i = bestS[i][l];
                            for (auto &itemS_i : bestS_i) {
                                int index_i = itemS_i.first;
                                State<T> stateS_i = itemS_i.second;
                                int k1, nuck1, nuci;
                                std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                                int k = k1 - 1;
                                if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i)){
                                    vector<int> nuck_list = Base_table.at(rna_seq[k]); //k list
                                    for(auto &nuck : nuck_list){
                                        
                                        if(_allowed_pairs[nuck][nucj] && IsLegal(con_seq, nuck, nuck1, k)){
                                            
                                            newscore = - v_score_single(k, j, i1, j_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score + stateS_i.score;

                                            if(isLastNuc(j) && (lambda))
                                                newscore += getCAI(ami_seq[j/3], nucj);
                                            if(isLastNuc(k) && (lambda))
                                                newscore += getCAI(ami_seq[k/3], nuck);
                                            
                                            //S+C->C
                                            index = GetIndex(k, nuck, nucj);
                                            update(bestC_j, index, newscore, index_i, index_j_1, MANNER_S_CtoC);

                            }}}}
                        }
                        
                    }
                    
            }} //Cj_1
            
        } //nucj
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestMulti_j_1 = bestMulti[j_1];
            for (auto &itemMulti_j_1 : bestMulti_j_1) { 
                int index_j_1 = itemMulti_j_1.first;
                State<T> stateMulti_j_1 = itemMulti_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    
                    int last_pair_pos_Multi = stateMulti_j_1.last_pair_pos;
                    if(j - last_pair_pos_Multi <= SINGLE_MAX_LEN){
                        
                        newscore = stateMulti_j_1.score;

                        if(isLastNuc(j) && lambda)
                            newscore += getCAI(ami_seq[j/3], nucj);

                        //Multi->Multi
                        index = GetIndex(i, nuci, nucj);
                        update(bestMulti_j, index, newscore, index_j_1, MANNER_Multi_EtoMulti, last_pair_pos_Multi);
                    }
                    
                    int i_1 = i - 1;
                    if(i_1 >= 0){
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                        if(IsLegal(con_seq, nuci_1, nuci, i_1) && _allowed_pairs[nuci_1][nucj]){
                            newscore = - v_score_multi(EnrBASE(nuci_1), EnrBASE(nucj)) + stateMulti_j_1.score;
                            
                            if(isLastNuc(j) && (lambda))
                                newscore += getCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i_1) && (lambda))
                                newscore += getCAI(ami_seq[i_1/3], nuci_1);

                            //Multi->C
                            index = GetIndex(i_1, nuci_1, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_MultitoC);
                        }}
                    }
                }
            }
        } //nucj
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);

        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = - v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score;

            //C->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_CtoM1, j);

            int i_1 = i - 1;
            if(i_1 >= HAIRPIN_GAP + 1){
                std::unordered_map<int, State<T>>& bestM1_i_1 = bestM1[i_1];
                for (auto &itemM1_i_1 : bestM1_i_1) {
                    int index_i_1 = itemM1_i_1.first;
                    State<T> stateM1_i_1 = itemM1_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = - v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score + stateM1_i_1.score;
                        // M1+C->M2
                        index = GetIndex(k, nuck, nucj);
                        update(bestM2_j, index, newscore, index_i_1 ,index_j, MANNER_M1_CtoM2, j);
                    }
                }
            }
        }
        
        BeamPrune(con_seq, rna_seq, bestM2_j, bestF, false);

        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestM1_j_1 = bestM1[j_1];
            for (auto &itemM1_j_1 : bestM1_j_1) {
                int index_j_1 = itemM1_j_1.first;
                State<T> stateM1_j_1 = itemM1_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);

                int last_pair_pos_M1 = stateM1_j_1.last_pair_pos;
                assert(last_pair_pos_M1 != -1);
                if(j - last_pair_pos_M1 <= SINGLE_MAX_LEN && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateM1_j_1.score;

                    if(isLastNuc(j) && (lambda))
                        newscore += getCAI(ami_seq[j/3], nucj);

                    //M1->M1
                    index = GetIndex(i, nuci, nucj);
                    update(bestM1_j, index, newscore, index_j_1, MANNER_M1_EtoM1, last_pair_pos_M1);
                }
            }
        } //nucj

        for (auto &itemM2_j : bestM2_j) {
            int index_j = itemM2_j.first;
            State<T> stateM2_j = itemM2_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = stateM2_j.score;
            int last_pair_pos_M2 = stateM2_j.last_pair_pos;

            //M2->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_M2toM1, last_pair_pos_M2);
            //M2->Multi
            update(bestMulti_j, index_j, newscore, index_j, MANNER_M2toMulti, last_pair_pos_M2);
            int i_1 = i - 1;
            if(i_1 >= 0){
                for(int l = 1; l <= min(SINGLE_MAX_LEN,i_1); l++){
                    std::unordered_map<int, State<T>>& bestS_i_1 = bestS[i_1][l];
                    for (auto &itemS_i_1 : bestS_i_1) {
                        int index_i_1 = itemS_i_1.first;
                        State<T> stateS_i_1 = itemS_i_1.second;
                        int k, nuck, nuci_1;
                        std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                        if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                            newscore = stateM2_j.score + stateS_i_1.score;
                            //S+M2->Multi
                            index = GetIndex(k, nuck, nucj);
                            update(bestMulti_j, index, newscore, index_i_1, index_j, MANNER_S_M2toMulti, last_pair_pos_M2);
                        }
                    }
                }            
            }
            
        }

        BeamPrune(con_seq, rna_seq, bestM1_j, bestF, false);
        BeamPrune(con_seq, rna_seq, bestMulti_j, bestF, false);

        {//conclusion!
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestF_j_1 = bestF[j_1];
            for (auto &itemF_j_1 : bestF_j_1) {
                int index_j_1 = itemF_j_1.first;
                State<T> stateF_j_1 = itemF_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateF_j_1.score;
                    if(isLastNuc(j) && (lambda))
                        newscore += getCAI(ami_seq[j/3], nucj);
                    //F->F
                    index = GetIndex(i, nuci, nucj);
                    update(bestF_j, index, newscore, index_j_1, MANNER_F_EtoF);
                }
            }
        }
        
        
        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);

            if(i == 0 ){
                newscore = stateC_j.score - v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                //C->F
                update(bestF_j, index_j, newscore, index_j, MANNER_CtoF);
            }

            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestF_i_1 = bestF[i_1];
                for (auto &itemF_i_1 : bestF_i_1) {
                    int index_i_1 = itemF_i_1.first;
                    State<T> stateF_i_1 = itemF_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if((k == 0 || j == seq_length-1) && IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateF_i_1.score + stateC_j.score - v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                        //F+C->F
                        index = GetIndex(k, nuck, nucj);
                        update(bestF_j, index, newscore, index_i_1, index_j, MANNER_F_CtoF);
            }}}
        }
        
        //BeamPrune(con_seq, rna_seq, bestF_j, bestF, false);
        }

    } //jend
}

template<typename T>
void LCDSfoldCAI_DN_beam(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq){
    //DERNA
    int seq_length = rna_seq.size();
    T newscore;
    int index;
    int preindex;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    vector<int> nuc0_list = Base_table.at(rna_seq[0]); //j list
    for(int c0 = 0; c0 < nuc0_list.size(); ++c0){
        int nuc0 = nuc0_list[c0];
        index = GetIndex(0, nuc0, nuc0);
        update(bestN[0], index, 0, -1, MANNER_NONEtoN);
        update(bestF[0], index, 0, -1, MANNER_NONEtoF);
        update(bestS[0][1], index, 0, -1, MANNER_NONEtoS);
    }
    
    for(int j = 1; j < seq_length; ++j){// 1 ~ seq_length
        Processing(seq_length);
        std::unordered_map<int, State<T>>& bestN_j = bestN[j];
        std::vector<std::unordered_map<int, State<T>>>& bestS_j = bestS[j];
        std::unordered_map<int, State<T>>& bestF_j = bestF[j];
        std::unordered_map<int, State<T>>& bestC_j = bestC[j];
        std::unordered_map<int, State<T>>& bestCS_j = bestCS[j];
        std::unordered_map<int, State<T>>& bestM1_j = bestM1[j];
        std::unordered_map<int, State<T>>& bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>>& bestMulti_j = bestMulti[j];
        
        vector<int> nucj_list = Base_table.at(rna_seq[j]); //j list
        int j_1 = j - 1;
        for (auto &nucj : nucj_list){// all possible nuc at position j 
            
            // NONE->S
            // NONE->N
            index = GetIndex(j, nucj, nucj);
            newscore = 0;

            if(isLastNuc(j) && (1-lambda))
                newscore += getCAI_DN(ami_seq[j/3], nucj);
            
            update(bestS_j[1], index, newscore, -1, MANNER_NONEtoS);// initialize bestS[j]
            update(bestN_j, index, newscore, -1, MANNER_NONEtoN);// initialize bestN[j] 

            /*all N -> another state*/
            std::unordered_map<int, State<T>>& bestN_j_1 = bestN[j_1];
            for (auto &itemN_j_1 : bestN_j_1) {// all different states -> N at position j-1
                int index_j_1 = itemN_j_1.first;// key
                State<T> stateN_j_1 = itemN_j_1.second;// state
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);// from which position i and nuc of i and j-1
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    //because N -> N means Non-pairing to Non-pairing, only need to consider CAI(at position j)
                    newscore = stateN_j_1.score;
                    
                    if(isLastNuc(j) && (1-lambda))
                        newscore += getCAI_DN(ami_seq[j/3], nucj);
                    // N -> N
                    index = GetIndex(i, nuci, nucj);
                    update(bestN_j, index, newscore, index_j_1, MANNER_N_EtoN);

                    int i_1 = i - 1;//the continuous Non-pairing is from i to j-1, so if we want to transfer N -> C ,the pairing must be i-1 and j (loop : i ~ j-1)
                    if((i_1 >= 0) && (j - i_1 > HAIRPIN_GAP)){//consider hairpin(gap > 3)
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                            if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                                if(_allowed_pairs[nuci_1][nucj]){
                                    //because N -> C means Non-pairing to Closing, need to consider hairpin structure
                                    newscore = - lambda*v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)) + stateN_j_1.score;
                                    if(isLastNuc(j) && (1-lambda))
                                        newscore += getCAI_DN(ami_seq[j/3], nucj);
                                    if(isLastNuc(i_1) && (1-lambda))
                                        newscore += getCAI_DN(ami_seq[i_1/3], nuci_1);
                                    //N->C
                                    index = GetIndex(i_1, nuci_1, nucj);
                                    update(bestC_j, index, newscore, index_j_1, MANNER_NtoC);
                            }}
                }}}
            }
        }
        BeamPrune(con_seq, rna_seq, bestN_j, bestF, !(1-lambda));
         /*all S -> another state*/
        for (auto &nucj : nucj_list){
            int counter = 0;
            for(int l = 1; l <= min(SINGLE_MAX_LEN-1,j) ; l++){
                std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1][l];
                for (auto &itemS_j_1 : bestS_j_1) {
                    int index_j_1 = itemS_j_1.first;
                    State<T> stateS_j_1 = itemS_j_1.second;
                    int i, nucj_1, nuci;
                    std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                    if(IsLegal(con_seq, nucj_1, nucj, j_1) && (j - i) < SINGLE_MAX_LEN){
                        //S means single length extend
                        newscore = stateS_j_1.score;
                        if(isLastNuc(j) && (1-lambda))
                            newscore += getCAI_DN(ami_seq[j/3], nucj);
                        // S -> S
                        index = GetIndex(i, nuci, nucj);
                        update(bestS_j[l+1], index, newscore, index_j_1, MANNER_S_EtoS);
                    }
                }
            }
        }
        
        // for(int l = 1; l <= (SINGLE_MAX_LEN,j) ; l++)
        //     BeamPrune(con_seq, rna_seq, bestS_j[l], bestF, !(1-lambda));


        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestCS_j_1 = bestCS[j_1];
            for (auto &itemCS_j_1 : bestCS_j_1) {
                int index_j_1 = itemCS_j_1.first;
                State<T> stateCS_j_1 = itemCS_j_1.second;
                int i1, nucj_1, nuci1, nuci1_pair, nucx,len;
                std::tie(i1, nuci1, nuci1_pair, nucj_1, nucx, len) = GetIndexTupleCS(index_j_1);//make_tuple(i, nuci, nucj, nucx, len)
                int i = i1 - 1;
                if(i >= 0 && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i list
                    for(auto &nuci : nuci_list){
                        if(IsLegal(con_seq, nuci, nuci1, i) && _allowed_pairs[nuci][nucj]){
                            //CS means the closing pair (nucleotide p,nucleotide q) plus all POSSIBLE single extend.
                            int pre_index_2 = stateCS_j_1.index_2;
                            int pre_index_1 = stateCS_j_1.index_1;
                            int s, nucs, nucl;
                            std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                            int s_1 = s - 1;
                            int a, nuca, nucs_1;
                            std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);
                            //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                            newscore = - lambda*v_score_single(i, j, i1, s_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                                EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score;
                            
                            if(isLastNuc(j) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[i/3], nuci);
                            
                            //CS->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CStoC);//left bulge
                        }
                    }
                    //cout<<"solve state S + CS to C"<<endl;
                    
                    for(int l = 1; l <= min(SINGLE_MAX_LEN-len,i) ; l++){
                        std::unordered_map<int, State<T>>& bestS_i = bestS[i][l];
                        for (auto &itemS_i : bestS_i) {
                            int index_i = itemS_i.first;
                            State<T> stateS_i = itemS_i.second;
                            int k1, nuck1, nuci;
                            std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                            int k = k1 - 1;
                            
                            if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i) && i - k + len <= SINGLE_MAX_LEN ){
                                vector<int> nuck_list = Base_table.at(rna_seq[k]); //i list
                                for(auto &nuck : nuck_list){
                                    if(IsLegal(con_seq, nuck, nuck1, k) && _allowed_pairs[nuck][nucj]){
                                        int pre_index_2 = stateCS_j_1.index_2;
                                        int pre_index_1 = stateCS_j_1.index_1;
                                        int s, nucs, nucl;
                                        std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                                        int s_1 = s - 1;
                                        int a, nuca, nucs_1;
                                        std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);
                                        //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                                        newscore = - lambda*v_score_single(k, j, i1, s_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score + stateS_i.score;
                                        
                                        if(isLastNuc(j) && (1-lambda))
                                            newscore += getCAI_DN(ami_seq[j/3], nucj);
                                        if(isLastNuc(k) && (1-lambda))
                                            newscore += getCAI_DN(ami_seq[k/3], nuck);
                                        
                                        //S+CS->C
                                        index = GetIndex(k, nuck, nucj);
                                        update(bestC_j, index, newscore, index_i ,index_j_1, MANNER_S_CStoC);
                                    }
                                }
                            }
                            
                            
                        }
                    }
                    
                }
            }
        }
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        for(int l = min(SINGLE_MAX_LEN,j-4);l >= 1; --l){
            for (auto &itemS_j : bestS_j[l]) {
                int index_j = itemS_j.first;
                State<T> stateS_j = itemS_j.second;
                int i, nucj, nuci;
                std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
                int i_1 = i - 1;
                if(i_1 >= 0){
                    std::unordered_map<int, State<T>>& bestC_i_1 = bestC[i_1];
                    
                    for (auto &itemC_i_1 : bestC_i_1) {
                        int index_i_1 = itemC_i_1.first;
                        State<T> stateC_i_1 = itemC_i_1.second;
                        int k, nuck, nuci_1;    
                        std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);

                        if(k >= 0 && IsLegal(con_seq, nuci_1, nuci, i_1)){
                            newscore = stateS_j.score + stateC_i_1.score;
                            // C+S->CS
                            //GetIndexCS(int i, int nuci, int nuci_pair,int nucj, int nucx ,int len)
                            index = GetIndexCS(k, nuck, nuci_1,nucj, nuci,l);
                            update(bestCS_j, index, newscore, index_i_1, index_j, MANNER_C_StoCS);
                        }
                    }
                }
        
            }
        }

        BeamPrune(con_seq, rna_seq, bestCS_j, bestF, false);
        
        for (auto &nucj : nucj_list){

            std::unordered_map<int, State<T>>& bestC_j_1 = bestC[j_1];
            for (auto &itemC_j_1 : bestC_j_1) {
                int index_j_1 = itemC_j_1.first;
                State<T> stateC_j_1 = itemC_j_1.second;
                int i1, nucj_1, nuci1;
                std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);
                int i = i1 - 1;

                
                
                if(IsLegal(con_seq, nucj_1, nucj, j_1) && (i >= 0) && (j - i > (HAIRPIN_GAP + 2))){
                    
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i-1 list
                    for(auto &nuci : nuci_list){
                        
                        if(_allowed_pairs[nuci][nucj] && IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j_1)){
                            
                            newscore = - lambda*v_score_single(i, j, i1, j_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score;

                            if(isLastNuc(j) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[i/3], nuci);
                            
                            //C->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CtoC);
                            
                            
                    }}

                    int i_1 = i - 1;
                    
                    if(i_1 >= 0){
                        for(int l = 1;l<=min(SINGLE_MAX_LEN,i);l++){
                            std::unordered_map<int, State<T>>& bestS_i = bestS[i][l];
                            for (auto &itemS_i : bestS_i) {
                                int index_i = itemS_i.first;
                                State<T> stateS_i = itemS_i.second;
                                int k1, nuck1, nuci;
                                std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                                int k = k1 - 1;
                                if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i)){
                                    vector<int> nuck_list = Base_table.at(rna_seq[k]); //k list
                                    for(auto &nuck : nuck_list){
                                        
                                        if(_allowed_pairs[nuck][nucj] && IsLegal(con_seq, nuck, nuck1, k)){
                                            
                                            newscore = - lambda*v_score_single(k, j, i1, j_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score + stateS_i.score;

                                            if(isLastNuc(j) && (1-lambda))
                                                newscore += getCAI_DN(ami_seq[j/3], nucj);
                                            if(isLastNuc(k) && (1-lambda))
                                                newscore += getCAI_DN(ami_seq[k/3], nuck);
                                            
                                            //S+C->C
                                            index = GetIndex(k, nuck, nucj);
                                            update(bestC_j, index, newscore, index_i, index_j_1, MANNER_S_CtoC);

                            }}}}
                        }
                        
                    }
                    
            }} //Cj_1
            
        } //nucj
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestMulti_j_1 = bestMulti[j_1];
            for (auto &itemMulti_j_1 : bestMulti_j_1) { 
                int index_j_1 = itemMulti_j_1.first;
                State<T> stateMulti_j_1 = itemMulti_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    
                    int last_pair_pos_Multi = stateMulti_j_1.last_pair_pos;
                    if(j - last_pair_pos_Multi <= SINGLE_MAX_LEN){
                        newscore = stateMulti_j_1.score;

                        if(isLastNuc(j) && (1-lambda))
                            newscore += getCAI_DN(ami_seq[j/3], nucj);

                        //Multi->Multi
                        index = GetIndex(i, nuci, nucj);
                        //newscore = INT32_MIN;//test
                        update(bestMulti_j, index, newscore, index_j_1, MANNER_Multi_EtoMulti, last_pair_pos_Multi);
                    }

                    int i_1 = i - 1;
                    if(i_1 >= 0){
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                        if(IsLegal(con_seq, nuci_1, nuci, i_1) && _allowed_pairs[nuci_1][nucj]){
                            newscore = - lambda*v_score_multi(EnrBASE(nuci_1), EnrBASE(nucj)) + stateMulti_j_1.score;
                            
                            if(isLastNuc(j) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[j/3], nucj);
                            if(isLastNuc(i_1) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[i_1/3], nuci_1);

                            //Multi->C
                            index = GetIndex(i_1, nuci_1, nucj);
                            //newscore = INT32_MIN;//test
                            update(bestC_j, index, newscore, index_j_1, MANNER_MultitoC);
                        }}
                    }
                }
            }
        } //nucj
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);

        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = - lambda*v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score;

            //C->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_CtoM1, j);

            int i_1 = i - 1;
            if(i_1 >= HAIRPIN_GAP + 1){
                std::unordered_map<int, State<T>>& bestM1_i_1 = bestM1[i_1];
                for (auto &itemM1_i_1 : bestM1_i_1) {
                    int index_i_1 = itemM1_i_1.first;
                    State<T> stateM1_i_1 = itemM1_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = - lambda*v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score + stateM1_i_1.score;
                        // M1+C->M2
                        index = GetIndex(k, nuck, nucj);
                        update(bestM2_j, index, newscore, index_i_1 ,index_j, MANNER_M1_CtoM2, j);
                    }
                }
            }
        }
        
        BeamPrune(con_seq, rna_seq, bestM2_j, bestF, false);

        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestM1_j_1 = bestM1[j_1];
            for (auto &itemM1_j_1 : bestM1_j_1) {
                int index_j_1 = itemM1_j_1.first;
                State<T> stateM1_j_1 = itemM1_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);

                int last_pair_pos_M1 = stateM1_j_1.last_pair_pos;
                assert(last_pair_pos_M1 != -1);
                if(j - last_pair_pos_M1 <= SINGLE_MAX_LEN && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateM1_j_1.score;

                    if(isLastNuc(j) && (1-lambda))
                        newscore += getCAI_DN(ami_seq[j/3], nucj);

                    //M1->M1
                    index = GetIndex(i, nuci, nucj);
                    update(bestM1_j, index, newscore, index_j_1, MANNER_M1_EtoM1, last_pair_pos_M1);
                }
            }
        } //nucj

        for (auto &itemM2_j : bestM2_j) {
            int index_j = itemM2_j.first;
            State<T> stateM2_j = itemM2_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = stateM2_j.score;
            int last_pair_pos_M2 = stateM2_j.last_pair_pos;

            //M2->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_M2toM1, last_pair_pos_M2);
            //M2->Multi
            update(bestMulti_j, index_j, newscore, index_j, MANNER_M2toMulti, last_pair_pos_M2);
            int i_1 = i - 1;
            if(i_1 >= 0){
                for(int l = 1; l <= min(SINGLE_MAX_LEN,i_1); l++){
                    std::unordered_map<int, State<T>>& bestS_i_1 = bestS[i_1][l];
                    for (auto &itemS_i_1 : bestS_i_1) {
                        int index_i_1 = itemS_i_1.first;
                        State<T> stateS_i_1 = itemS_i_1.second;
                        int k, nuck, nuci_1;
                        std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                        if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                            newscore = stateM2_j.score + stateS_i_1.score;
                            //S+M2->Multi
                            index = GetIndex(k, nuck, nucj);
                            update(bestMulti_j, index, newscore, index_i_1, index_j, MANNER_S_M2toMulti, last_pair_pos_M2);
                        }
                    }
                }
                
            }
            
        }

        BeamPrune(con_seq, rna_seq, bestM1_j, bestF, false);
        BeamPrune(con_seq, rna_seq, bestMulti_j, bestF, false);

        {//conclusion!
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestF_j_1 = bestF[j_1];
            for (auto &itemF_j_1 : bestF_j_1) {
                int index_j_1 = itemF_j_1.first;
                State<T> stateF_j_1 = itemF_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateF_j_1.score;
                    if(isLastNuc(j) && (1-lambda))
                        newscore += getCAI_DN(ami_seq[j/3], nucj);
                    //F->F
                    index = GetIndex(i, nuci, nucj);
                    update(bestF_j, index, newscore, index_j_1, MANNER_F_EtoF);
                }
            }
        }
        
        
        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);

            if(i == 0 ){
                newscore = stateC_j.score - lambda*v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                //C->F
                update(bestF_j, index_j, newscore, index_j, MANNER_CtoF);
            }

            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestF_i_1 = bestF[i_1];
                for (auto &itemF_i_1 : bestF_i_1) {
                    int index_i_1 = itemF_i_1.first;
                    State<T> stateF_i_1 = itemF_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if((k == 0 || j == seq_length-1) && IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateF_i_1.score + stateC_j.score - lambda*v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                        //F+C->F
                        index = GetIndex(k, nuck, nucj);
                        update(bestF_j, index, newscore, index_i_1, index_j, MANNER_F_CtoF);
            }}}
        }
        
        //BeamPrune(con_seq, rna_seq, bestF_j, bestF, false);
        }

    } //jend
}

template<typename T>
void LCDSfoldCAI_DN_exact(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq){
    //DERNA
    int seq_length = rna_seq.size();
    T newscore;
    int index;
    int preindex;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    vector<int> nuc0_list = Base_table.at(rna_seq[0]); //j list
    for(int c0 = 0; c0 < nuc0_list.size(); ++c0){
        int nuc0 = nuc0_list[c0];
        index = GetIndex(0, nuc0, nuc0);
        update(bestN[0], index, 0, -1, MANNER_NONEtoN);
        update(bestF[0], index, 0, -1, MANNER_NONEtoF);
        update(bestS[0][1], index, 0, -1, MANNER_NONEtoS);
    }
    
    for(int j = 1; j < seq_length; ++j){// 1 ~ seq_length
        Processing(seq_length);
        std::unordered_map<int, State<T>>& bestN_j = bestN[j];
        std::vector<std::unordered_map<int, State<T>>>& bestS_j = bestS[j];
        std::unordered_map<int, State<T>>& bestF_j = bestF[j];
        std::unordered_map<int, State<T>>& bestC_j = bestC[j];
        std::unordered_map<int, State<T>>& bestM1_j = bestM1[j];
        std::unordered_map<int, State<T>>& bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>>& bestMulti_j = bestMulti[j];
        
        vector<int> nucj_list = Base_table.at(rna_seq[j]); //j list
        int j_1 = j - 1;
        for (auto &nucj : nucj_list){// all possible nuc at position j 
            
            // NONE->S
            // NONE->N
            index = GetIndex(j, nucj, nucj);
            newscore = 0;

            if(isLastNuc(j) && (1-lambda))
                newscore += getCAI_DN(ami_seq[j/3], nucj);
            
            update(bestS_j[1], index, newscore, -1, MANNER_NONEtoS);// initialize bestS[j]
            update(bestN_j, index, newscore, -1, MANNER_NONEtoN);// initialize bestN[j] 

            /*all N -> another state*/
            std::unordered_map<int, State<T>>& bestN_j_1 = bestN[j_1];
            for (auto &itemN_j_1 : bestN_j_1) {// all different states -> N at position j-1
                int index_j_1 = itemN_j_1.first;// key
                State<T> stateN_j_1 = itemN_j_1.second;// state
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);// from which position i and nuc of i and j-1
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    //because N -> N means Non-pairing to Non-pairing, only need to consider CAI(at position j)
                    newscore = stateN_j_1.score;
                    
                    if(isLastNuc(j) && (1-lambda))
                        newscore += getCAI_DN(ami_seq[j/3], nucj);
                    // N -> N
                    index = GetIndex(i, nuci, nucj);
                    update(bestN_j, index, newscore, index_j_1, MANNER_N_EtoN);

                    int i_1 = i - 1;//the continuous Non-pairing is from i to j-1, so if we want to transfer N -> C ,the pairing must be i-1 and j (loop : i ~ j-1)
                    if((i_1 >= 0) && (j - i_1 > HAIRPIN_GAP)){//consider hairpin(gap > 3)
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                            if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                                if(_allowed_pairs[nuci_1][nucj]){
                                    //because N -> C means Non-pairing to Closing, need to consider hairpin structure
                                    newscore = - lambda*v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)) + stateN_j_1.score;
                                    if(isLastNuc(j) && (1-lambda))
                                        newscore += getCAI_DN(ami_seq[j/3], nucj);
                                    if(isLastNuc(i_1) && (1-lambda))
                                        newscore += getCAI_DN(ami_seq[i_1/3], nuci_1);
                                    //N->C
                                    index = GetIndex(i_1, nuci_1, nucj);
                                    update(bestC_j, index, newscore, index_j_1, MANNER_NtoC);
                            }}
                }}}
            }
        }

         /*all S -> another state*/
        for (auto &nucj : nucj_list){
            int counter = 0;
            for(int l = 1; l <= min(SINGLE_MAX_LEN-1,j) ; l++){
                std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1][l];
                for (auto &itemS_j_1 : bestS_j_1) {
                    int index_j_1 = itemS_j_1.first;
                    State<T> stateS_j_1 = itemS_j_1.second;
                    int i, nucj_1, nuci;
                    std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                    if(IsLegal(con_seq, nucj_1, nucj, j_1) && (j - i) < SINGLE_MAX_LEN){
                        //S means single length extend
                        newscore = stateS_j_1.score;
                        if(isLastNuc(j) && (1-lambda))
                            newscore += getCAI_DN(ami_seq[j/3], nucj);
                        // S -> S
                        index = GetIndex(i, nuci, nucj);
                        update(bestS_j[l+1], index, newscore, index_j_1, MANNER_S_EtoS);
                    }
                }
            }
        }

        for (auto &nucj : nucj_list){
            for(int l1 = min(SINGLE_MAX_LEN,j-4) ; l1 >= 1; --l1){
                std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1][l1];
                for (auto &itemS_j_1 : bestS_j_1) {
                    int index_j_1 = itemS_j_1.first;
                    State<T> stateS_j_1 = itemS_j_1.second;
                    int q1, nucq1, nucj_1;
                    std::tie(q1, nucq1, nucj_1) = GetIndexTuple(index_j_1);
                    int q = q1 - 1;
                    
                    if(q >= 0 &&  IsLegal(con_seq, nucj_1, nucj, j_1)){
                        std::unordered_map<int, State<T>>& bestC_q = bestC[q];    
                        for (auto &itemC_q : bestC_q) {
                            int index_q = itemC_q.first;
                            State<T> stateC_q = itemC_q.second;
                            int p, nucp, nucq;    
                            std::tie(p, nucp, nucq) = GetIndexTuple(index_q);
                            int p_1 = p - 1;

                            if(p_1 >= 0 && IsLegal(con_seq, nucq, nucq1, q)){
                                //C_StoS (right bulge)
                                int i = p_1;
                                vector<int> nuci_list = Base_table.at(rna_seq[i]); //i list
                                for(auto &nuci : nuci_list){
                                    if(IsLegal(con_seq, nuci, nucp, i) && _allowed_pairs[nuci][nucj]){
                                        //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                                        newscore = -lambda*v_score_single(i,j,p,q,EnrBASE(nuci),EnrBASE(nucp),EnrBASE(nucj_1),EnrBASE(nucj),EnrBASE(nuci),EnrBASE(nucp),EnrBASE(nucq),EnrBASE(nucq1))
                                        + stateS_j_1.score + stateC_q.score;  

                                        if(isLastNuc(j) && (1-lambda))
                                            newscore += getCAI_DN(ami_seq[j/3], nucj);
                                        if(isLastNuc(i) && (1-lambda))
                                            newscore += getCAI_DN(ami_seq[i/3], nuci);

                                        index = GetIndex(i, nuci, nucj);
                                        update(bestC_j, index, newscore, index_q, index_j_1, MANNER_C_StoC);
                                    }
                                }//end nuci list


                                //S_C_StoC (internal loop)
                                for(int l2 = 1; l2 <= min(SINGLE_MAX_LEN-l1 , p_1) ; ++l2){
                                    std::unordered_map<int, State<T>>& bestS_p_1 = bestS[p_1][l2];
                                    for (auto &itemS_p_1 : bestS_p_1) {
                                        int index_p_1 = itemS_p_1.first;
                                        State<T> stateS_p_1 = itemS_p_1.second;
                                        int i1, nuci1, nucp_1;
                                        std::tie(i1, nuci1, nucp_1) = GetIndexTuple(index_p_1);
                                        int i = i1 - 1;
                                        if(i >= 0 && IsLegal(con_seq, nucp_1, nucp, p_1)){
                                            vector<int> nuci_list = Base_table.at(rna_seq[i]); //i list
                                            for(auto &nuci : nuci_list){
                                                if(IsLegal(con_seq, nuci, nuci1, i) && _allowed_pairs[nuci][nucj]){
                                                    //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                                                    newscore = -lambda*v_score_single(i,j,p,q,EnrBASE(nuci),EnrBASE(nuci1),EnrBASE(nucj_1),EnrBASE(nucj),EnrBASE(nucp_1),EnrBASE(nucp),EnrBASE(nucq),EnrBASE(nucq1))
                                                    + stateS_j_1.score + stateC_q.score + stateS_p_1.score;  

                                                    if(isLastNuc(j) && (1-lambda))
                                                        newscore += getCAI_DN(ami_seq[j/3], nucj);
                                                    if(isLastNuc(i) && (1-lambda))
                                                        newscore += getCAI_DN(ami_seq[i/3], nuci);

                                                    index = GetIndex(i, nuci, nucj);
                                                    update(bestC_j, index, newscore, index_p_1, index_q, index_j_1, MANNER_S_C_StoC);
                                                }
                                            }//end nuci list
                                        
                                        }
                                        
                                    }//end bestS[p_1][l2]
                                }
                            }
                        }//end bestC_q
                    }
            
                }//end bestS[j_1]
            }
        }//end of nucj
        
        
        for (auto &nucj : nucj_list){

            std::unordered_map<int, State<T>>& bestC_j_1 = bestC[j_1];
            for (auto &itemC_j_1 : bestC_j_1) {
                int index_j_1 = itemC_j_1.first;
                State<T> stateC_j_1 = itemC_j_1.second;
                int i1, nucj_1, nuci1;
                std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);
                int i = i1 - 1;

                
                
                if(IsLegal(con_seq, nucj_1, nucj, j_1) && (i >= 0) && (j - i > (HAIRPIN_GAP + 2))){
                    
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i-1 list
                    for(auto &nuci : nuci_list){
                        
                        if(_allowed_pairs[nuci][nucj] && IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j_1)){
                            
                            newscore = - lambda*v_score_single(i, j, i1, j_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score;

                            if(isLastNuc(j) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[i/3], nuci);
                            
                            //C->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CtoC);
                            
                            
                    }}

                    int i_1 = i - 1;
                    
                    if(i_1 >= 0){
                        for(int l = 1;l<=min(SINGLE_MAX_LEN,i);l++){
                            std::unordered_map<int, State<T>>& bestS_i = bestS[i][l];
                            for (auto &itemS_i : bestS_i) {
                                int index_i = itemS_i.first;
                                State<T> stateS_i = itemS_i.second;
                                int k1, nuck1, nuci;
                                std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                                int k = k1 - 1;
                                if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i)){
                                    vector<int> nuck_list = Base_table.at(rna_seq[k]); //k list
                                    for(auto &nuck : nuck_list){
                                        
                                        if(_allowed_pairs[nuck][nucj] && IsLegal(con_seq, nuck, nuck1, k)){
                                            
                                            newscore = - lambda*v_score_single(k, j, i1, j_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score + stateS_i.score;

                                            if(isLastNuc(j) && (1-lambda))
                                                newscore += getCAI_DN(ami_seq[j/3], nucj);
                                            if(isLastNuc(k) && (1-lambda))
                                                newscore += getCAI_DN(ami_seq[k/3], nuck);
                                            
                                            //S+C->C
                                            index = GetIndex(k, nuck, nucj);
                                            update(bestC_j, index, newscore, index_i, index_j_1, MANNER_S_CtoC);

                            }}}}
                        }
                        
                    }
                    
            }} //Cj_1
            
        } //nucj
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestMulti_j_1 = bestMulti[j_1];
            for (auto &itemMulti_j_1 : bestMulti_j_1) { 
                int index_j_1 = itemMulti_j_1.first;
                State<T> stateMulti_j_1 = itemMulti_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    
                    int last_pair_pos_Multi = stateMulti_j_1.last_pair_pos;
                    if(j - last_pair_pos_Multi <= SINGLE_MAX_LEN){
                        newscore = stateMulti_j_1.score;

                        if(isLastNuc(j) && (1-lambda))
                            newscore += getCAI_DN(ami_seq[j/3], nucj);

                        //Multi->Multi
                        index = GetIndex(i, nuci, nucj);
                        //newscore = INT32_MIN;//test
                        update(bestMulti_j, index, newscore, index_j_1, MANNER_Multi_EtoMulti, last_pair_pos_Multi);
                    }

                    int i_1 = i - 1;
                    if(i_1 >= 0){
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                        if(IsLegal(con_seq, nuci_1, nuci, i_1) && _allowed_pairs[nuci_1][nucj]){
                            newscore = - lambda*v_score_multi(EnrBASE(nuci_1), EnrBASE(nucj)) + stateMulti_j_1.score;
                            
                            if(isLastNuc(j) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[j/3], nucj);
                            if(isLastNuc(i_1) && (1-lambda))
                                newscore += getCAI_DN(ami_seq[i_1/3], nuci_1);

                            //Multi->C
                            index = GetIndex(i_1, nuci_1, nucj);
                            //newscore = INT32_MIN;//test
                            update(bestC_j, index, newscore, index_j_1, MANNER_MultitoC);
                        }}
                    }
                }
            }
        } //nucj

        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = - lambda*v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score;

            //C->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_CtoM1, j);

            int i_1 = i - 1;
            if(i_1 >= HAIRPIN_GAP + 1){
                std::unordered_map<int, State<T>>& bestM1_i_1 = bestM1[i_1];
                for (auto &itemM1_i_1 : bestM1_i_1) {
                    int index_i_1 = itemM1_i_1.first;
                    State<T> stateM1_i_1 = itemM1_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = - lambda*v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score + stateM1_i_1.score;
                        // M1+C->M2
                        index = GetIndex(k, nuck, nucj);
                        update(bestM2_j, index, newscore, index_i_1 ,index_j, MANNER_M1_CtoM2, j);
                    }
                }
            }
        }
        
        BeamPrune(con_seq, rna_seq, bestM2_j, bestF, false);

        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestM1_j_1 = bestM1[j_1];
            for (auto &itemM1_j_1 : bestM1_j_1) {
                int index_j_1 = itemM1_j_1.first;
                State<T> stateM1_j_1 = itemM1_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);

                int last_pair_pos_M1 = stateM1_j_1.last_pair_pos;
                assert(last_pair_pos_M1 != -1);
                if(j - last_pair_pos_M1 <= SINGLE_MAX_LEN && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateM1_j_1.score;

                    if(isLastNuc(j) && (1-lambda))
                        newscore += getCAI_DN(ami_seq[j/3], nucj);

                    //M1->M1
                    index = GetIndex(i, nuci, nucj);
                    update(bestM1_j, index, newscore, index_j_1, MANNER_M1_EtoM1, last_pair_pos_M1);
                }
            }
        } //nucj

        for (auto &itemM2_j : bestM2_j) {
            int index_j = itemM2_j.first;
            State<T> stateM2_j = itemM2_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = stateM2_j.score;
            int last_pair_pos_M2 = stateM2_j.last_pair_pos;

            //M2->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_M2toM1, last_pair_pos_M2);
            //M2->Multi
            update(bestMulti_j, index_j, newscore, index_j, MANNER_M2toMulti, last_pair_pos_M2);
            int i_1 = i - 1;
            if(i_1 >= 0){
                for(int l = 1; l <= min(SINGLE_MAX_LEN,i_1); l++){
                    std::unordered_map<int, State<T>>& bestS_i_1 = bestS[i_1][l];
                    for (auto &itemS_i_1 : bestS_i_1) {
                        int index_i_1 = itemS_i_1.first;
                        State<T> stateS_i_1 = itemS_i_1.second;
                        int k, nuck, nuci_1;
                        std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                        if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                            newscore = stateM2_j.score + stateS_i_1.score;
                            //S+M2->Multi
                            index = GetIndex(k, nuck, nucj);
                            update(bestMulti_j, index, newscore, index_i_1, index_j, MANNER_S_M2toMulti, last_pair_pos_M2);
                        }
                    }
                }
            }
            
        }

        {//conclusion!
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestF_j_1 = bestF[j_1];
            for (auto &itemF_j_1 : bestF_j_1) {
                int index_j_1 = itemF_j_1.first;
                State<T> stateF_j_1 = itemF_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateF_j_1.score;
                    if(isLastNuc(j) && (1-lambda))
                        newscore += getCAI_DN(ami_seq[j/3], nucj);
                    //F->F
                    index = GetIndex(i, nuci, nucj);
                    update(bestF_j, index, newscore, index_j_1, MANNER_F_EtoF);
                }
            }
        }
        
        
        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);

            if(i == 0 ){
                newscore = stateC_j.score - lambda*v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                //C->F
                update(bestF_j, index_j, newscore, index_j, MANNER_CtoF);
            }

            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestF_i_1 = bestF[i_1];
                for (auto &itemF_i_1 : bestF_i_1) {
                    int index_i_1 = itemF_i_1.first;
                    State<T> stateF_i_1 = itemF_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if((k == 0 || j == seq_length-1) && IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateF_i_1.score + stateC_j.score - lambda*v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                        //F+C->F
                        index = GetIndex(k, nuck, nucj);
                        update(bestF_j, index, newscore, index_i_1, index_j, MANNER_F_CtoF);
            }}}
        }
        
        }

    } //jend
}

template<typename T>
std::vector<pair<T,T>> result_output(AllTables<T> &alltables,string& rna_seq, vector<int>& con_seq, string& ami_seq, string& output_txt, string& output_csv, double TimeSpend, bool show_score,bool is_DN){//return {CAI,MFE}
    std::ofstream outputfile(output_txt,std::ios::app);
    std::ofstream outputcsv(output_csv,std::ios::app);

    std::string rna_solution;
    std::string structure_solution;

    outputfile << "Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;
    T maxscore = BackTrack<T>(rna_solution, structure_solution, rna_seq, alltables, con_seq);
    T weighted_cai_score = GetCAIScore<T>(rna_solution,is_DN);
    // cout<<"weighted_CAI: "<<weighted_cai_score<<endl;
    double cai_value = GetUnweghtedCAIScore(rna_solution);

    cai_value = exp(cai_value/double(ami_seq.size()));
    //output score and results
    std::cout << "Coding sequence and its secondary structure:" << std::endl;
    outputfile << "Coding sequence and its secondary structure:" << std::endl;
    std::cout << rna_solution << std::endl;
    outputfile << rna_solution << std::endl;
    std::cout << structure_solution << std::endl;
    outputfile << structure_solution << std::endl;
    
    if (show_score){
        std::cerr << "Score: " << maxscore << std::endl;
        outputfile << "Score: " << maxscore << std::endl;
    }
           
    double mfe_value ;
    if(is_DN){
        mfe_value = -double(double(maxscore - weighted_cai_score)/(100.0*lambda));//*100
        std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
        outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    }
    else{
        if(lambda){
            mfe_value = -double(double(maxscore - weighted_cai_score)/(100.0));//*100
            std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
            outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;  
        }
        else{//only consider MFE
            mfe_value = -(maxscore/(100.0));//*100
            std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
            outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;  
        }
    }
    std::cout << "CAI: " << cai_value << std::endl;
    outputfile << "CAI: " <<  cai_value << std::endl;

    std::cout << "Total runtime: " << TimeSpend << " s" << std::endl;
    outputfile <<  "Total runtime: " << TimeSpend << " s" << std::endl;

    // std::cout << "======================================================================================================" << std::endl;
    // outputfile << "======================================================================================================" << std::endl;
    check_ami_solution(ami_seq, rna_solution);

    std::vector<pair<T,T>> result_container;
    // result_container.push_back({round_up(cai_value),round_up(mfe_value)});
    result_container.push_back({round_up(cai_value),round_up(mfe_value)});
    outputcsv<<std::to_string(lambda) +","+ std::to_string(mfe_value) +","+std::to_string(cai_value)<<std::endl;
    return result_container;

}

template<typename T>
void Pareto_solution(double threshold1, double threshold2, string& rna_seq, vector<int>& con_seq, string& ami_seq, vector<double>& cai_vector, string &output_txt, string &output_csv, bool show_score){
    double left_lambda,right_lambda;
    double left_cai,left_mfe,right_cai,right_mfe;
    timeval start,end;
    std::queue<pair<double,double>> lambda_vector;
    std::unordered_map<double,double> MFE_map,CAI_map;
    lambda_vector.push({0.00001,0.99999});
    while(true){
        std::cout<<"Remain: "<<lambda_vector.size()<<" different lambda"<<std::endl;
        if(!lambda_vector.empty()){ 
            left_lambda = lambda_vector.front().first;
            right_lambda = lambda_vector.front().second;
            lambda_vector.pop();
        }
        else break;// No lambda need to try
        
        if(MFE_map.find(left_lambda) != MFE_map.end()){//left lambda has tried before
            left_cai = CAI_map[left_lambda];
            left_mfe = MFE_map[left_lambda];
        }
        else{
            lambda = left_lambda;
            std::cout <<"Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;
            AllTables<T> alltables(rna_seq, rna_seq.size());
            std::string rna_solution;
            std::string structure_solution;
            initialize_CAI_table(cai_vector,true);

            gettimeofday(&start,0);
            initialize_Special_HP_DN<double>(alltables,rna_seq, con_seq, ami_seq);
            LCDSfoldCAI_DN_exact<double>(alltables, rna_seq, con_seq, ami_seq);
            gettimeofday(&end,0);
            double TimeSpend = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);

            std::vector<pair<double,double>> result = result_output<double>(alltables, rna_seq, con_seq, ami_seq,output_txt,output_csv,TimeSpend, show_score,true);
            CAI_map[left_lambda] = result.front().first;
            MFE_map[left_lambda] = result.front().second;
            left_cai = CAI_map[left_lambda];
            left_mfe = MFE_map[left_lambda];
                
        }

        if(MFE_map.find(right_lambda) != MFE_map.end()){//right lambda has tried before
            right_cai = CAI_map[right_lambda];
            right_mfe = MFE_map[right_lambda];
        }
        else{
            lambda = right_lambda;
            std::cout <<"Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;
            AllTables<T> alltables(rna_seq, rna_seq.size());
            std::string rna_solution;
            std::string structure_solution;
            initialize_CAI_table(cai_vector,true);

            gettimeofday(&start,0);
            initialize_Special_HP_DN<double>(alltables,rna_seq, con_seq, ami_seq);
            LCDSfoldCAI_DN_exact<double>(alltables, rna_seq, con_seq, ami_seq);
            gettimeofday(&end,0);
            double TimeSpend = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);

            std::vector<pair<double,double>> result = result_output<double>(alltables, rna_seq, con_seq, ami_seq,output_txt,output_csv,TimeSpend, show_score, true);
            CAI_map[right_lambda] = result.front().first;
            MFE_map[right_lambda] = result.front().second;
            right_cai = CAI_map[right_lambda];
            right_mfe = MFE_map[right_lambda];

        }
        if (!compare(left_cai, right_cai) && !compare(left_mfe, right_mfe)) {
            if (right_lambda < threshold1) {//0.0025
                if (!compare(left_lambda,right_lambda,threshold2)) {//0.00075
                    double m = (left_lambda + right_lambda) / 2;
                    lambda_vector.push({left_lambda, m});
                    lambda_vector.push({m, right_lambda});
                }
            } else {
                if (!compare(left_lambda,right_lambda,threshold1)) {//0.0025
                    double m = (left_lambda + right_lambda) / 2;
                    lambda_vector.push({left_lambda, m});
                    lambda_vector.push({m, right_lambda});
                }
            }

        }
        
    }
    
    
}

