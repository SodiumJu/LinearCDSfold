#pragma once


template<typename T>
unsigned long QuickselectPartition(std::vector<std::pair<T, int>>& scores, unsigned long lower, unsigned long upper) {
    T pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);

    }
    return upper;
}

template<typename T>
T QuickSelect(std::vector<std::pair<T, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = QuickselectPartition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return QuickSelect(scores, lower, split-1, k);
    else return QuickSelect(scores, split+1, upper, k - length);
}

template<typename T>
T BeamPrune(std::vector<int>& con_seq, std::string& rna_seq, std::unordered_map<int, State<T>> &bestMap, std::vector<std::unordered_map<int, State<T>>>& bestF, bool isN) {
    std::vector<std::pair<T, int>> scores;
    scores.clear();
    bool have_previous;
    T newscore, pre_newscore, best_pre_newscore;
    int seq_length = rna_seq.size();

    if (bestMap.size() <= beamsize || beamsize == 0) return VALUE_MIN<T>();

    for (auto &item : bestMap) {

        int tail_index = item.first;
        State<T> cand = item.second;
        int i, nuci, nuci_pair, nucj, nucx, len;
        if(cand.MANNER == MANNER_C_StoCS)
            std::tie(i, nuci, nuci_pair, nucj, nucx, len) = GetIndexTupleCS(tail_index);
        else
            std::tie(i, nuci, nucj) = GetIndexTuple(tail_index);

        
        int i_1 = i - 1;

        if(isN){
            newscore = i;
        }else if(i_1 >= 0){

            std::vector<int> i_1_nuclist = Base_table.at(rna_seq[i_1]);
            std::vector<int> f_nuclist = Base_table.at(rna_seq[0]);
            pre_newscore = VALUE_MIN<T>();
            best_pre_newscore = VALUE_MIN<T>();
            have_previous=false;

            for(auto& nucf : f_nuclist){
                for(auto& nuci_1 : i_1_nuclist){

                    int nucindex = GetIndex(0, nucf, nuci_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1) && bestF[i_1].count(nucindex)){
                        
                        pre_newscore = bestF[i_1][nucindex].score;
                        have_previous = true;
                        if(pre_newscore > best_pre_newscore)
                            best_pre_newscore = pre_newscore;
                        
                    }
            }}
            
            if(have_previous)
                newscore = best_pre_newscore + cand.score;
            else
                newscore = cand.score;

        }else{newscore = cand.score;}

        scores.push_back(std::make_pair(newscore, tail_index));
        
    }
    
    T threshold = QuickSelect(scores, 0, scores.size() - 1, scores.size() - beamsize);
    for (auto &p : scores) {
        if (p.first < threshold) bestMap.erase(p.second);
    }

    return threshold;
}