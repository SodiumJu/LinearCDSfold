void BackTrack(int tail, int maxindex, const State& maxstate,
std::string& rna_solution, std::string& structure_solution, std::string& rna_seq, const  AllTables& alltables){
    
    const std::vector<std::unordered_map<int, State>>& bestN = alltables.bestN;
    const std::vector<std::unordered_map<int, State>>& bestS = alltables.bestS;
    const std::vector<std::unordered_map<int, State>>& bestF = alltables.bestF;
    const std::vector<std::unordered_map<int, State>>& bestC = alltables.bestC;
    const std::vector<std::unordered_map<int, State>>& bestCS = alltables.bestCS;
    const std::vector<std::unordered_map<int, State>>& bestM1 = alltables.bestM1;
    const std::vector<std::unordered_map<int, State>>& bestM2 = alltables.bestM2;
    const std::vector<std::unordered_map<int, State>>& bestMulti = alltables.bestMulti;


    int seq_length = rna_seq.size();
    rna_solution.assign(seq_length, '.');
    structure_solution.assign(seq_length, '.');
    std::stack<std::tuple<int, int, State>> stk;

    int max_score = maxstate.score;
    stk.push(std::make_tuple(tail, maxindex, maxstate));
    while ( !stk.empty() ) {

        std::tuple<int, int, State> top = stk.top();
        int j = std::get<0>(top);
        int cindex = std::get<1>(top);
        State state = std::get<2>(top);
        int index1 = state.index_1;
        int index2 = state.index_2;
        int i, nuci, nucj;
        std::tie(i, nuci, nucj) = GetIndexTuple(cindex);

        Manner MANNER = state.MANNER;
        stk.pop();
        int k, nuck, nucl;
        std::cout << "score: " << -(max_score - state.score)/100.0 << std::endl;
        max_score = state.score;
        std::cout << mannerToString(MANNER) << std::endl;
        std::cout << "(" << i << " " << j << ")"  << std::endl;
        std::cout << "(" << reBASE(nuci) << " " << reBASE(nucj) << ")" << std::endl;
        std::cout << "---" << std::endl;
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
                stk.push(std::make_tuple(j-1 ,index1, bestN[j-1].find(index1)->second));
                break;
            case MANNER_NtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                stk.push(std::make_tuple(j-1 ,index1, bestN[j-1].find(index1)->second));
                break;
            case MANNER_S_EtoS:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                stk.push(std::make_tuple(j-1 ,index1, bestS[j-1].find(index1)->second));
                break;
            case MANNER_CtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                stk.push(std::make_tuple(j-1 ,index1, bestC[j-1].find(index1)->second));
                break;
            case MANNER_S_CtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                stk.push(std::make_tuple(k-1 ,index1, bestS[k-1].find(index1)->second));
                stk.push(std::make_tuple(j-1 ,index2, bestC[j-1].find(index2)->second));
                break;
            case MANNER_CStoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                stk.push(std::make_tuple(j-1 ,index1, bestCS[j-1].find(index1)->second));
                break;
            case MANNER_S_CStoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                stk.push(std::make_tuple(j-1 ,index2, bestCS[j-1].find(index2)->second));
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                stk.push(std::make_tuple(k-1 ,index1, bestS[k-1].find(index1)->second));
                break;
            case MANNER_C_StoCS:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                stk.push(std::make_tuple(k-1 ,index1, bestC[k-1].find(index1)->second));
                stk.push(std::make_tuple(j ,index2, bestS[j].find(index2)->second));
                break;
            case MANNER_CtoM1:
                stk.push(std::make_tuple(j ,index1, bestC[j].find(index1)->second));
                break;
            case MANNER_M1_CtoM2:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                stk.push(std::make_tuple(k-1 ,index1, bestM1[k-1].find(index1)->second));
                stk.push(std::make_tuple(j ,index2, bestC[j].find(index2)->second));
                break;
            case MANNER_M1_EtoM1:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                stk.push(std::make_tuple(j-1 ,index1, bestM1[j-1].find(index1)->second));
                break;
            case MANNER_M2toM1:
                stk.push(std::make_tuple(j ,index1, bestM2[j].find(index1)->second));
                break;
            case MANNER_M2toMulti:
                stk.push(std::make_tuple(j ,index1, bestM2[j].find(index1)->second));
                break;
            case MANNER_S_M2toMulti:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                stk.push(std::make_tuple(k-1 ,index1, bestS[k-1].find(index1)->second));
                stk.push(std::make_tuple(j ,index2, bestM2[j].find(index2)->second));
                break;
            case MANNER_Multi_EtoMulti:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                stk.push(std::make_tuple(j-1 ,index1, bestMulti[j-1].find(index1)->second));
                break;
            case MANNER_MultitoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                stk.push(std::make_tuple(j-1 ,index1, bestMulti[j-1].find(index1)->second));
                break;
            case MANNER_F_EtoF:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                stk.push(std::make_tuple(j-1 ,index1, bestF[j-1].find(index1)->second));
                break;
            case MANNER_CtoF:
                stk.push(std::make_tuple(j ,index1, bestC[j].find(index1)->second));
                break;
            case MANNER_F_CtoF:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                stk.push(std::make_tuple(k-1 ,index1, bestF[k-1].find(index1)->second));
                stk.push(std::make_tuple(j ,index2, bestC[j].find(index2)->second));
                break;
            default:
                std::cout << "Something error" << std::endl;
        }
    }


}