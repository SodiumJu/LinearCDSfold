#include <sstream>
#include <fstream>
template <typename T>
constexpr T VALUE_MIN() { return std::numeric_limits<T>::lowest(); }

#define SINGLE_MAX_LEN 20
#define MIN_STACKPAIR_GAP 5
#define MIN_MULTI_GAP 10
#define NOTON 16 // NUM_OF_TYPE_OF_NUCS
#define HAIRPIN_GAP 3 // NUM_OF_SMALLEST_DIS_CAN_PAIR

//Amino acids label
#define normal_ami -1
#define Arg 0
#define Leu 1
#define Ser 2

#define base_number 16
#define square_base_number 256
#define cube_base_number 4096
#define fourthPower_base_number 65536


//global variables
// double LDCONST = 100;
double LDCONST = 100;
double lambda;
int beamsize;
bool pareto;


bool _allowed_pairs[NOTON][NOTON];
bool xcon_table[3][2][NOTON][NOTON];
// std::unordered_map<char, std::unordered_map<int,int>> CodonSetCAIMap;
std::unordered_map<char, std::unordered_map<int, double>> CodonSetCAIMap;
std::unordered_map<char, std::unordered_map<int, double>> CodonSetCAIMap_LD;
std::unordered_map<char, std::unordered_map<int, double>> CodonSetCAIMap_DN;
// std::unordered_map<std::string, int> CalCAIMap; //for get weighted CAI score from codon string
std::unordered_map<std::string, double> CalCAIMap;
std::unordered_map<std::string, double> CalCAIMap_DN;
std::unordered_map<std::string, double> unweighted_CAIMap;
std::unordered_map<int, std::string> sp_backtrack;
std::unordered_map<std::string, float> cai_table_csv;
// declaration of functions
int BASE(const std::string& x);

std::unordered_map<std::string, std::string>  Amino_label = {
        {"A","GCN"},{"C","UGY"},{"D","GAY"},{"E","GAR"},{"F","UUY"},
        {"G","GGN"},{"H","CAY"},{"I","AUH"},{"K","AAR"},{"L","YVn"},
        {"M","AUG"},{"N","AAY"},{"P","CCN"},{"Q","CAR"},{"R","MOn"},
        {"S","WSz"},{"T","ACN"},{"V","GUN"},{"W","UGG"},{"Y","UAY"},
        {"*","NNN"}
    };

std::unordered_map<char, std::vector<std::string>>  Amino_to_nucs = {
    {'A',{"GCU","GCC","GCA","GCG"}},
    {'R',{"CGU","CGC","CGA","CGG","AGA","AGG"}},
    {'N',{"AAU","AAC"}},
    {'D',{"GAU","GAC"}},
    {'C',{"UGU","UGC"}},
    {'Q',{"CAA","CAG"}},
    {'E',{"GAA","GAG"}},
    {'G',{"GGU","GGC","GGA","GGG"}},
    {'H',{"CAU","CAC"}},
    {'M',{"AUG"}},
    {'I',{"AUU","AUC","AUA"}},
    {'L',{"CUU","CUC","CUA","CUG","UUA","UUG"}},
    {'K',{"AAA","AAG"}},
    {'F',{"UUU","UUC"}},
    {'P',{"CCU","CCC","CCA","CCG"}},
    {'S',{"UCU","UCC","UCA","UCG","AGU","AGC"}},
    {'T',{"ACU","ACC","ACA","ACG"}},
    {'W',{"UGG"}},
    {'Y',{"UAU","UAC"}},
    {'V',{"GUU","GUC","GUA","GUG"}},
};

std::vector<std::vector<std::string>> CodonSet=
{
    {"A","GCU"},{"A","GCC"},{"A","GCA"},{"A","GCG"},
    {"R","CGU"},{"R","CGC"},{"R","CGA"},{"R","CGG"},{"R","AGA"},{"R","AGG"},
    {"N","AAU"},{"N","AAC"},
    {"D","GAU"},{"D","GAC"},
    {"C","UGU"},{"C","UGC"},
    {"Q","CAA"},{"Q","CAG"},
    {"E","GAA"},{"E","GAG"},
    {"G","GGU"},{"G","GGC"},{"G","GGA"},{"G","GGG"},
    {"H","CAU"},{"H","CAC"},
    {"M","AUG"},
    {"I","AUU"},{"I","AUC"},{"I","AUA"},
    {"L","CUU"},{"L","CUC"},{"L","CUA"},{"L","CUG"},{"L","UUA"},{"L","UUG"},
    {"K","AAA"},{"K","AAG"},
    {"F","UUU"},{"F","UUC"},
    {"P","CCU"},{"P","CCC"},{"P","CCA"},{"P","CCG"},
    {"S","UCU"},{"S","UCC"},{"S","UCA"},{"S","UCG"},{"S","AGU"},{"S","AGC"},
    {"T","ACU"},{"T","ACC"},{"T","ACA"},{"T","ACG"},
    {"W","UGG"},
    {"Y","UAU"},{"Y","UAC"},
    {"V","GUU"},{"V","GUC"},{"V","GUA"},{"V","GUG"},
};

std::unordered_map<std::string,int> sp_loops = {
    {"CAACG",680},{"GUUAC",690},
    {"CAACGG",550},{"CCAAGG",330},{"CCACGG",370},{"CCCAGG",340},{"CCGAGG",350},{"CCGCGG",360},{"CCUAGG",370},
    {"CCUCGG",250},{"CUAAGG",360},{"CUACGG",280},{"CUCAGG",370},{"CUCCGG",270},{"CUGCGG",280},{"CUUAGG",350},
    {"CUUCGG",370},{"CUUUGG",370},
    {"ACAGUACU",280},{"ACAGUGAU",360},{"ACAGUGCU",290},{"ACAGUGUU",180}
};

std::vector<std::string> sp_loops_vec = {
    "CAACG","GUUAC",
    "CAACGG","CCAAGG","CCACGG","CCCAGG","CCGAGG","CCGCGG","CCUAGG","CCUCGG",
    "CUAAGG","CUACGG","CUCAGG","CUCCGG","CUGCGG","CUUAGG","CUUCGG","CUUUGG",
    "ACAGUACU","ACAGUGAU","ACAGUGCU","ACAGUGUU"
};

std::unordered_map<std::string, std::string>  re_Amino_label = 
{
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

std::unordered_map<char, std::vector<int>> Base_table={
    {'N',{BASE("A"),BASE("C"),BASE("G"),BASE("U")}},
    {'z',{BASE("A"),BASE("C"),BASE("G"),BASE("U"),BASE("U1"),BASE("C1")}},
    {'n',{BASE("A1"),BASE("G1"),BASE("A2"),BASE("G2"),BASE("C"),BASE("U")}},
    {'A',{BASE("A")}},
    {'U',{BASE("U")}},
    {'G',{BASE("G")}},
    {'C',{BASE("C")}},
    {'W',{BASE("A"),BASE("U")}},
    {'S',{BASE("G"),BASE("C")}},
    {'M',{BASE("A"),BASE("C")}},
    {'R',{BASE("A"),BASE("G")}},
    {'Y',{BASE("C"),BASE("U")}},
    {'H',{BASE("A"),BASE("C"),BASE("U")}},
    {'V',{BASE("U1"),BASE("U_AG2"),BASE("U_AG1")}},
    {'O',{BASE("G_CU"),BASE("G_AG2"),BASE("G_AG1")}}
};

class InputParser{
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }

        std::string getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            return "";
        }

        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

bool FindOption(int &argc, char **argv, std::string option){
    for (int i=1; i < argc; ++i){
        if(std::string(argv[i]) == option)
            return true;
    }
    return false;
}

void PrintHelp() {
    std::cout
          << "USAGE:\n"
          << "  LinearCDSfold [OPTIONS] <SEQUENCE_FILE>\n"
          << "\n"
          << "ARGUMENTS:\n"
          << "  <SEQUENCE_FILE>           Amino acid sequence in FASTA format\n"
          << "\n"
          << "OPTIONS:\n"
          << "  -c <CODON_USAGE_FILE>     Codon usage frequency table\n"
          << "  -O <OBJECTIVE_FUNCTION>   Objective function for MFE and CAI\n"
          << "  -l <LAMBDA>               Scaling weight between MFE and CAI\n"
          << "  -m <SEARCH_MODE>, --mode <SEARCH_MODE>\n"
          << "                            Search strategy:\n"
          << "                              exact   – Exact search (default)\n"
          << "                              beam    – Beam search\n"
          << "                              pareto  – Pareto-optimal search using DERNA objective\n"
          << "  -b <BEAM_SIZE>            Beam size for beam search\n"
          << "  -t <TAU1>, --tau1 <TAU1>  Termination threshold for Pareto-optimal search\n"
          << "  -u <TAU2>, --tau2 <TAU2>  Secondary threshold for Pareto-optimal search\n"
          << "  -o <FILE_NAME>            Output file for detailed results\n"
          << "  -f <FILE_NAME>            Output file for MFE/CAI results\n";
}

enum Manner {
    NONE = 0,
    MANNER_NONEtoN,
    MANNER_NONEtoF,
    MANNER_NONEtoS,
    MANNER_N_EtoN,
    MANNER_NtoC,
    MANNER_S_HP,
    MANNER_S_EtoS,
    MANNER_CtoM1,
    MANNER_M1_CtoM2,
    MANNER_M1_EtoM1,
    MANNER_M2toM1,
    MANNER_M2toMulti,
    MANNER_S_M2toMulti,
    MANNER_Multi_EtoMulti,
    MANNER_MultitoC,
    MANNER_F_EtoF,
    MANNER_CtoF,
    MANNER_F_CtoF,
    MANNER_C_StoCS,
    MANNER_S_CStoC,
    MANNER_S_C_StoC,
    MANNER_CStoC,
    MANNER_C_StoC,
    MANNER_CtoC,
    MANNER_S_CtoC
};

std::string mannerToString(Manner manner) {
    switch (manner) {
        case NONE:
            return "NONE";
        case MANNER_NONEtoN:
            return "MANNER_NONEtoN";
        case MANNER_NONEtoF:
            return "MANNER_NONEtoF";
        case MANNER_NONEtoS:
            return "MANNER_NONEtoS";
        case MANNER_N_EtoN:
            return "MANNER_N_EtoN";
        case MANNER_NtoC:
            return "MANNER_NtoC";
        case MANNER_S_HP:
            return "MANNER_S_HP";
        case MANNER_S_EtoS:
            return "MANNER_S_EtoS";
        case MANNER_CtoC:
            return "MANNER_CtoC";
        case MANNER_S_CtoC:
            return "MANNER_S_CtoC";
        case MANNER_CtoM1:
            return "MANNER_CtoM1";
        case MANNER_M1_CtoM2:
            return "MANNER_M1_CtoM2";
        case MANNER_M1_EtoM1:
            return "MANNER_M1_EtoM1";
        case MANNER_M2toM1:
            return "MANNER_M2toM1";
        case MANNER_M2toMulti:
            return "MANNER_M2toMulti";
        case MANNER_S_M2toMulti:
            return "MANNER_S_M2toMulti";
        case MANNER_Multi_EtoMulti:
            return "MANNER_Multi_EtoMulti";
        case MANNER_MultitoC:
            return "MANNER_MultitoC";
        case MANNER_F_EtoF:
            return "MANNER_F_EtoF";
        case MANNER_CtoF:
            return "MANNER_CtoF";
        case MANNER_F_CtoF:
            return "MANNER_F_CtoF";
        case MANNER_CStoC:
            return "MANNER_CStoC";
        case MANNER_C_StoC:
            return "MANNER_C_StoC";
        case MANNER_S_CStoC:
            return "MANNER_S_CStoC";
        case MANNER_S_C_StoC:
            return "MANNER_S_C_StoC";
        case MANNER_C_StoCS:
            return "MANNER_C_StoCS";
        default:
            return "UNKNOWN";
    }
}



// inline int GetIndexCS(int i, int nuci, int nuci_pair,int nucj, int nucx ,int len){//i : start posi, nuci :start nuc, nuci_pair : the nuc pair with nuci ,nucj:end nuc, nucx: single start, len : single length
//     return i*fourthPower_base_number*(SINGLE_MAX_LEN+1) + len * fourthPower_base_number+ nuci*cube_base_number+nuci_pair*square_base_number + nucj*base_number + nucx;
// }
// inline std::tuple<int, int, int, int, int, int> GetIndexTupleCS(int index_i_j) {
//     int i = index_i_j/(fourthPower_base_number*(SINGLE_MAX_LEN+1));
//     int temp = index_i_j % (fourthPower_base_number*(SINGLE_MAX_LEN+1));

//     int len = index_i_j/(fourthPower_base_number);
//     temp = index_i_j % (fourthPower_base_number);

//     int nuci = temp / (cube_base_number);
//     temp = temp % (cube_base_number);

//     int nuci_pair = temp  / square_base_number;
//     temp  = temp % square_base_number;

//     int nucj = temp  / base_number;
//     int nucx = temp  % base_number;

//     return std::make_tuple(i, nuci, nuci_pair,nucj, nucx, len);
// }

inline int GetIndex(int i, int nuci, int nucj){
    return i*square_base_number+nuci*base_number+nucj;
}

// speedup
inline int GetIndexSpeedup(int i, int nuci, int nuci1,int nucj_1, int nucj){//i : start posi, nuci :start nuc, nuci_pair : the nuc pair with nuci ,nucj:end nuc, nucx: single start, len : single length
    return i*fourthPower_base_number + nuci* cube_base_number+ nuci1*square_base_number + nucj_1*base_number + nucj;
}
inline std::tuple<int, int, int, int, int> GetIndexTupleSpeedup(int index_i_j) {
    int i = index_i_j/fourthPower_base_number;
    int temp = index_i_j % fourthPower_base_number;

    int nuci = temp/cube_base_number;
    temp = temp % cube_base_number;

    int nuci1 = temp / (square_base_number);
    temp = temp % (square_base_number);

    int nucj_1 = temp  / base_number;
    int nucj = temp  % base_number;

    return std::make_tuple(i, nuci,nuci1,nucj_1,nucj);
}

inline int GetIndexCS(int i, int nuci, int nuci_pair,int nucj, int nucx ,int len){//i : start posi, nuci :start nuc, nuci_pair : the nuc pair with nuci ,nucj:end nuc, nucx: single start, len : single length
    return i*cube_base_number*(SINGLE_MAX_LEN+1) + len * cube_base_number+ nuci*square_base_number + nucj*base_number + nucx;
}
inline std::tuple<int, int, int, int, int, int> GetIndexTupleCS(int index_i_j) {
    int i = index_i_j/(cube_base_number*(SINGLE_MAX_LEN+1));
    int temp = index_i_j % (cube_base_number*(SINGLE_MAX_LEN+1));

    int len = temp/(cube_base_number);
    temp = temp % (cube_base_number);

    int nuci = temp / (square_base_number);
    temp = temp % (square_base_number);

    int nucj = temp  / base_number;
    int nucx = temp  % base_number;

    return std::make_tuple(i, nuci, 0,nucj, nucx, len);
}

inline std::tuple<int, int, int> GetIndexTuple(int index_i_j) {
    int i = index_i_j / square_base_number;
    int nuc_index = index_i_j % square_base_number;
    int nuci = nuc_index / base_number;
    int nucj = nuc_index % base_number;

    return std::make_tuple(i, nuci, nucj);
}

template<typename T>
struct State {
    T score;
    int index_1;
    int index_2;
    int index_3;
    Manner MANNER;
    int last_pair_pos;

    State() : score(VALUE_MIN<T>()), index_1(-1), index_2(-1), index_3(-1),MANNER(NONE), last_pair_pos(-1){}
    State(double score_, const int& index_1_, const int& index_2_, const int& index_3_, Manner MANNER_)
        : score(score_), index_1(index_1_), index_2(index_2_), index_3(index_3_), MANNER(MANNER_){}
    State(double score_, const int& index_1_, const int& index_2_, const int& index_3_,Manner MANNER_, const int& last_pair_pos_)
        : score(score_), index_1(index_1_), index_2(index_2_), index_3(index_3_), MANNER(MANNER_), last_pair_pos(last_pair_pos_){}
};

template<typename T>
class AllTables {

public:
    std::vector<std::unordered_map<int, State<T>>> bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>> bestS;
    std::vector<std::unordered_map<int, State<T>>> bestF;
    std::vector<std::unordered_map<int, State<T>>> bestC;
    std::vector<std::unordered_map<int, State<T>>> bestCS;
    std::vector<std::unordered_map<int, State<T>>> bestM1;
    std::vector<std::unordered_map<int, State<T>>> bestM2;
    std::vector<std::unordered_map<int, State<T>>> bestMulti;
    std::vector<std::vector<std::unordered_map<int, State<T>>>> speedup;//speedup[j][l][index]

    AllTables(std::string rna_seq, int size_) {
        bestN = initi_vec(rna_seq, size_);
        bestS = initi_vec_S(rna_seq, size_, SINGLE_MAX_LEN+1);
        bestF = initi_vec(rna_seq, size_);
        bestC = initi_vec(rna_seq, size_);
        bestCS = initi_vec(rna_seq, size_);
        bestM1 = initi_vec(rna_seq, size_);
        bestM2 = initi_vec(rna_seq, size_);
        bestMulti = initi_vec(rna_seq, size_);
        speedup = initi_vec_S(rna_seq, size_, SINGLE_MAX_LEN+1);//speedup
    }

    std::vector<std::unordered_map<int, State<T>>> initi_vec(std::string& rna_seq, int size_){
        std::vector<std::unordered_map<int, State<T>>> seq_tables_(size_);
        return seq_tables_;
    }
    std::vector<std::vector<std::unordered_map<int, State<T>>>> initi_vec_S(std::string& rna_seq, int size_, int len){
        std::vector<std::vector<std::unordered_map<int, State<T>>>> seq_tables_(size_,std::vector<std::unordered_map<int, State<T>>>(len));
        return seq_tables_;
    }

};


int BASE(const std::string& x) {
    static const std::unordered_map<std::string, int> baseMap = {
        {"A", 1},
        {"C", 2},
        {"G", 3},
        {"U", 4},
        {"G_AG1", 5},
        {"G_AG2", 6},
        {"U_AG1", 7},
        {"U_AG2", 8},
        {"A1", 9},
        {"A2", 10},
        {"G1", 11},
        {"G2", 12},
        {"G_CU", 13},
        {"U1", 14},
        {"C1", 15}
    };

    auto it = baseMap.find(x);
    if (it != baseMap.end()) {
        return it->second;
    } else {
        return 0;
    }
}

static const int EnrBASEArr [] = { -1, 1, 2, 3, 4, 3, 3, 4, 4, 1, 1, 3, 3, 3, 4, 2 };

int EnrBASE(int x) {
    return EnrBASEArr[x];
}


std::string CheckReBASE(int x) {
    static const std::unordered_map<int, std::string> baseMap = {
        {1, "A"},
        {2, "C"},
        {3, "G"},
        {4, "U"},
        {5, "G_AG1"},
        {6, "G_AG2"},
        {7, "U_AG1"},
        {8, "U_AG2"},
        {9, "A1"},
        {10, "A2"},
        {11, "G1"},
        {12, "G2"},
        {13, "G_CU"},
        {14, "U1"},
        {15, "C1"}
    };

    auto it = baseMap.find(x);
    if (it != baseMap.end()) {
        return it->second;
    } else {
        return "*";
    }
}

char reBASE(int x) {
    static const std::unordered_map<int, char> baseMap = {
        {1, 'A'},
        {2, 'C'},
        {3, 'G'},
        {4, 'U'},
        {5, 'G'},
        {6, 'G'},
        {7, 'U'},
        {8, 'U'},
        {9, 'A'},
        {10, 'A'},
        {11, 'G'},
        {12, 'G'},
        {13, 'G'},
        {14, 'U'},
        {15, 'C'},
        {-1, '*'}
    };

    auto it = baseMap.find(x);
    if (it != baseMap.end()) {
        return it->second;
    } else {
        return '*';
    }
}

void initialize(){
    // AU
    // UA
    std::vector<int> base_A_vec = {BASE("A"), BASE("A1"), BASE("A2")};
    std::vector<int> base_U_vec = {BASE("U"), BASE("U_AG1"), BASE("U_AG2"), BASE("U1")};

    for (int i = 0; i < base_A_vec.size(); i++){
        for (int j = 0; j < base_U_vec.size(); j++){
            _allowed_pairs[base_A_vec[i]][base_U_vec[j]] = true;
            _allowed_pairs[base_U_vec[j]][base_A_vec[i]] = true;
        }
    }

    //CG
    //GC
    std::vector<int> base_C_vec = {BASE("C"), BASE("C1")};
    std::vector<int> base_G_vec = {BASE("G"), BASE("G_AG1"), BASE("G_AG2"), BASE("G1"), BASE("G2"), BASE("G_CU")};
    
    for (int i = 0; i < base_C_vec.size(); i++){
        for (int j = 0; j < base_G_vec.size(); j++){
            _allowed_pairs[base_C_vec[i]][base_G_vec[j]] = true;
            _allowed_pairs[base_G_vec[j]][base_C_vec[i]] = true;
        }
    }

    //UG
    //GU
    for (int i = 0; i < base_U_vec.size(); i++){
        for (int j = 0; j < base_G_vec.size(); j++){
            _allowed_pairs[base_U_vec[i]][base_G_vec[j]] = true;
            _allowed_pairs[base_G_vec[j]][base_U_vec[i]] = true;
        }
    }
     
    // Arginine
    xcon_table[Arg][0][BASE("A")][BASE("G_AG1")] = true; //R
    xcon_table[Arg][0][BASE("C")][BASE("G_AG2")] = true; //R
    xcon_table[Arg][0][BASE("C")][BASE("G_CU")] = true; //R
    xcon_table[Arg][1][BASE("G_AG1")][BASE("A1")] = true; //R
    xcon_table[Arg][1][BASE("G_AG1")][BASE("G1")] = true; //R
    xcon_table[Arg][1][BASE("G_AG2")][BASE("A2")] = true; //R 
    xcon_table[Arg][1][BASE("G_AG2")][BASE("G2")] = true; //R
    xcon_table[Arg][1][BASE("G_CU")][BASE("C")] = true; //R
    xcon_table[Arg][1][BASE("G_CU")][BASE("U")] = true; //R

    // Leucine
    xcon_table[Leu][0][BASE("U")][BASE("U_AG1")] = true; //L
    xcon_table[Leu][0][BASE("C")][BASE("U_AG2")] = true; //L
    xcon_table[Leu][0][BASE("C")][BASE("U1")] = true; //L
    xcon_table[Leu][1][BASE("U_AG1")][BASE("A1")] = true; //L
    xcon_table[Leu][1][BASE("U_AG1")][BASE("G1")] = true; //L
    xcon_table[Leu][1][BASE("U_AG2")][BASE("A2")] = true; //L
    xcon_table[Leu][1][BASE("U_AG2")][BASE("G2")] = true; //L
    xcon_table[Leu][1][BASE("U1")][BASE("C")] = true; //L
    xcon_table[Leu][1][BASE("U1")][BASE("U")] = true; //L

    //Serine
    xcon_table[Ser][0][BASE("U")][BASE("C")] = true; //S
    xcon_table[Ser][0][BASE("A")][BASE("G")] = true; //S
    xcon_table[Ser][1][BASE("C")][BASE("U")] = true; //S
    xcon_table[Ser][1][BASE("C")][BASE("C")] = true; //S
    xcon_table[Ser][1][BASE("C")][BASE("A")] = true; //S
    xcon_table[Ser][1][BASE("C")][BASE("G")] = true; //S
    xcon_table[Ser][1][BASE("G")][BASE("U1")] = true; //S
    xcon_table[Ser][1][BASE("G")][BASE("C1")] = true; //S

}

//...[i][j]...
inline bool IsLegal(std::vector<int>& con_seq, int base_i, int base_i1, int pos_i){
    //j must be i+1
    int amino_ = con_seq[pos_i];
    int pos_i1 = pos_i+1;
    if((pos_i/3) == (pos_i1/3)){
        if(amino_ != normal_ami){
            int amino_index = pos_i%3;
            if(xcon_table[amino_][amino_index][base_i][base_i1]){
                return true;}else{return false;}
        }else{return true;}   
    }else{return true;}
    
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void ami_to_rna(std::vector<std::string>& rna_seq_list, std::string& ami_seq){
    std::string rna_seq;
    std::string Amino1;
    std::string codon;
    for(int i = 0; i < ami_seq.size(); i++){
        Amino1 = ami_seq[i];
        codon = Amino_label.at(Amino1);
        for (int j = 0; j < codon.size(); j++){
            rna_seq.push_back(codon[j]);
        }       
    }
    rna_seq_list.push_back(rna_seq);
}

void add_con_seq(std::vector<std::vector<int>>& con_seq_list, std::string& ami_seq){
    std::vector<int> con_seq;
    char Amino1; 
    for(int i=0; i<ami_seq.size(); i++){
        Amino1 = ami_seq[i];
        if(Amino1=='R'){
            for(int j=0; j<3; j++)
                con_seq.push_back(Arg); //R==0
        }else if(Amino1=='L'){
            for(int j=0; j<3; j++)
                con_seq.push_back(Leu); //L==1
        }else if(Amino1=='S'){
            for(int j=0; j<3; j++)
                con_seq.push_back(Ser); //S==2
        }else{
            for(int j=0; j<3; j++)
                con_seq.push_back(normal_ami); //-1
        }
    }
    con_seq_list.push_back(con_seq);
}

void check_valid_ami(int seq_size, std::string result_seq, std::string& ami_seq){
    std::string sub_str;
    std::string index;
    for(int i=0; i<seq_size; i++){
        sub_str =  result_seq.substr(3*i,3);
        if(re_Amino_label.count(sub_str)){
            index = re_Amino_label.at(sub_str);    
        }else{
            std::cout << "error aminoacid: ";
            std::cout << sub_str << std::endl;
            std::cout << "error index: ";
            std::cout << i << std::endl;
            break;
        }
                
        ami_seq.push_back(index[0]);
        
    }
}

inline void log_CAI(std::vector<double>& cai_vector, std::vector<double>& log_cai_vector){

    for(int i = 0; i < cai_vector.size(); i++){
        if(cai_vector[i] == 0)
            log_cai_vector.push_back(-999999);
        else log_cai_vector.push_back(log(cai_vector[i]));// Round the number to three decimal places.
        //else log_cai_vector.push_back(logf(cai_vector[i]));
    }
}

inline int getLastExtendedNuc(char amino_ ,std::string codon_){
    int last_nuc = -1;

    if(amino_ == 'R'){
        if(codon_ == "AGA")
            last_nuc = BASE("A1");
        if(codon_ == "AGG")
            last_nuc = BASE("G1");
        
        if(codon_ == "CGA")
            last_nuc = BASE("A2");
        if(codon_ == "CGG")
            last_nuc = BASE("G2");
        if(codon_ == "CGC")
            last_nuc = BASE("C");
        if(codon_ == "CGU")
            last_nuc = BASE("U");
        // {"CGU","R"}, {"CGC","R"}, {"CGA","R"}, {"CGG","R"}, {"AGA","R"}, {"AGG","R"},
    }
    else if(amino_ == 'L'){
        if(codon_ == "CUU")
            last_nuc = BASE("U");
        if(codon_ == "CUC")
            last_nuc = BASE("C");
        if(codon_ == "CUA")
            last_nuc = BASE("A2");
        if(codon_ == "CUG")
            last_nuc = BASE("G2");
        
        if(codon_ == "UUA")
            last_nuc = BASE("A1");
        if(codon_ == "UUG")
            last_nuc = BASE("G1");
        // {"CUU","L"}, {"CUC","L"}, {"CUA","L"}, {"CUG","L"}, {"UUA","L"}, {"UUG","L"},
    }
    else if(amino_ == 'S'){
        if(codon_ == "UCU")
            last_nuc = BASE("U");
        if(codon_ == "UCC")
            last_nuc = BASE("C");
        if(codon_ == "UCA")
            last_nuc = BASE("A");
        if(codon_ == "UCG")
            last_nuc = BASE("G");
        
        if(codon_ == "AGU")
            last_nuc = BASE("U1");
        if(codon_ == "AGC")
            last_nuc = BASE("C1");
        // {"UCU","S"}, {"UCC","S"}, {"UCA","S"}, {"UCG","S"}, {"AGU","S"}, {"AGC","S"},
    }
    else{
        last_nuc = BASE(std::string(1, codon_[2]));
    }
    return last_nuc;
}

// inline void Read_CAI_file(std::string cai_file_path){
//     std::ifstream cai_file(cai_file_path);
//     if(!cai_file.is_open()){
//         std::cout<<"Can't open "<<cai_file_path<<"."<<std::endl;
//         return ;
//     }

//     std::string line;
//     bool first_line = true;
//     while(getline(cai_file, line)){
//         if(first_line){
//             first_line = false;
//             continue;
//         }

//         std::stringstream ss(line);
//         std::vector<std::string> row;
//         std::string cell;

//         while(getline(ss, cell ,',')){
//             row.push_back(cell);
//             std::cout<<cell<<",";
//         }
//         std::cout<<std::endl;

//         char amino_ = row[1][0];
//         std::string codon_ = row[0];
//         int last_nuc = BASE(std::string(1, codon_[2]));
//         double cai_value_ = std::stod(row[2]);
//         double cai_value = (cai_value_ == 0)? -999999 : log(cai_value_);

//         if(amino_ == 'R'){
//             if(codon_ == "AGA")
//                 last_nuc = BASE("A1");
//             if(codon_ == "AGG")
//                 last_nuc = BASE("G1");
            
//             if(codon_ == "CGA")
//                 last_nuc = BASE("A2");
//             if(codon_ == "CGG")
//                 last_nuc = BASE("G2");
//             if(codon_ == "CGC")
//                 last_nuc = BASE("C");
//             if(codon_ == "CGU")
//                 last_nuc = BASE("U");
//             // {"CGU","R"}, {"CGC","R"}, {"CGA","R"}, {"CGG","R"}, {"AGA","R"}, {"AGG","R"},
//         }
//         else if(amino_ == 'L'){
//             if(codon_ == "CUU")
//                 last_nuc = BASE("U");
//             if(codon_ == "CUC")
//                 last_nuc = BASE("C");
//             if(codon_ == "CUA")
//                 last_nuc = BASE("A2");
//             if(codon_ == "CUG")
//                 last_nuc = BASE("G2");
            
//             if(codon_ == "UUA")
//                 last_nuc = BASE("A1");
//             if(codon_ == "UUG")
//                 last_nuc = BASE("G1");
//             // {"CUU","L"}, {"CUC","L"}, {"CUA","L"}, {"CUG","L"}, {"UUA","L"}, {"UUG","L"},
//         }
//         else if(amino_ == 'S'){
//             if(codon_ == "UCU")
//                 last_nuc = BASE("U");
//             if(codon_ == "UCC")
//                 last_nuc = BASE("C");
//             if(codon_ == "UCA")
//                 last_nuc = BASE("A");
//             if(codon_ == "UCG")
//                 last_nuc = BASE("G");
            
//             if(codon_ == "AGU")
//                 last_nuc = BASE("U1");
//             if(codon_ == "AGC")
//                 last_nuc = BASE("C1"); 
//             // {"UCU","S"}, {"UCC","S"}, {"UCA","S"}, {"UCG","S"}, {"AGU","S"}, {"AGC","S"},
//         }
//         unweighted_CAIMap[codon_] = cai_value;
//         CodonSetCAIMap[amino_][last_nuc] = cai_value;
//     }
// }

// inline void initialize_CAI_table2(bool is_DN){
//     for (auto &item : CodonSetCAIMap) {
//         std::cout<<item.first<<",";
//         for(auto &item2 : item.second){
//             std::cout<<item2.first<<","<<item2.second<<std::endl;
//         }
//     }
// }

bool isCSV(const std::string& file_path) {
    // Check if file_path ends with ".csv"
    return file_path.size() >= 4 && file_path.substr(file_path.size() - 4) == ".csv";
}

inline void read_CAI_table_csv(std::string csv_file_name, std::vector<double>& cai_vector_){
    std::ifstream cai_csv_file(csv_file_name);

     if (!cai_csv_file.is_open()) {
        std::cerr << "[Error] Cannot open csv file: " << csv_file_name << std::endl;
        return;
    }

    std::unordered_map<std::string, float> max_amino_cai_map;
    std::string line;
    while (std::getline(cai_csv_file, line)) {
        std::stringstream ss(line);
        std::string key1, key2;
        float value;

        std::getline(ss, key1, ',');  // GCA
        std::getline(ss, key2, ',');  // A
        ss >> value;

        cai_table_csv[key1] = value;
        auto it = max_amino_cai_map.find(key2);
        if (it == max_amino_cai_map.end()) {
            max_amino_cai_map[key2] = value;
        } else {
            if (value > it->second) {
                it->second = value;
            }
        }
    }

    cai_csv_file.close();

    // elements in CodonSet
    for (auto& codon_info : CodonSet) {
        std::string codon_ = codon_info[1];
        std::string amino_ = codon_info[0];
        auto it = cai_table_csv.find(codon_);
        if (it == cai_table_csv.end()) {
            std::cout << "[Error] codon " << codon_ << " not found in file: " << csv_file_name << std::endl;
            continue;  // skip or handle error
        }
        double cai_value = ((it->second) / max_amino_cai_map[amino_]);
        cai_vector_.push_back(cai_value);
    }

    // // debug: check the items in cai_table_csv
    // for (const auto& pair : cai_table_csv) {
    //     std::cout << "Codon: " << pair.first << ", CAI Value: " << pair.second << std::endl;
    // }
}

inline void initialize_CAI_table(std::vector<double> cai_vector_, bool is_DN){
    std::vector<double> cai_vector;

    log_CAI(cai_vector_, cai_vector);


    for (int k = 0; k < CodonSet.size(); ++k) {
        char amino_ = CodonSet[k][0][0];
        if (is_DN) {
            if (CodonSetCAIMap_DN.find(amino_) == CodonSetCAIMap_DN.end())
                CodonSetCAIMap_DN[amino_] = std::unordered_map<int, double>();
        } else {
            if (CodonSetCAIMap.find(amino_) == CodonSetCAIMap.end())
                CodonSetCAIMap[amino_] = std::unordered_map<int, double>();
                // CodonSetCAIMap[amino_] = std::unordered_map<int, int>();
        }
    }


    for (int i = 0; i < CodonSet.size(); ++i){

        char amino_ = CodonSet[i][0][0];
        std::string codon_ = CodonSet[i][1];
        int last_nuc = BASE(std::string(1, codon_[2]));

        if(amino_ == 'R'){
            if(codon_ == "AGA")
                last_nuc = BASE("A1");
            if(codon_ == "AGG")
                last_nuc = BASE("G1");
            
            if(codon_ == "CGA")
                last_nuc = BASE("A2");
            if(codon_ == "CGG")
                last_nuc = BASE("G2");
            if(codon_ == "CGC")
                last_nuc = BASE("C");
            if(codon_ == "CGU")
                last_nuc = BASE("U");
            // {"CGU","R"}, {"CGC","R"}, {"CGA","R"}, {"CGG","R"}, {"AGA","R"}, {"AGG","R"},
        }
        else if(amino_ == 'L'){
            if(codon_ == "CUU")
                last_nuc = BASE("U");
            if(codon_ == "CUC")
                last_nuc = BASE("C");
            if(codon_ == "CUA")
                last_nuc = BASE("A2");
            if(codon_ == "CUG")
                last_nuc = BASE("G2");
            
            if(codon_ == "UUA")
                last_nuc = BASE("A1");
            if(codon_ == "UUG")
                last_nuc = BASE("G1");
            // {"CUU","L"}, {"CUC","L"}, {"CUA","L"}, {"CUG","L"}, {"UUA","L"}, {"UUG","L"},
        }
        else if(amino_ == 'S'){
            if(codon_ == "UCU")
                last_nuc = BASE("U");
            if(codon_ == "UCC")
                last_nuc = BASE("C");
            if(codon_ == "UCA")
                last_nuc = BASE("A");
            if(codon_ == "UCG")
                last_nuc = BASE("G");
            
            if(codon_ == "AGU")
                last_nuc = BASE("U1");
            if(codon_ == "AGC")
                last_nuc = BASE("C1"); 
            // {"UCU","S"}, {"UCC","S"}, {"UCA","S"}, {"UCG","S"}, {"AGU","S"}, {"AGC","S"},
        }
        if(is_DN){
            double weight_ = (1-lambda)*cai_vector[i]; //*100
            CodonSetCAIMap_DN[amino_][last_nuc] = weight_;
            CalCAIMap_DN[codon_] = weight_;
        }
        else{
            double weight_ = LDCONST*lambda * cai_vector[i];//*100
            // CodonSetCAIMap[amino_][last_nuc] = int(weight_);
            // CalCAIMap[codon_] = int(weight_);
            CodonSetCAIMap[amino_][last_nuc] = weight_;
            CalCAIMap[codon_] = weight_;

        }
        unweighted_CAIMap[codon_] = cai_vector[i];
        // std::cout << "[Debug] Codon: " << codon_ << ", CAI Value: " << cai_vector[i] << std::endl;
    }
}

void check_ami_solution(std::string ami_seq, std::string rna_solution){
    std::string check_ami_seq;
    check_valid_ami(ami_seq.size(), rna_solution, check_ami_seq);


    // //check correct
    if(ami_seq != check_ami_seq){
        std::cout << "Current backtracked amino acids sequence:" << std::endl;
        std::cout << check_ami_seq << std::endl;
        std::cout << "ERROR: Codon sequences cannot transfer into amino acids." << std::endl;
    }
}

template<typename T>
T GetCAIScore(const std::string& rna_solution,bool is_DN){
    T score = 0;
    if(is_DN){
        for (int i = 0; i <= rna_solution.size() - 3; i += 3){
            std::string codon = rna_solution.substr(i, 3);
            score = score + CalCAIMap_DN[codon];
        }
    }
    else{
        for (int i = 0; i <= rna_solution.size() - 3; i += 3){
            std::string codon = rna_solution.substr(i, 3);
            score = score + CalCAIMap[codon];
        }
    }
    
    return score;
}

double GetUnweghtedCAIScore(const std::string& rna_solution){
    double score = 0;
    for (int i = 0; i <= rna_solution.size() - 3; i += 3){
        std::string codon = rna_solution.substr(i, 3);
        score = score + unweighted_CAIMap[codon];
    }
    return score;
}
double round_up(double x){
    return round(x*1000000.0)/1000000.0;
}

int counter = 0;
double percent = 0;
std::string bar = "                                                 ";
void Processing(int size_){
    counter ++;
    double real_precent = double(counter)/double(size_);
    if(real_precent >= percent){
        percent += 0.02;
        bar = "=" + bar;
        bar = bar.substr(0,50);
    }
    std::cout << std::fixed<<std::setw(70)<< std::left<<("\rProcessing: [" + bar + "]  " + std::to_string(int(real_precent*100))) + "%"<<std::flush;
    
    if(size_-1 == counter){
        std::cout <<std::fixed<<std::setw(70)<< std::left<<("\rProcessing: [" + bar + "]  100%")<<std::endl;
         counter = 0;
         percent = 0;
         bar = "                                                                                                    ";
    }
       
}

inline void PrintInfo(std::string output_txt,std::string output_csv,std::string file, int beamsize, std::string cai_file_path, double lambda, std::string objective, bool pareto, bool is_rna_file,bool output){
    // std out information
    std::ofstream outputfile(output_txt);
    std::ofstream outputcsv(output_csv);

    /*AA file*/
    if(!is_rna_file){
        if(output) std::cout << "Amino acid file: " << file << std::endl;
        outputfile << "Amino acid file: " << file << std::endl;
    }
    else{
        if(output) std::cout << "RNA file: " << file << std::endl;
        outputfile << "Amino acid file: " << file << std::endl;
    }

    /*CAI file*/
    if(output) std::cout << "Codon usage table: " << cai_file_path << std::endl;
    outputfile << "Codon usage table: " << cai_file_path << std::endl;

    /*initialize CVS file*/
    outputcsv <<"lambda,MFE,CAI"<<std::endl;

    /*OJ function*/
    if (objective == "DN") {
        if(output)std::cout << "Objective function: DERNA" << std::endl;
        outputfile << "Objective function: DERNA" << std::endl;
    } else {
        if(output)std::cout << "Objective function: LinearDesign" << std::endl;
        outputfile << "Objective function: LinearDesign" << std::endl;
    }
    if(pareto){
        if(output)std::cout << "Search mode: Pareto-optimal search" << std::endl;
        outputfile << "Search mode: Pareto-optimal search" << std::endl;
    }
    else{
        if(output){
            if(beamsize != 0)
                std::cout << "Search mode: Beam search" << std::endl;
            else std::cout << "Search mode: Exact search" << std::endl;
        }
        if(beamsize != 0)
            outputfile << "Search mode: Beam search" << std::endl;
        else outputfile << "Search mode: Exact search" << std::endl;
    }
    if(output){
        if(beamsize != 0) 
            std::cout << "Beam size: " << beamsize << std::endl;
    } 
    if(beamsize != 0) 
        outputfile << "Beam size: " << beamsize << std::endl;
    // if(output)std::cout << "------------------------------------------ START ------------------------------------------" << std::endl;
    // outputfile << "------------------------------------------ START ------------------------------------------" << std::endl;

}