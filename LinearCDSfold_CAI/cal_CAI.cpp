#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <unordered_map>
#include <vector>


using namespace std;

vector<string> codonset={"AUU","AUA","AUC","CUA","CUC","CUG","CUU","UUA","UUG","GUU","GUA","GUC","GUG","UUU","UUC","AUG","UGU","UGC","GCA","GCC","GCG","GCU","GGU","GGC","GGA","GGG","CCU","CCC","CCA","CCG","ACU","ACC","ACA","ACG","UCU","UCC","UCA","UCG","AGU","AGC","UAU","UAC","UGG","CAA","CAG","AAU","AAC","CAU","CAC","GAA","GAG","GAU","GAC","AAA","AAG","CGU","CGC","CGA","CGG","AGA","AGG","UAA","UAG","UGA"};
std::unordered_map<int, float> codonset_CAI;

#define BASE(x) ((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: (x=='u'?5:(x=='g'?6:0)))))))

inline int cal_CAI_index(int base_0 , int base_1, int base_2){
    return base_0*100 + base_1*10 + base_2;
}

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

inline int cal_base_index(int i, int base_i, int base_j){
    return i*100+base_i*10+base_j;
}

inline float cal_CAI(string result_seq){
	int size = result_seq.size()/3;
    float CAI=0;
    int cai_index;
    for(int i = 0; i < size; i++){
        cai_index = cal_base_index(BASE(result_seq[3*i]),BASE(result_seq[3*i+1]),BASE(result_seq[3*i+2]));
        CAI = CAI + codonset_CAI[cai_index];
    }
    return pow(10.0, CAI/size);
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

inline void log_CAI(vector<float>& cai_vector, vector<float>& log_cai_vector){

    for(int i = 0; i < cai_vector.size(); i++){
        log_cai_vector.push_back(log10f(cai_vector[i]));
    }
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}


int main(int argc, char** argv){

	ifstream fasta_file;
    ifstream cai_file;

    std::vector<float> cai_vector;
    std::vector<float> log_cai_vector;

    string rna_seq, seq="";

	InputParser input(argc, argv);

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

    const std::string &file = input.getCmdOption("-f");
    if (!file.empty()){
        fasta_file.open(file);

        while (fasta_file >> seq){
        	rtrim(seq);
        	rna_seq += seq;
        	break;
            // cout << seq << endl;
            // if (seq.empty()) continue;
            // else if (seq[0] == '>' or seq[0] == ';'){

            //     continue;
            // }else{
            //     rtrim(seq);
            //     rna_seq += seq;
            // }
        }
    }

    cout << "RNA: " << rna_seq << endl;
    cout << "log(CAI): " << cal_CAI(rna_seq) << endl;

}
