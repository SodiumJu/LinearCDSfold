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
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <limits>
#include <cctype>
#include <chrono>
#include <ctime>
#include <math.h>


#include "LCDSfold.h"

using namespace std;

int main(int argc, char** argv){

    std::string seq="";
    std::ifstream fasta_file;
    std::ifstream cai_file;

    std::string rna_seq, ami_seq;
    std::vector<std::string> rna_seq_list, inseq_list;
    std::vector<std::vector<int>> con_seq_list;
    std::vector<float> cai_vector;

    bool is_rna_file = false;
    bool show_score = false;
    initialize();
    beamsize = 500;
    lambda = 3;

    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        PrintHelp();
        std::cout << std::endl;
        return 0;
    }else if(input.cmdOptionExists("-dev")){
        PrintHelp();
        std::cout
            << "  -rna:\n"
            << "         For checking if the calculation of RNA MFE is correct. Use this option if\n"
            << "         <SEQFILE> is RNA sequence fasta file.\n"
            << "  -score:\n"
            << "         Print the score of optimization.\n"
            << std::endl;
        return 0;
    }

    std::string file;
    if(argv[1][0] != '-' ){
        file = argv[1];
    }else{
        for (int i = 2; i < argc; ++i) { // Start from 1 to skip program name
            std::string arg_ = argv[i];
            std::string previous_arg_ = argv[i - 1];

            // Check if the argument does not start with "-"
            if (previous_arg_[0] != '-' && arg_[0] != '-') {
                file = arg_;
                break;
            }
            if ((previous_arg_ == "-rna" || previous_arg_ == "-score") && arg_[0] != '-') {
                file = arg_;
                break;
            }
        }
    }

    if (FindOption(argc, argv, "-rna")) {
        std::cout << "RNA MODE:  <SEQFILE> is RNA sequence fasta file." << std::endl;
        std::cout << "Validating MFE calculation of a RNA file." << std::endl;
        std::cout << std::endl;
        is_rna_file = true;
    }


    show_score = FindOption(argc, argv, "-score");

    // Read sequence file 
    if (!file.empty()) {
        bool start_sign_ = false;
        fasta_file.open(file);
        if (fasta_file.is_open()) {
            while (fasta_file >> seq) {
                if (seq.empty()) {
                    start_sign_ = false;
                    continue;
                }else if (seq[0] == '>' or seq[0] == ';'){
                    start_sign_ = true;
                    if (!ami_seq.empty())
                        inseq_list.push_back(ami_seq);
                    ami_seq.clear();
                    continue;
                }else if(start_sign_) {
                    rtrim(seq);
                    ami_seq += seq;
                }
            }
            if (!ami_seq.empty())
                inseq_list.push_back(ami_seq);
            fasta_file.close();
        }else{
            std::cerr << "Cannot open file:" << file << std::endl;
            return 0;
        }
    }else{
        for (seq; getline(cin, seq);) {
            if (seq.empty()) continue;
            if (!isalpha(seq[0])) {
                std::cerr << "Unrecognized sequence: " << seq << std::endl;
                continue;}
            inseq_list.push_back(seq);
        }
    }

    const std::string &beamsize_ = input.getCmdOption("-b");
    if (!beamsize_.empty())
        beamsize = stoi(beamsize_);

    std::string objective = "LD";
    const std::string &objective_ = input.getCmdOption("-o");
    if (!objective_.empty()) {
        objective = objective_;
        if (objective == "DN") {
            lambda = 1; // set default lambda = 1 when the objective is DN
        }
        if (objective != "LD" && objective != "DN") {
            std::cerr << "Error: objective is not LD or DN." << std::endl;
            return 0;
        }

    }
    
    const std::string &lambda_ = input.getCmdOption("-l");
    if (!lambda_.empty())  
        lambda = stof(lambda_);
    
    std::string cai_file_path = "yeast_relative_adaptiveness.txt";
    const std::string &cai_file_str = input.getCmdOption("-cai");
    if (!cai_file_str.empty())
        cai_file_path = cai_file_str;
    
    cai_file.open(cai_file_path);
    if (cai_file.is_open()) {
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
        if (objective == "DN" && lambda != 1) {
            initialize_CAI_table(cai_vector, true);
        } else if ((objective == "LD") || (objective == "DN" && lambda == 1)) {
            initialize_CAI_table(cai_vector);
        }
    }else{
        if(!cai_file_str.empty()) {
            std::cerr << "Codon usage table: cannot open <CAIFILE>: " << cai_file_path << std::endl;
            return 0;
        }
    }

    if(!is_rna_file) { // NORMAL MODE
        if ((objective == "LD") || (objective == "DN" && lambda == 1)) {
            // std out information
            PrintInfo(file, beamsize, cai_file_path, lambda, objective);

            if (objective == "DN" && lambda == 1) // if lambda == 1, then we are not maximizing CAI, we use DL, with lambda = 0 to maximize MFE
                lambda = 0;

            for(int i = 0; i < inseq_list.size(); i++) {
                ami_seq = inseq_list[i];
                transform(ami_seq.begin(), ami_seq.end(), ami_seq.begin(), ::toupper);
                ami_to_rna(rna_seq_list, ami_seq);
                add_con_seq(con_seq_list, ami_seq);

                rna_seq = rna_seq_list[i];
                std::vector<int> con_seq = con_seq_list[i];

                AllTables<int> alltables(rna_seq, rna_seq.size());
                LCDSfoldCAI<int>(alltables, rna_seq, con_seq, ami_seq);
                
                std::string rna_solution;
                std::string structure_solution;
                int maxscore = BackTrack<int>(rna_solution, structure_solution, rna_seq, alltables, con_seq);

                int weighted_cai_score = GetCAIScore(rna_solution);
                float cai_value = -1;
                if(lambda)
                    cai_value = expf(float(weighted_cai_score) / float(LDCONST * lambda * ami_seq.size()));
                else
                    cai_value = expf(float(weighted_cai_score) / float(LDCONST * ami_seq.size()));
                //output score and results
                std::cout << "Coding sequence and its secondary structure:" << std::endl;
                std::cout << rna_solution << std::endl;
                std::cout << structure_solution << std::endl;
                
                if (show_score)
                    std::cerr << "Score: " << maxscore << std::endl;

                if(!lambda)
                    std::cout << "Folding free energy: " << -(maxscore / 100.0) << " kcal/mol" << std::endl;
                else
                    std::cout << "Folding free energy: " << -float(maxscore - weighted_cai_score) / 100.0 << " kcal/mol" << std::endl;

                std::cout << "CAI: " << std::fixed << std::setprecision(3) << cai_value << std::endl;
                check_ami_solution(ami_seq, rna_solution);

            }
        } else { // DERNA MODE
            // std out information
            if (lambda <= 0 || lambda > 1){
                std::cerr << "Error: If objective is DN, lambda must be in (0, 1]." << std::endl;
                return 0;
            }
            PrintInfo(file, beamsize, cai_file_path, lambda, objective);

            for(int i = 0; i < inseq_list.size(); i++) {
                ami_seq = inseq_list[i];
                transform(ami_seq.begin(), ami_seq.end(), ami_seq.begin(), ::toupper);
                ami_to_rna(rna_seq_list, ami_seq);
                add_con_seq(con_seq_list, ami_seq);

                rna_seq = rna_seq_list[i];
                std::vector<int> con_seq = con_seq_list[i];

                AllTables<float> alltables(rna_seq, rna_seq.size());
                LCDSfoldCAI_derna<float>(alltables, rna_seq, con_seq, ami_seq);
                
                std::string rna_solution;
                std::string structure_solution;
                float maxscore = BackTrack<float>(rna_solution, structure_solution, rna_seq, alltables, con_seq);

                float weighted_cai_score = GetCAIScore_derna(rna_solution);
                float cai_value = -1;

                cai_value = expf(float(weighted_cai_score) / float((1.0 - lambda) * ami_seq.size()));

                //output score and results
                std::cout << "Coding sequence and its secondary structure:" << std::endl;
                std::cout << rna_solution << std::endl;
                std::cout << structure_solution << std::endl;
                
                if (show_score)
                    std::cerr << "Score: " << maxscore << std::endl;

                std::cout << "Folding free energy: " << -float(maxscore - weighted_cai_score) / (lambda * 100.0) << " kcal/mol" << std::endl;

                std::cout << "CAI: " << std::fixed << std::setprecision(3) << cai_value << std::endl;
                check_ami_solution(ami_seq, rna_solution);
            }
        }
    } else { //RNA MODE
        lambda = 0;
        std::cout << "RNA file: " << file << std::endl;
        std::cout << "Beam size: " << beamsize << std::endl;

        for(int i = 0; i < inseq_list.size(); i++) {
            rna_seq = inseq_list[i];
            std::vector<int> con_seq(rna_seq.size(), normal_ami);
            AllTables<int> alltables(rna_seq, rna_seq.size());
            LCDSfoldCAI(alltables, rna_seq, con_seq, ami_seq);
            
            std::string rna_solution;
            std::string structure_solution;
            int maxscore = BackTrack(rna_solution, structure_solution, rna_seq, alltables, con_seq);

            std::cout << "Coding sequence and its secondary structure:" << std::endl;
            std::cerr << rna_solution << std::endl;
            std::cerr << structure_solution << std::endl;
            std::cout << "Folding free energy: " << -(maxscore / 100.0) << " kcal/mol" << std::endl;
        }
    }
    
    return 0;
}
