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

#include "LCDSfold.h"
using namespace std;

int main(int argc, char** argv){
    std::string seq="";
    std::string objective ;
    std::ifstream fasta_file;
    std::ifstream cai_file;

    std::string rna_seq, ami_seq;
    std::vector<std::string> rna_seq_list, inseq_list;
    std::vector<std::vector<int>> con_seq_list;
    std::vector<double> cai_vector;

    bool is_rna_file = false;
    bool show_score = false;
    initialize();
    beamsize = 0;
    lambda = 0;
    pareto = false;
    objective = "LD";
    double threshold1 = 0.0025;
    double threshold2 = 0.00075;
    

    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--h") || input.cmdOptionExists("-help") || input.cmdOptionExists("--help")){
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
    
    /*Objective function*/
    std::string objective_ = input.getCmdOption("-O");
    if (!objective_.empty()) {
        objective = objective_;
        if(objective == "DN") lambda = 1;//change default lambda 

        if (objective != "LD" && objective != "DN") {
            std::cerr << "Error: objective is not LD or DN." << std::endl;
            return 0;
        }

    }

    // Read sequence file 
    if (!file.empty()){
        bool start_sign_ = false;
        fasta_file.open(file);
        if (fasta_file.is_open()){
            while (std::getline(fasta_file, seq)){
                if (seq.empty()){
                    start_sign_ = false;
                    continue;
                }else if (seq[0] == '>' or seq[0] == ';'){
                    start_sign_ = true;
                    if (!ami_seq.empty())
                        inseq_list.push_back(ami_seq);
                    ami_seq.clear();
                    continue;
                }else if(start_sign_){
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
        for (seq; getline(cin, seq);){
            if (seq.empty()) continue;
            if (!isalpha(seq[0])){
                std::cerr << "Unrecognized sequence: " << seq << std::endl;
                continue;}
            inseq_list.push_back(seq);
        }
    }

    /*switch mode*/
    bool use_beam_pruning = false;
    std::string mode_ = input.getCmdOption("-m");
    if (mode_.empty()) {
        mode_ = input.getCmdOption("--mode");
    }
    if (mode_.empty() || mode_ == "exact" || mode_ == "Exact") {
        beamsize = 0;
    } else if (mode_ == "beam" || mode_ == "Beam") {
        use_beam_pruning = true;
        beamsize = 500; // default beam size
        // beamsize will be set later based on -b option
    } else if (mode_ == "pareto" || mode_ == "Pareto") {
        pareto = true;
    } else {
        std::cerr << "Error: SEARCH_MODE " << mode_ << " is not recognized." << std::endl;
        return 0;
    }

    /*beam search*/
    std::string beamsize_ = input.getCmdOption("-b");
    if (!beamsize_.empty() && use_beam_pruning)
        beamsize = stoi(beamsize_);


    /*lambda*/
    std::string lambda_ = input.getCmdOption("-l");
    if (!lambda_.empty())  
        lambda = stof(lambda_);

    /*threshold*/
    std::string threshold1_ = input.getCmdOption("-t");
    if (threshold1_ == "") {
        threshold1_ = input.getCmdOption("--tau1");
    }
    if (!threshold1_.empty())
        threshold1 = stof(threshold1_);
    std::string threshold2_ = input.getCmdOption("-u");
    if (threshold2_ == "") {
        threshold2_ = input.getCmdOption("--tau2");
    }
    if (!threshold2_.empty())  
        threshold2 = stof(threshold2_);
    
    /*cai file*/
    std::string cai_file_path = "codon_usage_freq_table_human.csv";//default cai file (human)
    std::string cai_file_str = input.getCmdOption("-c");
    if (!cai_file_str.empty())
        cai_file_path = cai_file_str;

    if (isCSV(cai_file_path)){
        // std::cout << "[Debug] Reading CAI table from CSV file: " << cai_file_path << std::endl;
        read_CAI_table_csv(cai_file_path, cai_vector);
    } else {
        // in here, the input cai txt file must be relative cai values for codons in CodonSet order
        cai_file.open(cai_file_path);
        if (cai_file.is_open()){
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
        }else{
            if(!cai_file_str.empty()){
                std::cout << "[Error] Cannot open <CAIFILE>: " << cai_file_path << std::endl;
                return 0;
            }
        }

    }

    //Reading output file name
    std::string output_txt;
    std::string output_csv;
    if(input.cmdOptionExists("-o"))
        output_txt = input.getCmdOption("-o");//edited
    else output_txt = "result.txt";

    if(input.cmdOptionExists("-f"))
        output_csv = input.getCmdOption("-f");//edited
    else output_csv = "result.csv";

    if(pareto){
        objective = "DN";
        beamsize = 0;
    }
        
    PrintInfo(output_txt,output_csv,file,beamsize, cai_file_path, lambda, objective, pareto, is_rna_file,true);
    if(!is_rna_file){ // NORMAL MODE
                
        if((objective == "LD")){//LinearDesign mode
            for(int i = 0; i < inseq_list.size(); i++){
                timeval start,end;
                int sec,usec;
                ami_seq = inseq_list[i];
                transform(ami_seq.begin(), ami_seq.end(), ami_seq.begin(), ::toupper);
                ami_to_rna(rna_seq_list, ami_seq);
                add_con_seq(con_seq_list, ami_seq);
                rna_seq = rna_seq_list[i];
                std::vector<int> con_seq = con_seq_list[i];

                std::cout <<"Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;
                AllTables<double> alltables(rna_seq, rna_seq.size());
                
                initialize_CAI_table(cai_vector,false);
                gettimeofday(&start, 0);
                initialize_Special_HP_LD<double>(alltables,rna_seq, con_seq, ami_seq);
                if(beamsize)//beam search
                    LCDSfoldCAI_LD_beam<double>(alltables, rna_seq, con_seq, ami_seq);
                else
                    LCDSfoldCAI_LD_exact<double>(alltables, rna_seq, con_seq, ami_seq);

                gettimeofday(&end, 0);
                double TimeSpend = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);

                std::vector<pair<double,double>> result = result_output<double>(alltables, rna_seq, con_seq, ami_seq,output_txt,output_csv,TimeSpend, show_score, false);
            }
        }
        else{//DERNA mode
            for(int i = 0; i < inseq_list.size(); i++){
                timeval start,end;
                int sec,usec;
                ami_seq = inseq_list[i];
                transform(ami_seq.begin(), ami_seq.end(), ami_seq.begin(), ::toupper);
                ami_to_rna(rna_seq_list, ami_seq);
                add_con_seq(con_seq_list, ami_seq);

                
                rna_seq = rna_seq_list[i];
                std::vector<int> con_seq = con_seq_list[i];
                if(pareto)
                    Pareto_solution<double>(threshold1,threshold2,rna_seq, con_seq, ami_seq, cai_vector, output_txt, output_csv, show_score);

                else{
                    std::cout <<"Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;
                    if(lambda == 0 )//prevent MFE not predictable
                        lambda += 0.00001;
                    AllTables<double> alltables(rna_seq, rna_seq.size());
                    initialize_CAI_table(cai_vector,true);
                    gettimeofday(&start, 0);
                    initialize_Special_HP_DN<double>(alltables,rna_seq, con_seq, ami_seq);
                    if(beamsize)
                        LCDSfoldCAI_DN_beam<double>(alltables, rna_seq, con_seq, ami_seq);
                    else
                        LCDSfoldCAI_DN_exact<double>(alltables, rna_seq, con_seq, ami_seq);
                    gettimeofday(&end, 0);

                    double TimeSpend = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
                    std::vector<pair<double,double>> result = result_output<double>(alltables, rna_seq, con_seq, ami_seq,output_txt,output_csv,TimeSpend, show_score, true);
                    
                }
            }    
        }
    
    }else{ //RNA MODE
        lambda = 0;
        std::cout << "RNA file: " << file << std::endl;
        std::cout << "Beam size: " << beamsize << std::endl;
        
        for(int i = 0; i < inseq_list.size(); i++){
            timeval start,end;
            int sec,usec;
            rna_seq = inseq_list[i];
            std::vector<int> con_seq(rna_seq.size(), normal_ami);
            AllTables<double> alltables(rna_seq, rna_seq.size());
            gettimeofday(&start, 0);
            initialize_Special_HP_LD<double>(alltables,rna_seq, con_seq, ami_seq);
            LCDSfoldCAI_LD_exact<double>(alltables, rna_seq, con_seq, ami_seq);
            gettimeofday(&end, 0);
            
            double TimeSpend = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
            std::vector<pair<double,double>> result = result_output<double>(alltables, rna_seq, con_seq, ami_seq,output_txt,output_csv,TimeSpend, show_score, true);

        }

    }

    // cout<<"FINISH!"<<endl;
    
    return 0;
}
