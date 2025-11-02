# LinearCDSfold (version 1.5)

LinearCDSfold is a tool for designing coding sequences by jointly optimizing their secondary structure stability and codon usage preferences.

## To Compile:
Run the following command to compile LinearCDSfold.

```
make
```

## To Run:
To run **LinearCDSfold**, use the following command:

```
./LinearCDSfold [OPTIONS] <SEQUENCE_FILE>
```

### ARGUMENTS: 

`SEQUENCE_FILE` is an amino acid sequence file in FASTA format. The following is an example of `SEQUENCE_FILE`.

```
> P15421.fasta
MYGKIIFVLLLSGIVSISASSTTGVAMHTSTSSSVTKSYISSQTNGITLINWWAMARVIFEVMLVVVGMIILISYCIR
```

### OPTIONS:

```
-c <CODON_USAGE_FILE> 
```

`CODON_USAGE_FILE` is a CSV file containing codon usage frequencies, where each row specifies a codon triplet, the amino acid it encodes, and its relative usage frequency among synonymous codons. 

Format of `CODON_USAGE_FILE`: `[codon triplet], [amino acid], [relative frequency]`

The following is an example of the content in `CODON_USAGE_FILE`. If the first line of `CODON_USAGE_FILE` starts with `#`, it is treated as a comment.

```
# triplet, amino acid, frequency
GCU,A,0.26
GCC,A,0.4
GCA,A,0.23
GCG,A,0.11
UGU,C,0.45
UGC,C,0.55
GAU,D,0.46
GAC,D,0.54
...
```

`CODON_USAGE_FILE` defaults to `codon_usage_freq_table_human.csv`.

```
-O <OBJECTIVE_FUNCTION>
```

The `-O` option (uppercase letter O) specifies the objective function for joint optimization of MFE (Minimum Free Energy) and CAI (Codon Adaptation Index). The `OBJECTIVE_FUNCTION` parameter must be set to either `LD` (default) or `DN`.

When `OBJECTIVE_FUNCTION` is `LD`, the objective function defined by LinearDesign is applied:

MFE − `LAMBDA` × _l_ × log(CAI),

where _l_ is the length of the input amino acid sequence, and `LAMBDA` is a real value ranging from 0 to ∞.

Conversely, if `OBJECTIVE_FUNCTION` is `DN`, the objective function defined by DERNA is utilized:

`LAMBDA` × MFE − (1 − `LAMBDA`) × _l_ × log(CAI),

where `LAMBDA` is a real value ranging from 0 to 1.

```
-l <LAMBDA>
```

`LAMBDA` is a non-negative real-valued scaling parameter used to balance the contributions of the MFE of a coding sequence and its CAI value in the joint optimization objective.

When the objective function defined by LinearDesign is applied (i.e., `OBJECTIVE_FUNCTION` is set to `LD`), setting `LAMBDA` to `0` restricts the optimization to consider only MFE. Otherwise, both MFE and CAI are included, with higher values of `LAMBDA` increasing the weight assigned to CAI. If the LinearDesign objective function is used, the default value of `LAMBDA` is `0`.

Conversely, when the objective function defined by DERNA is utilized (i.e., `OBJECTIVE_FUNCTION` is set to `DN`), setting `LAMBDA` to `1` restricts the optimization to consider only MFE. Otherwise, both MFE and CAI are included, with larger values of `LAMBDA` increasing the contribution of MFE. If the DERNA objective function is used, the default value of `LAMBDA` is `1`.

```
-m <SEARCH_MODE> or --mode <SEARCH_MODE>
```

`SEARCH_MODE` is a string that selects the search strategy in LinearCDSfold. Available search modes are:

  - `exact` (default): Performs an exact search to obtain an optimal CDS.
  - `beam`: Uses beam search to quickly generate an approximate CDS.
  - `pareto`: Generates a set of Pareto-optimal CDSs using exact search and the DERNA objective function.

If no mode is specified, LinearCDSfold defaults to `exact` search mode.

```
-b <BEAM_SIZE>
```

`BEAM_SIZE` is a positive integer that specifies the beam size used by LinearCDSfold during beam search (i.e., `-m beam` is specified) and its default value is `500`.

```
-t <TAU1> or --tau1 <TAU1>
```

`TAU1` is a termination threshold used by LinearCDSfold when Pareto-optimal search is enabled (i.e., when `-m pareto` is specified). Its default value is `0.0025`. In general, smaller values of `TAU1` result in more Pareto-optimal CDSs being generated, at the cost of longer runtime.

```
-u <TAU2> or --tau2 <TAU2>
```

`TAU2` is another termination threshold used by LinearCDSfold when Pareto-optimal search is enabled (i.e., when `-m pareto` is specified). Essentially, `TAU2` is used to explore Pareto-optimal CDSs that are generated using `LAMBDA` values smaller than `TAU1`. Therefore, the value of `TAU2` should be smaller than that of `TAU1`. By default, it is set to `0.00075`. Smaller values of `TAU2` allow more Pareto-optimal CDSs to be generated, but increase runtime.

```
-o <FILE_NAME>
```

The `-o` option (lowercase letter o) specifies the name of the output file, which will contain detailed results from LinearCDSfold in plain text format. The file `FILE_NAME` should use the `.txt` extension. If `-o` is not specified, the default output file is `result.txt`.

```
-f <FILE_NAME>
```

The `-f` option specifies the name of an additional output file that contains only the MFE and CAI results returned by LinearCDSfold for each `LAMBDA` value in CSV format. The `FILE_NAME` file should use the `.csv` extension. If `-f` is not specified, the default output file is `result.csv`.

## Examples:

### Exact search using LinearDesign objective function

```
> ./LinearCDSfold -l 2 -o P15421.txt example/P15421.fasta
```

Output: `cat P15421.txt`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: LinearDesign
Search mode: Exact search
Lambda: 2.000
Processing: [==================================================]  100%
Coding sequence and its secondary structure:
AUGUAUGGCAAGAUCAUCUUUGUGCUGCUGCUGAGCGGGAUCGUGUCGAUCUCCGCCAGCAGCACCACAGGGGUGGCCAUGCAUACCAGCACCAGCAGUAGCGUGACCAAGAGCUACAUCUCCAGCCAGACCAACGGCAUCACCUUGAUCAACUGGUGGGCCAUGGCCCGCGUGAUCUUCGAGGUGAUGCUGGUGGUGGUGGGGAUGAUCAUCUUGAUCAGCUACUGCAUCCGC
(((((((((.....((((((((((.((((((((.(((((((((...))))).)))))))))))).)))))))))))))))))))....((....((((((((.(((.(((((....((((((((.((..((((.((((((((((((((.......((((((....)))))).......)))))))))))))).))))))))))))))....))))).)))))))))))....))
Folding free energy: -132.600 kcal/mol
CAI: 0.919
Total runtime: 3.296 s
```

### Exact search using DERNA objective function

```
> ./LinearCDSfold -O DN -l 2 -o P15421_DN.txt example/P15421.fasta
```

Output: `cat P15421_DN.txt`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: DERNA
Search mode: Exact search
Lambda: 2.000
Processing: [==================================================]  100%
Coding sequence and its secondary structure:
AUGUAUGGCAAAAUCAUCUUCGUCUUGCUGCUCUCCGGGAUCGUAUCGAUCUCGGCGAGCAGCACGACGGGGGUGGCCAUGCAUACGAGUACGAGCAGUAGCGUGACUAAGAGUUAUAUAUCCUCACAGACCAACGGCAUCACCUUGAUAAAUUGGUGGGCGAUGGCCCGCGUAAUUUUCGAGGUGAUGCUGGUGGUCGUGGGGAUGAUAAUUCUUAUCAGCUACUGCAUUCGU
(((((((((.....((((((((((.((((((((.(((((((((...))))))))).)))))))).)))))))))))))))))))((((((....((((((((.(((.((((((((((.(((((((((.(((((.((((((((((((((.((((((((((((....)))))).)))))))))))))))))))).)))))))))))))))))))))))))))))))))))))))))
Folding free energy: -148.700 kcal/mol
CAI: 0.697
Total runtime: 3.310 s
```

### Beam search using LinearDesign objective function

```
> ./LinearCDSfold -m beam -b 100 -l 2 -o P15421_beam.txt example/P15421.fasta
```

Output: `cat P15421_beam.txt`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: LinearDesign
Search mode: Beam search
Beam size: 100
Lambda: 2.000
Processing: [==================================================]  100%
Coding sequence and its secondary structure:
AUGUAUGGCAAGAUCAUCUUUGUGCUGCUGCUGAGCGGAAUUGUGAGCAUUUCCGCCAGCAGCACCACAGGGGUGGCCAUGCAUACCAGCACCAGCAGUAGCGUGACCAAGAGCUACAUCUCCAGCCAGACCAACGGCAUCACCUUGAUCAAUUGGUGGGCCAUGGCCCGCGUGAUUUUCGAGGUGAUGCUGGUGGUGGUGGGGAUGAUCAUCUUGAUCAGCUACUGCAUCCGC
(((((((((.....((((((((((.((((((((.((((((.((....)).)))))))))))))).)))))))))))))))))))....((....((((((((.(((.(((((....((((((((.((..((((.((((((((((((((..((((.((((((....))))))..)))).)))))))))))))).))))))))))))))....))))).)))))))))))....))
Folding free energy: -130.300 kcal/mol
CAI: 0.924
Total runtime: 0.258 s
```

### Beam search using DERNA objective function

```
> ./LinearCDSfold -m beam -b 100 -O DN -l 0.001 -o P15421_beam_DN.txt example/P15421.fasta
```

Output: `cat P15421_beam_DN.txt`

```
Amino acid file: example/P15421.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: DERNA
Search mode: Beam search
Beam size: 100
Lambda: 0.001
Processing: [==================================================]  100%
Coding sequence and its secondary structure:
AUGUACGGCAAGAUCAUCUUCGUGCUGCUGCUGAGCGGCAUCGUGUCCAUCAGCGCCAGCAGCACCACCGGCGUGGCCAUGCACACCUCCACCAGCAGCAGCGUGACCAAGAGCUACAUCAGCUCUCAGACCAAUGGCAUCACCCUGAUCAACUGGUGGGCCAUGGCCAGGGUGAUCUUCGAGGUGAUGCUGGUGGUGGUGGGCAUGAUCAUCCUGAUCAGCUACUGCAUCAGG
..((((((..(((...)))))))))((((((...))))))((((((((((((.((((((((.((((.(((.(.((((((((......((((((((...(((.(((((((((((((.....))))))........)))..)))).))).....)))))))).)))))))).).))......).)))).)))))))).))))))))))))....(((((((((...))).))))))
Folding free energy: -103.100 kcal/mol
CAI: 0.991
Total runtime: 0.263 s
```

### Pareto-optimal search using default termination thresholds

```
./LinearCDSfold -m pareto -o P15421_Pareto.txt -f P15421_Pareto.csv example/P15421.fasta
```

Output: `cat P15421_Pareto.csv`

```
lambda,MFE,CAI
0.000010,-88.000000,1.000000
0.999990,-148.700000,0.733673
0.500000,-148.700000,0.733673
0.250005,-148.700000,0.733673
0.125007,-148.700000,0.733673
0.062509,-148.300000,0.760887
0.031259,-148.300000,0.760887
0.093758,-148.700000,0.733673
0.015635,-146.100000,0.798641
0.078133,-148.700000,0.733673
0.007822,-137.700000,0.882926
0.023447,-148.300000,0.760887
0.070321,-148.700000,0.733673
0.003916,-128.400000,0.939532
0.011729,-146.100000,0.798641
0.019541,-148.300000,0.760887
0.066415,-148.700000,0.733673
0.001963,-124.300000,0.956380
0.005869,-134.700000,0.905781
0.009775,-143.400000,0.829243
0.017588,-147.500000,0.775186
0.064462,-148.300000,0.760887
0.000987,-103.100000,0.990816
0.000498,-97.200000,0.996524
0.001475,-118.800000,0.968451
```

## Contact Information:

Corresponding author: Prof. Chin Lung Lu (Email: cllu@cs.nthu.edu.tw)
