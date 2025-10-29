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

### SEQUENCE_FILE:

`SEQUENCE_FILE` is an amino acid sequence file in FASTA format. The following is an example of `SEQUENCE_FILE`.

```
> P15421.fasta
MYGKIIFVLLLSGIVSISASSTTGVAMHTSTSSSVTKSYISSQTNGITLINWWAMARVIFEVMLVVVGMIILISYCIR
```

### OPTIONS:

```
-b <BEAM_SIZE>
```

`BEAM_SIZE` is a non-negative integer parameter that specifies the beam size used by LinearCDSfold when running with beam pruning. By default, it is set to `0`, meaning that LinearCDSfold performs exact search, rather than beam search.

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

**Note:** The default value of `CODON_USAGE_FILE` is set to `codon_usage_freq_table_human.csv`.

```
-O <OBJECTIVE_FUNCTION>
```

The `-O` (uppercase letter O) parameter specifies the objective function. `OBJECTIVE_FUNCTION` defines the objective function for jointly optimizing both MFE (Minimum Free Energy) and CAI (Codon Adaptation Index), and must be set to either `LD` or `DN`.

When `OBJECTIVE_FUNCTION` is `LD`, the objective function defined by LinearDesign is applied:

MFE − `LAMBDA` × _l_ × log(CAI),

where _l_ is the length of the input amino acid sequence, and `LAMBDA` is a real value ranging from 0 to ∞.

Conversely, if `OBJECTIVE_FUNCTION` is `DN`, the objective function defined by DERNA is utilized:

`LAMBDA` × MFE − (1 − `LAMBDA`) × _l_ × log(CAI),

where `LAMBDA` is a real value ranging from 0 to 1.

The default value of `OBJECTIVE_FUNCTION` is `LD`.

```
-l <LAMBDA>
```

`LAMBDA` is a non-negative real-valued scaling parameter used to balance the contributions of the MFE of a coding sequence and its CAI value in the joint optimization objective.

When the objective function defined by LinearDesign is applied (i.e., `OBJECTIVE_FUNCTION` is set to `LD`), setting `LAMBDA` to `0` restricts the optimization to consider only MFE. Otherwise, both MFE and CAI are included, with higher values of `LAMBDA` increasing the weight assigned to CAI. If the LinearDesign objective function is used, the default value of `LAMBDA` is `0`.

Conversely, when the objective function defined by DERNA is utilized (i.e., `OBJECTIVE_FUNCTION` is set to `DN`), setting `LAMBDA` to `1` restricts the optimization to consider only MFE. Otherwise, both MFE and CAI are included, with larger values of `LAMBDA` increasing the contribution of MFE. If the DERNA objective function is used, the default value of `LAMBDA` is `1`.

```
-P 
```

`P` (uppercase letter P) is a flag that indicates whether Pareto-optimal search is enabled. When Pareto-optimal search is enabled (i.e., the `-P` option is used), LinearCDSfold automatically generates a set of Pareto-optimal CDSs using the DERNA objective function (instead of the LinearDesign objective function).

```
-t <tau1> or --tau1 <TAU1>
```

`TAU1` is a termination threshold used by LinearCDSfold when Pareto-optimal search is enabled (i.e., the `-P` option is specified). Its default value is `0.0025`. In principle, the smaller the value of `TAU1`, the more Pareto-optimal CDSs can be generated, and the longer the required runtime.

```
-u <TAU2> or --tau2 <TAU2>
```

`TAU2` is another termination threshold used by LinearCDSfold when Pareto-optimal search is enabled (i.e., the `-P` option is specified). Essentially, `TAU2` is used to explore Pareto-optimal CDSs that are generated using `LAMBDA` values smaller than `TAU1`. Therefore, the value of `TAU2` should be smaller than that of `TAU1`. By default, it is set to `0.00075`. The smaller the value of `TAU2`, the more Pareto-optimal CDSs can be generated, and the longer the required runtime.

```
-o <FILE_NAME>
```

The `-o` (lowercase letter o) parameter specifies the output file name. The specified file `FILE_NAME` will contain the detailed results from LinearCDSfold in plain text format. If not specified, the default is `result.txt`.


```
-f <FILE_NAME>
```

The `-f` parameter specifies the output file name. The specified file `FILE_NAME` will contain only the MFE and CAI results returned by LinearCDSfold for each `LAMBDA` value in CSV format. If not specified, the default is `result.csv`.

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
Total runtime: 2.898 s
```

### Beam search using LinearDesign objective function

```
> ./LinearCDSfold -b 100 -l 2 -o P15421_beam.txt example/P15421.fasta
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
total runtime: 0.287 s
```

### Beam search using DERNA objective function

```
> ./LinearCDSfold -O DN -b 100 -l 0.001 -o P15421_beam_DN.txt example/P15421.fasta
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
total runtime: 0.294 s
```

### Pareto-optimal search using default termination thresholds

```
./LinearCDSfold -O DN -P -o P15421_Pareto.txt -f P15421_Pareto.csv example/P15421.fasta
```

### Pareto-optimal search using modified termination thresholds

```
./LinearCDSfold -O DN -P -t 0.0001 -o P15421_Pareto.txt -f P15421_Pareto.csv example/P15421.fasta
```

## Contact Information:

Corresponding author: Prof. Chin Lung Lu (Email: cllu@cs.nthu.edu.tw)
