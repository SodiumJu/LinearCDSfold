# LinearCDSfold (version 1.5)

LinearCDSfold is a tool for designing a coding sequence (CDS) by jointly optimizing its secondary structure stability and codon usage.

To see LinearCDSfold (version 1.0), please go to
[LinearCDSfold v1](https://github.com/ablab-nthu/LinearCDSfold/tree/1dc4698e51576ecca6648d3f1d0449c3381a72da)
or [download this version as ZIP](https://github.com/ablab-nthu/LinearCDSfold/archive/1dc4698e51576ecca6648d3f1d0449c3381a72da.zip)


## To Compile:
Run the following command to compile LinearCDSfold.

```
make
```

## To Run:
To run **LinearCDSfold**, use the following command:

```bash
./LinearCDSfold [OPTIONS] SEQFILE
```

Example:

```
> ./LinearCDSfold -b 100 -l 2 -txt P15421_beam.txt -csv P15421_beam.csv example/P15421.fasta
```
Output:
```
Amino acid file: example/P15421.fasta
Codon usage table: cai_example/codon_usage_freq_table_human_LinearDesign.csv
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

### SEQFILE:

`SEQFILE` is an amino acid sequence file in FASTA format. The following is an example of `SEQFILE`.

```
> P15421.fasta
MYGKIIFVLLLSGIVSISASSTTGVAMHTSTSSSVTKSYISSQTNGITLINWWAMARVIFEVMLVVVGMIILISYCIR
```

### OPTIONS:

```
-b <BEAMSIZE>
```

`BEAMSIZE` is a non-negative integer parameter that specifies the beam size used by LinearCDSfold when running with beam pruning. By default, it is set to `0`, meaning that LinearCDSfold performs exact search rather than beam search.

```
-cai <CAIFILE>
```

`CAIFILE` is a CSV file containing codon usage frequencies.
Each entry lists a codon triplet, its corresponding amino acid, and the frequency value (as a fraction).

Example files:
- `cai_example/codon_usage_freq_table_human_LinearDesign.csv`
- `cai_example/codon_usage_freq_table_yeast_LinearDesign.csv`.

File format: `[triplet], [amino acid], [frequency fraction]`

Example:
```
#,,
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

**Note:** The default value of `CAIFILE` is set to `cai_example/codon_usage_freq_table_human_LinearDesign.csv`.

```
-o <OBJECTIVE>
```
`OBJECTIVE` represents the joint objective function for optimizing both MFE (Minimum Free Energy) and CAI (Codon Adaptation Index). It can be set to either `LD` or `DN`.

When `OBJECTIVE` is `LD`, the objective function defined by LinearDesign is applied:

MFE − `LAMBDA` × _l_ × log(CAI),

where _l_ is the length of the input amino acid sequence, and `LAMBDA` is a real value ranging from 0 to ∞.

Conversely, if `OBJECTIVE` is `DN`, the objective function defined by DERNA is utilized:

`LAMBDA` × MFE − (1 − `LAMBDA`) × _l_ × log(CAI),

where `LAMBDA` is a real value ranging from 0 to 1.

The default value of `OBJECTIVE` is `LD`.

```
-l <LAMBDA>
```

`LAMBDA` is a non-negative real-valued scaling parameter used to balance the contributions of the MFE of a coding sequence and its CAI value in the joint optimization objective.

When the objective function defined by LinearDesign is applied (i.e., `OBJECTIVE` is set to `LD`), setting `LAMBDA` to `0` restricts the optimization to consider only MFE. Otherwise, both MFE and CAI are included, with higher values of `LAMBDA` increasing the weight assigned to CAI. If the LinearDesign objective function is used, the default value of `LAMBDA` is `0`.

Conversely, when the objective function defined by DERNA is utilized (i.e., `OBJECTIVE` is set to `DN`), setting `LAMBDA` to `1` restricts the optimization to consider only MFE. Otherwise, both MFE and CAI are included, with larger values of `LAMBDA` increasing the contribution of MFE. If the DERNA objective function is used, the default value of `LAMBDA` is `1`.

```
-p 
```

`p` is a flag that indicates whether Pareto-optimal search is enabled. When Pareto-optimal search is enabled (i.e., the `-p` option is used), LinearCDSfold automatically generates a set of Pareto-optimal CDSs using the DERNA objective function (instead of the LinearDesign objective function).

```
-t1 <tau1>
```

`tau1` is a termination threshold used by LinearCDSfold when Pareto-optimal search is enabled (i.e., the `-p` option is specified). Its default value is `0.0025`. In principle, the smaller the value of `tau1`, the more Pareto-optimal CDSs can be generated, and the longer the required runtime.

```
-t2 <tau2>
```

`tau2` is another termination threshold used by LinearCDSfold when Pareto-optimal search is enabled (i.e., the `-p` option is specified). Essentially, `tau2` is used to explore Pareto-optimal CDSs that are generated using `LAMBDA` values smaller than `tau1`. Therefore, the value of `tau2` should be smaller than that of `tau1`. By default, it is set to `0.00075`. The smaller the value of `tau2`, the more Pareto-optimal CDSs can be generated, and the longer the required runtime.

```
-txt <FILENAME.txt>
```

`FILENAME.txt` specifies the output file that contains the detailed results from LinearCDSfold in plain text format. By default, it is set to `result.txt`.

```
-csv <FILENAME.csv>
```

`FILENAME.csv` specifies the output file that contains only the MFE and CAI results returned by LinearCDSfold for each `LAMBDA` value in CSV format. By default, it is set to `result.csv`.

## Examples:

### Exact search:

```
> ./LinearCDSfold -l 2 -txt P15421.txt -csv P15421.csv example/P15421.fasta
```
Output:
```
Amino acid file: example/P15421.fasta
Codon usage table: cai_example/codon_usage_freq_table_human_LinearDesign.csv
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
> ./LinearCDSfold -b 100 -l 2 -txt P15421_beam.txt -csv P15421_beam.csv example/P15421.fasta
```
Output:
```
Amino acid file: example/P15421.fasta
Codon usage table: cai_example/codon_usage_freq_table_human_LinearDesign.csv
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
> ./LinearCDSfold -o DN -b 100 -l 0.001 -txt P15421_beam_DN.txt -csv P15421_beam_DN.csv example/P15421.fasta
```
Output:
```
Amino acid file: example/P15421.fasta
Codon usage table: cai_example/codon_usage_freq_table_human_LinearDesign.csv
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
./LinearCDSfold -o DN -p -txt P15421_p.txt -csv P15421_p.csv example/P15421.fasta
```

### Pareto-optimal search using modified termination thresholds

```
./LinearCDSfold -o DN -t1 0.002 -t2 0.0007 -p -txt P15421_p.txt -csv P15421_p.csv example/P15421.fasta
```
## Contact Information:

Corresponding author: Prof. Chin Lung Lu (Email: cllu@cs.nthu.edu.tw)
