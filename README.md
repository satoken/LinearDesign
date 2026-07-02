<img src="/pic/baidu_research_logo.jpg"  width="40%" alt="Baidu Research Logo">

#  Algorithm for Optimized mRNA Design Improves Stability and Immunogenicity (LinearDesign)
![GitHub all releases](https://img.shields.io/github/downloads/LinearDesignSoftware/LinearDesign/total)


This repository contains the source code for the LinearDesign project.

He Zhang†, Liang Zhang†, Ang Lin†, Congcong Xu†, Ziyu Li, Kaibo Liu, Boxiang Liu, Xiaopin Ma, Fanfan Zhao, Huiling Jiang, Chunxiu Chen, Haifa Shen, Hangwen Li*, David H. Mathews*, Yujian Zhang*, Liang Huang†*<sup>#</sup>. Algorithm for Optimized mRNA Design Improves Stability and Immunogenicity. Nature [https://doi.org/10.1038/s41586-023-06127-z](https://doi.org/10.1038/s41586-023-06127-z) (2023)

† contributed equally, 
\* corresponding authors, 
<sup>#</sup> lead corresponding author

For questions, please contact the lead corresponding author at <liang.huang.sh@gmail.com>.

## Dependencies
Clang 11.0.0 (or above) or GCC 4.8.5 (or above)

python3

## To Compile
```
make
```

## To Run
The LinearDesign program can be run with:
```
echo SEQUENCE | ./lineardesign [OPTIONS]

OR

cat FASTA_FILE | ./lineardesign [OPTIONS]
```

OPTIONS:
```
--lambda LAMBDA or -l LAMBDA
```
Set LAMBDA, a hyperparameter balancing MFE and CAI. (default 0.0)
```
--codonusage FILE_NAME or -c FILE_NAME
```
Import a Codon Usage Frequency Table. See "codon_usage_freq_table_human.csv" for the format.
(default: using human codon usage frequency table)
```
--fixedprefix RNA_PREFIX or -p RNA_PREFIX
```
Fix the coding RNA prefix and design the RNA sequence for the input amino acids after it.
The fixed prefix must contain only A/C/G/U (T is accepted and converted to U), must be a
multiple of 3 nucleotides, and must not contain a stop codon.
```
--fixedutr RNA_PREFIX or -u RNA_PREFIX
```
Fix a 5' UTR RNA prefix before the coding sequence. The fixed UTR participates in
folding free-energy optimization, but it is excluded from translation and CAI scoring.
The fixed UTR may have any length and must contain only A/C/G/U (T is accepted and
converted to U).
```
--avoidpairrange START-END or -r START-END
--avoidpairpenalty KCAL
```
Add a pseudo-energy penalty for every base pair involving at least one nucleotide in
the 1-based inclusive range START-END. The range is counted on the final output RNA,
including any fixed 5' UTR. This can be used to discourage secondary structure around
the start codon. When this option is active, the reported folding free energy is the
pseudo-adjusted objective value.
```
--verbose or -v
```
Print out more details. (default False)

For Macbook, users may encounter a pop-up message at the first run.
For Mac-M1 system, the message is:
```
"LinearDesign_Mac_M1.so" can't be opened because Apple cannot check it for malicious software.
```
For Mac-Intel system, the message is:
```
"LinearDesign_Mac_Intel.so" cannot be opened because it is from an unidentified developer.
```
If so, please go to "System Preferences -> Security & Privacy -> General" to allow LinearDesign-Mac-M1.so (or LinearDesign-Mac-Intel.so) to open.

## Example: Single Sequence Design
```
echo MNDTEAI | ./lineardesign
mRNA sequence:  AUGAACGAUACGGAGGCGAUC
mRNA structure: ......(((.((....)))))
mRNA folding free energy: -1.10 kcal/mol; mRNA CAI: 0.695
```

## Example: Multiple Sequences Design with Option --lambda (-l)
```
cat testseq | ./lineardesign --lambda 3
>seq1
mRNA sequence:  AUGCCAAACACCCUGGCAUGCCCC
mRNA structure: ((((((.......)))))).....
mRNA folding free energy: -6.00 kcal/mol; mRNA CAI: 0.910

>seq2
mRNA sequence:  AUGCUGGAUCAGGUGAACAAGCUGAAGUACCCAGAGGUGAGCCUGACCUGA
mRNA structure: .....((.((((((..((...(((.......)))..))..))))))))...
mRNA folding free energy: -13.50 kcal/mol; mRNA CAI: 0.979
```

## Example: Option --codonusage (-c)
```
echo MNDTEAI | ./lineardesign -l 0.3 --codonusage codon_usage_freq_table_yeast.csv
mRNA sequence:  AUGAAUGAUACGGAAGCGAUC
mRNA structure: ......(((.((....)))))
mRNA folding free energy: -1.10 kcal/mol; mRNA CAI: 0.670
```

## Example: Fixed RNA Prefix
```
echo DTEAI | ./lineardesign --fixedprefix AUGAAC
mRNA sequence:  AUGAACGAUACGGAGGCGAUC
mRNA structure: ......(((.((....)))))
mRNA folding free energy: -1.10 kcal/mol; mRNA CAI: 0.695
```

## Example: Fixed 5' UTR
```
echo MNDTEAI | ./lineardesign --fixedutr GGAA
mRNA sequence:  GGAAAUGAACGACACCGAGGCCAUU
mRNA structure: ...((((..(........)..))))
mRNA folding free energy: -1.20 kcal/mol; mRNA CAI: 0.960
```

## Example: Discourage Base Pairing in a Range
```
echo MNDTEAI | ./lineardesign --avoidpairrange 7-20 --avoidpairpenalty 10
mRNA sequence:  AUGAACGACACGGAGGCGAUC
mRNA structure: .....................
mRNA pseudo-adjusted folding free energy: -0.00 kcal/mol; mRNA CAI: 0.711
```

## Example: Option --verbose (-v)
```
echo MNDTEAI | ./lineardesign --verbose
Input protein: MNDTEAI
Using lambda = 0; Using codon frequency table = codon_usage_freq_table_human.csv
mRNA sequence:  AUGAACGAUACGGAGGCGAUC
mRNA structure: ......(((.((....)))))
mRNA folding free energy: -1.10 kcal/mol; mRNA CAI: 0.695
Runtime: 0.002 seconds
```


## Declarations
Baidu Research has filed a patent for the LinearDesign algorithm that lists He Zhang, Liang Zhang, Ziyu Li, Kaibo Liu, Boxiang Liu, and Liang Huang as inventors.
