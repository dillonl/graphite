Graphite - Graph-based adjudication 
========================================

Overview
========================================
Graphite is a variant adjudication tool that aids in confirming or rejecting the presence of candidate alleles by constructing a graph using reference as well as alternate alleles.

Output
========================================
Graphite annotates existing VCFs by appending allele "counts" to the VCF that represents the number of reads supporting the reference/alternate alleles.

Installation
========================================
Building and installing Graphite is made simple with [CMAKE](https://cmake.org/). Additionally, Graphite takes advantage of the new and improved C11 C++ compiler. If you are using the GCC compiler please use 4.9.2 or greater.

To build Graphite, use the following commands (from the base Graphite directory):

```Shell
mkdir bin/
cd bin
cmake ../
make
make install #optional
```

Usage
========================================
Graphite takes in VCF(s), BAM file(s) and the FASTA file used to align the BAM(s). A variant graph representation is generated based on the reference and VCF variants. The reads from each sample within the BAM file(s) are then remapped to variant regions of the graph.

Graphite uses a modified Smith-Waterman algorithm for the read mapping. This read mapping sensitivity can be adjusted by modifying the match value, mismatch penalty, gap open penalty, gap extension penalty and Smith-Waterman percent match.

Use help for a complete list of Graphite's parameters and usage:
```Shell
./graphite -h
```
License
========================================
`Graphite` is freely available under the MIT [license](https://opensource.org/licenses/MIT).

