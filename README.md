# Presentation

This is a suite of badly written but useful programs to compute population genetics statistics using sequences data in fasta or phylip format.
These programs are written using the [**Bio ++** library ](https://biopp.github.io/) (Guéguen et al. 2013).



This include:
- `seq_stat_2pop` to compute statistics on two populations

- `seq_stat_2pop_2N` to compute statistics using only two diploid individuals. It was used in Allio et al. 2021.

- `seq_stat_coding` to compute statistics on only population using a alignment of coding sequences.

Last update of `seq_stat_2pop` (Vesion 2.0; source : seq_stat_2pop_bppV3.cpp) is written in [Bio++ V3](https://biopp.github.io/).

*Author:* Benoit Nabholz

*References:* 
- Guéguen L, Gaillard S, Boussau B, Gouy M, Groussin M, Rochette NC, Bigot T, Fournier D, Pouyet F, Cahais V, et al. 2013. Bio++: Efficient Extensible Libraries and Tools for Computational Molecular Evolution. Mol. Biol. Evol. 30:1745–1750.
- Allio R, Tilak M-K, Scornavacca C, Avenant NL, Kitchener AC, Corre E, Nabholz B, Delsuc F. 2021. High-quality carnivoran genomes from roadkill samples enable comparative species delineation in aardwolf and bat-eared fox. Elife 10:e63167.


--------
### Installation

You can use the static executable compiled for linux x64 computer (see [Release](https://github.com/benoitnabholz/seq_stat/releases/)). You can also compile the program assuming that you have [**Bio ++**](https://biopp.github.io/) installed (here the Bio++ library V2 is locally installed in `$HOME/local/bpp/dev/` directory):

```
g++ --static -std=c++14 -g  seq_stat_2pop_bppV3.cpp -o seq_stat_2pop \
 -Wall -lbpp-phyl3 -lbpp-popgen3 -lbpp-seq3 -lbpp-core3  \
 -I$HOME/local/include -L$HOME/local/lib
```

--------
## seq_stat_2pop

###  Usage:
```
seq_stat_2pop -seq [listSeq] -f [phylip or fasta] -coding [coding or non-coding] \
 -tvts [tv/ts ratio for computing NSS] -pop1 [prefix_pop1] -pop2 [prefix_pop2] \
 -outgroup [prefix_out] -o [out file]
```

## Options :

- seq : a text file with the list of the sequence to analysed.
- f : sequence format `fasta` ot `phylip`
- coding : `coding` if protein coding sequences (only Standard Genetic Code)
- tvts : the transition over transversion ratio used for the computation of the number of synonymous site 
- pop1 : the sequence name of the population 1 must include this prefix in their names (e.g. "PopA_ind1", "PopA_ind2" etc...).
- pop2 : the sequence name of the population 2 must include this prefix in their names (e.g. "PopB_ind1", "PopB_ind2" etc...).
- outgroup : the sequence name of the outgroup must include this prefix in their names.
- o : the name of the out file in csv format.

## Statistics :
- Size : Size of the alignment (bp)
- S_Pop : Number of polymorphic site
- Pi_Pop : Tajima's estimator of nucleotides diversity
- W_Pop : Watterson's estimator of nucleotides diversity
- D_Pop : Tajima's D
- Dxy or Pi between : Average number of pairwise differences between sequences from two populations, excluding all comparisons between sequences within populations (Nei 1987; Cruickshank & Hahn, 2014)
- FstHud : Fst computed as in Hudson et al. 1992 Genetics 132:153 eq. 3
- FstNei_w : Fst computed as 1 - mean_Pi_Intra_Pop / Pi_total with mean weigthed using sample size (Nei 1982) with w = 1.* n_x/(n_x + n_y) and meanPiIntra = w*p_x + (1-w)*p_y. With n_x and n_y is the sample size of the population x and y respectively.
- FstNei_uw : Fst computed as 1 - mean_Pi_Intra_Pop / Pi_total (Nei 1982)
- PiInter : Inter popultaion nucleotides diversity
- Ts_Pop : Number of transition
- gc : GC content
- sizeOut : Size of outgroup sequence excluding gap and unresolved site (bp)
- divOut : Mean divergence between outgroup sequence and population n°2

## Statistics for coding sequences :
- PS_Pop : Nucleotide diversity of synonymous sites
- PN_Pop : Nucleotide diversity of non-synonymous sites
- gc3 : GC content at the third codon position
- NSS_Pop : Number of synonymous site
- div_Syn_Out : Mean synonymous divergence between outgroup sequence and population n°2
- div_NonSyn_Out : Mean non-synonymous divergence between outgroup sequence and population n°2

-----
## Example file

Two fasta sequences are store in the directory `data/`.
There are two populations (`Tguttata` and `Pacuticauda`) but no outgroup.

To use the program:

``` 
# No outgroup are availbale in this example.

# store the sequence name in a list
ls data/*fasta >list_file

# compute sequence statistics according that the sequence are non-coding DNA 
seq_stat_2pop -seq list_file -f fasta -coding non-coding -tvts 1.0 \
 -pop1 Tguttata -pop2 Pacuticauda -outgroup NA -o out_non_coding.csv

# compute sequence statistics according that the sequence are coding DNA 
seq_stat_2pop -seq list_file -f fasta -coding coding -tvts 2.0 \
 -pop1 Tguttata -pop2 Pacuticauda -outgroup NA -o out_coding.csv

```

--------

## seq_stat_2pop_2N

###  Usage of :
```
seq_stat_2pop_2N -seq [listSeq] -f [phylip or fasta] -pop1 [prefix_pop1] -pop2 [prefix_pop2] -o [out file]
```

## Options :
- seq : a text file with the list of the sequence to analysed.
- f : sequence format `fasta` ot `phylip`
- pop1 : the sequence name of the population 1 must include this prefix in their names (e.g. "PopA_ind1", "PopA_ind2" etc...).
- pop2 : the sequence name of the population 2 must include this prefix in their names (e.g. "PopB_ind1", "PopB_ind2" etc...).
- o : the name of the out file in csv format.


## Statistics :
- Fixed : Substitution between individual 1 and 2
- PrivatePop1 &	PrivatePop2 : Heterozygous position unique to individual 1 and 2 respectively.
- Shared : Heterozygous position shared between individual 1 and 2
- Pi1 & Pi2 : Heterozygosity of individual 1 and 2
- PiTot : Tajima's estimator of nucleotides diversity of individual 1 and 2 combined

--------
## seq_stat_coding

###  Usage:
```
seq_stat_coding -seq [listSeq] -f [phylip or fasta] -tstv [ts/tv ratio for computing NSS] -code [univ or mtmam or mtinv or mtechi] -o [out file]
```
## Options :
- seq : a text file with the list of the sequence to analysed.
- f : sequence format `fasta` ot `phylip`
- tvts : the transition over transversion ratio used for the computation of the number of synonymous site 
- code : genetic code (univ = standard universal; mt for mitochondrial
- o : the name of the out file in csv format.



