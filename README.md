# The program compute population genetic statistics 

This is a badly written but useful program to compute population genetics statistics using sequences data in fasta or phylip format.
This program is written using the [**Bio ++** library ](http://biopp.univ-montp2.fr/wiki/index.php/Main_Page) (Guéguen et al. 2013).


*Author:* Benoit Nabholz

*Reference:* Guéguen L, Gaillard S, Boussau B, Gouy M, Groussin M, Rochette NC, Bigot T, Fournier D, Pouyet F, Cahais V, et al. 2013. Bio++: Efficient Extensible Libraries and Tools for Computational Molecular Evolution. Mol. Biol. Evol. 30:1745–1750.

--------
### Installation

You can use the static executable compiled for linux x64 computer (see [Release](https://github.com/benoitnabholz/seq_stat_2pop/releases/tag/v1)). You can also compile the program assuming that you have [**Bio ++**](http://biopp.univ-montp2.fr/wiki/index.php/Main_Page) installed (here the Bio++ library is locally installed in `$HOME/local/bpp/dev/` directory):

```
g++ -g seq_stat_2pop.cpp -o ./seq_stat_2pop \
 -DVIRTUAL_COV=yes -Wall \
 -I$HOME/local/bpp/dev/include  -L$HOME/local/bpp/dev/lib64 \
 -lbpp-popgen -lbpp-phyl -lbpp-seq -lbpp-core
```

--------

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
- tvts : the transition over transversion ratio used for the computation of the number of synonymous site ( see http://biopp.univ-montp2.fr/apidoc/bpp-seq/html/classbpp_1_1CodonSiteTools.html#a0fb6fb2b314dd557b52219dffac330e9 )
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
- Dxy: Mean pairwaise diverge between all ingroup sequences
- FstHud : Fst computed as in Hudson et al. 1992 Genetics 132:153 eq. 3
- FstNei_w : Fst computed as 1 - mean_Pi_Intra_Pop / Pi_total with mean weigthed using sample size (Nei 1982)
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

