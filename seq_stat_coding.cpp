/**

Compute basic sequence statistic without considering any sub-populations.

g++ -std=c++14 --static -g ~/pCloudDrive/Progr/C++/popgen/seq_stat_coding.cpp -o ~/bin/seq_stat_coding -I$HOME/local/bpp/dev/include/ -L$HOME/local/bpp/dev/lib/ -DVIRTUAL_COV=yes -Wall -lbpp-popgen -lbpp-seq -lbpp-core
strip ~/bin/seq_stat_coding
 
**/


#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/PopGen/SequenceStatistics.h>
#include <Bpp/PopGen/PolymorphismSequenceContainerTools.h>
#include <Bpp/PopGen/PolymorphismSequenceContainer.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Seq/GeneticCode/VertebrateMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/InvertebrateMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/EchinodermMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/AscidianMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/MoldMitochondrialGeneticCode.h>


#include <Bpp/Numeric/VectorTools.h>

#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <iostream>


using namespace std;
using namespace bpp;

/************************************************************************************************/
double varTajima83(const PolymorphismSequenceContainer& psc)
{
	unsigned int n = psc.getNumberOfSequences();
	vector <double> vpi;
	for(unsigned int i = 0; i < n; i++){
		for(unsigned int j = i; j < n; j++){
			if(i == j)
				continue;
			PolymorphismSequenceContainer * tmpPsc = new PolymorphismSequenceContainer(psc.getAlphabet());
			tmpPsc->addSequence(psc.getSequence(i));
			tmpPsc->addSequence(psc.getSequence(j));

			vpi.push_back(SequenceStatistics::tajima83( *tmpPsc, true ) / double (PolymorphismSequenceContainerTools::getNumberOfCompleteSites(*tmpPsc, false)));
			delete tmpPsc;
		}
	}
	return VectorTools::var<double, double>(vpi);
			
}
/******************************************************************************/
bool isSynonymousPolymorphic2(const Site & site, const GeneticCode & gc)
    throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  //Alphabet checking
  if(!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::isSynonymousPolymorphic: alphabet is not CodonAlphabet", site.getAlphabet());
  if(site.getAlphabet()->getAlphabetType() != gc.getSourceAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("CodonSiteTools::isSynonymousPolymorphic: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gc.getSourceAlphabet());
  //Empty site checking
  if(site.size() == 0)
    throw EmptySiteException("CodonSiteTools::isSynonymousPolymorphic: Incorrect specified site", &site);

  // Global polymorphism checking
  if (SiteTools::isConstant(site)) return false;

  // Synonymous polymorphism checking
  vector<int> prot;
  int first_aa = -1;
  for(unsigned int i = 0; i < site.size(); i++) {
		
		if(site.getAlphabet()->isUnresolved(site[i]) || site.getAlphabet()->isGap(site[i])){
			continue;
		}
		if(first_aa == -1){
			first_aa = gc.translate(site[i]);
			continue;
		}
		
    int aa = gc.translate(site[i]);
    if (aa != first_aa) return false;
  }
  return true;
}
/******************************************************************************/
/** Modification compare to bio++ : In case where gapFlag is "false", the function create for each site, 
 * a new site whithout unresolved or gap position and compute the mean number of synonymous position on that. **/
double meanSynonymousSitesNumber2(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double ratio, bool gapflag)
{
  double S = 0.;
  ConstSiteIterator* si = NULL;
  if (gapflag)
		si = new CompleteSiteContainerIterator(psc);
  else
		si = new SimpleSiteContainerIterator(psc);
  const Site* site = 0;
  while(si->hasMoreSites()) {
    site = si->nextSite();
    // creat a site without gap or unresolved position
	Site siteSg = Site(site->getAlphabet());
	for(unsigned int s = 0; s < site->size(); s++){
		if (!site->getAlphabet()->isGap(site->getValue(s)) && !site->getAlphabet()->isUnresolved(site->getValue(s)))
			siteSg.addElement(site->getValue(s));
	}
	S +=  CodonSiteTools::meanNumberOfSynonymousPositions(siteSg, gc, ratio);
  }
  delete si;
  return S;
}
/******************************************************************************/
/** in CodonSiteTools **/
double piSynonymous3(const Site & site, const GeneticCode & gc, bool minchange)
    throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  //Alphabet checking
  if(!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::piSynonymous: alphabet is not CodonAlphabet", site.getAlphabet());
  if(site.getAlphabet()->getAlphabetType() != gc.getSourceAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("CodonSiteTools::piSynonymous3: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gc.getSourceAlphabet());
  //Empty site checking
  if(site.size() == 0)
    throw EmptySiteException("CodonSiteTools::piSynonymous3: Incorrect specified site", &site);

  //General polymorphism checking
  if (SiteTools::isConstant(site)) return 0;
  //Computation
  map<int,double> freq;
  SiteTools::getFrequencies(site, freq);
  double pi = 0;
  for(map<int,double>::iterator it1 = freq.begin(); it1 != freq.end(); it1++)
  {
    for(map<int,double>::iterator it2 = freq.begin(); it2 != freq.end(); it2++) {
			
			//if(site.getAlphabet()->isUnresolved(it2 -> first) || site.getAlphabet()->isGap(it2 -> first) || site.getAlphabet()->isUnresolved(it1 -> first) || site.getAlphabet()->isGap(it1 -> first)){
				//continue;
			//}			
			pi += (it1 -> second) * (it2 -> second) * (CodonSiteTools::numberOfSynonymousDifferences(it1->first,it2->first,gc,minchange));
    }
  }
  unsigned int n = site.size();
  return pi * n / (n - 1);
}

/******************************************************************************/
/******************************************************************************/
/** in CodonSiteTools **/
double piNonSynonymous3(const Site & site, const GeneticCode & gc, bool minchange)
    throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  //Alphabet checking
  if(!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::piNonSynonymous3: alphabet is not CodonAlphabet", site.getAlphabet());
  if(site.getAlphabet()->getAlphabetType() != gc.getSourceAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("CodonSiteTools::piNonSynonymous3: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gc.getSourceAlphabet());
  //Empty site checking
  if(site.size() == 0)
    throw EmptySiteException("CodonSiteTools::piNonSynonymous3: Incorrect specified site", &site);

  //General polymorphism checking
  if(SiteTools::isConstant(site)) return 0;
  if(isSynonymousPolymorphic2(site,gc)) return 0;
  //Computation
  map<int,double> freq;
  CodonSiteTools::getFrequencies(site, freq);
  const CodonAlphabet * ca = dynamic_cast<const CodonAlphabet *>(site.getAlphabet());
  double pi = 0;
  for(map<int,double>::iterator it1 = freq.begin(); it1 != freq.end(); it1++) {
    for(map<int,double>::iterator it2 = freq.begin(); it2 != freq.end(); it2++) {
			
			/** manual modification here **/
			//if(site.getAlphabet()->isUnresolved(it2 -> first) || site.getAlphabet()->isGap(it2 -> first) || site.getAlphabet()->isUnresolved(it1 -> first) || site.getAlphabet()->isGap(it1 -> first)){
			//	continue;
			//}
			
      unsigned int nbtot = CodonSiteTools::numberOfDifferences(it1->first,it2->first, *ca);
      double nbsyn = CodonSiteTools::numberOfSynonymousDifferences(it1->first, it2 -> first, gc, minchange);
      pi += (it1 -> second) * (it2 -> second) * (nbtot - nbsyn);
    }
  }
  unsigned int n = site.size();
  return pi * n / (n - 1);
}

/******************************************************************************/
/** in SequenceStatistics **/
double piSynonymous2(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, bool minchange, bool gapflag)
{
  double S = 0.;
	ConstSiteIterator* si = 0;
  if (gapflag)
		si = new CompleteSiteContainerIterator(psc);
  else
		si = new SimpleSiteContainerIterator(psc);


  const Site* site = 0;
  while(si->hasMoreSites()) {
    site = si->nextSite();
    // creat a site without gap or unresolved site
	Site siteSg = Site(site->getAlphabet());
	for(unsigned int s = 0; s < site->size(); s++){
		if (!site->getAlphabet()->isGap(site->getValue(s)) && !site->getAlphabet()->isUnresolved(site->getValue(s)))
			siteSg.addElement(site->getValue(s));
	}
    S += piSynonymous3(siteSg, gc, minchange);
  }
  delete si;
  return S;
}
/******************************************************************************/
/** SequenceStatistics **/
double piNonSynonymous2(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, bool minchange, bool gapflag)
{
  double S = 0.;
	ConstSiteIterator* si = 0;
  if (gapflag)
		si = new CompleteSiteContainerIterator(psc);
  else
		si = new SimpleSiteContainerIterator(psc);
  const Site* site = 0;

  while(si->hasMoreSites()) {
    site = si->nextSite();
    // creat a site without gap or unresolved site
	Site siteSg = Site(site->getAlphabet());
	for(unsigned int s = 0; s < site->size(); s++){
		if (!site->getAlphabet()->isGap(site->getValue(s)) && !site->getAlphabet()->isUnresolved(site->getValue(s)))
			siteSg.addElement(site->getValue(s));
	}
    S += piNonSynonymous3(siteSg, gc, minchange);
  }
  delete si;
  return S;
}
/******************************************************************************/
bool isSingleton(const Site& site) {
  bool sing = false;
  map<int, unsigned long int> states_count;
  SymbolListTools::getCounts(site, states_count);
  for (map<int, unsigned long int>::iterator it = states_count.begin() ; it != states_count.end() ; it++)
	{
		
		if(site.getAlphabet()->isUnresolved(it->first) || site.getAlphabet()->isGap(it->first))
				continue;
		
		if (it->second == 1)
		{
      sing = true;
			break;
		}
	}
  return sing;
}
/******************************************************************************/
double getGCThirdCodonPosition(const Sequence * seq,  bool ignoreUnresolved=true, bool ignoreGap=true)
{
  const Alphabet * alphabet = seq->getAlphabet();
  if (!AlphabetTools::isNucleicAlphabet(alphabet)){
    cout << "getGCThirdCodonPosition. Method only works on nucleotides." << endl;
    return 0;
  }

  unsigned int l = seq->size();
  double gc = 0;
  double total = 0;
  unsigned int j = l/3;

  for( unsigned int k=1; k<=j; k ++){
    unsigned int entier = k * 3;
    unsigned int pos = entier-1;

    int state = seq->getValue(pos);
    if (state > -1) { // not a gap
      if (state == 1 || state == 2) { // G or C
        gc++;
        total++;
      } else if (state == 0 || state == 3) { // A, T or U
        total++;
      } else { // Unresolved character
        if (!ignoreUnresolved) {
          total++;
          switch(state) {
            case(7): gc++; break;// G or C
            case(4): gc+=0.5; break;// A or C
            case(5): gc+=0.5; break;// A or G
            case(6): gc+=0.5; break;// C or T
            case(9): gc+=0.5; break;// G or T
            case(10): gc+=2./3.; break;// A or C or G
            case(11): gc+=1./3.; break;// A or C or T
            case(12): gc+=1./3.; break;// A or G or T
            case(13): gc+=2./3.; break;// C or G or T
            case(14): gc+=0.5; break;// A or C or G or T
          }
        }
      }
    } else {
      if (!ignoreGap) total++;
    }
  }
  return total != 0 ? gc/total : 0;
}


int main (int argc, char* argv[]){
	
try{
 
 
if (argc == 1 || argc < 10)
{
	cout << "seq_stat_coding -seq [listSeq] -f [phylip or fasta] -tstv [ts/tv ratio for computing NSS] -code [univ or mtmam or mtinv or mtechi] -o [out file]" << endl;
	
	cout << "### Statistics :" << endl;
	cout << "\tSize : Size of the alignment (bp)" << endl;
	cout << "\tS : Number of polymorphic site" << endl;
	cout << "\tP : Tajima's estimator of nucleotides diversity" << endl;
	cout << "\tW : Watterson's estimator of nucleotides diversity" << endl;
	cout << "\tPs : Nucleotide diversity of synonymous sites" << endl;
	cout << "\tPn : Nucleotide diversity of non-synonymous sites" << endl;
	cout << "\tNSS : Number of synonymous site" << endl;
	cout << "\tD_Taj : Tajima's D" << endl;
	cout << "\tgc3 : GC content at the third codon position of the consensus sequence (ignoring gap and unknown characters)" << endl;
	return 0;
}

string listName, format, code, typeAlg, outF;
double tstv = -1.0;

/*************************************************************/
/** read command line **/
/*************************************************************/
unsigned int i = 1;
while (i < argc){
	string s = argv[i];
	if (s == "-seq")	{
		i++;
		if (i == argc) 	{
		cerr << "error in command: seq <listFile>\n";
		cerr << '\n';
		exit(1);
		}
		listName = argv[i];
	}
	if (s == "-f")	{
		i++;
		if (i == argc) 	{
			cerr << "error in command: -f <phylip or fasta>\n";
			cerr << '\n';
			exit(1);
		}
		format = argv[i];
	}
	if (s == "-tstv")	{
		i++;
		if (i == argc) 	{
			cerr << "error in command: -tstv <double tstv>\n";
			cerr << '\n';
			exit(1);
		}

		tstv = TextTools::toDouble(argv[i]);
	}
	if (s == "-code")	{
		i++;
		if (i == argc) 	{
		cerr << "error in command: -code <uni or mtmam or mtinv or mtechi or mtasci or mtmold or non-coding>\n";
		cerr << '\n';
		exit(1);
		}
		code = argv[i];
	}
	if (s == "-o")	{
		i++;
		if (i == argc) 	{
		cerr << "error in command: -o <outFile>\n";
		cerr << '\n';
		exit(1);
		}
		outF = argv[i];
	}
	i++;
}

if(tstv == -1){
  cerr << "ERROR!!! -tstv should be provided !!" << endl;
  cerr << "Command ligne should be like:" << endl;
  cerr << "seq_stat_coding -seq [listSeq] -f [phylip or fasta] -tstv [ts/tv ratio for computing NSS] -code [uni or mtmam or mtinv or mtechi or mtasci or mtmold or non-coding] -o [out file]" << endl;
  exit(-1);
}

if(!FileTools::fileExists( listName )) {
  cerr << "ERROR!!! File " << listName << " does not exists." << endl;
  exit(-1);
}

ifstream Filelist (listName.c_str(), ios::in);
ofstream Fileout (outF.c_str(), ios::out);

if(code == "noncoding")
	Fileout << "name\tsize\tN\tS\tPi\tW\tD_taj" <<  endl;
else
	Fileout << "name\tSize\tN\tS\tP\tW\tPs\tPn\tNSS\tD_Taj\tGC3" <<  endl;
	
while (!Filelist.eof ()){
  
	string nomfic = FileTools::getNextLine(Filelist);
	cout << nomfic << endl;
	if(TextTools::isEmpty(nomfic)){
		continue;
	}
	if(!FileTools::fileExists( nomfic )) {
		cerr << "ERROR!!! File " << nomfic << " does not exists." << endl;
		exit(-1);
	}
  
	const NucleicAlphabet * alpha = new DNA();
	const CodonAlphabet *codonAlpha = new CodonAlphabet(alpha);
	const GeneticCode *GC = NULL;
	if(code != "univ" && code != "mtmam" && code != "mtinv" && code != "mtechi" &&  code != "mtasci" &&  code != "mtmold" && code != "noncoding"){
		cerr << "error in command: -code <uni or mtmam or mtinv or mtechi or mtasci or mtmold or non-coding>\n";
		return 1;
	}
	
	
	if(code == "univ")
		GC = new StandardGeneticCode(alpha);
	if(code == "mtmam")
		GC = new VertebrateMitochondrialGeneticCode(alpha);
	if(code == "mtinv")
		GC = new InvertebrateMitochondrialGeneticCode(alpha);
	if(code == "mtechi")
		GC = new EchinodermMitochondrialGeneticCode(alpha);
	if(code == "mtasci")
		GC = new AscidianMitochondrialGeneticCode(alpha);
	if(code == "mtmold")
		GC = new MoldMitochondrialGeneticCode(alpha);
		
	VectorSequenceContainer * seqCont = NULL;

	Phylip * PhySeq = new Phylip(true, true);
	Fasta * FstSeq = new  Fasta;

	if(format == "phylip"){
		seqCont = PhySeq->readAlignment(nomfic , alpha);
	}
	if(format == "fasta"){
		seqCont = FstSeq->readAlignment(nomfic , alpha);
	}
  
	VectorSiteContainer *sitec = new VectorSiteContainer( *seqCont );
	PolymorphismSequenceContainer * psc1 = new PolymorphismSequenceContainer( *sitec);
	
	if(code == "noncoding"){
		unsigned int S = SequenceStatistics::numberOfPolymorphicSites(*psc1, false );
		Fileout << nomfic << "\t" << psc1->getNumberOfSequences() << "\t"<< psc1->getNumberOfSites() << "\t";
		Fileout << S << "\t" << SequenceStatistics::tajima83(*psc1, false) << "\t" << SequenceStatistics::watterson75(*psc1, false ) << "\t";
		if(S > 0)
			Fileout << SequenceStatistics::tajimaDss( *psc1, false ) << endl;
		else
			Fileout <<"-999" << endl;		
	}else{
	
		/*************************************************************/
		/** convert alphabet **/
		/*************************************************************/
	
		PolymorphismSequenceContainer * pscCodon = new PolymorphismSequenceContainer(codonAlpha);
		SequenceContainerTools::convertAlphabet(*psc1, *pscCodon);
	
		string sc = "no";
		for(unsigned int i = 0; i < pscCodon->getNumberOfSites(); i++){
			Site site = pscCodon->getSite(i);
			if(CodonSiteTools::hasStop(site, *GC)){
				sc = "yes";
				cout << nomfic << " has codon stop at codon " << i << endl;
				break;
			}
		}


		if(sc == "no"){
			unsigned int S = SequenceStatistics::numberOfPolymorphicSites(*psc1, false );
			Fileout << nomfic << "\t" << psc1->getNumberOfSites() << "\t" << psc1->getNumberOfSequences() << "\t";
			Fileout << S << "\t" << SequenceStatistics::tajima83(*psc1, false) << "\t" << SequenceStatistics::watterson75(*psc1, false ) << "\t";
			Fileout << piSynonymous2( *pscCodon,*GC, false, false) << "\t" << piNonSynonymous2( *pscCodon,*GC, false, false) << "\t";
			Fileout << meanSynonymousSitesNumber2( *pscCodon, *GC, tstv, false) <<"\t" ;
			if(S > 0)
				Fileout << SequenceStatistics::tajimaDss( *psc1, false ) << "\t";
			else
				Fileout <<"-999" << "\t";
				
			const Sequence * cons = SiteContainerTools::getConsensus(*psc1, "consensus", true, true );
			Fileout << getGCThirdCodonPosition(cons) << endl;
		}
		delete pscCodon;
	}
	

	delete FstSeq;
	delete PhySeq;
	delete alpha;
	delete sitec;
	delete seqCont;
	delete psc1;
	delete GC;
	
}

}
catch(exception & e){
  cout << e.what() << endl;
    }
return 0;
}
