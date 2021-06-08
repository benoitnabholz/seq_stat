/**

g++ -g /home/benoitnabholz/Progr/C++/popgen/seq_stat_2pop/seq_stat_2pop.cpp -o ./seq_stat_2pop \
 -DVIRTUAL_COV=yes -Wall \
 -I$HOME/local/bpp/dev/include  -L$HOME/local/bpp/dev/lib64 \
 -lbpp-popgen -lbpp-phyl -lbpp-seq -lbpp-core

g++ --static -std=c++14 -g ~/pCloudDrive/Progr/C++/popgen/seq_stat_2pop/seq_stat_2pop.cpp -o ~/bin/seq_stat_2pop -DVIRTUAL_COV=yes -Wall -I$HOME/local/bpp/dev/include  -L$HOME/local/bpp/dev/lib -lbpp-phyl -lbpp-popgen -lbpp-seq -lbpp-core -Wdeprecated 
strip ~/bin/seq_stat_2pop

* 
**/
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/StringSequenceTools.h>
#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include <Bpp/Seq/Container/SiteContainerIterator.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceTools.h>

#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Seq/Container/SiteContainerIterator.h>
#include <Bpp/PopGen/SequenceStatistics.h>
#include <Bpp/PopGen/PolymorphismSequenceContainerTools.h>
#include <Bpp/PopGen/PolymorphismSequenceContainer.h>

#include <Bpp/Numeric/VectorTools.h>

#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;
using namespace bpp;


/********************************************************************/
bool isSeqGapOrUnresolvedOnly(const Sequence& seq)
{
	// Main loop : for all characters in site
	for (unsigned int i = 0; i < seq.size(); i++)
	{
		if (!(seq.getAlphabet()->isGap(seq[i]) || seq.getAlphabet()->isUnresolved(seq[i]))) return false;
	}
	return true;
}
/********************************************************************/
bool hasGapOrUnresolvedOnly2(const Site& site)
{
  // Main loop : for all characters in site
  for (unsigned int i = 0; i < site.size(); i++)
  {
    if (site.getAlphabet()->isGap(site[i]) || site.getAlphabet()->isUnresolved(site[i])) return true;
  }
  return false;
}
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
double meanNumberOfSynonymousPositions2(const Site & site, const GeneticCode & gc, double ratio)
    throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  //Alphabet checking
  if(!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::meanNumberOfSynonymousPositions: alphabet is not CodonAlphabet", site.getAlphabet());
  if(site.getAlphabet()->getAlphabetType() != gc.getSourceAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("CodonSiteTools::meanNumberOfSynonymousPositions: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gc.getSourceAlphabet());
  //Empty site checking
  if(site.size() == 0)
    throw EmptySiteException("CodonSiteTools::meanNumberOfSynonymousPositions: Incorrect specified site", &site);

  //Computation
  double NbSyn=0;
  map<int,double> freq;
  SiteTools::getFrequencies(site, freq);
  for(map<int,double>::iterator it = freq.begin(); it != freq.end(); it++)
  {
    NbSyn += (it->second)*CodonSiteTools::numberOfSynonymousPositions(it->first, gc, ratio);
  }
  return NbSyn;
}

/******************************************************************************/
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
    // creat a site without gap or unresolved site
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
/******************  remove singleton *****************************************/
PolymorphismSequenceContainer * removeSingletion(const PolymorphismSequenceContainer& psc)
{
	PolymorphismSequenceContainer *NoSingPsc = new PolymorphismSequenceContainer(psc.getNumberOfSequences(), psc.getAlphabet());
	for(unsigned int i = 0; i < psc.getNumberOfSites(); i++)
	{
		Site site = psc.getSite(i);
		
		if(SiteTools::isConstant(site))
		{
			NoSingPsc->addSite(site);
			continue;
		}
		
		map< int, unsigned long int > counts;
		CodonSiteTools::getCounts(site, counts); 

		bool singleton = false;
		for(map<int, unsigned long int>::iterator it = counts.begin(); it != counts.end(); it++)
		{
			if(site.getAlphabet()->isUnresolved(it->first) || site.getAlphabet()->isGap(it->first))
				continue;
			
			if(it->second == 1)
			{
				singleton = true;				
			}
		}
		if(!singleton)
			NoSingPsc->addSite(site);
	}
	return NoSingPsc;
}
	
/******************************************************************************/
int getVariant(const Site& sitePol, const Site& siteMon){
	if(!SiteTools::isConstant(siteMon, true))
		throw ("Problem getVariant : siteMon is polymorphic");
	
	int alleleM = -1;
	int c = 0;
	while(alleleM < 0 || alleleM > 3){
		alleleM = siteMon.getValue(c);
		c++;
	}
	
	for(unsigned int i = 0; i < sitePol.size(); i++){
		int alleleP = sitePol.getValue(i);
		if(sitePol.getAlphabet()->isGap(alleleP) || sitePol.getAlphabet()->isUnresolved(alleleP))
			continue;
		
		if(alleleP != alleleM)
			return alleleP;
	}
	throw ("Problem getVariant");
	return -1; // should never arrived here
}
 
/*****************************************************************/
double getGCThirdCodonPosition(const Sequence & seq,  bool ignoreUnresolved=true, bool ignoreGap=true)
{
  const Alphabet * alphabet = seq.getAlphabet();
  if (!AlphabetTools::isNucleicAlphabet(alphabet)){
    cout << "getGCThirdCodonPosition. Method only works on nucleotides." << endl;
    return 0;
  }

  unsigned int l = seq.size();
  double gc = 0;
  double total = 0;
  unsigned int j = l/3;

  for( unsigned int k=1; k<=j; k ++){
    unsigned int entier = k * 3;
    unsigned int pos = entier-1;

    int state = seq.getValue(pos);
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
/********************************************************************/
/** Compute the mean number of differences between an outgroup sequence (seqOut) 
 *  and several ingroup sequences (PolymorphismSequenceContainer)
 * 
 *  This take into account the potentially shared polymorphism 
 *  
 **/
/********************************************************************/
double computeMeanNumberOfDifference(const PolymorphismSequenceContainer& psc, BasicSequence seqOut, unsigned int minNumSites){
	
	vector <double> diff;
	double meanDiff;
	
	for(unsigned int i = 0; i < psc.getNumberOfSequences(); i++){
		const Sequence& seq1 = psc.getSequence(i);
		
		/** check number of complete site before computing the divergence **/
		PolymorphismSequenceContainer * tmpSc = new PolymorphismSequenceContainer(seqOut.getAlphabet());
		tmpSc->addSequence(seq1);
		tmpSc->addSequence(seqOut);
		SiteContainer * tmpComp = SiteContainerTools::getCompleteSites(*tmpSc);
		
		if(tmpComp->getNumberOfSites() > minNumSites){
			/** compute the divergence **/
			double d = SiteContainerTools::computeSimilarity(seqOut, seq1, true, "no gap", true)*tmpComp->getNumberOfSites();
			diff.push_back(d);
// 			cout << d << endl;
		}
		delete tmpSc;
		delete tmpComp;
	}
	meanDiff = VectorTools::mean<double, double>(diff);
	
	return meanDiff;
}
/********************************************************************/
vector <double> computeMeanSynonymousAndNonSynonymousDifferences(const PolymorphismSequenceContainer& psc, 
  const Sequence& seqOut, const GeneticCode& gc, double tvts){
	
	vector <double> diffNS, diffN, diffS;
	double NSS;
	for(unsigned int i = 0; i < psc.getNumberOfSequences(); i++){
		const Sequence& seq1 = psc.getSequence(i);
		
		/** check number of complete site before computing the divergence **/
		PolymorphismSequenceContainer * tmpSc = new PolymorphismSequenceContainer(seqOut.getAlphabet());
		tmpSc->addSequence(seq1);
		tmpSc->addSequence(seqOut);
		SiteContainer * tmpComp = SiteContainerTools::getCompleteSites(*tmpSc);
		
		NSS = meanSynonymousSitesNumber2(*tmpComp,gc,tvts, false);

		//if(tmpComp->getNumberOfSites() > 15){
			
			double dS = SequenceStatistics::numberOfSynonymousSubstitutions(*tmpComp, gc, 0.0);
			diffS.push_back(dS);
			
			double dN = SequenceStatistics::numberOfSynonymousSubstitutions(*tmpComp, gc, 0.0);
			diffN.push_back(dN);;
		//}
		delete tmpSc;
		delete tmpComp;
	}
	diffNS.push_back(VectorTools::mean<double, double>(diffN));
	diffNS.push_back(VectorTools::mean<double, double>(diffS));
	diffNS.push_back(NSS);
	
	return diffNS;
}
/********************************************************************/
/** PiInter : compute the Tajima's Pi between population ************/
/********************************************************************/

double PiInter(const PolymorphismSequenceContainer& psc, unsigned int id1, unsigned int id2)
{
	vector<double> vdiff;
	double piInter;
	
	PolymorphismSequenceContainer *Pop1 = PolymorphismSequenceContainerTools::extractGroup(psc, id1);
	PolymorphismSequenceContainer *Pop2 = PolymorphismSequenceContainerTools::extractGroup(psc, id2);
	SiteContainerTools::changeUnresolvedCharactersToGaps (*Pop1);
	SiteContainerTools::changeUnresolvedCharactersToGaps (*Pop2);
	unsigned int N, totLength, meanLength = 0;
	N = totLength = meanLength = 0;
	for(unsigned int i = 0; i<Pop1->getNumberOfSequences(); i++){
		const Sequence &s1 = Pop1->getSequence(i);
		for(unsigned int j = 0; j<Pop2->getNumberOfSequences(); j++){
			N++;
			const Sequence &s2 = Pop2->getSequence(j);
			vdiff.push_back(SiteContainerTools::computeSimilarity(s1,s2,true,"no gap",true));
			totLength += SequenceTools::getNumberOfPositionsWithoutGap(s1,s2);
		}
	}
	piInter = (VectorTools::sum(vdiff) / N)*(totLength/N);
	
	delete Pop1;
	delete Pop2;
	
	return piInter;
	
}



/**********************************************************************************               ****************/
/** Fst from eq. 3 of  Hudson, Slatkin and Maddison 1992 Genetics 132:153 ****************************************/
/********************************************************************************               ******************/
/*
double FstHudson92(const PolymorphismSequenceContainer& psc, unsigned int id1, unsigned int id2, bool gapflag)
{
	vector<double> vdiff;
	double piIntra1, piIntra2, meanPiIntra, piInter, Fst;
	
	PolymorphismSequenceContainer *Pop1 = PolymorphismSequenceContainerTools::extractGroup(psc, id1);
	PolymorphismSequenceContainer *Pop2 = PolymorphismSequenceContainerTools::extractGroup(psc, id2);
	
  piIntra1 = SequenceStatistics::tajima83(*Pop1, gapflag);
	piIntra2 = SequenceStatistics::tajima83(*Pop2, gapflag);
	
	meanPiIntra = (piIntra1+piIntra2)/2;
	
	unsigned int N = 0;
	for(unsigned int i = 0; i<Pop1->getNumberOfSequences(); i++){
		const Sequence &s1 = Pop1->getSequence(i);
		for(unsigned int j = 0; j<Pop2->getNumberOfSequences(); j++){
			N++;
			const Sequence &s2 = Pop2->getSequence(j);
			vdiff.push_back(SiteContainerTools::computeSimilarity(s1,s2,true,"no gap",true));
		}
	}
	piInter = (VectorTools::sum(vdiff) / N)*psc.getNumberOfSites();
	
	Fst = 1.0 - meanPiIntra/piInter;
	
	delete Pop1;
	delete Pop2;
	
	return Fst;
	
}
* */
double FstHudson92(const PolymorphismSequenceContainer& psc, unsigned int id1, unsigned int id2, bool gapflag)
{
	vector<double> vdiff;
	double piIntra1, piIntra2, meanPiIntra, piInter, Fst;
	
	PolymorphismSequenceContainer *Pop1 = PolymorphismSequenceContainerTools::extractGroup(psc, id1);
	PolymorphismSequenceContainer *Pop2 = PolymorphismSequenceContainerTools::extractGroup(psc, id2);
	
	piIntra1 = SequenceStatistics::tajima83(*Pop1, false);
	piIntra2 = SequenceStatistics::tajima83(*Pop2, false);
	
	meanPiIntra = (piIntra1+piIntra2)/2;
	
	piInter = PiInter(psc, id1, id2);
	
	Fst = 1.0 - meanPiIntra/piInter;
	
	delete Pop1;
	delete Pop2;
	
	return Fst;
	
}
/**********************************************************************************               ****************/
/** Nei 1982 Fst =  (Pi Total - Pi intra) / Pi Total  ***************************************/
/** Pi intra could be weighted mean or not 
 * with w = 1.* n_x/(n_x + n_y) and meanPiIntra = w*p_x + (1-w)*p_y
 * no weighted is w = 0.5
/********************************************************************************               ******************/

double FstNei82(const PolymorphismSequenceContainer& psc, unsigned int id1, unsigned int id2, bool weightedPiIntra, bool gapflag)
{
	vector<double> vdiff;
	double piTotal, piIntra1, piIntra2, meanPiIntra, Fst;
	double w = 0.5;

	/** select ingroup sequences **/
	SequenceSelection ss;
	PolymorphismSequenceContainer* psci = dynamic_cast<PolymorphismSequenceContainer*>(psc.clone());
	for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++){
		if (psc.getGroupId(i) != id1 && psc.getGroupId(i) != id2)
			ss.push_back(i);
	}
	for (unsigned int i = ss.size(); i > 0; i--){
		psci->getGroupId(i-1);
		psci->deleteSequence(ss[i - 1]);
	}
	piTotal = SequenceStatistics::tajima83(*psci, gapflag);
	
	PolymorphismSequenceContainer *Pop1 = PolymorphismSequenceContainerTools::extractGroup(psc, id1);
	PolymorphismSequenceContainer *Pop2 = PolymorphismSequenceContainerTools::extractGroup(psc, id2);
	
	piIntra1 = SequenceStatistics::tajima83(*Pop1, gapflag);
	piIntra2 = SequenceStatistics::tajima83(*Pop2, gapflag);
	
	if(weightedPiIntra){
		int n_x = Pop1->getNumberOfSequences();
		int n_y = Pop2->getNumberOfSequences();
		w = 1.* n_x/(n_x + n_y);
	}
	//cout << w << endl;

	meanPiIntra = w * piIntra1 + w * piIntra2;
	
	Fst = 1 - meanPiIntra/piTotal;
	
	delete Pop1;
	delete Pop2;
	delete psci;
	return Fst;
	
}

/************************************************************************************************/
bool isGapOrUnresolvedOnly2(const Site& site)
{
  // Main loop : for all characters in site
  for (unsigned int i = 0; i < site.size(); i++)
  {
    if (!site.getAlphabet()->isGap(site[i]) && !site.getAlphabet()->isUnresolved(site[i])) return false;
  }
  return true;
}
/**************************************************************************************************/
size_t getNumberOfTransitions2(const PolymorphismSequenceContainer& psc)
{
	size_t nbT = 0;
	ConstSiteIterator* si = new SimpleSiteContainerIterator(psc);
	const Site* site = 0;
	while (si->hasMoreSites())
	{
		site = si->nextSite();
		// creat a site without gap or unresolved site
		Site siteSg = Site(site->getAlphabet());
		for(unsigned int s = 0; s < site->size(); s++){
			if (!site->getAlphabet()->isGap(site->getValue(s)) && !site->getAlphabet()->isUnresolved(site->getValue(s))){
				siteSg.addElement(site->getValue(s));
			}
		}
		
		// if (SiteTools::isConstant(*site) || SiteTools::isTriplet(*site)) continue;
		if (SiteTools::getNumberOfDistinctCharacters(siteSg) != 2)
			continue;
		
		int state1 = (*site)[0];
		int state2 = (*site)[0];
		for (size_t i = 1; i < site->size(); i++)
		{
			if (state1 != (*site)[i])
			{
				state2 = (*site)[i];
				break;
			}
		}
		if(((state1 == 0 && state2 == 2) || (state1 == 2 && state2 == 0)) ||
			((state1 == 1 && state2 == 3) || (state1 == 3 && state2 == 1)))
		{
			nbT++;
		}
	}
	delete si;
	return nbT;
}

/** **********************************************************************************************/


/********************************************************************/
/** Number of Fixed, Private and Shared polymorphism between population *******/ 
/** */
/********************************************************************/

vector<int> NumberOfDifferenceBetweenPopulations(const PolymorphismSequenceContainer& psc, unsigned int id1, unsigned int id2)
{
	vector<int> vdiff;
		
	PolymorphismSequenceContainer *Pop1 = PolymorphismSequenceContainerTools::extractGroup(psc, id1);
	PolymorphismSequenceContainer *Pop2 = PolymorphismSequenceContainerTools::extractGroup(psc, id2);
	
	SiteContainerTools::changeUnresolvedCharactersToGaps (*Pop1);
	SiteContainerTools::changeUnresolvedCharactersToGaps (*Pop2);
	
	int Fixed = 0; 
	int PrivatePop1 = 0;
	int PrivatePop2 = 0;
	int Shared = 0;
	
	for(unsigned int i = 0; i < psc.getNumberOfSites(); i++){
		
		const Site& site = psc.getSite(i);
		const Site& sitePop1 = Pop1->getSite(i);
		const Site& sitePop2 = Pop2->getSite(i);
		string ok = "no";
		if(!SiteTools::isConstant(site, true, true)){
			if(SiteTools::isConstant(sitePop1, true, true) && SiteTools::isConstant(sitePop2, true, true)){
				Fixed++;
				ok = "ok";
			}
			
			if(!SiteTools::isConstant(sitePop1, true, true) && SiteTools::isConstant(sitePop2, true, true)){
				PrivatePop1++;
				ok = "ok";
			}
				
			if(SiteTools::isConstant(sitePop1, true, true) && !SiteTools::isConstant(sitePop2, true, true)){
				PrivatePop2++;
				ok = "ok";
			}
			if(!SiteTools::isConstant(sitePop1, true, true) && !SiteTools::isConstant(sitePop2, true, true)){
				//cout << "Shared " << i << endl;
				Shared++;
				ok = "ok";
			}
			if( ok == "no" ){
				cout << site.toString() << endl;
			}
		}
		
		
	}
	
	vdiff.push_back(Fixed);
	vdiff.push_back(PrivatePop1);
	vdiff.push_back(PrivatePop2);
	vdiff.push_back(Shared);
	
	delete Pop1;
	delete Pop2;
	
	return vdiff;
	
}


/** **********************************************************************************************/


int main (int argc, char* argv[]){
	
try{
if (argc == 1 || argc < 13)
{
	cout << "\n#####################";
	cout << "\nVersion 1.1\n";
	cout << "\n##### USAGE :\n";
	cout << "\nseq_stat_2pop -seq [listSeq] -f [phylip or fasta] -coding [coding or non-coding] -tvts [tv/ts ratio for computing NSS] -pop1 [prefix_pop1] -pop2 [prefix_pop2] -outgroup [prefix_out] -o [out file]\n" << endl;
	
	cout << "### Statistics :" << endl;
	cout << "\tSize : Size of the alignment (bp)" << endl;
	cout << "\tS_Pop : Number of polymorphic site" << endl;
	cout << "\tPi_Pop : Tajima's estimator of nucleotides diversity" << endl;
	cout << "\tW_Pop : Watterson's estimator of nucleotides diversity" << endl;
	cout << "\tD_Pop : Tajima's D" << endl;
	cout << "\tDxy or Pi between : Average number of pairwise differences between sequences from two populations, excluding all comparisons between sequences within populations (Nei 1987; Cruickshank & Hahn, 2014)" << endl;
	cout << "\tFstHud : Fst computed as 1 - mean_Pi_Intra_Pop / Pi Between (or Dxy) (Hudson et al. 1992 Genetics 132:153 eq. 3)" << endl;
	cout << "\tFstNei : Fst computed as 1 - mean_Pi_Intra_Pop / Pi Total (Nei 1982)" << endl;
	cout << "\tFstNei_weighted : weigth by the sample size of each population with w = 1.* n_x/(n_x + n_y) and meanPiIntra = w*p_x + (1-w)*p_y" << endl;
	cout << "\tFstNei_unweighted : with w = 0.5 and meanPiIntra = w*p_x + (1-w)*p_y" << endl;
	cout << "\tPiTotal : Nucleotides diversity" << endl;
	cout << "\tTs_Pop : Number of transition" << endl;
	cout << "\tgc : GC content" << endl;
	cout << "\tsizeOut : Size of outgroup sequence excluding gap and unresolved site (bp)" << endl;
	cout << "\tdivOut : Mean divergence between outgroup sequence and population n°2" << endl;

	cout << "### Statistics for coding sequences only :" << endl;
	cout << "\tPS_Pop : Nucleotide diversity of synonymous sites" << endl;
	cout << "\tPN_Pop : Nucleotide diversity of non-synonymous sites" << endl;
	cout << "\tgc3 : GC content at the third codon position" << endl;
	cout << "\tNSS_Pop : Number of synonymous site" << endl;
	cout << "\tdiv_Syn_Out : Mean synonymous divergence between outgroup sequence and population n°2" << endl;
	cout << "\tdiv_NonSyn_Out : Mean non-synonymous divergence between outgroup sequence and population n°2\n\n" << endl;

	return 0;
}
string listName, format, typeAlg, outF, coding, pop1, pop2, outgroup;
double tvts;

/*************************************************************/
/** read command line **/
/*************************************************************/
int i = 1;
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
	if (s == "-coding")	{
		i++;
		if (i == argc) 	{
			cerr << "error in command: -coding <coding or non-coding>\n";
			cerr << '\n';
			exit(1);
		}
		coding = argv[i];
	}
	if (s == "-tvts")	{
		i++;
		if (i == argc) 	{
			cerr << "error in command: -tvts <double tvts>\n";
			cerr << '\n';
			exit(1);
		}
		tvts = TextTools::toDouble(argv[i]);
	}
	if (s == "-pop1")	{
		i++;
		if (i == argc) 	{
			cerr << "error in command: -pop1 <prefix pop1>\n";
			cerr << '\n';
			exit(1);
		}
		pop1 = argv[i];
	}
	if (s == "-pop2")	{
		i++;
		if (i == argc) 	{
			cerr << "error in command: -pop2 <prefix pop2>\n";
			cerr << '\n';
			exit(1);
		}
		pop2 = argv[i];
	}
	if (s == "-outgroup")	{
		i++;
		if (i == argc) 	{
			cerr << "error in command: -outgroup <prefix outgroup>\n";
			cerr << '\n';
			exit(1);
		}
		outgroup = argv[i];
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

if(!FileTools::fileExists( listName )) {
	cerr << "ERROR!!! File " << listName << " does not exists." << endl;
	exit(-1);
}

ifstream Filelist (listName.c_str(), ios::in);
ofstream Fileout (outF.c_str(), ios::out);

if(coding == "coding"){
	Fileout << "name\tsize\t";
	Fileout << "S_Pop1\tS_Pop2\tPi_Pop1\tPi_Pop2\tW_Pop1\tW_Pop2\tD_Pop1\tD_Pop2\tDxy\tPiTotal\tFstHud\tFstNei_weighted\tFstNei_unweighted\t";
	Fileout << "Ts_Pop1\tTs_Pop2\t";
	Fileout << "PS_Pop1\tPN_Pop1\tNSS_Pop1\tPS_Pop2\tPN_Pop2\tNSS_Pop2\t";
	Fileout << "gc3\tNumberOfStop\tInFrameStop\tsizeOut1\tsizeOut2\tdivOut1\tdivOut2\tNSS_Out1\tdiv_NonSyn_Out1\tdiv_Syn_Out1\tNSS_Out2\tdiv_NonSyn_Out2\tdiv_Syn_Out2" <<  endl;


}
if(coding == "non-coding"){
	Fileout << "name\tsize\t";
	Fileout << "S_Pop1\tS_Pop2\tPi_Pop1\tPi_Pop2\tW_Pop1\tW_Pop2\tD_Pop1\tD_Pop2\tDxy\tPiTotal\tFstHud\tFstNei_weighted\tFstNei_unweighted\t";
	Fileout << "gc\tsizeOut1\tsizeOut2\tdivOut1\tdivOut2\tFixed\tPrivatePop1\tPrivatePop2\tShared" <<  endl;
}

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
  	double FstHud, FstNei_w, FstNei_uw;
	const NucleicAlphabet * alpha = new DNA();
	const CodonAlphabet *codonAlpha = new CodonAlphabet(alpha);
	const GeneticCode *GC = new StandardGeneticCode(alpha);
	VectorSequenceContainer * seqCont = NULL;

	Phylip * PhySeq = new Phylip(true, true);
	Fasta * FstSeq = new  Fasta(10000000);
	
	if(format == "phylip"){
		seqCont = PhySeq->readAlignment(nomfic , alpha);
	}
	if(format == "fasta"){
		seqCont = FstSeq->readAlignment(nomfic , alpha);
	}

	VectorSiteContainer *sitec = new VectorSiteContainer( *seqCont );
	PolymorphismSequenceContainer * psc1 = new PolymorphismSequenceContainer( *sitec);
	
	/** set variables for statistic **/
	double diffOut1, diffOut2;
	unsigned int Size, S_Pop2, S_Pop1;
	unsigned int sizeSeqOut1 = 0;
	unsigned int sizeSeqOut2 = 0;
	double PiPop1, PiPop2, WPop1, WPop2, TajimaDPop1, TajimaDPop2, Dxy, piTotal;

	if(coding == "coding"){
	
		/*************************************************************/
		/** convert alphabet and remove stop codon **/
		/*************************************************************/
		PolymorphismSequenceContainer * pscCodon = new PolymorphismSequenceContainer(codonAlpha);
		SequenceContainerTools::convertAlphabet(*psc1, *pscCodon);	
		PolymorphismSequenceContainer * pscFinal = new PolymorphismSequenceContainer(pscCodon->getNumberOfSequences(), codonAlpha);
		pscFinal->setSequencesNames(psc1->getSequencesNames(), false);
	
		/** record stop codon **/
		unsigned int numberOfStopCodon = 0;
		string InFrameStop ="no";
		for(unsigned int i = 0; i < pscCodon->getNumberOfSites(); i++) {
		
			Site site = pscCodon->getSite(i);
		
			if(CodonSiteTools::hasStop(site, *GC)){
				numberOfStopCodon++;
				if(i < pscCodon->getNumberOfSites()-1){
					InFrameStop = "yes";
				}
				continue;
			}		
			pscFinal->addSite(site);
		}
	
		/** Split dataset in population **/
		string Pop1, Pop2, Out;
		Pop1 = Pop2 = Out = "no";

		for(unsigned int i = 0; i < pscFinal->getNumberOfSequences(); i ++){
		   
			if(TextTools::hasSubstring(pscFinal->getSequence(i).getName(), pop1)){
				Pop1 = "yes";
				pscFinal->setGroupId(i, 1);
				continue;
			}
			if(TextTools::hasSubstring(pscFinal->getSequence(i).getName(), pop2)){
				Pop2 = "yes";
				pscFinal->setGroupId(i, 2);
				continue;
			}
			if(TextTools::hasSubstring(pscFinal->getSequence(i).getName(), outgroup)){
				Out = "yes";
				pscFinal->setGroupId(i, 3);
				continue;
			}
		}
  
		PolymorphismSequenceContainer * pscPop1 = PolymorphismSequenceContainerTools::extractGroup (*pscFinal, 1);
		PolymorphismSequenceContainer * pscPop2 = PolymorphismSequenceContainerTools::extractGroup (*pscFinal, 2);

		/** create nucleotides alignment **/
		PolymorphismSequenceContainer * pscPop1Nuc = new PolymorphismSequenceContainer(alpha);
		SequenceContainerTools::convertAlphabet(*pscPop1, *pscPop1Nuc);
	
		PolymorphismSequenceContainer * pscPop2Nuc = new PolymorphismSequenceContainer(alpha);
		SequenceContainerTools::convertAlphabet(*pscPop2, *pscPop2Nuc);	
	
		PolymorphismSequenceContainer * pscFinalNuc = new PolymorphismSequenceContainer(alpha);
	
		for(unsigned int i = 0; i < pscFinal->getNumberOfSequences(); i ++){
			BasicSequence seq = BasicSequence(pscFinal->getSequence(i).getName(), pscFinal->getSequence(i).toString(), alpha);
			pscFinalNuc->addSequence(seq);
			pscFinalNuc->setGroupId(seq.getName(), pscFinal->getGroupId(i));
		}

		/******************************************************************/
		/** compute sequence statistic coding                            **/
		/******************************************************************/
  
		/** compute GC **/
		const Sequence * cons = SiteContainerTools::getConsensus(*pscPop2);

		BasicSequence seqPop2 = BasicSequence(cons->getName(), cons->toString(), alpha);
	
		double gc3 = getGCThirdCodonPosition(seqPop2);
	
		/** compute divergence **/
		PolymorphismSequenceContainer * pscDiv = new PolymorphismSequenceContainer(alpha);
		SequenceContainerTools::convertAlphabet(*pscPop2, *pscDiv);
		vector<double> diffNSOut1, diffNSOut2;
	
		/** compute number of complete site out1 and then divergence **/
		if(Out == "yes"){
			PolymorphismSequenceContainer * pscOut = PolymorphismSequenceContainerTools::extractGroup (*pscFinal, 3);
			PolymorphismSequenceContainer * pscOutNuc = new PolymorphismSequenceContainer(alpha);
			SequenceContainerTools::convertAlphabet(*pscOut, *pscOutNuc);
			const Sequence & seqOut1Codon = pscOut->getSequence(0);
			const Sequence & seqOut2Codon = pscOut->getSequence(1);

			BasicSequence seqOut1 = BasicSequence(seqOut1Codon.getName(), seqOut1Codon.toString(), alpha);
			BasicSequence seqOut2 = BasicSequence(seqOut2Codon.getName(), seqOut2Codon.toString(), alpha);
			
			PolymorphismSequenceContainer * tmpSc = new PolymorphismSequenceContainer(seqOut1.getAlphabet());
			tmpSc->addSequence(seqOut1);
			SiteContainer * tmpComp = SiteContainerTools::getCompleteSites(*tmpSc);
			sizeSeqOut1 = tmpComp->getNumberOfSites();
	
			if( !isSeqGapOrUnresolvedOnly(seqOut1) && sizeSeqOut1 > 1 ){
				diffOut1 = computeMeanNumberOfDifference(*pscDiv, seqOut1, 1);
			}else{
				diffOut1 = -999;
			}

			/** compute number of complete site out2 and then divergence **/
			tmpSc->clear();
			tmpComp->clear();
			tmpSc->addSequence(seqOut2);
			tmpComp = SiteContainerTools::getCompleteSites(*tmpSc);
			sizeSeqOut2 = tmpComp->getNumberOfSites();
			delete tmpComp;
	
	
			if( !isSeqGapOrUnresolvedOnly(seqOut2) && sizeSeqOut2 > 1){
				diffOut2 = computeMeanNumberOfDifference(*pscDiv, seqOut2, 1);
			}else{
				diffOut2 = -999;
			}
	
			/** compute number of complete codon in out1 and then synonymous and non-synonymous divergence **/
			PolymorphismSequenceContainer * tmpScCodon1 = new PolymorphismSequenceContainer(seqOut1Codon.getAlphabet());
			tmpScCodon1->addSequence(seqOut1Codon);
			SiteContainer * tmpCompCodon1 = SiteContainerTools::getCompleteSites(*tmpScCodon1);
			unsigned int sizeSeqOutCodon1 = tmpCompCodon1->getNumberOfSites();
			delete tmpScCodon1;
			delete tmpCompCodon1;
	
			if( !isSeqGapOrUnresolvedOnly(seqOut1) && sizeSeqOutCodon1 > 1 ){
				diffNSOut1 = computeMeanSynonymousAndNonSynonymousDifferences(*pscPop2, seqOut1Codon, *GC, tvts);
			}else{
				diffNSOut1.push_back(-999);
				diffNSOut1.push_back(-999);
				diffNSOut1.push_back(-999);
			}
		
			/** compute number of complete codon in out2 and then synonymous and non-synonymous divergence **/
			PolymorphismSequenceContainer * tmpScCodon = new PolymorphismSequenceContainer(seqOut2Codon.getAlphabet());
			tmpScCodon->addSequence(seqOut2Codon);
			SiteContainer * tmpCompCodon = SiteContainerTools::getCompleteSites(*tmpScCodon);
			unsigned int sizeSeqOutCodon2 = tmpCompCodon->getNumberOfSites();
			delete tmpScCodon;
			delete tmpCompCodon;
	
			if( !isSeqGapOrUnresolvedOnly(seqOut2) && sizeSeqOutCodon2 > 1 ){
				diffNSOut2 = computeMeanSynonymousAndNonSynonymousDifferences(*pscPop2, seqOut2Codon, *GC, tvts);
			}else{
				diffNSOut2.push_back(-999);
				diffNSOut2.push_back(-999);
				diffNSOut2.push_back(-999);
			}
			delete pscOut;
			delete pscOutNuc;
		}else{ // if no outgroup
			diffNSOut1.push_back(-999);
			diffNSOut1.push_back(-999);
			diffNSOut1.push_back(-999);
			diffNSOut2.push_back(-999);
			diffNSOut2.push_back(-999);
			diffNSOut2.push_back(-999);
		}
		/** ************************************************ **/
		/** ***** Compute Pop Gen Statistic *******************/
		/** ************************************************ **/
		Size = pscPop2Nuc->getNumberOfSites();
		S_Pop1 = SequenceStatistics::numberOfPolymorphicSites( *pscPop1Nuc, false );
		S_Pop2 = SequenceStatistics::numberOfPolymorphicSites( *pscPop2Nuc, false );
		
		PiPop1 = SequenceStatistics::tajima83(*pscPop1Nuc, false);
		PiPop2 = SequenceStatistics::tajima83(*pscPop2Nuc, false);
		
		WPop1 = SequenceStatistics::watterson75(*pscPop1Nuc, false );
		WPop2 = SequenceStatistics::watterson75(*pscPop2Nuc, false );

		if(S_Pop1 > 0)
			TajimaDPop1 =  SequenceStatistics::tajimaDss( *pscPop1Nuc, false );
		else
			TajimaDPop1 = -999.0;
		if(S_Pop2 > 0)
			TajimaDPop2 =  SequenceStatistics::tajimaDss( *pscPop2Nuc, false );
		else
			TajimaDPop2 = -999.0;
	

	
		if(S_Pop1 == 0 && S_Pop2 == 0){
			FstHud = -999.0;
			FstNei_w = -999.0;
			FstNei_uw = -999.0;
			Dxy = 0.0;
			piTotal = 0.0;
		}else{
			FstHud = FstHudson92(*psc1, 1, 2, false);
			FstNei_w = FstNei82(*psc1, 1, 2, true, false);
			FstNei_uw = FstNei82(*psc1, 1, 2, false, false);
			Dxy = PiInter(*pscFinalNuc, 1, 2);
			piTotal = SequenceStatistics::tajima83(*psc1, false);
		}

	
		Fileout << nomfic << "\t" << Size << "\t" << S_Pop1 << "\t" << S_Pop2 << "\t" << PiPop1 << "\t" << PiPop2 << "\t" << WPop1 << "\t" << WPop2 << "\t"; 
		Fileout << TajimaDPop1 << "\t" << TajimaDPop2 << "\t";
		Fileout <<  Dxy << "\t" << FstHud << "\t" << FstNei_w << "\t" << FstNei_uw << "\t";
		Fileout << getNumberOfTransitions2(*pscPop1Nuc) << "\t" << getNumberOfTransitions2(*pscPop2Nuc) << "\t" ;
		Fileout << piSynonymous2( *pscPop1,*GC, false, false)  << "\t" << piNonSynonymous2( *pscPop1, *GC, false, false) << "\t" << meanSynonymousSitesNumber2( *pscPop1, *GC, tvts, false) << "\t";
		Fileout << piSynonymous2( *pscPop2,*GC, false, false)  << "\t" << piNonSynonymous2( *pscPop2, *GC, false, false) << "\t" << meanSynonymousSitesNumber2( *pscPop2, *GC, tvts, false) << "\t";
		Fileout << gc3 << "\t" << numberOfStopCodon  << "\t" << InFrameStop << "\t" << sizeSeqOut1 << "\t" << sizeSeqOut2 << "\t"<< diffOut1 << "\t" << diffOut2 << "\t";
		Fileout << diffNSOut1[2] << "\t" << diffNSOut1[0] << "\t" << diffNSOut1[1] << "\t";
		Fileout << diffNSOut2[2] << "\t" << diffNSOut2[0] << "\t" << diffNSOut2[1] << endl;
		
		delete pscCodon;
		delete pscPop2;
		delete pscPop1;
		delete pscDiv;
		delete pscFinal;
		delete pscPop1Nuc;
		delete pscPop2Nuc;
		delete pscFinalNuc;
	
	} // end of coding statistic
	
	if(coding == "non-coding"){
		cout << "Non-coding DNA" << endl;
		/** Split dataset in population **/
		string Pop1, Pop2, Out;
		Pop1 = Pop2 = Out = "no";

		for(unsigned int i = 0; i < psc1->getNumberOfSequences(); i ++){
		   
			if(TextTools::hasSubstring(psc1->getSequence(i).getName(), pop1)){
				Pop1 = "yes";
				psc1->setGroupId(i, 1);
				continue;
			}
			if(TextTools::hasSubstring(psc1->getSequence(i).getName(), pop2)){
				Pop2 = "yes";
				psc1->setGroupId(i, 2);
				continue;
			}
			if(TextTools::hasSubstring(psc1->getSequence(i).getName(), outgroup)){
				Out = "yes";
				psc1->setGroupId(i,3);
				continue;
			}
		}
		
		PolymorphismSequenceContainer * pscPop1 = PolymorphismSequenceContainerTools::extractGroup (*psc1, 1);
		PolymorphismSequenceContainer * pscPop2 = PolymorphismSequenceContainerTools::extractGroup (*psc1, 2);

		PolymorphismSequenceContainer * pscIngroup = new PolymorphismSequenceContainer(alpha);
	
		for(unsigned int i = 0; i < psc1->getNumberOfSequences(); i ++){
			if(psc1->getGroupId(i) == 1 || psc1->getGroupId(i) == 2){
				BasicSequence seq = BasicSequence(psc1->getSequence(i).getName(), psc1->getSequence(i).toString(), alpha);
				pscIngroup->addSequence(seq);
				pscIngroup->setGroupId(seq.getName(), psc1->getGroupId(i));
			}
		}

		/** ************************************************ **/
		/** ***** Compute Pop Gen Statistic *******************/
		/** ************************************************ **/
		
		double gc = SequenceStatistics::gcContent(*pscPop2);
		Size = pscPop2->getNumberOfSites();
		S_Pop2 = SequenceStatistics::numberOfPolymorphicSites( *pscPop2, false );
		S_Pop1 = SequenceStatistics::numberOfPolymorphicSites( *pscPop1, false );
		PiPop1 = SequenceStatistics::tajima83(*pscPop1, false);
		PiPop2 = SequenceStatistics::tajima83(*pscPop2, false);
		WPop1 = SequenceStatistics::watterson75(*pscPop1, false );
		WPop2 = SequenceStatistics::watterson75(*pscPop2, false );
		Dxy = SequenceStatistics::tajima83(*pscIngroup, false);
	

		if(S_Pop1 > 0)
			TajimaDPop1 =  SequenceStatistics::tajimaDss( *pscPop1, false );
		else
			TajimaDPop1 = -999.0;
		if(S_Pop2 > 0)
			TajimaDPop2 =  SequenceStatistics::tajimaDss( *pscPop2, false );
		else
			TajimaDPop2 = -999.0;

		if(S_Pop1 == 0 && S_Pop2 == 0){
			FstHud = -999.0;
			FstNei_w = -999.0;
			FstNei_uw = -999.0;
			Dxy = 0.0;
			piTotal = 0.0;
		}else{
			FstHud = FstHudson92(*psc1, 1, 2, false);
			FstNei_w = FstNei82(*psc1, 1, 2, true, false);
			FstNei_uw = FstNei82(*psc1, 1, 2, false, false);
			Dxy = PiInter(*psc1, 1, 2);
			piTotal = SequenceStatistics::tajima83(*psc1, false);
		}
		
		/** Outgroup **/
		diffOut1 = diffOut2 = -999;
		if(Out == "yes"){
			PolymorphismSequenceContainer * pscOut = PolymorphismSequenceContainerTools::extractGroup (*psc1, 3);
			const Sequence & seqOut1 = pscOut->getSequence(0);
			const Sequence & seqOut2 = pscOut->getSequence(1);
	
			sizeSeqOut1 = SequenceTools::getNumberOfCompleteSites(seqOut1);
			sizeSeqOut2 = SequenceTools::getNumberOfCompleteSites(seqOut2);
	
			if( !isSeqGapOrUnresolvedOnly(seqOut1) && sizeSeqOut1 > 1 ){
				diffOut1 = computeMeanNumberOfDifference(*pscPop2, seqOut1, 1);
			}else{
				diffOut1 = -999;
			}
	
			if( !isSeqGapOrUnresolvedOnly(seqOut2) && sizeSeqOut2 > 1 ){
				diffOut2 = computeMeanNumberOfDifference(*pscPop2, seqOut2, 1);
			}else{
				diffOut2 = -999;
			}
			delete pscOut;
		}
		
		vector<int> SharedFixed = NumberOfDifferenceBetweenPopulations(*psc1, 1, 2);
		
		Fileout << nomfic << "\t" << Size << "\t" << S_Pop1 << "\t" << S_Pop2 << "\t";
		Fileout << PiPop1 << "\t" << PiPop2 << "\t" << WPop1 << "\t" << WPop2 << "\t"; 
		Fileout << TajimaDPop1 << "\t" << TajimaDPop2 << "\t";
		Fileout << Dxy << "\t" << piTotal << "\t";
		Fileout << FstHud << "\t" << FstNei_w  << "\t" << FstNei_uw << "\t";
		Fileout << gc << "\t" << sizeSeqOut1 << "\t" << sizeSeqOut2 << "\t" <<  diffOut1 << "\t" << diffOut2 << "\t";
		Fileout << SharedFixed[0] << "\t" << SharedFixed[1] << "\t"  <<SharedFixed[2] << "\t" << SharedFixed[3] << endl;
		
		delete pscPop1;
		delete pscPop2;

	}
	
	delete FstSeq;
	delete PhySeq;
	delete GC;	
	delete alpha;
	delete codonAlpha;
	delete sitec;
	delete seqCont;
	delete psc1;
	
}

}
catch(exception & e){
  cout << e.what() << endl;
    }
return 0;
}
