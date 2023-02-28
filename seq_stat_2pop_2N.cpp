/**

g++ --static -std=c++14 -g ~/pCloudDrive/Progr/C++/popgen/seq_stat_2pop/seq_stat_2pop_2N.cpp -o ~/bin/seq_stat_2pop_2N
 -DVIRTUAL_COV=yes -Wall \
 -I$HOME/local/bpp/dev/include  -L$HOME/local/bpp/dev/lib \
 -lbpp-phyl -lbpp-popgen -lbpp-seq -lbpp-core 
 
strip ~/bin/seq_stat_2pop_2N


## bio++ V3
export LIBRARY_PATH=$HOME/local/bpp/dev/lib
g++ -static -static-libgcc -std=c++14 -Wall -Wshadow -Wconversion  -g ~/pCloudDrive/Progr/C++/popgen/seq_stat_2pop/seq_stat_2pop_2N.cpp -o ~/bin/seq_stat_2pop_2N  -I$HOME/local/bpp/dev/include -L$HOME/local/bpp/dev/lib  -DVIRTUAL_COV=yes  -lbpp-popgen3  -lbpp-phyl3 -lbpp-seq3 -lbpp-core3

Population genetics statistic computed with single diploid genome : For Roadkills genome project

* 
**/


#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

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
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Seq/Container/SiteContainerIterator.h>
#include <Bpp/PopGen/SequenceStatistics.h>
#include <Bpp/PopGen/PolymorphismSequenceContainerTools.h>
#include <Bpp/PopGen/PolymorphismSequenceContainer.h>

#include <Bpp/Numeric/VectorTools.h>


using namespace bpp;


/********************************************************************/
/** Number of Fixed, Private and Shared polymorphism between population *******/ 
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
if (argc == 1 || argc < 8)
{
	cout << "\n#####################\n### Usage : ";
	cout << "\nStat_2Pop_2N -seq [listSeq] -f [phylip or fasta]-pop1 [prefix_pop1] -pop2 [prefix_pop2] -o [out file]\n" << endl;
	

	return 0;
}
string listName, format, typeAlg, outF, coding, pop1, pop2, outgroup;
//double tvts;

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

/** set variables for statistic **/
unsigned int Size, S_Pop2, S_Pop1;
double PiTot, PiPop1, PiPop2;

Fileout << "name	size	Fixed	PrivatePop1	PrivatePop2	Shared	Pi1	Pi2	PiTot"<<  endl;

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
	const CodonAlphabet *GC = new StandardGeneticCode(alpha);
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
	
	/** exclude all 'N' **/
	SiteContainer * siteContComplete = SiteContainerTools::getCompleteSites(*sitec);
	
	PolymorphismSequenceContainer * psc1 = new PolymorphismSequenceContainer( *siteContComplete);

	/** Split dataset in population **/
	string Pop1, Pop2;
	Pop1 = Pop2 ="no";
	vector<string> seq2delete;
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
		// remove unidentified sequenced
		seq2delete.push_back(psc1->getSequence(i).getName());
	}
	
	for(unsigned int i = 0; i < seq2delete.size(); i ++){
		psc1->deleteSequence(seq2delete[i]);
	}
	
	PolymorphismSequenceContainer * pscPop1 = PolymorphismSequenceContainerTools::extractGroup (*psc1, 1);
	PolymorphismSequenceContainer * pscPop2 = PolymorphismSequenceContainerTools::extractGroup (*psc1, 2);

	/** ************************************************ **/
	/** ***** Compute Pop Gen Statistic *******************/
	/** ************************************************ **/
	
	
	Size = psc1->getNumberOfSites();
	// S_Pop1 = SequenceStatistics::numberOfPolymorphicSites( *pscPop1, true, true );
	// S_Pop2 = SequenceStatistics::numberOfPolymorphicSites( *pscPop2, true, true );
	vector<int> SharedFixed = NumberOfDifferenceBetweenPopulations(*psc1, 1, 2);
	
	PiPop1 = SequenceStatistics::tajima83(*pscPop1, false);
	PiPop2 = SequenceStatistics::tajima83(*pscPop2, false);
	PiTot = SequenceStatistics::tajima83(*psc1, false);
	
	Fileout << nomfic << "\t" << Size <<  "\t";
	Fileout << SharedFixed[0] << "\t" << SharedFixed[1] << "\t"  <<SharedFixed[2] << "\t" << SharedFixed[3] << "\t" << PiPop1 << "\t" << PiPop2 << "\t" << PiTot <<  endl;
	
	delete pscPop1;
	delete pscPop2;
	delete FstSeq;
	delete PhySeq;
	delete alpha;
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
