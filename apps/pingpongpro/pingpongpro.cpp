// ==========================================================================
//				pingpongpro
// ==========================================================================
// todo: copyright

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/arg_parse.h>

#include <map>
#include <vector>
#include <list>
#include <ctime>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>

using namespace std;
using namespace seqan;

// ==========================================================================
// Types & Classes
// ==========================================================================

#if defined(WIN32) || defined(_WIN32) 
#define PATH_delimiter '\\' 
#else 
#define PATH_delimiter '/'
#endif

// type to hold the list of input files given as arguments to the program
typedef vector<CharString> TInputFiles;

// todo: description
enum TCountMultiHits { multiHitsWeighted, multiHitsDiscard, multiHitsUnique };

// todo: description
enum TFileFormat { fileFormatBED, fileFormatCSV, fileFormatGFF, fileFormatGTF, fileFormatTSV };

// struct to store the options from the command line
struct AppOptions
{
	bool bedGraph;
	TInputFiles inputFiles;
	unsigned int minAlignmentLength;
	unsigned int maxAlignmentLength;
	unsigned int minStackHeight;
	TCountMultiHits countMultiHits;
	CharString output;
	bool plot;
	TInputFiles transposonFiles;
	unsigned int verbosity;
};

// types to store @SQ header lines of BAM/SAM files
typedef StringSet<CharString> TNameStore;
typedef Iterator<TNameStore>::Type TNameStoreIterator;

// constants to refer to + and - strands throughout the program
const unsigned int STRAND_PLUS = 0;
const unsigned int STRAND_MINUS = 1;

// for every locus (position) on the genome the following attributes are calculated:
//  - reads: the number of reads which begin at this position
//  - UAt5PrimeEnd: whether the reads of the stack have a U at the 5' end
struct TCountsPosition
{
	float reads;
	bool UAt5PrimeEnd;
	TCountsPosition():
		reads(0),
		UAt5PrimeEnd(false)
	{}
};

// The following types define nested arrays to store the above stats for every position in the genome.
// The stats are grouped by strand and contig/chromosome.
typedef map< unsigned int, TCountsPosition > TCountsContig;
typedef map< unsigned int, TCountsContig > TCountsStrand;
typedef TCountsStrand TCountsGenome[2];

// true ping-pong stacks overlap by this many nt
#define PING_PONG_OVERLAP 10

// stacks with overlaps between MIN_ARBITRARY_OVERLAP and MAX_ARBITRARY_OVERLAP (except for PING_PONG_OVERLAP)
// are used to estimate what is background noise
#define MIN_ARBITRARY_OVERLAP 3
#define MAX_ARBITRARY_OVERLAP 23

// type to store a score for each stack height found in the input file
typedef map< unsigned int, float > THeightScoreMap;

// todo: description
#define IS_URIDINE 0
#define IS_NOT_URIDINE 1
#define IS_ABOVE_COVERAGE 0
#define IS_BELOW_COVERAGE 1
#define HEIGHT_SCORE_BINS 1000
typedef vector< vector< vector< vector< float > > > > TGroupedStackCounts;
typedef vector< TGroupedStackCounts > TGroupedStackCountsByOverlap;

// todo: description
struct TPingPongSignature
{
	unsigned int position: 32; // position on contig where the ping-pong overlap is located
	unsigned int heightScoreBin: 29; // must be big enough to hold <HEIGHT_SCORE_BINS>
	unsigned int localHeightScoreBin: 1; // holds either <IS_ABOVE_COVERAGE> or <IS_BELOW_COVERAGE>
	unsigned int UAt5PrimeEndOnPlusStrandBin: 1; // holds either <IS_URIDINE> or <IS_NOT_URIDINE>
	unsigned int UAt5PrimeEndOnMinusStrandBin: 1; // holds either <IS_URIDINE> or <IS_NOT_URIDINE>
	unsigned int readsOnPlusStrand: 32; // stack height on + strand
	unsigned int readsOnMinusStrand: 32; // stack height on - strand
	float fdr; // chances of this putative ping-pong overlap being a false discovery

	// constructor to initialize with values
	TPingPongSignature(unsigned int position, unsigned int heightScoreBin, unsigned int localHeightScoreBin, unsigned int UAt5PrimeEndOnPlusStrandBin, unsigned int UAt5PrimeEndOnMinusStrandBin, unsigned int readsOnPlusStrand, unsigned int readsOnMinusStrand):
		position(position), heightScoreBin(heightScoreBin), localHeightScoreBin(localHeightScoreBin), UAt5PrimeEndOnPlusStrandBin(UAt5PrimeEndOnPlusStrandBin), UAt5PrimeEndOnMinusStrandBin(UAt5PrimeEndOnMinusStrandBin), readsOnPlusStrand(readsOnPlusStrand), readsOnMinusStrand(readsOnMinusStrand), fdr(0)
	{
	}

	// copy constructor (needed to add instance to a STL list)
	TPingPongSignature(const TPingPongSignature &pingPongSignature)
	{
		position = pingPongSignature.position;
		heightScoreBin = pingPongSignature.heightScoreBin;
		localHeightScoreBin = pingPongSignature.localHeightScoreBin;
		UAt5PrimeEndOnPlusStrandBin = pingPongSignature.UAt5PrimeEndOnPlusStrandBin;
		UAt5PrimeEndOnMinusStrandBin = pingPongSignature.UAt5PrimeEndOnMinusStrandBin;
		readsOnPlusStrand = pingPongSignature.readsOnPlusStrand;
		readsOnMinusStrand = pingPongSignature.readsOnMinusStrand;
		fdr = pingPongSignature.fdr;
	}
};
// type to store all ping-pong signatures of a single contig
typedef list< TPingPongSignature > TPingPongSignaturesPerContig;
// type to store all ping-pong signatures of the entire genome
typedef map< unsigned int, TPingPongSignaturesPerContig > TPingPongSignaturesPerGenome;
// type to store all ping-pong signatures with a certain overlap (including overlaps other than 10 nt, i.e., no real ping-pong signatures)
typedef vector< TPingPongSignaturesPerGenome > TPingPongSignaturesByOverlap;

// type to store the region of a single transposon
struct TTransposon
{
	string identifier;
	unsigned int strand;
	unsigned int start;
	unsigned int end;
	float pValue;
	float qValue;
	float normalizedSignatureCount;

	// constructor to initialize with values
	TTransposon(string identifier, unsigned int strand, unsigned int start, unsigned int end):
		identifier(identifier), strand(strand), start(start), end(end), pValue(1), qValue(1), normalizedSignatureCount(0)
	{
	}

	// copy constructor (needed to add instance to a STL list
	TTransposon(const TTransposon &transposon)
	{
		identifier = transposon.identifier;
		strand = transposon.strand;
		start = transposon.start;
		end = transposon.end;
		pValue = transposon.pValue;
		qValue = transposon.qValue;
		normalizedSignatureCount = transposon.normalizedSignatureCount;
	}

	// needed for sorting of the list by genomic position
	inline bool operator< (const TTransposon &transposon)
	{
		if (start == transposon.start)
			return end < transposon.end;
		else
			return start < transposon.start;
	}
};
// type to store all transposons of a single contig
typedef list< TTransposon > TTransposonsPerContig;
// type to store all transposons of the entire genome
typedef map< unsigned int, TTransposonsPerContig > TTransposonsPerGenome;

// ==========================================================================
// Functions
// ==========================================================================

// function to parse command-line arguments
ArgumentParser::ParseResult parseCommandLine(AppOptions &options, int argc, char const ** argv)
{
	ArgumentParser parser("pingpongpro");

	// define usage and description
	addUsageLine(parser, "[\\fIOPTIONS\\fP] [-i \\fISAM_INPUT_FILE\\fP [-i ...]] [-o \\fIOUTPUT_DIRECTORY\\fP]");
	setShortDescription(parser, "Find ping-pong signatures like a pro");
	// todo: define long description
	addDescription(parser, "PingPongPro scans piRNA-Seq data for signs of ping-pong cycle activity. The ping-pong cycle produces piRNA molecules with complementary 5'-ends. These molecules appear as stacks of aligned reads whose 5'-ends overlap with the 5'-ends of reads on the opposite strand by exactly 10 bases.");
	setVersion(parser, "0.1");
	setDate(parser, "Mar 2014");

	addOption(parser, ArgParseOption("b", "bedgraph", "Output loci with ping-pong signature in bedGraph format. Default: off."));

	addOption(parser, ArgParseOption("s", "min-stack-height", "Omit stacks with fewer than the specified number of reads from the output.", ArgParseArgument::INTEGER, "NUMBER_OF_READS", true));
	setDefaultValue(parser, "min-stack-height", 0);
	setMinValue(parser, "min-stack-height", "0");

	addOption(parser, ArgParseOption("i", "input", "Input file(s) in SAM/BAM format. \"-\" means stdin.", ArgParseArgument::INPUTFILE, "PATH", true));
	setDefaultValue(parser, "input", "-");
	setValidValues(parser, "input", ".bam .sam -");

	addOption(parser, ArgParseOption("l", "min-alignment-length", "Ignore alignments in the input file that are shorter than the specified length.", ArgParseArgument::INTEGER, "LENGTH", true));
	setDefaultValue(parser, "min-alignment-length", 24);
	setMinValue(parser, "min-alignment-length", "1");
	addOption(parser, ArgParseOption("L", "max-alignment-length", "Ignore alignments in the input file that are longer than the specified length.", ArgParseArgument::INTEGER, "LENGTH", true));
	setDefaultValue(parser, "max-alignment-length", 32);
	setMinValue(parser, "max-alignment-length", "1");

	addOption(parser, ArgParseOption("m", "multi-hits", "How to count multi-mapping reads.", ArgParseArgument::STRING, "METHOD", true));
	setDefaultValue(parser, "multi-hits", "weighted");
	setValidValues(parser, "multi-hits", "weighted discard unique");

	addOption(parser, ArgParseOption("o", "output", "Write output to specified directory. Default: current working directory.", ArgParseArgument::OUTPUTFILE, "PATH", true));

	addOption(parser, ArgParseOption("p", "plot", "Generate R plots on background noise estimation. Requires Rscript. Default: off."));

	addOption(parser, ArgParseOption("t", "transposons", "Check if the transposons given in the file are suppressed through ping-pong activity. \"-\" means stdin.", ArgParseArgument::INPUTFILE, "PATH", true));
	setValidValues(parser, "transposons", ".bed .csv .gff .gtf .tsv");

	addOption(parser, ArgParseOption("v", "verbose", "Print messages to stderr about the current progress. Default: off."));

	// parse command line
	ArgumentParser::ParseResult parserResult = parse(parser, argc, argv);
	if (parserResult != ArgumentParser::PARSE_OK)
		return parserResult;

	// extract options, if parsing was successful
	options.bedGraph = isSet(parser, "bedgraph");
	options.inputFiles.resize(getOptionValueCount(parser, "input")); // store input files in vector
	if (options.inputFiles.size() == 0)
	{
		options.inputFiles.push_back("/dev/stdin"); // read from stdin, if no input file is given
	}
	else
	{
		for (vector< string >::size_type i = 0; i < options.inputFiles.size(); i++)
		{
			getOptionValue(options.inputFiles[i], parser, "input", i);
			if (options.inputFiles[i] == "-")
				options.inputFiles[i] = "/dev/stdin";
		}
	}
	string countMultiHits;
	getOptionValue(countMultiHits, parser, "multi-hits");
	if (countMultiHits == "unique")
	{
		options.countMultiHits = multiHitsUnique;
	}
	else if (countMultiHits == "discard")
	{
		options.countMultiHits = multiHitsDiscard;
	}
	else
	{
		options.countMultiHits = multiHitsWeighted;
	}
	getOptionValue(options.minStackHeight, parser, "min-stack-height");
	getOptionValue(options.minAlignmentLength, parser, "min-alignment-length");
	getOptionValue(options.maxAlignmentLength, parser, "max-alignment-length");
	if (options.minAlignmentLength > options.maxAlignmentLength)
	{
		cerr << getAppName(parser) << ": maximum alignment length (" << options.maxAlignmentLength << ") must not be lower than minimum alignment length (" << options.minAlignmentLength << ")" << endl;
		return ArgumentParser::PARSE_ERROR;
	}
	getOptionValue(options.output, parser, "output");
	if ((length(options.output) > 0) && (options.output[length(options.output)-1] != PATH_delimiter))
		options.output += PATH_delimiter; // append slash to output path, if missing
	options.plot = isSet(parser, "plot");
	options.transposonFiles.resize(getOptionValueCount(parser, "transposons")); // store input files in vector
	for (vector< string >::size_type i = 0; i < options.transposonFiles.size(); i++)
		getOptionValue(options.transposonFiles[i], parser, "transposons", i);
	if (isSet(parser, "verbose"))
	{
		options.verbosity = 3;
	}
	else
	{
		options.verbosity = 0;
	}

	return parserResult;
}

// function to measure time between the first and second invocation of the function
unsigned int stopwatch(const string &operation, unsigned int verbosity)
{
	static time_t start = 0;
	unsigned int elapsedSeconds = 0;
	if (start != 0)
	{
		elapsedSeconds = time(NULL) - start;
		start = 0;
	}
	else
	{
		start = time(NULL);
	}
	if (verbosity >= 3)
	{
		if (operation != "")
			cerr << operation << " ... ";
		else
			cerr << "done (" << elapsedSeconds << " seconds)" << endl;
	}
	cerr.flush();
	return elapsedSeconds;
}
unsigned int stopwatch(unsigned int verbosity)
{
	return stopwatch("", verbosity);
}

// Function which sums up the number of reads that start at a given position in the genome.
// Additionally, it counts the number of reads with uridine at the 5' end.
// Parameters:
//   bamFile: the BAM/SAM file from where to load the reads
//   readCounts: stats for positions were reads on the minus strand overlap the 5' ends of reads on the plus strand
// todo: document parameters
int countReadsInBamFile(BamStream &bamFile, TCountsGenome &readCounts, const unsigned int minAlignmentLength, const unsigned int maxAlignmentLength, TCountMultiHits countMultiHits, float &totalReadCount)
{
	TCountsPosition *position;

	BamAlignmentRecord record;
	while (!atEnd(bamFile))
	{
		if (readRecord(record, bamFile) != 0)
		{
			cerr << "Failed to read record" << endl;
			return 1;
		}

		if ((record.beginPos != BamAlignmentRecord::INVALID_POS) && (record.beginPos != -1)) // skip unmapped reads
		{
			// calculate length of alignment using CIGAR string
			size_t alignmentLength = 0;
			for (unsigned int cigarIndex = 0; cigarIndex < length(record.cigar); ++cigarIndex)
			{
				if ((record.cigar[cigarIndex].operation == 'M') || (record.cigar[cigarIndex].operation == 'N') || (record.cigar[cigarIndex].operation == 'D') || (record.cigar[cigarIndex].operation == '=') || (record.cigar[cigarIndex].operation == 'X')) // these CIGAR elements indicate alignment
					alignmentLength += record.cigar[cigarIndex].count;
			}

			// skip read, if alignment is too long or too short
			if ((alignmentLength < minAlignmentLength) || (alignmentLength > maxAlignmentLength))
				continue;

			// the stack height is increased by the value of this variable (depends on how multi-hits are handled)
			float readWeight;

			if (countMultiHits == multiHitsUnique) // we do not distinguish between multi-hits and unique hits
			{
				// increase stack height by 1, regardless of whether it is a multi-hit or unique hit
				readWeight = 1;
			}
			else
			{
				// check if the record is a multi-hit by examining the optional tags
				BamTagsDict tagsDictionary(record.tags);
				unsigned int tagIndex;
				unsigned int multiHits = 1;
				if (findTagKey(tagIndex, tagsDictionary, "NH"))
					extractTagValue(multiHits, tagsDictionary, tagIndex);

				if (countMultiHits == multiHitsWeighted)
				{
					// increase stack height by fraction
					readWeight = 1.0 / multiHits;
				}
				else /*if (countMultiHits == multiHitsDiscard)*/
				{
					// discard read (i.e., set readWeight to 0), if there is more than 1 instance in the SAM file
					readWeight = (multiHits == 1) ? 1 : 0;
				}
			}

			if (readWeight <= 0)
				continue; // skip to next read, if read is to be discarded

			if (hasFlagRC(record)) // read maps to minus strand
			{
				// get a pointer to counter of the position of the read
				position = &(readCounts[STRAND_MINUS][record.rID][record.beginPos+alignmentLength]);

				// check if base at 5' end is uridine
				size_t clippedBasesAt5PrimeEnd = 0;
				if ((length(record.cigar) > 1) && (record.cigar[length(record.cigar)-1].operation == 'S'))
					clippedBasesAt5PrimeEnd = record.cigar[length(record.cigar)-1].count;
				if ((record.seq[length(record.seq)-clippedBasesAt5PrimeEnd-1] == 'A') || (record.seq[length(record.seq)-clippedBasesAt5PrimeEnd-1] == 'a')) // check if last base is uridine (we check for adenine, because reads on the - strand are stored as the complement in SAM files
					position->UAt5PrimeEnd = true;
			}
			else // read maps to plus strand
			{
				// get a pointer to counter of the position of the read
				position = &(readCounts[STRAND_PLUS][record.rID][record.beginPos]);

				// check if base at 5' end is uridine
				size_t clippedBasesAt5PrimeEnd = 0;
				if (record.cigar[0].operation == 'S')
					clippedBasesAt5PrimeEnd = record.cigar[0].count;
				if ((record.seq[clippedBasesAt5PrimeEnd] == 'T') || (record.seq[clippedBasesAt5PrimeEnd] == 't'))
					position->UAt5PrimeEnd = true;
			}

			// increase stack height
			position->reads += readWeight;

			// increase read counter
			totalReadCount += readWeight;
		}
	}

	return 0;
}

// convert stack heights to scores
// todo: document parameters
void mapHeightsToScores(TCountsGenome &readCounts, THeightScoreMap &heightScoreMap)
{
	// iterate through all strands, contigs and positions to count how many stacks there are of any given height
	for (unsigned int strand = STRAND_PLUS; strand <= STRAND_MINUS; ++strand)
		for (TCountsStrand::iterator contig = readCounts[strand].begin(); contig != readCounts[strand].end(); ++contig)
			for (TCountsContig::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
				heightScoreMap[0.5 + position->second.reads] += 1;
}

// todo: description
void countStacksByGroup(TCountsGenome &readCounts, THeightScoreMap &heightScoreMap, TGroupedStackCountsByOverlap &groupedStackCountsByOverlap, TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap)
{
	// the following loop initializes a multi-dimensional array of stack counts with the following boundaries:
	// MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1 (one for each possible overlap)
	// HEIGHT_SCORE_BINS (one of each bin of the height scores)
	// 2 (one for reads with uridine at the 5' end of reads on the + strand and one for those with a different base)
	// 2 (one for reads with uridine at the 5' end of reads on the - strand and one for those with a different base)
	// 2 (one for stack heights below the local coverage and one for stack heights above)
	groupedStackCountsByOverlap.resize(MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1);
	for (TGroupedStackCountsByOverlap::iterator i = groupedStackCountsByOverlap.begin(); i != groupedStackCountsByOverlap.end(); ++i)
	{
		i->resize(HEIGHT_SCORE_BINS);
		for (TGroupedStackCounts::iterator j = i->begin(); j != i->end(); ++j)
		{
			j->resize(2);
			for (vector< vector< vector< float > > >::iterator k = j->begin(); k != j->end(); ++k)
			{
				k->resize(2);
				for (vector< vector< float > >::iterator l = k->begin(); l != k->end(); ++l)
					l->resize(2);
			}
		}
	}

	pingPongSignaturesByOverlap.resize(MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1);

	// find the highest possible score that two overlapping stacks can get
	float maxHeightScore = 0;
	for (THeightScoreMap::iterator heightScore = heightScoreMap.begin(); heightScore != heightScoreMap.end(); ++heightScore)
		if (heightScore->second > maxHeightScore)
			maxHeightScore = heightScore->second;
	maxHeightScore = log10(maxHeightScore * maxHeightScore);

	// iterate through all strands, contigs and positions to find those positions where a stack on the plus strand overlaps a stack on the minus strand by <overlap> nt
	for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
	{
		TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
		if (contigMinusStrand != readCounts[STRAND_MINUS].end())
		{
			for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
			{
				vector< TCountsContig::iterator > stacksOnMinusStrand(MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1, contigMinusStrand->second.end());
				float meanStackHeightInVicinity = 0;
				float maxStackHeightInVicinity = 0;
				for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
				{
					TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
					if (positionMinusStrand != contigMinusStrand->second.end())
					{
						// calculate mean of stack heights in the vicinity
						meanStackHeightInVicinity += positionMinusStrand->second.reads;
						// find highest stacks height in the vicinity
						if (positionMinusStrand->second.reads > maxStackHeightInVicinity)
							maxStackHeightInVicinity = positionMinusStrand->second.reads;
						// remember the stacks that we found, so we do not have to search them again
						stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP] = positionMinusStrand;
					}
				}
				meanStackHeightInVicinity /= stacksOnMinusStrand.size();

				if (maxStackHeightInVicinity > 0) // only continue, if there are any stacks in the vicinity at all
				{
					float heightScorePlus = heightScoreMap[0.5 + positionPlusStrand->second.reads];

					for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
					{
						if (stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP] != contigMinusStrand->second.end())
						{
							// calculate score based on heights of overlapping stacks
							float heightScore = heightScorePlus * heightScoreMap[0.5 + stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP]->second.reads];
							// find the bin for the score
							unsigned int heightScoreBin =
								static_cast<int>(0.5 // add 0.5 for arithmetic rounding when casting float to int
								+ log10(heightScore) // take logarithm of score
								/ maxHeightScore * (groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP].size() - 1)); // assign every score to a bin

							// calculate score based on how much higher the stack is compared to the stacks in the vicinity
							float localHeightScore = (stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP]->second.reads - (meanStackHeightInVicinity - stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP]->second.reads/stacksOnMinusStrand.size())) / maxStackHeightInVicinity;
							// 0.2 seems to be the magical threshold that best segregates ping-pong overlaps from arbitrary overlaps
							unsigned int localHeightScoreBin = (localHeightScore < 0.2) ? IS_BELOW_COVERAGE : IS_ABOVE_COVERAGE;

							// calculate score based on whether the stack on the + strand has Uridine at the 5' end
							unsigned int uridinePlusBin = (positionPlusStrand->second.UAt5PrimeEnd) ? IS_URIDINE : IS_NOT_URIDINE;
							// calculate score based on whether the stack on the - strand has Uridine at the 5' end
							unsigned int uridineMinusBin = (stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP]->second.UAt5PrimeEnd) ? IS_URIDINE : IS_NOT_URIDINE;

							// increase bin counter
							groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][heightScoreBin][uridinePlusBin][uridineMinusBin][localHeightScoreBin]++;

							// keep a list of putative ping-pong signatures, so we can analyze later, which of them are (likely) true
							pingPongSignaturesByOverlap[overlap - MIN_ARBITRARY_OVERLAP][contigPlusStrand->first].push_back(TPingPongSignature(positionPlusStrand->first, heightScoreBin, localHeightScoreBin, uridinePlusBin, uridineMinusBin, positionPlusStrand->second.reads, stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP]->second.reads));
						}
					}
				}
			}

			// free memory of contig on - strand
			contigMinusStrand->second.clear();
		}
		// free memory of contig on + strand
		contigPlusStrand->second.clear();
	}

	// free the rest of memory that might potentially not have been freed yet
	for (TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].begin(); contigMinusStrand != readCounts[STRAND_MINUS].end(); ++contigMinusStrand)
		contigMinusStrand->second.clear();
}

//todo: description
void collapseBins(TGroupedStackCountsByOverlap &groupedStackCountsByOverlap, TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap)
{
	// create a new container to hold the collapsed bin counts
	TGroupedStackCountsByOverlap collapsed = groupedStackCountsByOverlap;

	// keep track of which bins are mapped to which collapsed bins
	vector< unsigned int > oldBinCollapsedBinMap(groupedStackCountsByOverlap.begin()->size());

	unsigned int collapsedBin = 0;
	unsigned int bin = 0;
	while (bin < groupedStackCountsByOverlap.begin()->size()){

		// initialize collapsed bin with 0
		for (TGroupedStackCountsByOverlap::iterator i = collapsed.begin(); i != collapsed.end(); ++i)
			for (vector< vector< vector< float > > >::iterator j = (*i)[collapsedBin].begin(); j != (*i)[collapsedBin].end(); ++j)
				for (vector< vector< float > >::iterator k = j->begin(); k != j->end(); ++k)
					for (vector< float >::iterator l = k->begin(); l != k->end(); ++l)
						*l = 0;

		unsigned int emptyBins;
		do {
			emptyBins = 0;

			// collapse bins
			for (unsigned int overlap = 0; overlap < groupedStackCountsByOverlap.size(); overlap++)
				for (unsigned int i = 0; i < groupedStackCountsByOverlap[overlap][bin].size(); i++)
					for (unsigned int j = 0; j < groupedStackCountsByOverlap[overlap][bin][i].size(); j++)
						for (unsigned int k = 0; k < groupedStackCountsByOverlap[overlap][bin][i][j].size(); k++)
						{
							collapsed[overlap][collapsedBin][i][j][k] += groupedStackCountsByOverlap[overlap][bin][i][j][k];

							// check if collapsedBin is still empty to decide whether to collapse even more
							if (collapsed[overlap][collapsedBin][i][j][k] <= 0)
								emptyBins++;
						}

			// remember with which new bin the current bin was merged
			oldBinCollapsedBinMap[bin] = collapsedBin;

			bin++;
		// stop collapsing bins, when there are no empty bins or when all bins have been merged
		} while ((emptyBins > 0) && (bin < groupedStackCountsByOverlap.begin()->size()));

		collapsedBin++;
	}

	// shrink container to new number of bins
	for (TGroupedStackCountsByOverlap::iterator i = collapsed.begin(); i != collapsed.end(); ++i)
		i->resize(collapsedBin);

	for (TPingPongSignaturesByOverlap::iterator pingPongSignaturesPerGenome = pingPongSignaturesByOverlap.begin(); pingPongSignaturesPerGenome != pingPongSignaturesByOverlap.end(); ++pingPongSignaturesPerGenome)
		for (TPingPongSignaturesPerGenome::iterator contig = pingPongSignaturesPerGenome->begin(); contig != pingPongSignaturesPerGenome->end(); ++contig)
			for (TPingPongSignaturesPerContig::iterator pingPongSignature = contig->second.begin(); pingPongSignature != contig->second.end(); ++pingPongSignature)
				pingPongSignature->heightScoreBin = oldBinCollapsedBinMap[pingPongSignature->heightScoreBin];

	// return collapsed bins as result
	groupedStackCountsByOverlap = collapsed;
}

//todo: description
void calculateFDRs(TGroupedStackCountsByOverlap &groupedStackCountsByOverlap, TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap)
{
	TGroupedStackCountsByOverlap FDRs = groupedStackCountsByOverlap; // the assignment shall only ensure that <FDRs> has the same dimensions as <groupedStackCountsByOverlap>

	// estimate how many of the putative ping-pong signatures are random noise
	for (unsigned int i = 0; i < groupedStackCountsByOverlap.begin()->size(); i++)
		for (unsigned int j = 0; j < (*groupedStackCountsByOverlap.begin())[i].size(); j++)
			for (unsigned int k = 0; k < (*groupedStackCountsByOverlap.begin())[i][j].size(); k++)
				for (unsigned int l = 0; l < (*groupedStackCountsByOverlap.begin())[i][j][k].size(); l++)
				{
					// calculate mean for all the bins of arbitrary overlaps
					float meanOfArbitraryOverlaps = 0;
					for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
						if (overlap != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
							meanOfArbitraryOverlaps += groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k][l];
					meanOfArbitraryOverlaps = meanOfArbitraryOverlaps / (groupedStackCountsByOverlap.size() - 1 /* minus the one bin for ping-pong overlaps */);

					// calculate standard deviation for all the bins of arbitrary overlaps
					float stdDevOfArbitraryOverlaps = 0;
					for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
						if (overlap != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
							stdDevOfArbitraryOverlaps += pow(groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k][l] - meanOfArbitraryOverlaps, 2);
					stdDevOfArbitraryOverlaps = sqrt(1.0 / (groupedStackCountsByOverlap.size() - 1 - 1 /* minus 1 for corrected sample STDDEV */) * stdDevOfArbitraryOverlaps);
					if (stdDevOfArbitraryOverlaps <= 1E-10)
						stdDevOfArbitraryOverlaps = 1E-10; // prevent division by 0, in case the STDDEV is 0

					// calculate expected number of false positives among the putative ping-pong signatures
					for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
					{
						const double approximationAccuracy = 0.01;
						const double approximationRange = 5;
						double fdr = 0;
						for (double x = -approximationRange; x <= +approximationRange; x += approximationAccuracy)
						{
							if (x >= (groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k][l] - meanOfArbitraryOverlaps) / stdDevOfArbitraryOverlaps)
							{
								fdr += approximationAccuracy * 1/sqrt(2*M_PI)*exp(-0.5*x*x) * 1;
							}
							else
							{
								fdr += approximationAccuracy * 1/sqrt(2*M_PI)*exp(-0.5*x*x) * (x * stdDevOfArbitraryOverlaps + meanOfArbitraryOverlaps) / groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k][l];
							}
						}

						if (fdr < 0)
						{
							fdr = 0;
						}
						else if (fdr > 1)
						{
							fdr = 1;
						}

						FDRs[overlap - MIN_ARBITRARY_OVERLAP][i][j][k][l] = fdr;
					}
				}

	// assign a FDR to every putative ping-pong signature
	for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
		for (TPingPongSignaturesPerGenome::iterator contig = pingPongSignaturesByOverlap[overlap - MIN_ARBITRARY_OVERLAP].begin(); contig != pingPongSignaturesByOverlap[overlap - MIN_ARBITRARY_OVERLAP].end(); ++contig)
			for (TPingPongSignaturesPerContig::iterator pingPongSignature = contig->second.begin(); pingPongSignature != contig->second.end(); ++pingPongSignature)
				pingPongSignature->fdr = FDRs[overlap - MIN_ARBITRARY_OVERLAP][pingPongSignature->heightScoreBin][pingPongSignature->UAt5PrimeEndOnPlusStrandBin][pingPongSignature->UAt5PrimeEndOnMinusStrandBin][pingPongSignature->localHeightScoreBin];
}

// todo: description
void plotHistograms(TGroupedStackCountsByOverlap &groupedStackCountsByOverlap, unsigned int dimension, const string &title, vector< string > xAxisLabels, bool logScale = false)
{
	// determine size of given dimension and give histogram plot as many bars
	unsigned int histogramBars = 0;
	switch (dimension)
	{
		case 0:
			histogramBars = groupedStackCountsByOverlap.begin()->size();
			break;
		case 1:
			histogramBars = groupedStackCountsByOverlap.begin()->begin()->size();
			break;
		case 2:
			histogramBars = groupedStackCountsByOverlap.begin()->begin()->begin()->size();
			break;
		case 3:
			histogramBars = groupedStackCountsByOverlap.begin()->begin()->begin()->begin()->size();
			break;
	}
	vector< vector< float > > histograms(groupedStackCountsByOverlap.size(), vector< float >(histogramBars, 0));

	// sum up bins grouped by the given dimension
	for (unsigned int overlap = 0; overlap < groupedStackCountsByOverlap.size(); overlap++)
		for (unsigned int i = 0; i < groupedStackCountsByOverlap[overlap].size(); i++)
			for (unsigned int j = 0; j < groupedStackCountsByOverlap[overlap][i].size(); j++)
				for (unsigned int k = 0; k < groupedStackCountsByOverlap[overlap][i][j].size(); ++k)
					for (unsigned int l = 0; l < groupedStackCountsByOverlap[overlap][i][j][k].size(); ++l)
						switch (dimension)
						{
							case 0:
								histograms[overlap][i] += groupedStackCountsByOverlap[overlap][i][j][k][l];
								break;
							case 1:
								histograms[overlap][j] += groupedStackCountsByOverlap[overlap][i][j][k][l];
								break;
							case 2:
								histograms[overlap][k] += groupedStackCountsByOverlap[overlap][i][j][k][l];
								break;
							case 3:
								histograms[overlap][l] += groupedStackCountsByOverlap[overlap][i][j][k][l];
								break;
						}

	// generate an R script that produces a histogram plot
	string fileName = title;
	replace(fileName.begin(), fileName.end(), ' ', '_'); // replace blanks in file name
	ofstream rScript(toCString(fileName + ".R"));

	rScript << "histograms <- data.frame("; // store histogram counts in a data frame
	for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
	{
			// print height of bars
			rScript << endl << "overlap_";
			if (overlap < 0)
				rScript << "minus_";
			rScript << abs(overlap) << "=c(";
			for (unsigned int bar = 0; bar < histograms[overlap-MIN_ARBITRARY_OVERLAP].size(); bar++)
			{
				if (bar % 10 == 0)
					rScript << endl; // insert a line-break every once in a while, because R cannot parse very long lines
				rScript << histograms[overlap-MIN_ARBITRARY_OVERLAP][bar];
				if (bar < histograms[overlap-MIN_ARBITRARY_OVERLAP].size() - 1)
					rScript << ", "; // separate values by comma, unless it is the last one
			}
			rScript << ")" << endl; // close column of data frame

			if (overlap < MAX_ARBITRARY_OVERLAP)
				rScript << ", "; // separate columns of data frame by comma, unless it is the last column
	}
	rScript << ")" << endl; // close data frame


	// save plot as PNG
	rScript
		<< "options(bitmapType='cairo')" << endl
		<< "png('" << fileName << ".png')" << endl;

	// wrap histogram counts in "log10()", if y-axis should be log-scaled
	if (logScale)
		rScript << "histograms <- log10(histograms)" << endl;

	rScript	<< "plot(0, 0, xlim=c(0," << histograms[0].size() << "), type='n', xlab='" << title << "'";
	
	if (logScale)
	{
		rScript << ", ylim=c(0,max(histograms,0)), ylab='log10(frequency)'";
	}
	else
	{
		rScript << ", ylim=c(0,max(histograms)), ylab='frequency'";
	}

	rScript << ", xaxt='n')" << endl;

	// draw x-axis
	if (xAxisLabels.size() == 0)
	{
		// auto-generate x-axis based on quantiles, if no custom x-axis labels are given
		rScript << "axis(1, at=quantile(c(0," << histograms[0].size() << "), probs = seq(0, 1, 0.2))+0.5, labels=quantile(c(0," << histograms[0].size() << "), probs = seq(0, 1, 0.2)))" << endl;
	}
	else
	{
		// use custom x-axis labels to generate x-axis
		rScript << "axis(1, at=0:" << (xAxisLabels.size() - 1) << "+0.5, labels=c(";
		for (unsigned int i = 0; i < xAxisLabels.size() - 1; i++)
			rScript << "'" << xAxisLabels[i] << "', ";
		rScript << "'" << xAxisLabels[xAxisLabels.size() - 1] << "'))" << endl;
	}

	rScript
		// draw bars for arbitrary overlaps
		<< "for (overlap in " << MIN_ARBITRARY_OVERLAP << ":" << MAX_ARBITRARY_OVERLAP << ")" << endl
		<< "	if (overlap != " << PING_PONG_OVERLAP << ")" << endl
		<< "		barplot(histograms[,gsub('-', 'minus_', paste('overlap_', overlap, sep=''))], col=rgb(0,0,0,alpha=0.1), border=NA, axes=FALSE, add=TRUE, width=1, space=0)" << endl
		// draw a red line for ping-pong overlaps
		<< "for (bin in 1:" << histograms[0].size() << ")" << endl
		<< "	lines(c(bin-1, bin), c(histograms[bin, 'overlap_10'], histograms[bin, 'overlap_10']), type='l', col='red', lwd=2)" << endl
		// draw legend
		<< "legend(x='top', c('" << PING_PONG_OVERLAP << " nt overlap', 'arbitrary overlaps'), col=c('red', 'black'), ncol=2, lwd=c(3,3), xpd=TRUE, inset=-0.1)" << endl
		<< "garbage <- dev.off()" << endl;

	// close R script
	rScript.close();

	// execute R script with "Rscript"
	string RCommand = "Rscript '" + fileName + ".R'";
	system(toCString(RCommand));
}
void plotHistograms(TGroupedStackCountsByOverlap &groupedStackCountsByOverlap, unsigned int dimension, const string &title, bool logScale = false)
{
	vector< string > xAxisLabels;
	plotHistograms(groupedStackCountsByOverlap, dimension, title, xAxisLabels, logScale);
}

// todo: document parameters
//void generateBedGraph(TPingPongSignaturesPerGenome &pingPongSignaturesPerGenome, const TNameStore &bamNameStore, unsigned int minStackHeight)
void generateBedGraph(TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap, const TNameStore &bamNameStore, unsigned int minStackHeight)
{
for (unsigned int i = 0; i < pingPongSignaturesByOverlap.size(); i++)
{
TPingPongSignaturesPerGenome &pingPongSignaturesPerGenome = pingPongSignaturesByOverlap[i];
	// open bedGraph files
	stringstream ss;
	ss << "reads_on_plus_strand_" << i << ".bedGraph";
	ofstream readsOnPlusStrandBedGraph(ss.str().c_str());
	ss.str(""); ss << "reads_on_minus_strand_" << i << ".bedGraph";
	ofstream readsOnMinusStrandBedGraph(ss.str().c_str());
	ss.str(""); ss << "score_" << i << ".bedGraph";
	ofstream scoreBedGraph(ss.str().c_str());
	if (readsOnPlusStrandBedGraph.fail() || readsOnMinusStrandBedGraph.fail() || scoreBedGraph.fail())
	{
		cerr << "Failed to create bedGraph files" << endl;
		return;
	}

	// use scientific formatting for floating point numbers
	readsOnPlusStrandBedGraph.setf(ios::scientific, ios::floatfield);
	readsOnMinusStrandBedGraph.setf(ios::scientific, ios::floatfield);
	scoreBedGraph.setf(ios::scientific, ios::floatfield);

	// write track headers
	readsOnPlusStrandBedGraph << "track type=\"bedGraph\" name=\"_" << i << "read stacks on + strand\" description=\"height of ping-pong stacks on the + strand\" visibility=full color=0,0,0 altColor=0,0,0 priority=20" << endl;
	readsOnMinusStrandBedGraph << "track type=\"bedGraph\" name=\"_" << i << "read stacks on - strand\" description=\"height of ping-pong stacks on the - strand\" visibility=full color=0,0,0 altColor=0,0,0 priority=20" << endl;
	scoreBedGraph << "track type=\"bedGraph\" name=\"_" << i << "scores\" description=\"scores of ping-pong stacks (1 - FDR)\" visibility=full color=0,0,0 altColor=0,0,0 priority=20 viewLimits=0.0:1.0 autoScale=off" << endl;

	// write a line for each ping-pong signature
	for (TPingPongSignaturesPerGenome::iterator contig = pingPongSignaturesPerGenome.begin(); contig != pingPongSignaturesPerGenome.end(); ++contig)
		for (TPingPongSignaturesPerContig::iterator pingPongSignature = contig->second.begin(); pingPongSignature != contig->second.end(); ++pingPongSignature)
//			if ((pingPongSignature->readsOnPlusStrand >= minStackHeight) && (pingPongSignature->readsOnMinusStrand >= minStackHeight))
			{
				readsOnPlusStrandBedGraph << bamNameStore[contig->first] << '\t' << pingPongSignature->position << '\t' << (pingPongSignature->position+1) << '\t' << pingPongSignature->readsOnPlusStrand << endl;
				readsOnMinusStrandBedGraph << bamNameStore[contig->first] << '\t' << pingPongSignature->position << '\t' << (pingPongSignature->position+1) << '\t' << pingPongSignature->readsOnMinusStrand << endl;
				scoreBedGraph << bamNameStore[contig->first] << '\t' << pingPongSignature->position << '\t' << (pingPongSignature->position+1) << '\t' << (1-pingPongSignature->fdr) << endl;
			}

	// close bedGraph files
	readsOnPlusStrandBedGraph.close();
	readsOnMinusStrandBedGraph.close();
	scoreBedGraph.close();
}
}

void readTransposonsFromFile(ifstream &transposonFile, TFileFormat fileFormat, TTransposonsPerGenome &transposons, TNameStore &nameStore)
{
	// the following variables store the numbers of the columns of the respective fields
	unsigned int identifierField, strandField, contigField, startField, endField;
	int offset;
	// character that separates columns
	char delimiter;

	switch (fileFormat)
	{
		case fileFormatBED:
			identifierField = 4;
			strandField = 6;
			contigField = 1;
			startField = 2;
			endField = 3;
			offset = 0;
			delimiter = ' ';
			break;
		case fileFormatCSV:
			identifierField = 1;
			strandField = 2;
			contigField = 3;
			startField = 4;
			endField = 5;
			offset = 1;
			delimiter = ',';
			break;
		case fileFormatGFF:
		case fileFormatGTF:
			identifierField = 9;
			strandField = 7;
			contigField = 2;
			startField = 4;
			endField = 5;
			offset = 1;
			delimiter = '\t';
			break;
		case fileFormatTSV:
		default:
			identifierField = 1;
			strandField = 2;
			contigField = 3;
			startField = 4;
			endField = 5;
			offset = 1;
			delimiter = '\t';
			break;
	}

	unsigned int fieldNumber = 1;
	string fieldValue = "";
	string transposonIdentifier = "";
	int transposonStrand = -1;
	int transposonContig = -1;
	int transposonStart = -1;
	int transposonEnd = -1;
	bool newLine = false; // set to true, when a line-feed is read
	bool quotesOpen = false; // in CSV files, keeps track of opening and closing double-quotes

	// read the transposon file character by character
	while (transposonFile.good())
	{
		char character = transposonFile.get();
		
		if (!transposonFile.good())
			character = '\n'; // make sure the last line is processed, even if it does not end on a line feed

		if (character == '\n')
		{
			newLine = true;
			character = delimiter; // when we encounter a new line, we basically do the same as when we encounter a delimiter
		}
		else if ((delimiter == ' ') && (character == '\t')) // when delimiter is set to a blank, then treat any white-space (i.e., ' ' and '\t') as a delimiter
		{
			character = ' ';
		}

		if ((character == delimiter) && !((fileFormat = fileFormatCSV) && quotesOpen))
		{
			if (!fieldValue.empty() || newLine)
			{
				if (fieldNumber == identifierField)
				{
					transposonIdentifier = fieldValue;
				}
				else if (fieldNumber == strandField)
				{
					if (fieldValue == "+")
						transposonStrand = STRAND_PLUS;
					else if (fieldValue == "-")
						transposonStrand = STRAND_MINUS;
				}
				else if (fieldNumber == contigField)
				{
					for (unsigned int i = 0; (transposonContig == -1) && (i < length(nameStore)); i++)
						if (nameStore[i] == fieldValue)
							transposonContig = i;

					if (transposonContig == -1) // the contig was not found in the name store
					{
						// add a new element to the name store
						appendValue(nameStore, fieldValue);
						transposonContig = length(nameStore) - 1;
					}
				}
				else if (fieldNumber == startField)
				{
					transposonStart = atoi(fieldValue.c_str());
					if ((transposonStart == 0) && (fieldValue != "0")) // a conversion error occurred
					{
						transposonStart = -1;
					}
					else
					{
						transposonStart -= offset; // some file formats have 0-based coordinates (offset=0), others 1-based (offset=1)
					}
				}
				else if (fieldNumber == endField)
				{
					transposonEnd = atoi(fieldValue.c_str());
					if ((transposonEnd == 0) && (fieldValue != "0")) // a conversion error occurred
					{
						transposonEnd = -1;
					}
					else
					{
						transposonEnd -= offset; // some file formats have 0-based coordinates (offset=0), others 1-based (offset=1)
					}
				}

				if (newLine)
				{
					// skip lines that could not be parsed
					if ((transposonStrand >= 0) && (transposonContig >= 0) && (transposonStart >= 0) && (transposonEnd >= 0) && !transposonIdentifier.empty())
						transposons[transposonContig].push_back(TTransposon(transposonIdentifier, transposonStrand, transposonStart, transposonEnd));

					// reset fields
					newLine = false;
					transposonStrand = transposonContig = transposonStart = transposonEnd = -1;
					transposonIdentifier.clear();
					fieldNumber = 1;
				}
				else
				{
					fieldNumber++;
				}

				fieldValue = "";
			}
		}
		else 
		{
			if (!((character == '\r') && (transposonFile.peek() == '\n'))) // skip carriage-returns, if the next character is a line feed
			{
				if ((fileFormat == fileFormatCSV) && (character == '"')) // in the CSV format, fields may be enclosed by double-quotes
				                                                         // we need to keep track of whether all opening quotes have a closing counterpart
				{
					if (!fieldValue.empty() && !quotesOpen)
						fieldValue += character; // two consecutive double-quotes are collapsed to one
					quotesOpen = !quotesOpen;
				}
				else
				{
					fieldValue += character;
				}
			}
		}
	}

	// sort transposons by genomic position
	for (TTransposonsPerGenome::iterator contig = transposons.begin(); contig != transposons.end(); ++contig)
		contig->second.sort();
}

// Function to sort list of Transposons by p-value for multiple testing correction using Benjamini-Hochberg procedure.
// The function is used by list::sort to compare which of two TTransposon objects is lower, based on their pValue attribute.
// Parameters:
// 	transposon1, transposon2: the transposons to compare
inline bool compareTransposonsByPValue(const TTransposonsPerContig::iterator &transposon1, const TTransposonsPerContig::iterator &transposon2)
{
	return transposon1->pValue < transposon2->pValue;
}

// todo: description
void findSuppressedTransposons(TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap, TTransposonsPerGenome &transposons, const float totalReadCount)
{
	// slide over the genome and calculate a z-score for every transposon overlapping the current position
	for (TTransposonsPerGenome::iterator contig = transposons.begin(); contig != transposons.end(); ++contig)
	{
		// for every overlap we need to keep track of the iterator that points to the ping-pong signature where we are currently at
		vector< TPingPongSignaturesPerContig::iterator > positionByOverlap(pingPongSignaturesByOverlap.size());
		for (unsigned int overlap = 0; overlap < pingPongSignaturesByOverlap.size(); overlap++)
			if (pingPongSignaturesByOverlap[overlap].find(contig->first) != pingPongSignaturesByOverlap[overlap].end())
				positionByOverlap[overlap] = pingPongSignaturesByOverlap[overlap][contig->first].begin();
			else
				positionByOverlap[overlap] = pingPongSignaturesByOverlap[overlap][contig->first].end(); // there are no ping-pong stacks for the given contig and overlap

		for (TTransposonsPerContig::iterator transposon = contig->second.begin(); transposon != contig->second.end(); ++transposon)
		{
			float scoreOfPingPongOverlap = 0;
			// calculate mean transposon score of all arbitrary overlaps
			float meanOfArbitraryOverlaps = 0;
			for (unsigned int overlap = 0; overlap < positionByOverlap.size(); overlap++)
			{
				// move iterator of ping-pong signature to start of current transposon
				while ((positionByOverlap[overlap] != pingPongSignaturesByOverlap[overlap][contig->first].begin()) && (positionByOverlap[overlap]->position > transposon->start))
					--(positionByOverlap[overlap]);
				while ((positionByOverlap[overlap] != pingPongSignaturesByOverlap[overlap][contig->first].end()) && (positionByOverlap[overlap]->position < transposon->start))
					++(positionByOverlap[overlap]);

				// sum up the scores of all ping-pong signatures within the transposon region
				float sumOfScores = 0;
				if (positionByOverlap[overlap]->position >= transposon->start)
				{
					while ((positionByOverlap[overlap]->position <= transposon->end) && (positionByOverlap[overlap] != pingPongSignaturesByOverlap[overlap][contig->first].end()))
					{
						sumOfScores += (1 - positionByOverlap[overlap]->fdr);
						++(positionByOverlap[overlap]);
					}
				}

				if (overlap + MIN_ARBITRARY_OVERLAP != PING_PONG_OVERLAP) // ignore ping-pong overlaps in the mean calculation, since they would skew the result
					meanOfArbitraryOverlaps += sumOfScores;
				else
					scoreOfPingPongOverlap = sumOfScores; // remember the transposon score for ping-pong overlap
			}
			meanOfArbitraryOverlaps = meanOfArbitraryOverlaps / (positionByOverlap.size() - 1 /* minus the one bin for ping-pong overlaps */);
			if ((meanOfArbitraryOverlaps == 0) && (scoreOfPingPongOverlap == 0)) // there are no ping-pong signatures in the region of the transposon
			{
				transposon->pValue = 1;
				transposon->qValue = 1;
				transposon->normalizedSignatureCount = 0;
			}
			else
			{
				// calculate standard deviation of transposon score of all arbitrary overlaps
				float stdDevOfArbitraryOverlaps = 0;
				for (unsigned int overlap = 0; overlap < positionByOverlap.size(); overlap++)
					if (overlap + MIN_ARBITRARY_OVERLAP != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
					{
						// move iterator of ping-pong signature to start of current transposon
						while ((positionByOverlap[overlap] != pingPongSignaturesByOverlap[overlap][contig->first].begin()) && (positionByOverlap[overlap]->position > transposon->start))
							--(positionByOverlap[overlap]);
						while ((positionByOverlap[overlap] != pingPongSignaturesByOverlap[overlap][contig->first].end()) && (positionByOverlap[overlap]->position < transposon->start))
							++(positionByOverlap[overlap]);

						// sum up the scores of all ping-pong signatures within the transposon region
						float sumOfScores = 0;
						if (positionByOverlap[overlap]->position >= transposon->start)
						{
							while ((positionByOverlap[overlap]->position <= transposon->end) && (positionByOverlap[overlap] != pingPongSignaturesByOverlap[overlap][contig->first].end()))
							{
								sumOfScores += (1 - positionByOverlap[overlap]->fdr);
								++(positionByOverlap[overlap]);
							}
						}

						stdDevOfArbitraryOverlaps += pow(sumOfScores - meanOfArbitraryOverlaps, 2);
					}
				stdDevOfArbitraryOverlaps = sqrt(1.0 / (positionByOverlap.size() - 1 - 1 /* minus 1 for corrected sample STDDEV */) * stdDevOfArbitraryOverlaps);
				if (stdDevOfArbitraryOverlaps <= 1E-10)
					stdDevOfArbitraryOverlaps = 1E-10; // prevent division by 0, in case the STDDEV is 0

				// calculate significance of transposon score of ping-pong overlap vs. arbitrary overlaps
				const double approximationAccuracy = 0.01;
				const double approximationRange = 5;
				double zValue = (scoreOfPingPongOverlap - meanOfArbitraryOverlaps) / stdDevOfArbitraryOverlaps;
				double pValue = 0;
				for (double x = zValue; x <= zValue + approximationRange; x += approximationAccuracy)
					pValue += approximationAccuracy * 1/sqrt(2*M_PI)*exp(-0.5*x*x);

				transposon->pValue = pValue;
				// the normalized signature count is the number of ping-pong signatures per kilobase per million mapped reads
				transposon->normalizedSignatureCount = scoreOfPingPongOverlap / ((static_cast<float>(transposon->end) - transposon->start)/1000) / (totalReadCount/1000000);
			}
		}
	}

	// sort transposons by p-value for multiple testing-correction with Benjamini-Hochberg procedure (FDR)
	list< TTransposonsPerContig::iterator > transposonsSortedByPValue;
	for (TTransposonsPerGenome::iterator contig = transposons.begin(); contig != transposons.end(); ++contig)
		for (TTransposonsPerContig::iterator transposon = contig->second.begin(); transposon != contig->second.end(); ++transposon)
			transposonsSortedByPValue.push_back(transposon);
	transposonsSortedByPValue.sort(compareTransposonsByPValue);
	
	unsigned i = 1;
	float previousPValue = 1;
	float previousQValue = 1;
	list< TTransposonsPerContig::iterator >::reverse_iterator transposon = transposonsSortedByPValue.rbegin();
	while (transposon != transposonsSortedByPValue.rend())
	{
		if ((*transposon)->pValue == previousPValue)
		{
			// If two transposons have the same p-value, re-use the previously calculated q-value.
			// This is pretty unlikely, but ensures that two transposons with the same p-value also
			// get the same q-value.
			(*transposon)->qValue = previousQValue;
		}
		else
		{
			(*transposon)->qValue = (*transposon)->pValue * transposonsSortedByPValue.size() / i;
			if ((*transposon)->qValue > 1)
				(*transposon)->qValue = 1;

			previousPValue = (*transposon)->pValue;
			previousQValue = (*transposon)->qValue;
		}

		++transposon;
		i++;
	}
}

void writeTransposonsToTSV(TTransposonsPerGenome &transposons, TNameStore &bamNameStore)
{
	ofstream transposonsTSV("transposons.tsv");
	if (transposonsTSV.fail())
	{
		cerr << "Failed to create transposon file" << endl;
		return;
	}

	// use scientific formatting for floating point numbers in the output file
	transposonsTSV.setf(ios::scientific, ios::floatfield);

	transposonsTSV << "identifier\tstrand\tcontig\tstart\tend\tpValue\tqValue\tnormalizedSignatureCount" << endl;
	for (TTransposonsPerGenome::iterator contig = transposons.begin(); contig != transposons.end(); ++contig)
		for (TTransposonsPerContig::iterator transposon = contig->second.begin(); transposon != contig->second.end(); ++transposon)
			transposonsTSV << transposon->identifier << '\t' << ((transposon->strand == STRAND_PLUS) ? '+' : '-') << '\t' << bamNameStore[contig->first] << '\t' << transposon->start << '\t' << transposon->end << '\t' << transposon->pValue << '\t' << transposon->qValue << '\t' << transposon->normalizedSignatureCount << endl;

	transposonsTSV.close();
}

// program entry point
int main(int argc, char const ** argv)
{
	// parse the command line options
	AppOptions options;
	if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
		return 1;

	TCountsGenome readCounts; // stats about positions where reads on the minus strand overlap with the 5' ends of reads on the plus strand

	float totalReadCount = 0;

	TNameStore bamNameStore; // structure to store contig names

	// read all BAM/SAM files
	if (options.verbosity >= 3)
		cerr << "Counting reads in SAM/BAM files" << endl;
	for(TInputFiles::iterator inputFile = options.inputFiles.begin(); inputFile != options.inputFiles.end(); ++inputFile)
	{
		stopwatch(toCString(*inputFile), options.verbosity);

		// open SAM/BAM file
		BamStream bamFile(toCString(*inputFile));
		if (!isGood(bamFile))
		{
			cerr << "Failed to open input file: " << *inputFile << endl;
			return 1;
		}

		// for every position in the genome, count the number of reads that start at a given position
		if (countReadsInBamFile(bamFile, readCounts, options.minAlignmentLength, options.maxAlignmentLength, options.countMultiHits, totalReadCount) != 0)
			return 1;

		// remember @SQ header lines from BAM file for mapping of contig IDs to human-readable names
		if (length(bamNameStore) == 0)
		{
			bamNameStore = nameStore(bamFile.bamIOContext);
		}
		else
		{
			// if multiple BAM files are given, check if headers are identical
			TNameStore tempBamNameStore = nameStore(bamFile.bamIOContext);
			bool nameStoresDiffer = false;
			TNameStoreIterator bamNameStoreIterator = begin(bamNameStore);
			TNameStoreIterator tempBamNameStoreIterator = begin(tempBamNameStore);
			while ((bamNameStoreIterator != end(bamNameStore)) && (tempBamNameStoreIterator != end(tempBamNameStore)) && !nameStoresDiffer)
			{
				if (value(bamNameStoreIterator) != value(tempBamNameStoreIterator))
					nameStoresDiffer = true;
				++bamNameStoreIterator;
				++tempBamNameStoreIterator;
			}
			if (nameStoresDiffer || (length(bamNameStore) != length(tempBamNameStore)))
			{
				cerr << "@SQ header lines of '" << *inputFile << "' differ from those of previous input files" << endl;
				return 1;
			}
		}

		// close SAM/BAM file
		close(bamFile);

		stopwatch(options.verbosity);
	}

	// read transposons, if files are given
	TTransposonsPerGenome transposons;
	if (options.transposonFiles.size() > 0)
	{
		if (options.verbosity >= 3)
			cerr << "Loading transposon coordinates" << endl;
		for (TInputFiles::iterator transposonFile = options.transposonFiles.begin(); transposonFile != options.transposonFiles.end(); ++transposonFile)
		{
			stopwatch(toCString(*transposonFile), options.verbosity);

			// try to open file
			ifstream fileStream(toCString(*transposonFile));
			if (fileStream.fail())
			{
				cerr << "Failed to open transposon file \"" << (*transposonFile) << "\"." << endl;
				return 1;
			}

			// determine type of input file
			TFileFormat fileFormat;
			if (_compareExtension(toCString(*transposonFile), ".bed"))
				fileFormat = fileFormatBED;
			else if (_compareExtension(toCString(*transposonFile), ".csv"))
				fileFormat = fileFormatCSV;
			else if (_compareExtension(toCString(*transposonFile), ".gff"))
				fileFormat = fileFormatGFF;
			else if (_compareExtension(toCString(*transposonFile), ".gtf"))
				fileFormat = fileFormatGTF;
			else
				fileFormat = fileFormatTSV;

			readTransposonsFromFile(fileStream, fileFormat, transposons, bamNameStore);
			fileStream.close();
			stopwatch(options.verbosity);
		}
	}

	// go to output directory
	if (length(options.output) > 0)
	{
		mkdir(toCString(options.output), 0777);
		if (chdir(toCString(options.output)) != 0)
		{
			cerr << "Failed to open output directory: " << options.output;
			return 1;
		}
	}

	stopwatch("Binning stacks", options.verbosity);
	THeightScoreMap heightScoreMap;
	mapHeightsToScores(readCounts, heightScoreMap);
	TGroupedStackCountsByOverlap groupedStackCountsByOverlap;
	TPingPongSignaturesByOverlap pingPongSignaturesByOverlap;
	countStacksByGroup(readCounts, heightScoreMap, groupedStackCountsByOverlap, pingPongSignaturesByOverlap);
	stopwatch(options.verbosity);

	stopwatch("Collapsing bins", options.verbosity);
	collapseBins(groupedStackCountsByOverlap, pingPongSignaturesByOverlap);
	stopwatch(options.verbosity);

	stopwatch("Calculating FDR for putative ping-pong signatures", options.verbosity);
	calculateFDRs(groupedStackCountsByOverlap, pingPongSignaturesByOverlap);
	stopwatch(options.verbosity);

	if (options.plot)
	{
		stopwatch("Generating R plots", options.verbosity);
		plotHistograms(groupedStackCountsByOverlap, 0, "height score", true);
		vector< string > xAxisLabels;
		xAxisLabels.push_back("uridine"); xAxisLabels.push_back("not uridine");
		plotHistograms(groupedStackCountsByOverlap, 1, "base content at 5-prime end on forward strand", xAxisLabels);
		xAxisLabels.resize(0);
		xAxisLabels.push_back("uridine"); xAxisLabels.push_back("not uridine");
		plotHistograms(groupedStackCountsByOverlap, 2, "base content at 5-prime end on reverse strand", xAxisLabels);
		xAxisLabels.resize(0);
		xAxisLabels.push_back("average"); xAxisLabels.push_back("above average");
		plotHistograms(groupedStackCountsByOverlap, 3, "local height score", xAxisLabels);
		stopwatch(options.verbosity);
	}

	if (options.bedGraph)
	{
		stopwatch("Generating bedGraph files", options.verbosity);
//		generateBedGraph(pingPongSignaturesByOverlap[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP], bamNameStore, options.minStackHeight);
		generateBedGraph(pingPongSignaturesByOverlap, bamNameStore, options.minStackHeight);
		stopwatch(options.verbosity);
	}

	if (options.transposonFiles.size() > 0)
	{
		stopwatch("Finding suppressed transposons", options.verbosity);
		findSuppressedTransposons(pingPongSignaturesByOverlap, transposons, totalReadCount);
		writeTransposonsToTSV(transposons, bamNameStore);
		stopwatch(options.verbosity);
	}

	return 0;
}

