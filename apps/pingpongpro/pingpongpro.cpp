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
	bool browserTracks;
	TInputFiles inputFiles;
	unsigned int minAlignmentLength;
	unsigned int maxAlignmentLength;
	unsigned int minStackHeight;
	TCountMultiHits countMultiHits;
	CharString output;
	bool plot;
	TInputFiles transposonFiles;
	bool predictTransposons;
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

// types to draw a histogram
typedef vector< float > THistogram;
typedef vector< THistogram > THistograms;

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
	THistogram histogram;

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
		histogram = transposon.histogram;
	}

	// needed for sorting of the list by genomic position
	inline bool operator<(const TTransposon &transposon)
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

#define PREDICT_TRANSPOSONS_SLIDING_WINDOW 1000
#define PREDICT_TRANSPOSONS_LENGTH 30

#define APPROXIMATION_ACCURACY 0.01
#define APPROXIMATION_RANGE 5

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

	addOption(parser, ArgParseOption("b", "browserTracks", "Generate genome browser tracks for loci with ping-pong signature and (if -t is specified) for transposons with ping-pong activity. Default: off."));

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

	addOption(parser, ArgParseOption("t", "transposons", "Check if the transposons given in the file are suppressed through ping-pong activity.", ArgParseArgument::INPUTFILE, "PATH", true));
	setValidValues(parser, "transposons", ".bed .csv .gff .gtf .tsv");

	addOption(parser, ArgParseOption("T", "predict-transposons", "Predict the location of suppressed transposons based on regions with high ping-pong activity. Default: off."));

	addOption(parser, ArgParseOption("v", "verbose", "Print messages to stderr about the current progress. Default: off."));

	// parse command line
	ArgumentParser::ParseResult parserResult = parse(parser, argc, argv);
	if (parserResult != ArgumentParser::PARSE_OK)
		return parserResult;

	// extract options, if parsing was successful
	options.browserTracks = isSet(parser, "browserTracks");
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
	options.predictTransposons = isSet(parser, "predict-transposons");
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
						double fdr = 0;
						double degreesOfFreedom = MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP;
						for (double x = -APPROXIMATION_RANGE; x <= +APPROXIMATION_RANGE; x += APPROXIMATION_ACCURACY)
						{
							if (x >= (groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k][l] - meanOfArbitraryOverlaps) / stdDevOfArbitraryOverlaps)
							{
								fdr += APPROXIMATION_ACCURACY * tgamma((degreesOfFreedom+1)/2) / sqrt(degreesOfFreedom * M_PI) / tgamma(degreesOfFreedom/2) * pow(1 + x * x / degreesOfFreedom, -(degreesOfFreedom + 1)/2) * 1;
							}
							else
							{
								fdr += APPROXIMATION_ACCURACY * tgamma((degreesOfFreedom+1)/2) / sqrt(degreesOfFreedom * M_PI) / tgamma(degreesOfFreedom/2) * pow(1 + x * x / degreesOfFreedom, -(degreesOfFreedom + 1)/2) * (x * stdDevOfArbitraryOverlaps + meanOfArbitraryOverlaps) / groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k][l];
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

void stringReplace(string &subjectString, const string &searchString, const string &replaceString) {
	size_t i = 0;
	while((i = subjectString.find(searchString, i)) != string::npos)
	{
	        subjectString.replace(i, searchString.length(), replaceString);
        	i += replaceString.length(); // skip the string that we just inserted
	}
}

// todo: description
void plotHistogram(const string &fileName, const vector< string > &titles, const THistograms &histograms)
{
	// generate an R script that produces a histogram plot
	ofstream rScript(toCString(fileName + ".R"), ios_base::out);
	if (rScript.fail())
	{
		cerr << "Failed to create R script file" << endl;
		return;
	}

	rScript << "histograms = data.frame(" << endl; // store histograms in a data frame
	rScript << "plotTitle = c("; // store plot titles in vector
	for (unsigned int i = 0; i < titles.size(); i++)
	{
		// escape all single-quotes (') and escape slashes (\) in the title
		string title = titles[i];
		stringReplace(title, "\\", "\\\\");
		stringReplace(title, "'", "\\'");
		rScript << "'" << title << "'";

		if (i < titles.size() - 1)
			rScript << "," << endl; // separate titles by a comma, unless it is the last one
	}
	rScript << ")," << endl; // close vector of plot titles
	for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
	{
		// print one column in the data frame for every overlap
		rScript << endl << "overlap_";
		if (overlap < 0)
			rScript << "minus_";
		rScript << abs(overlap) << "=c(";

		for (unsigned int i = 0; i < histograms.size(); i++)
		{
			if (i % 100 == 0)
				rScript << endl; // insert a line-break every once in a while, because R cannot parse long lines
			rScript << histograms[i][overlap-MIN_ARBITRARY_OVERLAP];
			if (i < histograms.size() - 1)
				rScript << ","; // separate histogram values by comma, unless it is the last one
		}

		rScript << ")" << endl; // close column of data frame
		if (overlap < MAX_ARBITRARY_OVERLAP)
			rScript << ","; // separate columns of data frame by comma, unless it is the last column
	}
	rScript << ")" << endl; // close data frame of histograms

	rScript
	 	<< "# convert absolute values to z-scores" << endl
		<< "means <- apply(histograms[,!colnames(histograms) %in% c('plotTitle', 'overlap_" << PING_PONG_OVERLAP << "')], 1, mean)" << endl
		<< "sds <- apply(histograms[,!colnames(histograms) %in% c('plotTitle', 'overlap_" << PING_PONG_OVERLAP << "')], 1, sd)" << endl
		<< "sds <- ifelse(sds < 1e-10, 1e-10, sds)" << endl
		<< "for (column in colnames(histograms[,colnames(histograms) != 'plotTitle'])) {" << endl
		<< "	histograms[,column] = (histograms[,column] - means) / sds" << endl
		<< "}" << endl
		<< "# save plots to a single PDF" << endl
		<< "pdf('" << fileName << ".pdf', onefile=TRUE)" << endl
		<< "par(font.lab=2, mar=c(5.1, 5.1, 5.1, 2.1))" << endl
		<< "# draw a red bar for ping-pong signatures and a grey bar for arbitrary overlaps" << endl
		<< "barColors <- ifelse(colnames(histograms[,colnames(histograms) != 'plotTitle']) != 'overlap_" << PING_PONG_OVERLAP << "', rgb(0.7,0.7,0.7), rgb(0.8,0.4,0.4))" << endl
		<< "# every row of the data frame <histograms> is rendered as a barplot" << endl
		<< "for (i in 1:nrow(histograms)) {" << endl
		<< "	histogram <- histograms[i,]" << endl
		<< "	plotTitle <- histogram$plotTitle" << endl
		<< "	histogram <- c(t(histogram[, colnames(histogram) != 'plotTitle']))" << endl
		<< "	plot(0, 0, xlim=c(0,length(histogram)), type='n', xlab='overlap', ylim=c(min(c(0, histogram)), max(histogram)), ylab='z-score', xaxt='n', main=plotTitle, cex.axis=0.75)" << endl
		<< "	axis(1, at=1:length(histogram)-0.25, labels=" << MIN_ARBITRARY_OVERLAP << ":" << MAX_ARBITRARY_OVERLAP << ", cex.axis=0.75)" << endl
		<< "	barplot(histogram, col=barColors, border=NA, axes=FALSE, add=TRUE, width=0.5, space=1)" << endl
		<< "}" << endl
		<< "garbage <- dev.off()" << endl;

	// close R script
	rScript.close();

	// execute R script with "Rscript"
	// running the script via the source command is faster that running it directly
	string RCommand = "Rscript -e 'source(\"" + fileName + ".R\")'";
	system(toCString(RCommand));
}

// todo: document parameters
void writePingPongSignaturesToFile(TPingPongSignaturesPerGenome &pingPongSignaturesPerGenome, const TNameStore &bamNameStore, unsigned int minStackHeight, bool browserTracks)
{
	// open files to write ping-pong signatures to
	ofstream signaturesTSV("ping-pong_signatures.tsv", ios_base::out);
	ofstream readsOnPlusStrandBedGraph;
	ofstream readsOnMinusStrandBedGraph;
	ofstream scoresBedGraph;
	if (browserTracks)
	{
		readsOnPlusStrandBedGraph.open("ping-pong_signatures_read_stacks_on_plus_strand.bedGraph", ios_base::out);
		readsOnMinusStrandBedGraph.open("ping-pong_signatures_read_stacks_on_minus_strand.bedGraph", ios_base::out);
		scoresBedGraph.open("ping-pong_signatures_scores.bedGraph", ios_base::out);
		if (readsOnPlusStrandBedGraph.fail() || readsOnMinusStrandBedGraph.fail() || scoresBedGraph.fail())
		{
			cerr << "Failed to create browser track files for ping-pong signatures" << endl;
			return;
		}
	}

	// use scientific formatting for floating point numbers
	signaturesTSV.setf(ios::scientific, ios::floatfield);
	if (browserTracks)
	{
		readsOnPlusStrandBedGraph.setf(ios::scientific, ios::floatfield);
		readsOnMinusStrandBedGraph.setf(ios::scientific, ios::floatfield);
		scoresBedGraph.setf(ios::scientific, ios::floatfield);
	}

	// write track headers
	signaturesTSV << "contig\tposition\tFDR\tstackHeightOnPlusStrand\tstackHeightOnMinusStrand" << endl;
	if (browserTracks)
	{
		readsOnPlusStrandBedGraph << "track type=bedGraph name=\"read stacks on + strand\" description=\"height of read stacks on the + strand\" visibility=full" << endl;
		readsOnMinusStrandBedGraph << "track type=bedGraph name=\"read stacks on - strand\" description=\"height of read stacks on the - strand\" visibility=full" << endl;
		scoresBedGraph << "track type=bedGraph name=\"scores\" description=\"scores of ping-pong stacks (1 - FDR)\" visibility=full viewLimits=0.0:1.0 autoScale=off" << endl;
	}

	// write a line for each ping-pong signature
	for (TPingPongSignaturesPerGenome::iterator contig = pingPongSignaturesPerGenome.begin(); contig != pingPongSignaturesPerGenome.end(); ++contig)
		for (TPingPongSignaturesPerContig::iterator pingPongSignature = contig->second.begin(); pingPongSignature != contig->second.end(); ++pingPongSignature)
			if ((pingPongSignature->readsOnPlusStrand >= minStackHeight) && (pingPongSignature->readsOnMinusStrand >= minStackHeight))
			{
				signaturesTSV << bamNameStore[contig->first] << '\t' << pingPongSignature->position << '\t' << pingPongSignature->fdr << '\t' << pingPongSignature->readsOnPlusStrand << '\t' << pingPongSignature->readsOnMinusStrand << endl;
				if (browserTracks)
				{
					readsOnPlusStrandBedGraph << bamNameStore[contig->first] << '\t' << pingPongSignature->position << '\t' << (pingPongSignature->position+1) << '\t' << pingPongSignature->readsOnPlusStrand << endl;
					readsOnMinusStrandBedGraph << bamNameStore[contig->first] << '\t' << pingPongSignature->position << '\t' << (pingPongSignature->position+1) << '\t' << pingPongSignature->readsOnMinusStrand << endl;
					scoresBedGraph << bamNameStore[contig->first] << '\t' << pingPongSignature->position << '\t' << (pingPongSignature->position+1) << '\t' << (1-pingPongSignature->fdr) << endl;
				}
			}

	// close files
	signaturesTSV.close();
	if (browserTracks)
	{
		readsOnPlusStrandBedGraph.close();
		readsOnMinusStrandBedGraph.close();
		scoresBedGraph.close();
	}
}

void generateGroupedStackCountsPlot(TGroupedStackCountsByOverlap &groupedStackCountsByOverlap)
{
	int binCount =
		groupedStackCountsByOverlap.begin()->size() *
		groupedStackCountsByOverlap.begin()->begin()->size() *
		groupedStackCountsByOverlap.begin()->begin()->begin()->size() *
		groupedStackCountsByOverlap.begin()->begin()->begin()->begin()->size();
	THistograms histograms(binCount);
	vector< string > plotTitles(binCount);

	unsigned int x = 0;
	stringstream ss;
	for (int i = groupedStackCountsByOverlap.begin()->size() - 1; i >= 0; i--)
		for (unsigned int j = 0; j < (*groupedStackCountsByOverlap.begin())[i].size(); j++)
			for (unsigned int k = 0; k < (*groupedStackCountsByOverlap.begin())[i][j].size(); k++)
				for (unsigned int l = 0; l < (*groupedStackCountsByOverlap.begin())[i][j][k].size(); l++)
				{
					histograms[x].resize(groupedStackCountsByOverlap.size());
					for (unsigned int overlap = 0; overlap < groupedStackCountsByOverlap.size(); overlap++)
					{
						histograms[x][overlap] = groupedStackCountsByOverlap[overlap][i][j][k][l];

						ss << "z-scores of signatures with the following properties:" << endl;
						ss << "stack height score of " << (groupedStackCountsByOverlap.begin()->size() - 1 - i) << endl;
						if ((j == IS_URIDINE) && (k == IS_URIDINE))
							ss << "uridine at the 5'-end on both strands" << endl;
						else if ((j == IS_URIDINE) && (k != IS_URIDINE))
							ss << "uridine at the 5'-end on the + strand" << endl;
						else if ((j != IS_URIDINE) && (k == IS_URIDINE))
							ss << "uridine at the 5'-end on the - strand" << endl;
						else if ((j != IS_URIDINE) && (k != IS_URIDINE))
							ss << "no uridine at the 5'-end on either strand" << endl;
						ss << "stack height " << ((l == IS_ABOVE_COVERAGE) ? "above" : "below") << " the local coverage";
						plotTitles[x] = ss.str();
						ss.str("");
					}
					x++;
				}
	plotHistogram("ping-pong_signature_z-scores", plotTitles, histograms);
}

void readTransposonsFromFile(ifstream &transposonFile, TFileFormat fileFormat, TTransposonsPerGenome &transposons, TNameStore &nameStore)
{
	// the following variables store the numbers of the columns of the respective fields
	unsigned int identifierField, strandField, contigField, startField, endField;
	// some file formats address the first position on the contig with 0, others with 1
	// some file formats include the last position in regions (half-open), others not (closed)
	unsigned int offset, halfOpenOrClosed;
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
			halfOpenOrClosed = 0;
			delimiter = ' ';
			break;
		case fileFormatCSV:
			identifierField = 1;
			strandField = 2;
			contigField = 3;
			startField = 4;
			endField = 5;
			offset = 1;
			halfOpenOrClosed = 0;
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
			halfOpenOrClosed = 1;
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
			halfOpenOrClosed = 0;
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

		if ((character == delimiter) && !((fileFormat == fileFormatCSV) && quotesOpen))
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
						transposonStart -= offset; // shift coordinate, in case the file format starts counting at 1 instead of 0
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
						transposonEnd = transposonEnd - offset + halfOpenOrClosed; // shift coordinate, in case the file format starts counting at 1 instead of 0 or if the format is closed instead of half open
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
			// calculate score for transposon for each overlap
			transposon->histogram.resize(positionByOverlap.size());
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

				transposon->histogram[overlap] = sumOfScores;
			}

			// calculate mean transposon score of all arbitrary overlaps
			float meanOfArbitraryOverlaps = 0;
			for (unsigned int overlap = 0; overlap < transposon->histogram.size(); overlap++)
				if (overlap + MIN_ARBITRARY_OVERLAP != PING_PONG_OVERLAP) // ignore ping-pong overlaps in the mean calculation, since they would skew the result
					meanOfArbitraryOverlaps += transposon->histogram[overlap];
			meanOfArbitraryOverlaps = meanOfArbitraryOverlaps / (positionByOverlap.size() - 1 /* minus the one bin for ping-pong overlaps */);

			if ((meanOfArbitraryOverlaps == 0) && (transposon->histogram[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP] == 0)) // there are no ping-pong signatures in the region of the transposon
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
						stdDevOfArbitraryOverlaps += pow(transposon->histogram[overlap] - meanOfArbitraryOverlaps, 2);
				stdDevOfArbitraryOverlaps = sqrt(1.0 / (transposon->histogram.size() - 1 - 1 /* minus 1 for corrected sample STDDEV */) * stdDevOfArbitraryOverlaps);
				if (stdDevOfArbitraryOverlaps <= 1E-10)
					stdDevOfArbitraryOverlaps = 1E-10; // prevent division by 0, in case the STDDEV is 0

				// calculate significance of transposon score of ping-pong overlap vs. arbitrary overlaps
				double zValue = (transposon->histogram[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP] - meanOfArbitraryOverlaps) / stdDevOfArbitraryOverlaps;
				double pValue = 0;
				double degreesOfFreedom = MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP;
				for (double x = zValue; x <= zValue + APPROXIMATION_RANGE; x += APPROXIMATION_ACCURACY)
					pValue += APPROXIMATION_ACCURACY * tgamma((degreesOfFreedom+1)/2) / sqrt(degreesOfFreedom * M_PI) / tgamma(degreesOfFreedom/2) * pow(1 + x * x / degreesOfFreedom, -(degreesOfFreedom + 1)/2);

				transposon->pValue = pValue;
				// the normalized signature count is the number of ping-pong signatures per kilobase per million mapped reads
				transposon->normalizedSignatureCount = transposon->histogram[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP] / ((static_cast<float>(transposon->end) - transposon->start)/1000) / (totalReadCount/1000000);
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

// todo: description
void predictSuppressedTransposons(TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap, TTransposonsPerGenome &putativeTransposons, const float totalReadCount, TNameStore &bamNameStore)
{
	// define a putative transposon around every ping-pong signature
	for (TPingPongSignaturesPerGenome::iterator contig = pingPongSignaturesByOverlap[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP].begin(); contig != pingPongSignaturesByOverlap[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP].end() && contig->second.size() > 0; ++contig)
	{
		unsigned int putativeTransposonStart = contig->second.begin()->position;
		unsigned int putativeTransposonEnd = putativeTransposonStart + 1;
		for (TPingPongSignaturesPerContig::iterator pingPongSignature = contig->second.begin(); pingPongSignature != contig->second.end(); ++pingPongSignature)
		{
			if (putativeTransposonEnd + PREDICT_TRANSPOSONS_SLIDING_WINDOW >= pingPongSignature->position) // merge close-by windows
			{
				putativeTransposonEnd = pingPongSignature->position + 1;
			}
			else // the ping-pong signatures are so far apart, that they likely do not belong to the same transposon
			{
				if (putativeTransposonStart + PREDICT_TRANSPOSONS_LENGTH <= putativeTransposonEnd) // skip regions that are too short
				{
					// name putative transposon after genomic location
					stringstream putativeTransposonIdentifier;
					putativeTransposonIdentifier << bamNameStore[contig->first] << ":" << putativeTransposonStart << "-" << putativeTransposonEnd;

					// add transposon to list
					putativeTransposons[contig->first].push_back(TTransposon(putativeTransposonIdentifier.str(), STRAND_PLUS, putativeTransposonStart, putativeTransposonEnd));
				}

				// start a new region
				putativeTransposonStart = pingPongSignature->position;
				putativeTransposonEnd = putativeTransposonStart + 1;
			}
		}
	}
	// check putative transposons for ping-pong activity
	findSuppressedTransposons(pingPongSignaturesByOverlap, putativeTransposons, totalReadCount);
}

void writeTransposonsToFile(TTransposonsPerGenome &transposons, TNameStore &bamNameStore, bool browserTracks, string fileName)
{
	// open files to write transposon data to
	ofstream transposonsTSV((fileName + ".tsv").c_str(), ios_base::out);
	if (transposonsTSV.fail())
	{
		cerr << "Failed to create TSV file for transposons" << endl;
		return;
	}
	ofstream transposonsBED;
	if (browserTracks)
	{
		transposonsBED.open((fileName + ".bed").c_str(), ios_base::out);
		if (transposonsBED.fail())
		{
			cerr << "Failed to create browser track file for transposons" << endl;
			return;
		}
	}

	// use scientific formatting for floating point numbers in the output files
	transposonsTSV.setf(ios::scientific, ios::floatfield);
	if (browserTracks)
		transposonsBED.setf(ios::scientific, ios::floatfield);

	// write file headers
	transposonsTSV << "identifier\tstrand\tcontig\tstart\tend\tpValue\tqValue\tnormalizedSignatureCount" << endl;
	if (browserTracks)
	{
		// remove underscores (_) from fileName for the track name
		stringReplace(fileName, "_", " ");
		transposonsBED << "track name=\"" << fileName << "\" description=\"" << fileName << " shaded by ping-pong activity (1000 * (1 - corrected p-value))\" useScore=1 visibility=dense" << endl;
	}

	// write transposon data in TSV/BED format
	for (TTransposonsPerGenome::iterator contig = transposons.begin(); contig != transposons.end(); ++contig)
		for (TTransposonsPerContig::iterator transposon = contig->second.begin(); transposon != contig->second.end(); ++transposon)
		{
			transposonsTSV << transposon->identifier << '\t' << ((transposon->strand == STRAND_PLUS) ? '+' : '-') << '\t' << bamNameStore[contig->first] << '\t' << transposon->start << '\t' << transposon->end << '\t' << transposon->pValue << '\t' << transposon->qValue << '\t' << transposon->normalizedSignatureCount << endl;
			if (browserTracks)
				transposonsBED << bamNameStore[contig->first] << '\t' << transposon->start << '\t' << transposon->end << '\t' << transposon->identifier << '\t' << static_cast<int>(round((1 - transposon->qValue) * 1000)) << '\t' << ((transposon->strand == STRAND_PLUS) ? '+' : '-') << endl;
		}

	// close output files
	transposonsTSV.close();
	if (browserTracks)
		transposonsBED.close();
}

void generateTransposonsPlot(TTransposonsPerGenome &transposons, const string &fileName)
{
	int transposonCount = 0;
	for (TTransposonsPerGenome::iterator contig = transposons.begin(); contig != transposons.end(); ++contig)
		transposonCount += contig->second.size();

	THistograms histograms(transposonCount);
	vector< string > plotTitles(transposonCount);
	stringstream ss;
	unsigned int i = 0;
	for (TTransposonsPerGenome::iterator contig = transposons.begin(); contig != transposons.end(); ++contig)
		for (TTransposonsPerContig::iterator transposon = contig->second.begin(); transposon != contig->second.end(); ++transposon)
		{
			histograms[i] = transposon->histogram;
			ss
				<< "z-scores of transposon " << transposon->identifier << endl
				<< "(p-value for overlap of " << PING_PONG_OVERLAP << " nt = " << transposon->pValue << ")";
			plotTitles[i] = ss.str();
			ss.str("");
			i++;
		}
	plotHistogram(fileName, plotTitles, histograms);
}

// program entry point
int main(int argc, char const ** argv)
{
	// parse the command line options
	AppOptions options;
	if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
		return 1;

	TCountsGenome readCounts; // stats about positions where reads on the minus strand overlap with the 5' ends of reads on the plus strand

	float totalReadCount = 0; // number of reads read from the SAM/BAM file (depending on the command-line arguments, multi-hits may count less than 1)

	TNameStore bamNameStore; // structure to store contig names

	// read all BAM/SAM files
	if (options.verbosity >= 3)
		cerr << "Counting reads in SAM/BAM files" << endl;
	for(TInputFiles::iterator inputFile = options.inputFiles.begin(); inputFile != options.inputFiles.end(); ++inputFile)
	{
		stopwatch((string("  ") + toCString(*inputFile)).c_str(), options.verbosity);

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
			stopwatch((string("  ") + toCString(*transposonFile)).c_str(), options.verbosity);

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
		stopwatch("Rendering plots for z-scores of ping-pong signatures", options.verbosity);
		generateGroupedStackCountsPlot(groupedStackCountsByOverlap);
		stopwatch(options.verbosity);
	}
	groupedStackCountsByOverlap.clear();

	stopwatch("Writing ping-pong signatures to file", options.verbosity);
	writePingPongSignaturesToFile(pingPongSignaturesByOverlap[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP], bamNameStore, options.minStackHeight, options.browserTracks);
	stopwatch(options.verbosity);

	if (options.transposonFiles.size() > 0)
	{
		stopwatch("Checking input transposons for ping-pong activity", options.verbosity);
		findSuppressedTransposons(pingPongSignaturesByOverlap, transposons, totalReadCount);
		stopwatch(options.verbosity);
		stopwatch("Writing input transposons to file", options.verbosity);
		writeTransposonsToFile(transposons, bamNameStore, options.browserTracks, "transposons");
		stopwatch(options.verbosity);
		if (options.plot)
		{
			stopwatch("Rendering plots for z-scores of input transposons", options.verbosity);
			generateTransposonsPlot(transposons, "transposons_z-scores");
			stopwatch(options.verbosity);
		}
	}

	if (options.predictTransposons)
	{
		stopwatch("Predicting transposons based on ping-pong activity", options.verbosity);
		TTransposonsPerGenome putativeTransposons;
		predictSuppressedTransposons(pingPongSignaturesByOverlap, putativeTransposons, totalReadCount, bamNameStore);
		stopwatch(options.verbosity);
		stopwatch("Writing predicted transposons to file", options.verbosity);
		writeTransposonsToFile(putativeTransposons, bamNameStore, options.browserTracks, "predicted_transposons");
		stopwatch(options.verbosity);
		if (options.plot)
		{
			stopwatch("Rendering plots for z-scores of predicted transposons", options.verbosity);
			generateTransposonsPlot(putativeTransposons, "predicted_transposons_z-scores");
			stopwatch(options.verbosity);
		}

	}

	return 0;
}

