// ==========================================================================
//				pingpongpro
// ==========================================================================
/*

Copyright (c) 2014, Sebastian Uhrig

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

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

// options how to handle multi-mapped reads for calculation of read stack heights
// multiHitsWeighted = multi-hits are counted as 1 / number of hits
// multiHitsDiscard = multi-hits are counted as 0
// multiHitsUnique = multi-hits are counted as 1 (i.e., no distinction is made between unique hits and multi-hits)
enum TCountMultiHits { multiHitsWeighted, multiHitsDiscard, multiHitsUnique };

// constants for various output and input file formats
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
	unsigned int predictTransposonsRange;
	unsigned int verbosity;
};

// types to store @SQ header lines of BAM/SAM files
typedef StringSet<CharString> TNameStore;
typedef Iterator<TNameStore>::Type TNameStoreIterator;

// constants to refer to + and - strands throughout the program
const unsigned int STRAND_PLUS = 0;
const unsigned int STRAND_MINUS = 1;

// for every position on the genome the following attributes are calculated:
//  - reads: the number of reads which begin at this position
//  - UAt5PrimeEnd: whether the reads of the stack have adenine at position 10
struct TReadStack
{
	float reads;
	bool AAtPosition10;
	TReadStack():
		reads(0),
		AAtPosition10(false)
	{}
};

// The following types define nested arrays to store the above stats for every position in the genome.
// The stats are grouped by strand and contig/chromosome.
typedef map< unsigned int, TReadStack > TReadStacksPerContig;
typedef map< unsigned int, TReadStacksPerContig > TReadStacksPerStrand;
typedef TReadStacksPerStrand TReadStacksPerGenome[2];

// true ping-pong stacks overlap by this many nt
const int PING_PONG_OVERLAP = 10;

// stacks with overlaps between MIN_ARBITRARY_OVERLAP and MAX_ARBITRARY_OVERLAP (except for PING_PONG_OVERLAP)
// are used to estimate what is background noise
const int MIN_ARBITRARY_OVERLAP = 3;
const int MAX_ARBITRARY_OVERLAP = 23;

// type to store a score for each stack height found in the input file
typedef map< unsigned int, float > THeightScoreMap;

// a ping-pong signature is defined as two stacks of reads, which are on opposite strands
// and overlap by 10 nt at their 5' ends
// the following type stores various attributes for every found signature
struct TPingPongSignature
{
	unsigned int position: 32; // position on contig where the ping-pong overlap is located
	unsigned int heightScoreBin: 30; // must be big enough to hold <HEIGHT_SCORE_BINS>
	unsigned int localHeightScoreBin: 1; // holds either <IS_ABOVE_COVERAGE> or <IS_BELOW_COVERAGE>
	unsigned int baseBiasBin: 1; // holds either <HAS_BASE_BIAS> or <HAS_NO_BASE_BIAS>
	float readsOnPlusStrand; // stack height on + strand
	float readsOnMinusStrand; // stack height on - strand
	float fdr; // chances of this putative ping-pong overlap being a false discovery

	// constructor to initialize with values
	TPingPongSignature(unsigned int position, unsigned int heightScoreBin, unsigned int localHeightScoreBin, unsigned int baseBiasBin, float readsOnPlusStrand, float readsOnMinusStrand):
		position(position), heightScoreBin(heightScoreBin), localHeightScoreBin(localHeightScoreBin), baseBiasBin(baseBiasBin), readsOnPlusStrand(readsOnPlusStrand), readsOnMinusStrand(readsOnMinusStrand), fdr(0)
	{
	}

	// copy constructor (needed to add instance to a STL list)
	TPingPongSignature(const TPingPongSignature &pingPongSignature)
	{
		position = pingPongSignature.position;
		heightScoreBin = pingPongSignature.heightScoreBin;
		localHeightScoreBin = pingPongSignature.localHeightScoreBin;
		baseBiasBin = pingPongSignature.baseBiasBin;
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

// ping-pong signatures are grouped and counted by the following criteria
// - the height of the overlapping stacks
// - whether the reads have adenine at position 10
// - whether the height of the stacks are above or below the local coverage
const unsigned int HAS_BASE_BIAS = 0; // bin for signatures with adenine at position 10
const unsigned int HAS_NO_BASE_BIAS = 1; // bin for signatures with another base than adenine at position 10
const unsigned int IS_ABOVE_COVERAGE = 0; // bin for signatures with stacks that are higher than the local coverage
const unsigned int IS_BELOW_COVERAGE = 1; // bin for signatures with stacks that are lower than the local coverage
const unsigned int HEIGHT_SCORE_BINS = 1000; // divide signatures by height into this many bins
typedef vector< vector< vector< float > > > TGroupedStackCounts;
typedef vector< TGroupedStackCounts > TGroupedStackCountsByOverlap;

// types to plot histograms
typedef vector< float > THistogram; // every vector element represents the height of a bar
typedef vector< THistogram > THistograms; // a collection of histograms, which are printed into a single PDF

// type to store attributes about transposons
struct TTransposon
{
	string identifier; // name of the transposon
	unsigned int strand;
	unsigned int start;
	unsigned int end;
	float pValue; // probability that there is ping-pong activity in the region of the transposon
	float qValue; // multiple-testing corrected <pValue>
	float readsOnPlusStrand;
	float readsOnMinusStrand;
	THistogram histogram; // number of ping-pong signatures within the transposon region for every overlap between <MIN_ARBITRARY_OVERLAP> and <MAX_ARBITRARY_OVERLAP>

	// constructor to initialize with values
	TTransposon(string identifier, unsigned int strand, unsigned int start, unsigned int end):
		identifier(identifier), strand(strand), start(start), end(end), pValue(1), qValue(1), readsOnPlusStrand(0), readsOnMinusStrand(0)
	{
	}

	// copy constructor (needed to add instance to a STL list)
	TTransposon(const TTransposon &transposon)
	{
		identifier = transposon.identifier;
		strand = transposon.strand;
		start = transposon.start;
		end = transposon.end;
		pValue = transposon.pValue;
		qValue = transposon.qValue;
		histogram = transposon.histogram;
		readsOnPlusStrand = transposon.readsOnPlusStrand;
		readsOnMinusStrand = transposon.readsOnMinusStrand;
	}

	// operator for sorting of a list of transposons by genomic position
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

// parameters for transposon prediction based on ping-pong activity
const unsigned int PREDICT_TRANSPOSONS_MIN_LENGTH = 30; // predicted transposons shorter than this are discarded

// the calculation of p-values is based on integrals
// the following constants are parameters for the precision of p-values calculation
const double APPROXIMATION_ACCURACY = 0.01; // step size with which integrals are calculated; smaller means more accurate
const double APPROXIMATION_RANGE = 5; // span of integral calculation; wider means more accurate
const double MIN_STANDARD_DEVIATION = 1E-10; // if the STDDEV is smaller than this, assume this fixed value to avoid division by 0

// ==========================================================================
// Functions
// ==========================================================================

// function to parse command-line arguments
// Input parameters:
//	argc: number of command-line arguments as passed to the function <main>
//	argv: array of command-line arguments as passed to the function <main>
// Output parameters:
//	options: parsed options
// Return value: status code about whether the command-line could be parsed
ArgumentParser::ParseResult parseCommandLine(AppOptions &options, int argc, char const ** argv)
{
	ArgumentParser parser("pingpongpro");

	// define usage and description
	addUsageLine(parser, "[\\fIOPTIONS\\fP] [-i \\fISAM_INPUT_FILE\\fP [-i ...]] [-o \\fIOUTPUT_DIRECTORY\\fP]");
	setShortDescription(parser, "Find ping-pong signatures like a pro");
	addDescription(parser, "PingPongPro scans piRNA-Seq data for signs of ping-pong cycle activity. The ping-pong cycle produces piRNA molecules with complementary 5'-ends. These molecules appear as stacks of aligned reads whose 5'-ends overlap with the 5'-ends of reads on the opposite strand by exactly 10 bases.");
	setVersion(parser, "1.0");
	setDate(parser, "Apr 2014");

	// define parameters
	addOption(parser, ArgParseOption("b", "browserTracks", "Generate genome browser tracks for loci with ping-pong signature and (if -t or -T is specified) for transposons with ping-pong activity. Default: \\fIoff\\fP."));

	addOption(parser, ArgParseOption("s", "min-stack-height", "Omit stacks with fewer than the specified number of reads from the output.", ArgParseArgument::INTEGER, "NUMBER_OF_READS"));
	setDefaultValue(parser, "min-stack-height", 0);
	setMinValue(parser, "min-stack-height", "0");

	addOption(parser, ArgParseOption("i", "input", "Input file(s) in SAM/BAM format. \"-\" means stdin.", ArgParseArgument::INPUTFILE, "PATH", true));
	setDefaultValue(parser, "input", "-");
	setValidValues(parser, "input", ".bam .sam -");

	addOption(parser, ArgParseOption("l", "min-alignment-length", "Ignore alignments in the input file that are shorter than the specified length.", ArgParseArgument::INTEGER, "LENGTH"));
	setDefaultValue(parser, "min-alignment-length", 24);
	setMinValue(parser, "min-alignment-length", "1");
	addOption(parser, ArgParseOption("L", "max-alignment-length", "Ignore alignments in the input file that are longer than the specified length.", ArgParseArgument::INTEGER, "LENGTH"));
	setDefaultValue(parser, "max-alignment-length", 32);
	setMinValue(parser, "max-alignment-length", "1");

	addOption(parser, ArgParseOption("m", "multi-hits", "How to count multi-mapping reads.", ArgParseArgument::STRING, "METHOD"));
	setDefaultValue(parser, "multi-hits", "weighted");
	setValidValues(parser, "multi-hits", "weighted discard unique");

	addOption(parser, ArgParseOption("o", "output", "Write output to specified directory. Default: current working directory.", ArgParseArgument::OUTPUTFILE, "PATH"));

	addOption(parser, ArgParseOption("p", "plot", "Generate R plots on how z-scores are calculated for ping-pong signatures and (if -t or -T is specified) for transposons. Requires Rscript. Default: \\fIoff\\fP."));

	addOption(parser, ArgParseOption("t", "transposons", "Check if the transposons given in the file \\fIPATH\\fP are suppressed through ping-pong activity.", ArgParseArgument::INPUTFILE, "PATH", true));
	setValidValues(parser, "transposons", ".bed .csv .gff .gtf .tsv");

	addOption(parser, ArgParseOption("T", "predict-transposons", "Predict the location of suppressed transposons based on regions with high ping-pong activity. Consider adjacent ping-pong signatures within a range of \\fIRANGE\\fP to belong to the same transposon. Default: \\fIoff\\fP.", ArgParseArgument::INTEGER, "RANGE"));
	stringstream ss;
	ss << PREDICT_TRANSPOSONS_MIN_LENGTH;
	setMinValue(parser, "predict-transposons", ss.str());

	addOption(parser, ArgParseOption("v", "verbose", "Print messages about the current progress to stderr. Default: \\fIoff\\fP."));

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

	if (isSet(parser, "predict-transposons"))
	{
		getOptionValue(options.predictTransposonsRange, parser, "predict-transposons");
	}
	else
	{
		options.predictTransposonsRange = 0;
	}

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

// Function to measure time between the first and second invocation of the function.
// Input parameters:
//	operation: a description of the task being measured
//	verbosity: if <operation> is not empty, the description is printed to stderr, if the verbosity level is >= INFO
//	           if <operation is empty, the function prints the number of seconds since its last invocation
// Return value: the number of seconds since the last invocation of the function
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

// Function which finds stacks of reads in a BAM file.
// Input parameters:
//	bamFile: the BAM/SAM file from where to load the reads
//	minAlignmentLength: reads in the <bamFile> which are shorter than this are ignored
//	maxAlignmentLength: reads in the <bamFile> which are longer than this are ignored
//	countMultiHits: how to count multi-mapped reads (see declaration of TCountMultiHits)
// Output parameters:
//	readStacks: stacks of reads that were found by the function
// Return value: 1, if the <bamFile> could not be read; 0 otherwise
int countReadsInBamFile(BamStream &bamFile, TReadStacksPerGenome &readStacks, const unsigned int minAlignmentLength, const unsigned int maxAlignmentLength, TCountMultiHits countMultiHits)
{
	TReadStack *position;

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
				position = &(readStacks[STRAND_MINUS][record.rID][record.beginPos+alignmentLength]);

				// check if base at position 10 is adenine
				size_t clippedBasesAt5PrimeEnd = 0;
				if ((length(record.cigar) > 1) && (record.cigar[length(record.cigar)-1].operation == 'S'))
					clippedBasesAt5PrimeEnd = record.cigar[length(record.cigar)-1].count;
				if ((record.seq[length(record.seq)-clippedBasesAt5PrimeEnd-1-9] == 'T') || (record.seq[length(record.seq)-clippedBasesAt5PrimeEnd-1-9] == 't')) // check if 10th base is adenine (we check for uracil, because reads on the - strand are stored as the complement in SAM files
					position->AAtPosition10 = true;
			}
			else // read maps to plus strand
			{
				// get a pointer to counter of the position of the read
				position = &(readStacks[STRAND_PLUS][record.rID][record.beginPos]);

				// check if base at position 10 is adenine
				size_t clippedBasesAt5PrimeEnd = 0;
				if (record.cigar[0].operation == 'S')
					clippedBasesAt5PrimeEnd = record.cigar[0].count;
				if ((record.seq[clippedBasesAt5PrimeEnd+9] == 'A') || (record.seq[clippedBasesAt5PrimeEnd+9] == 'a'))
					position->AAtPosition10 = true;
			}

			// increase stack height
			position->reads += readWeight;
		}
	}

	return 0;
}

// Function, which converts every stack height into a score.
// The score is directly based on how often a stack with a certain height is found in the input files of the program.
// Therefore, the score maps every stack height to the empirical frequency of encountering such a stack in the input dataset.
// Input parameters:
//	readStacks: the read stacks that were found by the function <countReadsInBamFile>
// Output parameters:
//	heightScoreMap: output of the function
//	                a mapping of [stack height -> empirical frequency of stacks with this height]
void mapHeightsToScores(TReadStacksPerGenome &readStacks, THeightScoreMap &heightScoreMap)
{
	// iterate through all strands, contigs and positions to count how many stacks there are of any given height
	for (unsigned int strand = STRAND_PLUS; strand <= STRAND_MINUS; ++strand)
		for (TReadStacksPerStrand::iterator contig = readStacks[strand].begin(); contig != readStacks[strand].end(); ++contig)
			for (TReadStacksPerContig::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
				heightScoreMap[0.5 + position->second.reads] += 1;
}

// Function, which groups read stacks by all possible combinations of the following criteria:
// - the height of the overlapping stacks
// - whether the reads have adenine at position 10
// - whether the height of the stacks are above or below the local coverage
// For every group, the number of stacks falling into that particular group is counted.
// Input parameters:
//	readStacks: the read stacks that were found by the function <countReadsInBamFile>
//	            the variable is emptied by the function to conserve memory
//	heightScoreMap: a mapping of [stack height -> empirical frequency of stacks with this height] as produced by the function <mapHeightsToScores>
// Output parameters:
//	groupedStackCountsByOverlap: for every overlap between <MIN_ARBITRARY_OVERLAP> and <MAX_ARBITRARY_OVERLAP>, the number of read stacks falling into all possible groups
//	pingPongSignaturesByOverlap: for every overlap between <MIN_ARBITRARY_OVERLAP> and <MAX_ARBITRARY_OVERLAP>, the ping-pong signatures that were found
void countStacksByGroup(TReadStacksPerGenome &readStacks, THeightScoreMap &heightScoreMap, TGroupedStackCountsByOverlap &groupedStackCountsByOverlap, TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap)
{
	// the following loop initializes a multi-dimensional array of stack counts with the following boundaries:
	// MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1 (one for each possible overlap)
	// HEIGHT_SCORE_BINS (one of each bin of the height scores)
	// 2 (one for reads with adenine at position 10 and one for those with a different base)
	// 2 (one for stack heights below the local coverage and one for stack heights above)
	groupedStackCountsByOverlap.resize(MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1);
	for (TGroupedStackCountsByOverlap::iterator i = groupedStackCountsByOverlap.begin(); i != groupedStackCountsByOverlap.end(); ++i)
	{
		i->resize(HEIGHT_SCORE_BINS);
		for (TGroupedStackCounts::iterator j = i->begin(); j != i->end(); ++j)
		{
			j->resize(2);
			for (vector< vector< float > >::iterator k = j->begin(); k != j->end(); ++k)
			{
				k->resize(2);
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
	for (TReadStacksPerStrand::iterator contigPlusStrand = readStacks[STRAND_PLUS].begin(); contigPlusStrand != readStacks[STRAND_PLUS].end(); ++contigPlusStrand)
	{
		TReadStacksPerStrand::iterator contigMinusStrand = readStacks[STRAND_MINUS].find(contigPlusStrand->first);
		if (contigMinusStrand != readStacks[STRAND_MINUS].end())
		{
			for (TReadStacksPerContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
			{
				vector< TReadStacksPerContig::iterator > stacksOnMinusStrand(MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1, contigMinusStrand->second.end());
				float meanStackHeightInVicinity = 0;
				float maxStackHeightInVicinity = 0;
				for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
				{
					TReadStacksPerContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
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

							// calculate score based on whether the stack has adenine at position 10
							unsigned int baseBiasBin = (positionPlusStrand->second.AAtPosition10 || stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP]->second.AAtPosition10) ? HAS_BASE_BIAS : HAS_NO_BASE_BIAS;

							// increase bin counter
							groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][heightScoreBin][baseBiasBin][localHeightScoreBin]++;

							// keep a list of putative ping-pong signatures, so we can analyze later, which of them are (likely) true
							pingPongSignaturesByOverlap[overlap - MIN_ARBITRARY_OVERLAP][contigPlusStrand->first].push_back(TPingPongSignature(positionPlusStrand->first, heightScoreBin, localHeightScoreBin, baseBiasBin, positionPlusStrand->second.reads, stacksOnMinusStrand[overlap - MIN_ARBITRARY_OVERLAP]->second.reads));
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
	for (TReadStacksPerStrand::iterator contigMinusStrand = readStacks[STRAND_MINUS].begin(); contigMinusStrand != readStacks[STRAND_MINUS].end(); ++contigMinusStrand)
		contigMinusStrand->second.clear();
}

// The groups of stacks as produced by the function <countStacksByGroup> may be empty.
// This function merges adjacent groups until there are no empty groups left.
// Input/output parameters:
// 	groupedStackCountsByOverlap: the grouped stack counts as produced by the function <countStacksByGroup>
//	pingPongSignaturesByOverlap: the ping-pong signatures found by the function <countStacksByGroup>
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
			for (vector< vector< float > >::iterator j = (*i)[collapsedBin].begin(); j != (*i)[collapsedBin].end(); ++j)
				for (vector< float >::iterator k = j->begin(); k != j->end(); ++k)
					*k = 0;

		unsigned int emptyBins;
		do {
			emptyBins = 0;

			// collapse bins
			for (unsigned int overlap = 0; overlap < groupedStackCountsByOverlap.size(); overlap++)
				for (unsigned int i = 0; i < groupedStackCountsByOverlap[overlap][bin].size(); i++)
					for (unsigned int j = 0; j < groupedStackCountsByOverlap[overlap][bin][i].size(); j++)
					{
						collapsed[overlap][collapsedBin][i][j] += groupedStackCountsByOverlap[overlap][bin][i][j];

						// check if collapsedBin is still empty to decide whether to collapse even more
						if (collapsed[overlap][collapsedBin][i][j] <= 0)
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

// This function assigns a FDR to every ping-pong stacks based on how often ping-pong stacks
// with the properties of a given ping-pong stack occur by chance compared to how often they occur
// when looking at overlaps of 10 nt.
// Input paramters:
// 	groupedStackCountsByOverlap: the collapsed grouped stack counts as modified by the function <collapseBins>
// Input/output parameters:
//	pingPongSignaturesByOverlap: the ping-pong signatures as modified by the function <collapseBins>
void calculateFDRs(TGroupedStackCountsByOverlap &groupedStackCountsByOverlap, TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap)
{
	TGroupedStackCountsByOverlap FDRs = groupedStackCountsByOverlap; // the assignment shall only ensure that <FDRs> has the same dimensions as <groupedStackCountsByOverlap>

	// estimate how many of the putative ping-pong signatures are random noise
	for (unsigned int i = 0; i < groupedStackCountsByOverlap.begin()->size(); i++)
		for (unsigned int j = 0; j < (*groupedStackCountsByOverlap.begin())[i].size(); j++)
			for (unsigned int k = 0; k < (*groupedStackCountsByOverlap.begin())[i][j].size(); k++)
			{
				// calculate mean for all the bins of arbitrary overlaps
				float meanOfArbitraryOverlaps = 0;
				for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
					if (overlap != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
						meanOfArbitraryOverlaps += groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k];
				meanOfArbitraryOverlaps = meanOfArbitraryOverlaps / (groupedStackCountsByOverlap.size() - 1 /* minus the one bin for ping-pong overlaps */);

				// calculate standard deviation for all the bins of arbitrary overlaps
				float stdDevOfArbitraryOverlaps = 0;
				for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
					if (overlap != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
						stdDevOfArbitraryOverlaps += pow(groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k] - meanOfArbitraryOverlaps, 2);
				stdDevOfArbitraryOverlaps = sqrt(1.0 / (groupedStackCountsByOverlap.size() - 1 - 1 /* minus 1 for corrected sample STDDEV */) * stdDevOfArbitraryOverlaps);
				if (stdDevOfArbitraryOverlaps <= MIN_STANDARD_DEVIATION)
					stdDevOfArbitraryOverlaps = MIN_STANDARD_DEVIATION; // prevent division by 0, in case the STDDEV is 0

				// calculate expected number of false positives among the putative ping-pong signatures
				for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
				{
					double fdr = 0;
					for (double x = -APPROXIMATION_RANGE; x <= +APPROXIMATION_RANGE; x += APPROXIMATION_ACCURACY)
					{
						if (x >= (groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k] - meanOfArbitraryOverlaps) / stdDevOfArbitraryOverlaps)
						{
							fdr += APPROXIMATION_ACCURACY * 1/sqrt(2*M_PI)*exp(-0.5*x*x) * 1;
						}
						else
						{
							fdr += APPROXIMATION_ACCURACY * 1/sqrt(2*M_PI)*exp(-0.5*x*x) * (x * stdDevOfArbitraryOverlaps + meanOfArbitraryOverlaps) / groupedStackCountsByOverlap[overlap - MIN_ARBITRARY_OVERLAP][i][j][k];
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

					FDRs[overlap - MIN_ARBITRARY_OVERLAP][i][j][k] = fdr;
				}
			}

	// assign a FDR to every putative ping-pong signature
	for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
		for (TPingPongSignaturesPerGenome::iterator contig = pingPongSignaturesByOverlap[overlap - MIN_ARBITRARY_OVERLAP].begin(); contig != pingPongSignaturesByOverlap[overlap - MIN_ARBITRARY_OVERLAP].end(); ++contig)
			for (TPingPongSignaturesPerContig::iterator pingPongSignature = contig->second.begin(); pingPongSignature != contig->second.end(); ++pingPongSignature)
				pingPongSignature->fdr = FDRs[overlap - MIN_ARBITRARY_OVERLAP][pingPongSignature->heightScoreBin][pingPongSignature->baseBiasBin][pingPongSignature->localHeightScoreBin];
}

// function to replace all occurrences of a string within a string for another string
// Input/output paramters:
//	subjectString: the string to modify
// Input parameters:
//	searchString: the string to search for
//	replaceString: the string that replaces <searchString>
void stringReplace(string &subjectString, const string &searchString, const string &replaceString) {
	size_t i = 0;
	while((i = subjectString.find(searchString, i)) != string::npos)
	{
	        subjectString.replace(i, searchString.length(), replaceString);
        	i += replaceString.length(); // skip the string that we just inserted
	}
}

// This function uses Rscript to generate histogram plots.
// Multiple plots are written to a single PDF.
// Input parameters:
//	fileName: the name of the R script and PDF file to be generated, without the file extension
//	titles: the titles of all histogram plots
//	histograms: a collection of histograms to plot
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
		<< "sds <- ifelse(sds < " << MIN_STANDARD_DEVIATION << "," << MIN_STANDARD_DEVIATION << ", sds)" << endl
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

// function to write ping-pong signatures found by the function <countStacksByGroup> to a TSV file
// Input parameters:
//	pingPongSignaturesPerGenome: the ping-pong signatures to write to a file as found by the function <countStacksByGroup>
//	bamNameStore: mapping of numeric contig IDs to human-readable names
//	minStackHeight: ping-pong signatures with a smaller stack height than this are omitted from the output
//	browserTracks: if set to true, then a bedGraph file is generated in addition to the TSV file
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
		scoresBedGraph << "track type=bedGraph name=\"scores\" description=\"scores of ping-pong signatures (1 - FDR)\" visibility=full viewLimits=0.0:1.0 autoScale=off" << endl;
	}

	// write a line for each ping-pong signature
	for (TPingPongSignaturesPerGenome::iterator contig = pingPongSignaturesPerGenome.begin(); contig != pingPongSignaturesPerGenome.end(); ++contig)
		for (TPingPongSignaturesPerContig::iterator pingPongSignature = contig->second.begin(); pingPongSignature != contig->second.end(); ++pingPongSignature)
			if ((pingPongSignature->readsOnPlusStrand >= minStackHeight) && (pingPongSignature->readsOnMinusStrand >= minStackHeight))
			{
				signaturesTSV
					<< bamNameStore[contig->first] << '\t'
					<< pingPongSignature->position << '\t'
					<< pingPongSignature->fdr << '\t'
					<< pingPongSignature->readsOnPlusStrand << '\t'
					<< pingPongSignature->readsOnMinusStrand << endl;
				if (browserTracks)
				{
					readsOnPlusStrandBedGraph
						<< bamNameStore[contig->first] << '\t'
						<< pingPongSignature->position << '\t'
						<< (pingPongSignature->position+1) << '\t'
						<< pingPongSignature->readsOnPlusStrand << endl;
					readsOnMinusStrandBedGraph
						<< bamNameStore[contig->first] << '\t'
						<< pingPongSignature->position << '\t'
						<< (pingPongSignature->position+1) << '\t'
						<< pingPongSignature->readsOnMinusStrand << endl;
					scoresBedGraph
						<< bamNameStore[contig->first] << '\t'
						<< pingPongSignature->position << '\t'
						<< (pingPongSignature->position+1) << '\t'
						<< (1-pingPongSignature->fdr) << endl;
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

// generate plots that illustrate the difference in stack counts by overlap
// Input parameters:
//	groupedStackCountsByOverlap: grouped stack counts as processed by the function <collapseBins>
void generateGroupedStackCountsPlot(TGroupedStackCountsByOverlap &groupedStackCountsByOverlap)
{
	// find out in how many bins the stacks were grouped,
	// because we need to generate a histogram for every bin
	int binCount =
		groupedStackCountsByOverlap.begin()->size() *
		groupedStackCountsByOverlap.begin()->begin()->size() *
		groupedStackCountsByOverlap.begin()->begin()->begin()->size();
	THistograms histograms(binCount);
	vector< string > plotTitles(binCount);

	// collect histogram values and generate plot names
	unsigned int x = 0;
	stringstream ss;
	for (int i = groupedStackCountsByOverlap.begin()->size() - 1; i >= 0; i--)
		for (unsigned int j = 0; j < (*groupedStackCountsByOverlap.begin())[i].size(); j++)
			for (unsigned int k = 0; k < (*groupedStackCountsByOverlap.begin())[i][j].size(); k++)
			{
				histograms[x].resize(groupedStackCountsByOverlap.size());
				for (unsigned int overlap = 0; overlap < groupedStackCountsByOverlap.size(); overlap++)
				{
					histograms[x][overlap] = groupedStackCountsByOverlap[overlap][i][j][k];

					ss << "z-scores of signatures with the following properties:" << endl;
					ss << "stack height score of " << (groupedStackCountsByOverlap.begin()->size() - 1 - i) << endl;
					if (j != HAS_BASE_BIAS)
						ss << "no ";
					ss << "adenine at position 10" << endl;
					ss << "stack height " << ((k == IS_ABOVE_COVERAGE) ? "above" : "below") << " the local coverage";
					plotTitles[x] = ss.str();
					ss.str("");
				}
				x++;
			}

	// render histograms
	plotHistogram("ping-pong_signature_z-scores", plotTitles, histograms);
}

// this functions reads genomic regions of transposons from a file
// the transposons are checked for ping-pong activity by the function <findSuppressedTransposons>
// Input parameters:
//	transposonFile: the file from where to read the transposons
//	fileFormat: the format of the file
// Input/output parameters:
//	bamNameStore: a mapping of numeric contig IDs to human readable names
//	              the name store is extended by names that it does not contain, but that are used in the transposon file
// Output parameters:
//	transposons: the transposons read from the file
void readTransposonsFromFile(ifstream &transposonFile, TFileFormat fileFormat, TTransposonsPerGenome &transposons, TNameStore &bamNameStore)
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

	unsigned int fieldNumber = 1; // index of the column that is currently being read from the input file
	string fieldValue = ""; // every column in the input file is first read into this variable

	// when a field has been fully read, it is assigned to one of the following transposon attributes:
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

		if ((character == delimiter) && !((fileFormat == fileFormatCSV) && quotesOpen)) // encountered delimiter character => assign current field to appropriate transposon attribute and then, begin a new field
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
					for (unsigned int i = 0; (transposonContig == -1) && (i < length(bamNameStore)); i++)
						if (bamNameStore[i] == fieldValue)
							transposonContig = i;

					if (transposonContig == -1) // the contig was not found in the name store
					{
						// add a new element to the name store
						appendValue(bamNameStore, fieldValue);
						transposonContig = length(bamNameStore) - 1;
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
// Input parameters:
// 	transposon1, transposon2: the transposons to compare
// Return value: true, if transposon1 has a lower genomic coordinate than transposon2; false otherwise
inline bool compareTransposonsByPValue(const TTransposonsPerContig::iterator &transposon1, const TTransposonsPerContig::iterator &transposon2)
{
	return transposon1->pValue < transposon2->pValue;
}

// This function checks given transposons for ping-pong activity.
// A transposon is assumed to be suppressed by ping-pong activity, if there are significantly more ping-pong signatures within its region
// than there are arbitrary signatures.
// Input parameters:
//	pingPongSignaturesByOverlap: the ping-pong signtures found by function <countStacksByGroup>
// Input/output parameters
//	transposons: a list of transposons to check for ping-pong activity
//	             every transposon is assigned a p-value and a q-value indicating the statistical significance of ping-pong activity
void findSuppressedTransposons(TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap, TTransposonsPerGenome &transposons)
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

				// sum up the scores of all signatures (ping-pong or arbitrary) within the transposon region
				float sumOfScores = 0;
				if (positionByOverlap[overlap]->position >= transposon->start)
				{
					while ((positionByOverlap[overlap]->position <= transposon->end) && (positionByOverlap[overlap] != pingPongSignaturesByOverlap[overlap][contig->first].end()))
					{
						// sum up scores of all signatures (ping-pong or arbitrary) within the transposon region
						sumOfScores += (1 - positionByOverlap[overlap]->fdr);
						++(positionByOverlap[overlap]);

						// sum up the number of reads on each strand (for ping-pong overlaps only)
						if (static_cast<int>(overlap) + MIN_ARBITRARY_OVERLAP == PING_PONG_OVERLAP)
						{
							transposon->readsOnPlusStrand += positionByOverlap[overlap]->readsOnPlusStrand;
							transposon->readsOnMinusStrand += positionByOverlap[overlap]->readsOnMinusStrand;
						}
					}
				}

				transposon->histogram[overlap] = sumOfScores;
			}

			// calculate mean transposon score of all arbitrary overlaps
			float meanOfArbitraryOverlaps = 0;
			for (unsigned int overlap = 0; overlap < transposon->histogram.size(); overlap++)
				if (static_cast<int>(overlap) + MIN_ARBITRARY_OVERLAP != PING_PONG_OVERLAP) // ignore ping-pong overlaps in the mean calculation, since they would skew the result
					meanOfArbitraryOverlaps += transposon->histogram[overlap];
			meanOfArbitraryOverlaps = meanOfArbitraryOverlaps / (positionByOverlap.size() - 1 /* minus the one bin for ping-pong overlaps */);

			if ((meanOfArbitraryOverlaps == 0) && (transposon->histogram[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP] == 0)) // there are no ping-pong signatures in the region of the transposon
			{
				transposon->pValue = 1;
				transposon->qValue = 1;
			}
			else
			{
				// calculate standard deviation of transposon score of all arbitrary overlaps
				float stdDevOfArbitraryOverlaps = 0;
				for (unsigned int overlap = 0; overlap < positionByOverlap.size(); overlap++)
					if (static_cast<int>(overlap) + MIN_ARBITRARY_OVERLAP != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
						stdDevOfArbitraryOverlaps += pow(transposon->histogram[overlap] - meanOfArbitraryOverlaps, 2);
				stdDevOfArbitraryOverlaps = sqrt(1.0 / (transposon->histogram.size() - 1 - 1 /* minus 1 for corrected sample STDDEV */) * stdDevOfArbitraryOverlaps);
				if (stdDevOfArbitraryOverlaps <= MIN_STANDARD_DEVIATION)
					stdDevOfArbitraryOverlaps = MIN_STANDARD_DEVIATION; // prevent division by 0, in case the STDDEV is 0

				// calculate significance of transposon score of ping-pong overlap vs. arbitrary overlaps
				double zValue = (transposon->histogram[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP] - meanOfArbitraryOverlaps) / stdDevOfArbitraryOverlaps;
				double pValue = 0;
				for (double x = zValue; x <= zValue + APPROXIMATION_RANGE; x += APPROXIMATION_ACCURACY)
					pValue += APPROXIMATION_ACCURACY * 1/sqrt(2*M_PI)*exp(-0.5*x*x);

				transposon->pValue = pValue;
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

// Similar to the function <findSuppressedTransposons>, this function transposons for ping-pong activity.
// In contrast to <findSuppressedTransposons>, this function does not take transposons as an input argument,
// but tries to find transposons automatically based on where there is a lot of ping-pong activity.
// Input parameters:
//	pingPongSignaturesByOverlap: the ping-pong signtures found by function <countStacksByGroup>
//      bamNameStore: a mapping of numeric contig IDs to human readable names
//	range: ping-pong signatures that are this close to one another are considered to belong to the same transposon
// Output parameters:
//	putativeTransposons: putative transposons that were found by the function, with p- and q-values
void predictSuppressedTransposons(TPingPongSignaturesByOverlap &pingPongSignaturesByOverlap, TTransposonsPerGenome &putativeTransposons, TNameStore &bamNameStore, unsigned int range)
{
	// define a putative transposon around every ping-pong signature
	for (TPingPongSignaturesPerGenome::iterator contig = pingPongSignaturesByOverlap[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP].begin(); contig != pingPongSignaturesByOverlap[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP].end() && contig->second.size() > 0; ++contig)
	{
		unsigned int putativeTransposonStart = contig->second.begin()->position;
		unsigned int putativeTransposonEnd = putativeTransposonStart + 1;
		for (TPingPongSignaturesPerContig::iterator pingPongSignature = contig->second.begin(); pingPongSignature != contig->second.end(); ++pingPongSignature)
		{
			if (putativeTransposonEnd + range >= pingPongSignature->position) // merge close-by windows
			{
				putativeTransposonEnd = pingPongSignature->position + 1;
			}
			else // the ping-pong signatures are so far apart, that they likely do not belong to the same transposon
			{
				if (putativeTransposonStart + PREDICT_TRANSPOSONS_MIN_LENGTH < putativeTransposonEnd) // skip regions that are too short
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
	findSuppressedTransposons(pingPongSignaturesByOverlap, putativeTransposons);
}

// Function to write transposons to a TSV file.
// Input paramters:
//	transposons: the transposons to write to the file
//      bamNameStore: a mapping of numeric contig IDs to human readable names
//	browserTracks: if set to true, then a BED file is generated in addition to the TSV file
//	fileName: the name of the file that the transposons are written to, without the file extension
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
	transposonsTSV << "identifier\tstrand\tcontig\tstart\tend\tpValue\tqValue\tsignatureCount\tsignatureCountPerKB\tstrandRatio" << endl;
	if (browserTracks)
	{
		// remove underscores (_) from fileName for the track name
		stringReplace(fileName, "_", " ");
		transposonsBED << "track name=\"" << fileName << "\" description=\"" << fileName << " shaded by ping-pong activity (1000 * (1 - q-value))\" useScore=1 visibility=dense" << endl;
	}

	// write transposon data in TSV/BED format
	for (TTransposonsPerGenome::iterator contig = transposons.begin(); contig != transposons.end(); ++contig)
		for (TTransposonsPerContig::iterator transposon = contig->second.begin(); transposon != contig->second.end(); ++transposon)
		{
			transposonsTSV
				<< transposon->identifier << '\t'
				<< ((transposon->strand == STRAND_PLUS) ? '+' : '-') << '\t'
				<< bamNameStore[contig->first] << '\t'
				<< transposon->start << '\t'
				<< transposon->end << '\t'
				<< transposon->pValue << '\t'
				<< transposon->qValue << '\t'
				<< transposon->histogram[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP] << '\t'
				<< transposon->histogram[PING_PONG_OVERLAP - MIN_ARBITRARY_OVERLAP] / ((static_cast<float>(transposon->end) - transposon->start)/1000) << '\t'
				<< ((transposon->readsOnMinusStrand > 0) ? transposon->readsOnPlusStrand/transposon->readsOnMinusStrand : 1) << endl;
			if (browserTracks)
				transposonsBED
					<< bamNameStore[contig->first] << '\t'
					<< transposon->start << '\t'
					<< transposon->end << '\t'
					<< transposon->identifier << '\t'
					<< static_cast<int>(round((1 - transposon->qValue) * 1000)) << '\t'
					<< ((transposon->strand == STRAND_PLUS) ? '+' : '-') << endl;
		}

	// close output files
	transposonsTSV.close();
	if (browserTracks)
		transposonsBED.close();
}

// generate plots that illustrate the statistical significance of ping-pong activity for a list of transposons
// Input paramters:
//	transposons: a list of transposons; a plot is generated for each of them
//	fileName: name of the file that the plots are written to
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

	TReadStacksPerGenome readStacks; // stats about positions where reads on the minus strand overlap with the 5' ends of reads on the plus strand

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
		if (countReadsInBamFile(bamFile, readStacks, options.minAlignmentLength, options.maxAlignmentLength, options.countMultiHits) != 0)
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
	mapHeightsToScores(readStacks, heightScoreMap);
	TGroupedStackCountsByOverlap groupedStackCountsByOverlap;
	TPingPongSignaturesByOverlap pingPongSignaturesByOverlap;
	countStacksByGroup(readStacks, heightScoreMap, groupedStackCountsByOverlap, pingPongSignaturesByOverlap);
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
		findSuppressedTransposons(pingPongSignaturesByOverlap, transposons);
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

	if (options.predictTransposonsRange > 0)
	{
		stopwatch("Predicting transposons based on ping-pong activity", options.verbosity);
		TTransposonsPerGenome putativeTransposons;
		predictSuppressedTransposons(pingPongSignaturesByOverlap, putativeTransposons, bamNameStore, options.predictTransposonsRange);
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

