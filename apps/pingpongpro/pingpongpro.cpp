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
#include <ctime>
#include <fstream>
#include <cmath>

using namespace std;
using namespace seqan;

// ==========================================================================
// Types & Classes
// ==========================================================================

// type to hold the list of input files given as arguments to the program
typedef vector<CharString> TInputFiles;

// struct to store the options from the command line
struct AppOptions
{
	TInputFiles inputFiles;
	CharString outputBedGraph;
	unsigned int verbosity;
	unsigned int minReadLength;
	unsigned int maxReadLength;
	unsigned int minStackHeight;

	AppOptions():
		verbosity(0),
		minReadLength(1),
		maxReadLength(1000),
		minStackHeight(100)
	{}
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
	uint32_t reads: 31;
	bool UAt5PrimeEnd: 1;
	TCountsPosition():
		reads(0),
		UAt5PrimeEnd(true)
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
#define MIN_ARBITRARY_OVERLAP 0
#define MAX_ARBITRARY_OVERLAP 20

// todo: description
#define HEIGHT_SCORE_HISTOGRAM_BINS 1000
#define LOCAL_HEIGHT_SCORE_HISTOGRAM_BINS 100
typedef map< unsigned int, unsigned int > TScoreHistogram;
typedef map< unsigned int, TScoreHistogram > TScoreHistogramsByOverlap;

typedef map< unsigned int, double > TScoresByBin;

// type to store a score for each stack height found in the input file
typedef map< unsigned int, double > THeightScoreMap;

// ==========================================================================
// Functions
// ==========================================================================

// function to parse command-line arguments
ArgumentParser::ParseResult parseCommandLine(AppOptions &options, int argc, char const ** argv)
{
	ArgumentParser parser("pingpongpro");

	// define usage and description
	addUsageLine(parser, "[\\fIOPTIONS\\fP] [-i \\fISAM_INPUT_FILE\\fP [-i ...]] -b \\fIBEDGRAPH_OUTPUT_FILE\\fP");
	setShortDescription(parser, "Find ping-pong signatures like a pro");
	// todo: define long description
	addDescription(parser, "pingpongpro scans piRNA-Seq data for signs of ping-pong cycle activity. The ping-pong cycle produces piRNA molecules with complementary 5'-ends. These molecules appear as stacks of aligned reads whose 5'-ends overlap with the 5'-ends of reads on the opposite strand by exactly 10 bases.");
	setVersion(parser, "0.1");
	setDate(parser, "Jan 2014");

	addOption(parser, ArgParseOption("b", "output-bedgraph", "Output loci with ping-pong signature to specified file in bedGraph format.", ArgParseArgument::OUTPUTFILE, "PATH", true));
	setRequired(parser, "output-bedgraph");

	addOption(parser, ArgParseOption("s", "min-stack-height", "Omit stacks with fewer than the specified number of reads from the output.", ArgParseArgument::INTEGER, "NUMBER_OF_READS", true));
	setDefaultValue(parser, "min-stack-height", options.minStackHeight);
	setMinValue(parser, "min-stack-height", "1");

	addOption(parser, ArgParseOption("i", "input", "Input file(s) in SAM/BAM format.", ArgParseArgument::INPUTFILE, "PATH", true));
	setDefaultValue(parser, "input", "-");
	vector< string > acceptedInputFormats;
	acceptedInputFormats.push_back(".sam");
	acceptedInputFormats.push_back(".bam");
	setValidValues(parser, "input", acceptedInputFormats);

	addOption(parser, ArgParseOption("l", "min-read-length", "Ignore reads in the input file that are shorter than the specified length.", ArgParseArgument::INTEGER, "LENGTH", true));
	setDefaultValue(parser, "min-read-length", options.minReadLength);
	setMinValue(parser, "min-read-length", "1");
	addOption(parser, ArgParseOption("L", "max-read-length", "Ignore reads in the input file that are longer than the specified length.", ArgParseArgument::INTEGER, "LENGTH", true));
	setDefaultValue(parser, "max-read-length", options.maxReadLength);
	setMinValue(parser, "max-read-length", "1");

	addOption(parser, ArgParseOption("v", "verbose", "Print messages to stderr  about the current progress. Default: off."));

	// parse command line
	ArgumentParser::ParseResult parserResult = parse(parser, argc, argv);
	if (parserResult != ArgumentParser::PARSE_OK)
		return parserResult;

	// extract options, if parsing was successful
	getOptionValue(options.outputBedGraph, parser, "output-bedgraph");
	getOptionValue(options.minStackHeight, parser, "min-stack-height");
	getOptionValue(options.minReadLength, parser, "min-read-length");
	getOptionValue(options.maxReadLength, parser, "max-read-length");
	if (options.minReadLength > options.maxReadLength)
	{
		cerr << getAppName(parser) << ": maximum read length (" << options.maxReadLength << ") must not be lower than minimum read length (" << options.minReadLength << ")" << endl;
		return ArgumentParser::PARSE_ERROR;
	}
	options.inputFiles.resize(getOptionValueCount(parser, "input"));
	for (vector< string >::size_type i = 0; i < options.inputFiles.size(); i++)
		getOptionValue(options.inputFiles[i], parser, "input", i);
	if (isSet(parser, "verbose"))
		options.verbosity = 3;

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
		if (elapsedSeconds > 0)
			cerr << "done (" << elapsedSeconds << " seconds)" << endl;
		else
			cerr << operation << " ... ";
	}
	return elapsedSeconds;
}

// Function which sums up the number of reads that start at a given position in the genome.
// Additionally, it counts the number of reads with uridine at the 5' end.
// Parameters:
//   bamFile: the BAM/SAM file from where to load the reads
//   readCounts: stats for positions were reads on the minus strand overlap the 5' ends of reads on the plus strand
// todo: document parameters
int countReadsInBamFile(BamStream &bamFile, TCountsGenome &readCounts, unsigned int minReadLength, unsigned int maxReadLength)
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

		//todo: weighted counting of multimapped reads
		if ((record.beginPos != BamAlignmentRecord::INVALID_POS) && (record.beginPos != -1) && // skip unmapped reads
		    (length(record.seq) >= minReadLength) && (length(record.seq) <= maxReadLength)) // skip reads which are not within the specified length range

		{
			// calculate start of alignment using CIGAR string
			size_t clippedBases = 0;
			if (record.cigar[0].operation == 'S')
				clippedBases = record.cigar[0].count;

			// calculate length of alignment using CIGAR string
			size_t alignmentLength = 0;
			for (unsigned int cigarIndex = 0; cigarIndex < length(record.cigar); ++cigarIndex)
			{
				if ((record.cigar[cigarIndex].operation == 'M') || (record.cigar[cigarIndex].operation == 'N') || (record.cigar[cigarIndex].operation == 'D') || (record.cigar[cigarIndex].operation == '=') || (record.cigar[cigarIndex].operation == 'X')) // these CIGAR elements indicate alignment
					alignmentLength += record.cigar[cigarIndex].count;
			}

			if (hasFlagRC(record)) // read maps to minus strand
			{
				position = &(readCounts[STRAND_MINUS][record.rID][record.beginPos+alignmentLength]); // get a pointer to counter of the position of the read
				if ((record.seq[clippedBases+alignmentLength] == 'A') || (record.seq[clippedBases+alignmentLength] == 'a')) // check if last base is an A
					position->UAt5PrimeEnd = true;
			}
			else // read maps to plus strand
			{
				position = &(readCounts[STRAND_PLUS][record.rID][record.beginPos]); // get a pointer to the counter of the position of the read
				if ((record.seq[clippedBases] == 'T') || (record.seq[clippedBases] == 't')) // check if first base is a U
					position->UAt5PrimeEnd = true;
			}

			// increase read counter for the given position on the genome
			position->reads++;
		}
	}

	return 0;
}

// todo: description
/*void calculateFDR(TScoreHistogramsByOverlap &scoreHistogramsByOverlap, TScoresByBin &scoresByBin, unsigned int bins, unsigned int testOverlap)
{
	map< unsigned int, double> means;
	map< unsigned int, double> stddevs;

	// calculate mean for every bin
	for (unsigned int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
		if (overlap != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
			if (overlap != testOverlap) // skip stacks with the overlap that we are testing
				for (unsigned int bin = 0; bin < bins; bin++)
					means[bin] += scoreHistogramsByOverlap[overlap][bin];
	for (unsigned int bin = 0; bin < bins; bin++)
		means[bin] /= MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1;

	// calculate standard deviation for every bin
	for (unsigned int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
		if (overlap != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
			if (overlap != testOverlap) // skip stacks with the overlap that we are testing
				for (unsigned int bin = 0; bin < bins; bin++)
					stddevs[bin] += pow(scoreHistogramsByOverlap[overlap][bin] - means[bin], 2);
	for (unsigned int bin = 0; bin < bins; bin++)
		stddevs[bin] = sqrt(1 / (MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP) * stddevs[bin]);

	// calculate score for each bin, based on divergence from mean and variance
	for (unsigned int bin = 0; bin < bins; bin++)
	{
		double mean = means[bin];
		double stddev = stddevs[bin];
			unsigned int value = scoreHistogramsByOverlap[testOverlap][bin];
		scoresByBin[bin] = statisticalSignificance(mean, stddev, value) * value / (value + mean);
	}
}*/

// convert stack heights to scores
// todo: document parameters
void mapHeightsToScores(TCountsGenome &readCounts, THeightScoreMap &heightScoreMap)
{
	// iterate through all strands, contigs and positions to count how many stacks there are of any given height
	for (unsigned int strand = STRAND_PLUS; strand <= STRAND_MINUS; ++strand)
		for (TCountsStrand::iterator contig = readCounts[strand].begin(); contig != readCounts[strand].end(); ++contig)
			for (TCountsContig::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
				heightScoreMap[position->second.reads] += 1;
}

// todo: description
void calculateHeightScoreHistograms(TCountsGenome &readCounts, THeightScoreMap &heightScoreMap, TScoreHistogramsByOverlap &scoreHistogramsByOverlap)
{
	// the highest possible score is that of two overlapping stacks which both have a height of just 1 read
	const double maxHeightScore = log10(heightScoreMap[1] * heightScoreMap[1]);

	// iterate through all strands, contigs and positions to find those positions where a stack on the plus strand overlaps a stack on the minus strand by <overlap> nt
	for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
	{
		TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
		if (contigMinusStrand != readCounts[STRAND_MINUS].end())
		{
			for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
			{
				for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
				{
					TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
					if (positionMinusStrand != contigMinusStrand->second.end())
					{
						// calculate score based on heights of overlapping stacks
						double heightScore = heightScoreMap[positionPlusStrand->second.reads] * heightScoreMap[positionMinusStrand->second.reads];
						// find the bin for the score
						unsigned int bin =
							static_cast<int>(0.5 // add 0.5 for arithmetic rounding when casting double to int
							+ log10(heightScore) // take logarithm of score
							/ maxHeightScore * (HEIGHT_SCORE_HISTOGRAM_BINS - 1)); // assign every score to a bin
						scoreHistogramsByOverlap[overlap][bin]++; // increase bin counter
					}
				}
			}
		}
	}
}

// todo: description
void calculateLocalHeightScoreHistograms(TCountsGenome &readCounts, TScoreHistogramsByOverlap &scoreHistogramsByOverlap)
{
	// iterate through all strands, contigs and positions to find those positions where a stack on the plus strand overlaps a stack on the minus strand by <overlap> nt
	for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
	{
		TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
		if (contigMinusStrand != readCounts[STRAND_MINUS].end())
		{
			for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
			{
				map< int, unsigned int > stackHeightsInVicinity;
				double mean = 0;
				unsigned int max = 0;
				for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
				{
					TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
					if (positionMinusStrand != contigMinusStrand->second.end())
					{
						stackHeightsInVicinity[overlap] = positionMinusStrand->second.reads;
						mean += positionMinusStrand->second.reads;
						if (positionMinusStrand->second.reads > max)
							max = positionMinusStrand->second.reads;
					}
				}
				mean /= (MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1);

				for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
				{
					// calculate score based on how much higher the stack is compared to the stacks in the vicinity
					double localHeightScore = (stackHeightsInVicinity[overlap] - (mean - stackHeightsInVicinity[overlap]/(MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1))) / max;
					// find bin for score
					unsigned int bin =
						static_cast<int>(0.5 // add 0.5 for arithmetic rounding when casting double to int
						+ (localHeightScore + 1) / 2 // score can take values between -1 and +1 => normalize score to values between 0 and 1
						* (LOCAL_HEIGHT_SCORE_HISTOGRAM_BINS - 1) // assign every normalized score to a bin
						);
					// increase the counter for stacks falling into the bin returned by getHeightScoreBin()
					scoreHistogramsByOverlap[overlap][bin]++;
				}
			}
		}
	}
}

/*int findOverlappingReads(ofstream &bedGraphFile, TCountsGenome &readCounts, const TNameStore &bamNameStore, unsigned int minStackHeight, THeightScoreMap &heightScoreMap, TCombinedStackScoreMap &combinedStackScoreMap)
{
	for (unsigned int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
	{

		// iterate through all strands, contigs and positions to find those positions where more than <minStackHeight> reads on both strands overlap by 10 nucleotides
		for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
		{
			TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
			if (contigMinusStrand != readCounts[STRAND_MINUS].end())
			{
			const unsigned int slidingWindowSize = 1000;
			bool slidingWindow[slidingWindowSize] = {false};
			unsigned int stacksWithinWindow = 0;
			unsigned int slidingWindowPosition = 0;

				for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
				{
					TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
					if (positionMinusStrand != contigMinusStrand->second.end())
					{
					
					while ((slidingWindowPosition < positionPlusStrand->first) && (stacksWithinWindow > 0))
					{
						++slidingWindowPosition; // move sliding window by 1 nt
						if (slidingWindow[slidingWindowPosition % slidingWindowSize] == true) // a stack fell out of the sliding window
						{
							slidingWindow[slidingWindowPosition % slidingWindowSize] = false; // there is no stack at the position that entered the sliding window
							--stacksWithinWindow;
						}
						if (slidingWindow[(slidingWindowPosition + slidingWindowSize/2) % slidingWindowSize] == true) // there is a stack in the center of the sliding window
							bedGraphFile << slidingWindowPosition - slidingWindowSize/2 - 1 << " " << stacksWithinWindow << endl; // output the number of stacks surrounding the center stack
					}
					slidingWindowPosition = positionPlusStrand->first;
					slidingWindow[slidingWindowPosition % slidingWindowSize] = true;
					++stacksWithinWindow;
					//todo: stacks at contig ends are ignored

					}
				}
			}
		}
	}

	return 0;
}*/

/*unsigned int getTotalScoreBin(double heightScore)
{
	double bin;
	return static_cast<int>(bin + 0.5); // +0.5 ensures that float variable is rounded correctly when casting
}*/

// find those positions on the genome where reads on opposite strands overlap by 10 nucleotides
// todo: document parameters
/*int findOverlappingReads(ofstream &bedGraphFile, TCountsGenome &readCounts, const TNameStore &bamNameStore, unsigned int minStackHeight, THeightScoreMap &heightScoreMap, TCombinedStackScoreMap &combinedStackScoreMap)
{
	// find those positions where reads on both strands overlap by <overlap> nucleotides
	for (unsigned int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
	{
		if (overlap == PING_PONG_OVERLAP)
			continue;

		// iterate through all contigs on the plus strand
		for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
		{
			// check if there are any stacks for the given contig on the minus strand
			TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
			if (contigMinusStrand != readCounts[STRAND_MINUS].end())
			{
				// iterate through all stacks on the plus strand
				for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
				{
					// check if there is a stack on the minus strand which overlaps with the stack on the plus strand by <overlap> nucleotides
					TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
					if (positionMinusStrand != contigMinusStrand->second.end())
					{
					}
				}
			}
		}
	}

	return 0;
}*/

//todo: support log of y-axis
//todo: x-axis labels
void plotHistograms(TScoreHistogramsByOverlap &scoreHistogramsByOverlap, unsigned int binCount, string fileName)
{
	ofstream rScript(toCString(fileName + ".R"));
	for (unsigned int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
	{
		rScript << "png(\"" << fileName << "_" << overlap << ".png\")" << endl;
		rScript << "barplot(c(";
		//if (overlap != PING_PONG_OVERLAP) // ignore ping-pong stacks, since they would skew the result
			for (unsigned int bin = 0; bin < binCount - 1; bin++)
			{
				if (bin % 100 == 0)
					rScript << endl; // insert a line-break every once in a while, because R cannot parse very long lines
				rScript << scoreHistogramsByOverlap[overlap][bin] << ", ";
			}
		rScript << scoreHistogramsByOverlap[overlap][binCount - 1] << "))" << endl;
		rScript << "dev.off()" << endl;
	}
	rScript.close();
}

// program entry point
int main(int argc, char const ** argv)
{
	// parse the command line
	AppOptions options;
	if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
		return 1;

	// read from stdin, if no input file is given
	if (options.inputFiles.size() == 0)
		options.inputFiles.push_back("/dev/stdin");
	//todo: scan input files for "-" and replace with "/dev/stdin"

	if (options.outputBedGraph == '-')
		options.outputBedGraph = "/dev/stdout";

	TCountsGenome readCounts; // stats about positions where reads on the minus strand overlap with the 5' ends of reads on the plus strand

	TNameStore bamNameStore;

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
		if (countReadsInBamFile(bamFile, readCounts, options.minReadLength, options.maxReadLength) != 0)
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

		stopwatch("", options.verbosity);
	}

	stopwatch("Calculating scores for stack heights", options.verbosity);
	THeightScoreMap heightScoreMap;
	mapHeightsToScores(readCounts, heightScoreMap);
	stopwatch("", options.verbosity);

	stopwatch("Generating histograms of stack heights", options.verbosity);
	TScoreHistogramsByOverlap heightScoreHistogramsByOverlap;
	calculateHeightScoreHistograms(readCounts, heightScoreMap, heightScoreHistogramsByOverlap);
	plotHistograms(heightScoreHistogramsByOverlap, HEIGHT_SCORE_HISTOGRAM_BINS, "height_score");
	stopwatch("", options.verbosity);

	stopwatch("Generating histograms of local stack heights", options.verbosity);
	TScoreHistogramsByOverlap localHeightScoreHistogramsByOverlap;
	calculateLocalHeightScoreHistograms(readCounts, localHeightScoreHistogramsByOverlap);
	plotHistograms(localHeightScoreHistogramsByOverlap, LOCAL_HEIGHT_SCORE_HISTOGRAM_BINS, "local_height_score");
	stopwatch("", options.verbosity);

/*	stopwatch("Scanning for overlapping reads", options.verbosity);
	ofstream bedGraphFile(toCString(options.outputBedGraph));
	if (bedGraphFile.fail())
	{
		cerr << "Failed to open bedGraph file: " << toCString(options.outputBedGraph) << endl;
		return 1;
	}
	// find positions where reads on the minus strand overlap with the 5' ends of reads on the plus strand
	if (findOverlappingReads(bedGraphFile, readCounts, bamNameStore, options.minStackHeight, heightScoreMap, combinedStackScoreMap) != 0)
		return 1;
	bedGraphFile.close();
	stopwatch("", options.verbosity);
*/
	return 0;
}


