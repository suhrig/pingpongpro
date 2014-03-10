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

#if defined(WIN32) || defined(_WIN32) 
#define PATH_SEPARATOR '\\' 
#else 
#define PATH_SEPARATOR '/'
#endif

// type to hold the list of input files given as arguments to the program
typedef vector<CharString> TInputFiles;

// struct to store the options from the command line
struct AppOptions
{
	bool bedGraph;
	TInputFiles inputFiles;
	unsigned int minReadLength;
	unsigned int maxReadLength;
	unsigned int minStackHeight;
	CharString output;
	bool plot;
	unsigned int verbosity;

	AppOptions():
		bedGraph(false),
		minReadLength(1),
		maxReadLength(1000),
		minStackHeight(100),
		plot(false),
		verbosity(0)
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
	unsigned int reads: 31;
	bool UAt5PrimeEnd: 1;
	unsigned int heightScoreBin: 10; // must be big enough to hold HEIGHT_SCORE_HISTOGRAM_BINS
	unsigned int localHeightScoreBin: 7; // must be big enough to hold LOCAL_HEIGHT_SCORE_HISTOGRAM_BINS
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
#define MIN_ARBITRARY_OVERLAP 0
#define MAX_ARBITRARY_OVERLAP 20

// todo: description
#define HEIGHT_SCORE_HISTOGRAM_BINS 1000
#define LOCAL_HEIGHT_SCORE_HISTOGRAM_BINS 100
#define CLUSTER_SCORE_HISTOGRAM_BINS 100
typedef map< unsigned int, double > TScoreHistogram;
typedef map< unsigned int, TScoreHistogram > TScoreHistogramsByOverlap;

#define BASE_SCORE_PLUS_STRAND_URIDINE 0
#define BASE_SCORE_PLUS_STRAND_NOT_URIDINE 1
#define BASE_SCORE_MINUS_STRAND_URIDINE 2
#define BASE_SCORE_MINUS_STRAND_NOT_URIDINE 3

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
	addUsageLine(parser, "[\\fIOPTIONS\\fP] [-i \\fISAM_INPUT_FILE\\fP [-i ...]] [-o \\fIOUTPUT_DIRECTORY\\fP]");
	setShortDescription(parser, "Find ping-pong signatures like a pro");
	// todo: define long description
	addDescription(parser, "PingPongPro scans piRNA-Seq data for signs of ping-pong cycle activity. The ping-pong cycle produces piRNA molecules with complementary 5'-ends. These molecules appear as stacks of aligned reads whose 5'-ends overlap with the 5'-ends of reads on the opposite strand by exactly 10 bases.");
	setVersion(parser, "0.1");
	setDate(parser, "Mar 2014");

	addOption(parser, ArgParseOption("b", "bedgraph", "Output loci with ping-pong signature in bedGraph format. Default: off."));

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

	addOption(parser, ArgParseOption("o", "output", "Write output to specified directory. Default: current working directory.", ArgParseArgument::OUTPUTFILE, "PATH", true));

	addOption(parser, ArgParseOption("p", "plot", "Generate R plots on background noise estimation. Requires Rscript. Default: off."));

	addOption(parser, ArgParseOption("v", "verbose", "Print messages to stderr about the current progress. Default: off."));

	// parse command line
	ArgumentParser::ParseResult parserResult = parse(parser, argc, argv);
	if (parserResult != ArgumentParser::PARSE_OK)
		return parserResult;

	// extract options, if parsing was successful
	options.bedGraph = isSet(parser, "bedgraph");
	options.inputFiles.resize(getOptionValueCount(parser, "input")); // store input files in vector
	for (vector< string >::size_type i = 0; i < options.inputFiles.size(); i++)
	{
		getOptionValue(options.inputFiles[i], parser, "input", i);
		if (options.inputFiles[i] == "-")
			options.inputFiles[i] = "/dev/stdin";
	}
	if (options.inputFiles.size() == 0)
		options.inputFiles.push_back("/dev/stdin"); // read from stdin, if no input file is given
	getOptionValue(options.minStackHeight, parser, "min-stack-height");
	getOptionValue(options.minReadLength, parser, "min-read-length");
	getOptionValue(options.maxReadLength, parser, "max-read-length");
	if (options.minReadLength > options.maxReadLength)
	{
		cerr << getAppName(parser) << ": maximum read length (" << options.maxReadLength << ") must not be lower than minimum read length (" << options.minReadLength << ")" << endl;
		return ArgumentParser::PARSE_ERROR;
	}
	getOptionValue(options.output, parser, "output");
	if ((length(options.output) > 0) && (options.output[length(options.output)-1] != PATH_SEPARATOR))
		options.output += PATH_SEPARATOR; // append slash to output path, if missing
	options.plot = isSet(parser, "plot");
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

/*{
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

// convert absolute counts to relative frequencies
// todo: parameters
void getRelativeFrequencies(TScoreHistogramsByOverlap &scoreHistogramsByOverlap)
{
	for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
	{
		double total = 0;
		for (TScoreHistogram::iterator bar = scoreHistogramsByOverlap[overlap].begin(); bar != scoreHistogramsByOverlap[overlap].end(); ++bar)
			total += bar->second;
		for (TScoreHistogram::iterator bar = scoreHistogramsByOverlap[overlap].begin(); bar != scoreHistogramsByOverlap[overlap].end(); ++bar)
			bar->second /= total;
	}
}

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
							HEIGHT_SCORE_HISTOGRAM_BINS - 1 // invert scale
							- static_cast<int>(0.5 // add 0.5 for arithmetic rounding when casting double to int
							+ log10(heightScore) // take logarithm of score
							/ maxHeightScore * (HEIGHT_SCORE_HISTOGRAM_BINS - 1)); // assign every score to a bin
						scoreHistogramsByOverlap[overlap][bin]++; // increase bin counter
					}
				}
			}
		}
	}
	getRelativeFrequencies(scoreHistogramsByOverlap);
}


// todo: description
void calculateLocalHeightScoreHistograms(TCountsGenome &readCounts, TScoreHistogramsByOverlap &scoreHistogramsByOverlap)
{
	// iterate through all stacks
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

				if (max > 0)
				{
					mean /= (MAX_ARBITRARY_OVERLAP - MIN_ARBITRARY_OVERLAP + 1);

					for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
					{
						if (stackHeightsInVicinity[overlap] > 0)
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
	}
	getRelativeFrequencies(scoreHistogramsByOverlap);
}

// todo: description
void calculateBaseScoreHistograms(TCountsGenome &readCounts, TScoreHistogramsByOverlap &scoreHistogramsByOverlap)
{
	// iterate through all stacks
	for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
	{
		TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
		if (contigMinusStrand != readCounts[STRAND_MINUS].end())
		{
			for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
			{
				TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + PING_PONG_OVERLAP);
				if (positionMinusStrand != contigMinusStrand->second.end())
				{
					// count how many stacks on the plus strand have Uridine at the 5' end
					if (positionPlusStrand->second.UAt5PrimeEnd)
					{
						scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_PLUS_STRAND_URIDINE]++;
					}
					else
					{
						scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_PLUS_STRAND_NOT_URIDINE]++;
					}
					// count how many stacks on the minus strand have Uridine at the 5' end
					if (positionMinusStrand->second.UAt5PrimeEnd)
					{
						scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_MINUS_STRAND_URIDINE]++;
					}
					else
					{
						scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_MINUS_STRAND_NOT_URIDINE]++;
					}
				}
			}
		}
	}

	// convert absolute counts to relative counts
	scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_PLUS_STRAND_URIDINE] /= (scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_PLUS_STRAND_URIDINE] + scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_PLUS_STRAND_NOT_URIDINE]);
	scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_PLUS_STRAND_NOT_URIDINE] = 1 - scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_PLUS_STRAND_URIDINE];
	scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_MINUS_STRAND_URIDINE] /= (scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_MINUS_STRAND_URIDINE] + scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_MINUS_STRAND_NOT_URIDINE]);
	scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_MINUS_STRAND_NOT_URIDINE] = 1 - scoreHistogramsByOverlap[PING_PONG_OVERLAP][BASE_SCORE_MINUS_STRAND_URIDINE];

	// for arbitrary overlaps, assume a fixed frequency of 25% of having Uridine at the 5' end of reads
	for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
	{
		if (overlap == PING_PONG_OVERLAP)
			continue;

		scoreHistogramsByOverlap[overlap][BASE_SCORE_PLUS_STRAND_URIDINE] = 0.25;
		scoreHistogramsByOverlap[overlap][BASE_SCORE_PLUS_STRAND_NOT_URIDINE] = 1 - scoreHistogramsByOverlap[overlap][BASE_SCORE_PLUS_STRAND_URIDINE];
		scoreHistogramsByOverlap[overlap][BASE_SCORE_MINUS_STRAND_URIDINE] = 0.25;
		scoreHistogramsByOverlap[overlap][BASE_SCORE_MINUS_STRAND_NOT_URIDINE] = 1 - scoreHistogramsByOverlap[overlap][BASE_SCORE_MINUS_STRAND_URIDINE];
	}
}

// todo: description
void calculateClusterScoreHistograms(TCountsGenome &readCounts, TScoreHistogramsByOverlap &scoreHistogramsByOverlap)
{
	// iterate through all stacks
	for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
	{
		TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
		if (contigMinusStrand != readCounts[STRAND_MINUS].end())
		{
			for (int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
			{
				// the following loop moves a sliding window over the contig
				// and counts the number of stacks within the window surrounding
				// the stack in the center of the window
				const unsigned int slidingWindowSize = 1000;
				bool slidingWindow[slidingWindowSize] = {false}; // initialize sliding window with "false", i.e., set all positions in the window to "no stack"
				unsigned int stacksWithinWindow = 0; // keeps track of how many stacks there are in the window
				unsigned int slidingWindowPosition = 0; // current position on the genome
				for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
				{
					TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
					if (positionMinusStrand != contigMinusStrand->second.end())
					{
						while ((slidingWindowPosition < positionPlusStrand->first + slidingWindowSize/2) && (stacksWithinWindow > 0))
						{
							++slidingWindowPosition; // move sliding window by 1 nt
							if (slidingWindow[slidingWindowPosition % slidingWindowSize]) // a stack fell out of the sliding window
							{
								slidingWindow[slidingWindowPosition % slidingWindowSize] = false; // there is no stack at the position that entered the sliding window
								--stacksWithinWindow;
							}

							if (slidingWindow[(slidingWindowPosition + slidingWindowSize/2) % slidingWindowSize]) // there is a stack at the center of the sliding window
							{
								// find bin for score
								unsigned int bin =
									static_cast<int>(0.5 // add 0.5 for arithmetic rounding when casting double to int
									+ static_cast<double>(stacksWithinWindow) / slidingWindowSize // normalize cluster density to values between 0 and 1
									* (CLUSTER_SCORE_HISTOGRAM_BINS - 1) // assign every normalized score to a bin
									);
								// increase the counter for stacks falling into the bin returned by getHeightScoreBin()
								scoreHistogramsByOverlap[overlap][bin]++;
							}
						}
						slidingWindowPosition = positionPlusStrand->first;
						slidingWindow[slidingWindowPosition % slidingWindowSize] = true;
						++stacksWithinWindow;
					}
				}
			}
		}
	}
	getRelativeFrequencies(scoreHistogramsByOverlap);
}

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

//todo: x-axis labels
void plotHistograms(TScoreHistogramsByOverlap &scoreHistogramsByOverlap, const string &title, unsigned int binCount, vector< string > xAxisLabels, bool logScale = false)
{
	// generate an R script that produces a histogram plot
	string fileName = title;
	replace(fileName.begin(), fileName.end(), ' ', '_'); // replace blanks in file name
	ofstream rScript(toCString(fileName + ".R"));

	rScript << "histograms <- data.frame("; // store histogram counts in a data frame
	for (unsigned int overlap = MIN_ARBITRARY_OVERLAP; overlap <= MAX_ARBITRARY_OVERLAP; overlap++)
	{
			// print height of bars
			rScript << endl << "overlap" << overlap << "=c(";
			for (unsigned int bin = 0; bin < binCount; bin++)
			{
				if (bin % 10 == 0)
					rScript << endl; // insert a line-break every once in a while, because R cannot parse very long lines
				rScript << scoreHistogramsByOverlap[overlap][bin];
				if (bin < binCount - 1)
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
	{
		rScript
			<< "histograms <- log10(histograms)" << endl
			<< "lowestValue <- floor(min(histograms[histograms != -Inf]))" << endl
			<< "histograms <- histograms - lowestValue" << endl
			<< "plot(0, 0, xlim=c(0," << binCount << "), ylim=c(0, -lowestValue), type='n', xlab='" << title << "', ylab='log10(relative frequency)', xaxt='n', yaxt='n')" << endl
			<< "axis(2, at=0:-lowestValue, labels=lowestValue:0)" << endl;
	}
	else
	{
		rScript	<< "plot(0, 0, xlim=c(0," << binCount << "), ylim=c(0,max(histograms)), type='n', xlab='" << title << "', ylab='relative frequency', xaxt='n')" << endl;
	}

	// draw x-axis
	if (xAxisLabels.size() == 0)
	{
		// auto-generate x-axis based on quantiles, if no custom x-axis labels are given
		rScript << "axis(1, at=quantile(c(0," << binCount << "), probs = seq(0, 1, 0.2))+0.5, labels=quantile(c(0," << binCount << "), probs = seq(0, 1, 0.2)))" << endl;
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
		<< "	if (overlap != 10)" << endl
		<< "		barplot(histograms[,paste('overlap', overlap, sep='')], col=rgb(0,0,0,alpha=0.1), border=NA, axes=FALSE, add=TRUE, width=1, space=0)" << endl
		// draw a red line for ping-pong overlaps
		<< "for (bin in 1:" << binCount << ")" << endl
		<< "	lines(c(bin-1, bin), c(histograms$overlap10[bin], histograms$overlap10[bin]), type='l', col='red', lwd=2)" << endl
		// draw legend
		<< "legend(x='top', c('10 nt overlap', 'arbitrary overlaps'), col=c('red', 'black'), ncol=2, lwd=c(3,3), xpd=TRUE, inset=-0.1)" << endl
		<< "garbage <- dev.off()" << endl;

	// close R script
	rScript.close();

	// execute R script with "Rscript"
	string RCommand = "Rscript '" + fileName + ".R'";
	system(toCString(RCommand));
}
void plotHistograms(TScoreHistogramsByOverlap &scoreHistogramsByOverlap, const string &title, unsigned int binCount, bool logScale = false)
{
	vector< string > xAxisLabels;
	plotHistograms(scoreHistogramsByOverlap, title, binCount, xAxisLabels, logScale);
}
void plotHistograms(TScoreHistogramsByOverlap &scoreHistogramsByOverlap, const string &title, vector< string > xAxisLabels, bool logScale = false)
{
	plotHistograms(scoreHistogramsByOverlap, title, xAxisLabels.size(), xAxisLabels, logScale);
}


// program entry point
int main(int argc, char const ** argv)
{
	// parse the command line options
	AppOptions options;
	if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
		return 1;

	TCountsGenome readCounts; // stats about positions where reads on the minus strand overlap with the 5' ends of reads on the plus strand

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

		stopwatch(options.verbosity);
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

	stopwatch("Generating histograms of stack heights", options.verbosity);
	THeightScoreMap heightScoreMap;
	mapHeightsToScores(readCounts, heightScoreMap);
	TScoreHistogramsByOverlap heightScoreHistogramsByOverlap;
	calculateHeightScoreHistograms(readCounts, heightScoreMap, heightScoreHistogramsByOverlap);
	stopwatch(options.verbosity);

	stopwatch("Generating histograms of local stack heights", options.verbosity);
	TScoreHistogramsByOverlap localHeightScoreHistogramsByOverlap;
	calculateLocalHeightScoreHistograms(readCounts, localHeightScoreHistogramsByOverlap);
	stopwatch(options.verbosity);

	stopwatch("Generating histograms of base content of the 5' end of reads", options.verbosity);
	TScoreHistogramsByOverlap baseScoreHistogramsByOverlap;
	calculateBaseScoreHistograms(readCounts, baseScoreHistogramsByOverlap);
	stopwatch(options.verbosity);

	stopwatch("Generating histograms of clustering densities", options.verbosity);
	TScoreHistogramsByOverlap clusterScoreHistogramsByOverlap;
	calculateClusterScoreHistograms(readCounts, clusterScoreHistogramsByOverlap);
	stopwatch(options.verbosity);

	if (options.plot)
	{
		stopwatch("Generating R plots", options.verbosity);
		plotHistograms(heightScoreHistogramsByOverlap, "height score", HEIGHT_SCORE_HISTOGRAM_BINS, true);
		plotHistograms(localHeightScoreHistogramsByOverlap, "local height score", LOCAL_HEIGHT_SCORE_HISTOGRAM_BINS);
		vector< string > xAxisLabels;
		xAxisLabels.push_back("+ strand, U"); xAxisLabels.push_back("+ strand, not U"); xAxisLabels.push_back("- strand, U"); xAxisLabels.push_back("- strand, not U");
		plotHistograms(baseScoreHistogramsByOverlap, "base content at 5-prime end", xAxisLabels);
		plotHistograms(clusterScoreHistogramsByOverlap, "clustering density score", CLUSTER_SCORE_HISTOGRAM_BINS);
		stopwatch(options.verbosity);
	}

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
	stopwatch(options.verbosity);
*/
	return 0;
}


