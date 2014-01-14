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
	unsigned int minCoverage;

	AppOptions():
		verbosity(0),
		minReadLength(1),
		maxReadLength(1000),
		minCoverage(100)
	{}
};

// constants to refer to + and - strands throughout the program
const unsigned int STRAND_PLUS = 0;
const unsigned int STRAND_MINUS = 1;

// for every locus (position) on the genome the following attributes are calculated:
//  - reads: the number of reads which begin at this position
//  - readsWithAAtBase10: how many of these reads have an Adenine or a Uracil at base 10 (which is typical for piRNAs)
//  - readsWithNonTemplateBase: how many of these reads have terminal base which is different from the reference genome (also typical for piRNAs)
struct TCountsPosition
{
	unsigned int reads;
	unsigned int readsWithAAtBase10;
	unsigned int readsWithUAtBase10FromEnd;
	unsigned int readsWithNonTemplateBase;
	TCountsPosition():
		reads(0),
		readsWithAAtBase10(0),
		readsWithUAtBase10FromEnd(0),
		readsWithNonTemplateBase(0)
	{}
};

// The following types define nested arrays to store the above stats for every position in the genome.
// The stats are grouped by strand and contig/chromosome.
typedef map< unsigned int, TCountsPosition > TCountsContig;
typedef map< unsigned int, TCountsContig > TCountsStrand;
typedef TCountsStrand TCountsGenome[2];

// types to store @SQ header lines of BAM/SAM files
typedef seqan::StringSet<seqan::CharString> TNameStore;
typedef Iterator<TNameStore>::Type TNameStoreIterator;

// type to store a score for each stack height found in the input file
typedef map< unsigned int, int64_t > TStackScoreMap;

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

	addOption(parser, ArgParseOption("c", "min-coverage", "Omit loci with fewer than the specified number of mapped reads from the output.", ArgParseArgument::INTEGER, "NUMBER_OF_READS", true));
	setDefaultValue(parser, "min-coverage", options.minCoverage);
	setMinValue(parser, "min-coverage", "1");

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
	getOptionValue(options.minCoverage, parser, "min-coverage");
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
unsigned int stopwatch()
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
	return elapsedSeconds;
}

// Function which sums up the number of reads that start at a given position in the genome.
// Additionally, it counts the number of reads with an A/U at base 10 as well as
// the number of reads with a terminal base different from the reference genome.
// The latter two figures help identify the significance of a putative piRNA signature.
// Parameters:
//   bamFile: the BAM/SAM file from where to load the reads
//   readCountsUpstream: stats for positions were reads on the minus strand overlap the upstream (5') ends of reads on the plus strand
int countReadsInBamFile(BamStream &bamFile, TCountsGenome &readCountsUpstream, unsigned int minReadLength, unsigned int maxReadLength)
{
	TCountsPosition *positionUpstream;

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
			// calculate alignment length using CIGAR string, i.e., distance from first matching base on the reference to the last matching base on the reference
			size_t alignmentLength = 0;
			for (unsigned int cigarIndex = 0; cigarIndex < length(record.cigar); ++cigarIndex)
				if ((record.cigar[cigarIndex].operation == 'M') || (record.cigar[cigarIndex].operation == 'N') || (record.cigar[cigarIndex].operation == 'D') || (record.cigar[cigarIndex].operation == '=') || (record.cigar[cigarIndex].operation == 'X'))
					alignmentLength += record.cigar[cigarIndex].count;

			if (hasFlagRC(record)) // read maps to minus strand
			{
				positionUpstream = &(readCountsUpstream[STRAND_MINUS][record.rID][record.beginPos+alignmentLength]);
				if ((record.cigar[0].operation == 'S') && (record.cigar[0].count == 1))
					positionUpstream->readsWithNonTemplateBase++;
			}
			else // read maps to plus strand
			{
				positionUpstream = &(readCountsUpstream[STRAND_PLUS][record.rID][record.beginPos]);
				if ((record.cigar[length(record.cigar)-1].operation == 'S') && (record.cigar[length(record.cigar)-1].count == 1))
					positionUpstream->readsWithNonTemplateBase++;
			}

			// increase read counters for the given position on the genome
			positionUpstream->reads++;

			// check if 10th base is an Adenine or if the 10th base from the end of the read is a Uracil
			if (alignmentLength >= 10)
			{
				if ((record.seq[9] == 'A') || (record.seq[9] == 'a'))
					positionUpstream->readsWithAAtBase10++;
				if ((record.seq[length(record.seq)-10] == 'T') || (record.seq[length(record.seq)-10] == 't'))
					positionUpstream->readsWithUAtBase10FromEnd++;
			}
		}
	}

	return 0;
}

void calculateStackScoreMap(TCountsGenome &readCounts, TStackScoreMap &stackScoreMap)
{
	// iterate through all strands, contigs and positions to count how many stacks there are of a given height
	for (unsigned int strand = STRAND_PLUS; strand <= STRAND_MINUS; ++strand)
		for (TCountsStrand::iterator contig = readCounts[strand].begin(); contig != readCounts[strand].end(); ++contig)
			for (TCountsContig::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
				stackScoreMap[position->second.reads] += 1;
}

// find those positions on the genome where reads on opposite strands overlap by 10 nucleotides
int findOverlappingReads(ofstream &bedGraphFile, TCountsGenome &readCounts, const TNameStore &bamNameStore, unsigned int coverage, TStackScoreMap &stackScoreMap)
{
	const int overlap = 10;

	// write one bedGraph track per strand
	for (unsigned int strand = STRAND_PLUS; strand <= STRAND_MINUS; ++strand)
	{
		bedGraphFile << "track type=bedGraph name=\"stacks on + strand\" description=\"height of ping-pong stacks on the + strand\" visibility=full color=0,0,0 altColor=0,0,0 priority=20" << endl;

		// iterate through all strands, contigs and positions to find those positions where more than <coverage> reads on both strands overlap by 10 nucleotides
		for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
		{
			TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
			if (contigMinusStrand != readCounts[STRAND_MINUS].end())
			{
				for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
				{
					TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
					if (positionMinusStrand != contigMinusStrand->second.end())
					{
						if ((positionPlusStrand->second.reads >= coverage) && (positionMinusStrand->second.reads >= coverage))
						{
							/*cout
								<< bamNameStore[contigPlusStrand->first] << "\t"
								<< (positionPlusStrand->first+1) << "\t"
								<< STRAND_PLUS << "\t"
								<< positionPlusStrand->second.reads << "\t"
								<< positionPlusStrand->second.readsWithAAtBase10 << "\t"
								<< positionPlusStrand->second.readsWithUAtBase10FromEnd << "\t"
								<< positionPlusStrand->second.readsWithNonTemplateBase << "\t"
								<< STRAND_MINUS << "\t"
								<< positionMinusStrand->second.reads << "\t"
								<< positionMinusStrand->second.readsWithAAtBase10 << "\t"
								<< positionMinusStrand->second.readsWithUAtBase10FromEnd << "\t"
								<< positionMinusStrand->second.readsWithNonTemplateBase
								<< endl;*/
							bedGraphFile
								<< bamNameStore[contigPlusStrand->first] << " "
								<< positionPlusStrand->first << " "
								<< positionPlusStrand->first+1 << " "
								<< ((strand == STRAND_PLUS) ? positionPlusStrand->second.reads : positionMinusStrand->second.reads) << endl;
							//todo: test error handling
							if (bedGraphFile.bad())
							{
								cerr << "Failed to write bedGraph record" << endl;
								return 1;
							}
						}
					}
				}
			}
		}
	}

	bedGraphFile << "track type=bedGraph name=\"scores for ping-pong stacks\" description=\"scores for ping-pong stacks\" visibility=full color=0,0,0 altColor=0,0,0 priority=20" << endl;

	// iterate through all strands, contigs and positions to find those positions where more than <coverage> reads on both strands overlap by 10 nucleotides
	for (TCountsStrand::iterator contigPlusStrand = readCounts[STRAND_PLUS].begin(); contigPlusStrand != readCounts[STRAND_PLUS].end(); ++contigPlusStrand)
	{
		TCountsStrand::iterator contigMinusStrand = readCounts[STRAND_MINUS].find(contigPlusStrand->first);
		if (contigMinusStrand != readCounts[STRAND_MINUS].end())
		{
			for (TCountsContig::iterator positionPlusStrand = contigPlusStrand->second.begin(); positionPlusStrand != contigPlusStrand->second.end(); ++positionPlusStrand)
			{
				TCountsContig::iterator positionMinusStrand = contigMinusStrand->second.find(positionPlusStrand->first + overlap);
				if (positionMinusStrand != contigMinusStrand->second.end())
				{
					if ((positionPlusStrand->second.reads >= coverage) && (positionMinusStrand->second.reads >= coverage))
					{
						bedGraphFile
							<< bamNameStore[contigPlusStrand->first] << " "
							<< positionPlusStrand->first << " "
							<< positionPlusStrand->first+1 << " "
							<< stackScoreMap[positionPlusStrand->second.reads] * stackScoreMap[positionMinusStrand->second.reads] << endl;
						//todo: test error handling
						if (bedGraphFile.bad())
						{
							cerr << "Failed to write bedGraph record" << endl;
							return 1;
						}
					}
				}
			}
		}
	}

	return 0;
}

struct TAggregatedCountsPosition
{
	double numberOfPositions;
	double readsWithAAtBase10;
	double readsWithUAtBase10FromEnd;
	double readsWithNonTemplateBase;
	TAggregatedCountsPosition():
		numberOfPositions(0),
		readsWithAAtBase10(0),
		readsWithUAtBase10FromEnd(0),
		readsWithNonTemplateBase(0)
	{}
};

// find those positions on the genome where reads on opposite strands overlap by 10 nucleotides
void aggregateOverlappingReadsByCoverage(TCountsGenome &readCounts, const unsigned int upstreamStrand)
{
	// set downstream strand to the opposite of upstream strand
	unsigned int downstreamStrand = (upstreamStrand == STRAND_PLUS) ? STRAND_MINUS : STRAND_PLUS;

	unsigned int coverage = 1;
	
	TAggregatedCountsPosition aggregatedCounts[2];
	do
	{
		// reset counters
		aggregatedCounts[upstreamStrand] = aggregatedCounts[downstreamStrand] = TAggregatedCountsPosition();

		// iterate through all strands, contigs and positions to find those positions where more than <coverage> reads on opposite strands overlap by 10 nucleotides
		for (TCountsStrand::iterator contig = readCounts[downstreamStrand].begin(); contig != readCounts[downstreamStrand].end(); ++contig)
		{
			TCountsStrand::iterator contigOnOppositeStrand = readCounts[upstreamStrand].find(contig->first);
			if (contigOnOppositeStrand != readCounts[upstreamStrand].end())
			{
				for (TCountsContig::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
				{
					TCountsContig::iterator positionOnOppositeStrand = contigOnOppositeStrand->second.find(position->first + 10);
					if (positionOnOppositeStrand != contigOnOppositeStrand->second.end())
					{
						if ((position->second.reads >= coverage) && (positionOnOppositeStrand->second.reads >= coverage))
						{
//							if (position->second.reads < coverage+10)
//							{
								aggregatedCounts[downstreamStrand].numberOfPositions++;
								aggregatedCounts[downstreamStrand].readsWithAAtBase10 += (double) position->second.readsWithAAtBase10 / position->second.reads;
								aggregatedCounts[downstreamStrand].readsWithUAtBase10FromEnd += (double) position->second.readsWithUAtBase10FromEnd / position->second.reads;
								aggregatedCounts[downstreamStrand].readsWithNonTemplateBase += (double) position->second.readsWithNonTemplateBase / position->second.reads;
								aggregatedCounts[upstreamStrand].numberOfPositions++;
								aggregatedCounts[upstreamStrand].readsWithAAtBase10 += (double) positionOnOppositeStrand->second.readsWithAAtBase10 / positionOnOppositeStrand->second.reads;
								aggregatedCounts[upstreamStrand].readsWithUAtBase10FromEnd += (double) positionOnOppositeStrand->second.readsWithUAtBase10FromEnd / positionOnOppositeStrand->second.reads;
								aggregatedCounts[upstreamStrand].readsWithNonTemplateBase += (double) positionOnOppositeStrand->second.readsWithNonTemplateBase / positionOnOppositeStrand->second.reads;
//							}
						}
						else
						{
							// free memory of positions that we won't need again
							contig->second.erase(position);
							contigOnOppositeStrand->second.erase(positionOnOppositeStrand);
						}
					}
				}
			}
		}

		cout
			<< coverage << "\t"
			<< downstreamStrand << "\t"
			<< aggregatedCounts[downstreamStrand].numberOfPositions << "\t"
			<< aggregatedCounts[downstreamStrand].readsWithAAtBase10/aggregatedCounts[downstreamStrand].numberOfPositions << "\t"
			<< aggregatedCounts[downstreamStrand].readsWithUAtBase10FromEnd/aggregatedCounts[downstreamStrand].numberOfPositions << "\t"
			<< aggregatedCounts[downstreamStrand].readsWithNonTemplateBase/aggregatedCounts[downstreamStrand].numberOfPositions << "\t"
			<< upstreamStrand << "\t"
			<< aggregatedCounts[upstreamStrand].numberOfPositions << "\t"
			<< aggregatedCounts[upstreamStrand].readsWithAAtBase10/aggregatedCounts[upstreamStrand].numberOfPositions << "\t"
			<< aggregatedCounts[upstreamStrand].readsWithUAtBase10FromEnd/aggregatedCounts[upstreamStrand].numberOfPositions << "\t"
			<< aggregatedCounts[upstreamStrand].readsWithNonTemplateBase/aggregatedCounts[upstreamStrand].numberOfPositions
			<< endl;
		coverage += 10;
	} while (aggregatedCounts[downstreamStrand].numberOfPositions > 0);
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

	TCountsGenome readCountsUpstream; // stats about positions where reads on the minus strand overlap with the upstream ends of reads on the plus strand

	TNameStore bamNameStore;

	// read all BAM/SAM files
	for(TInputFiles::iterator inputFile = options.inputFiles.begin(); inputFile != options.inputFiles.end(); ++inputFile)
	{
		if (options.verbosity >= 3)
		{
			cerr << "Counting reads in SAM/BAM file " << toCString(*inputFile) << " ... ";
			stopwatch();
		}

		// open SAM/BAM file
		BamStream bamFile(toCString(*inputFile));
		if (!isGood(bamFile))
		{
			cerr << "Failed to open input file: " << *inputFile << endl;
			return 1;
		}

		// for every position in the genome, count the number of reads that start at a given position
		if (countReadsInBamFile(bamFile, readCountsUpstream, options.minReadLength, options.maxReadLength) != 0)
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

		if (options.verbosity >= 3)
			cerr << " done (" << stopwatch() << " seconds)" << endl;
	}

	if (options.verbosity >= 3)
	{
		cerr << "Calculating stack scores ... ";
		stopwatch();
	}

	TStackScoreMap stackScoreMap;
	calculateStackScoreMap(readCountsUpstream, stackScoreMap);

	if (options.verbosity >= 3)
		cerr << " done (" << stopwatch() << " seconds)" << endl;

	if (options.verbosity >= 3)
	{
		cerr << "Scanning for overlapping reads ... ";
		stopwatch();
	}

	ofstream bedGraphFile(toCString(options.outputBedGraph));
	if (bedGraphFile.fail())
	{
		cerr << "Failed to open bedGraph file: " << toCString(options.outputBedGraph) << endl;
		return 1;
	}
	// find positions where reads on the minus strand overlap with the upstream ends of reads on the plus strand
	if (findOverlappingReads(bedGraphFile, readCountsUpstream, bamNameStore, options.minCoverage, stackScoreMap) != 0)
		return 1;
	bedGraphFile.close();

/*	if (options.verbosity >= 3)
		cerr << " done (" << stopwatch() << " seconds)" << endl;

	if (options.verbosity >= 3)
	{
		cerr << "Aggregating overlapping reads by coverage ... ";
		stopwatch();
	}

	aggregateOverlappingReadsByCoverage(readCountsUpstream, STRAND_MINUS);

	if (options.verbosity >= 3)
		cerr << " done (" << stopwatch() << " seconds)" << endl;
*/
	return 0;
}


