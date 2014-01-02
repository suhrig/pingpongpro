// ==========================================================================
//				pingpongpal
// ==========================================================================
// todo: copyright

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/arg_parse.h>

#include <map>
#include <vector>
#include <ctime>

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

	AppOptions()
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
	int reads;
	int readsWithAAtBase10;
	int readsWithUAtBase10FromEnd;
	int readsWithNonTemplateBase;
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

// ==========================================================================
// Functions
// ==========================================================================

// function to parse command-line arguments
ArgumentParser::ParseResult parseCommandLine(AppOptions &options, int argc, char const ** argv)
{
	ArgumentParser parser("pingpongpal");

	// define usage and description
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
	setShortDescription(parser, "Find ping-pong signatures with your ping-pong pal.");
	// todo: define long description
	addDescription(parser, "This is the application skelleton and you should modify this string.");
	setVersion(parser, "0.1");
	setDate(parser, "Jan 2014");
	// todo: Add Examples Section.
	addTextSection(parser, "Examples");
	addListItem(parser, "\\fBpingpongpal\\fP \\fB-v\\fP \\fItext\\fP", "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

	addOption(parser, ArgParseOption("i", "input", "Path to the input file(s)", ArgParseArgument::INPUTFILE, "IN", true));
	vector< string > acceptedInputFormats;
	acceptedInputFormats.push_back(".sam");
	acceptedInputFormats.push_back(".bam");
	setValidValues(parser, "input", acceptedInputFormats);

	// parse command line
	ArgumentParser::ParseResult parserResult = parse(parser, argc, argv);
	if (parserResult != ArgumentParser::PARSE_OK)
		return parserResult;

	// extract options, if parsing was successful
	options.inputFiles.resize(getOptionValueCount(parser, "input"));
	for (vector< string >::size_type i = 0; i < options.inputFiles.size(); i++)
		getOptionValue(options.inputFiles[i], parser, "input", i);

	return parserResult;
}

// function to measure time between the first and second invocation of the function
void tic()
{
	static time_t start = 0;
	if (start != 0)
	{
		cout << "Elapsed: " << (time(NULL) - start) << endl;
		start = 0;
	}
	else
	{
		start = time(NULL);
	}
}

// Function which sums up the number of reads that start at a given position in the genome.
// Additionally, it counts the number of reads with an A/U at base 10 as well as
// the number of reads with a terminal base different from the reference genome.
// The latter two figures help identify the significance of a putative piRNA signature.
// Parameters:
//   bamFile: the BAM/SAM file from where to load the reads
//   readCountsUpstream: stats for positions were reads on the minus strand overlap the upstream (5') ends of reads on the plus strand
//   readCountsDownstream: stats for positions were reads on the minus strand overlap the downstream (3') ends of reads on the plus strand
int countReadsInBamFile(BamStream &bamFile, TCountsGenome &readCountsUpstream, TCountsGenome &readCountsDownstream)
{
	TCountsPosition *positionUpstream;
	TCountsPosition *positionDownstream;

	size_t alignmentLength;
	BamAlignmentRecord record;
	while (!atEnd(bamFile))
	{
		if (readRecord(record, bamFile) != 0)
		{
			cerr << "Failed to read record" << endl;
			return 1;
		}

		//todo: weighted counting of multimapped reads
		if (record.beginPos != BamAlignmentRecord::INVALID_POS && record.beginPos != -1)
		{
			// calculate alignment length using CIGAR string, i.e., distance from first matching base on the reference to the last matching base on the reference
			alignmentLength = 0;
			for (unsigned int cigarIndex = 0; cigarIndex < length(record.cigar); ++cigarIndex)
				if ((record.cigar[cigarIndex].operation == 'M') || (record.cigar[cigarIndex].operation == 'N') || (record.cigar[cigarIndex].operation == 'D') || (record.cigar[cigarIndex].operation == '=') || (record.cigar[cigarIndex].operation == 'X'))
					alignmentLength += record.cigar[cigarIndex].count;

			if (hasFlagRC(record)) // read maps to minus strand
			{
				positionUpstream = &(readCountsUpstream[STRAND_MINUS][record.rID][record.beginPos+alignmentLength]);
				positionDownstream = &(readCountsDownstream[STRAND_MINUS][record.rID][record.beginPos]);
				if ((record.cigar[length(record.cigar)-1].operation == 'S') && (record.cigar[length(record.cigar)-1].count == 1))
					positionDownstream->readsWithNonTemplateBase++;
				if ((record.cigar[0].operation == 'S') && (record.cigar[0].count == 1))
					positionUpstream->readsWithNonTemplateBase++;
			}
			else // read maps to plus strand
			{
				positionUpstream = &(readCountsUpstream[STRAND_PLUS][record.rID][record.beginPos]);
				positionDownstream = &(readCountsDownstream[STRAND_PLUS][record.rID][record.beginPos+alignmentLength]);
				if ((record.cigar[length(record.cigar)-1].operation == 'S') && (record.cigar[length(record.cigar)-1].count == 1))
					positionUpstream->readsWithNonTemplateBase++;
				if ((record.cigar[0].operation == 'S') && (record.cigar[0].count == 1))
					positionDownstream->readsWithNonTemplateBase++;
			}

			// increase read counters for the given position on the genome
			positionUpstream->reads++;
			positionDownstream->reads++;

			// check if 10th base is an Adenine or if the 10th base from the end of the read is a Uracil
			if (alignmentLength >= 10)
			{
				if ((record.seq[9] == 'A') || (record.seq[9] == 'a'))
				{
					positionUpstream->readsWithAAtBase10++;
					positionDownstream->readsWithAAtBase10++;
				}
				if ((record.seq[length(record.seq)-10] == 'T') || (record.seq[length(record.seq)-10] == 't'))
				{
					positionUpstream->readsWithUAtBase10FromEnd++;
					positionDownstream->readsWithUAtBase10FromEnd++;
				}
			}
		}
	}

	return 0;
}

// find those positions on the genome, where reads on opposite strands overlap by 10 nucleotides
void findOverlappingReads(TCountsGenome &readCounts, const unsigned int upstreamStrand, const TNameStore &bamNameStore)
{
	unsigned int downstreamStrand = (upstreamStrand == STRAND_PLUS) ? STRAND_MINUS : STRAND_PLUS;

	TCountsStrand::iterator contigOnOppositeStrand;
	TCountsContig::iterator positionOnOppositeStrand;
	for(TCountsStrand::iterator contig = readCounts[downstreamStrand].begin(); contig != readCounts[downstreamStrand].end(); ++contig)
	{
		contigOnOppositeStrand = readCounts[upstreamStrand].find(contig->first);
		if (contigOnOppositeStrand != readCounts[upstreamStrand].end())
		{
			for(TCountsContig::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
			{
				positionOnOppositeStrand = contigOnOppositeStrand->second.find(position->first + 10);
				if (positionOnOppositeStrand != contigOnOppositeStrand->second.end())
				{
//					if ((position->second.reads >= 100) && (positionOnOppositeStrand->second.reads >= 100))
						cout
							<< bamNameStore[contig->first] << "\t"
							<< (position->first+1) << "\t"
							<< downstreamStrand << "\t"
							<< position->second.reads << "\t"
							<< position->second.readsWithAAtBase10 << "\t"
							<< position->second.readsWithUAtBase10FromEnd << "\t"
							<< position->second.readsWithNonTemplateBase << "\t"
							<< upstreamStrand << "\t"
							<< positionOnOppositeStrand->second.reads << "\t"
							<< positionOnOppositeStrand->second.readsWithAAtBase10 << "\t"
							<< positionOnOppositeStrand->second.readsWithUAtBase10FromEnd << "\t"
							<< positionOnOppositeStrand->second.readsWithNonTemplateBase
							<< endl;
				}
			}
		}
	}
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

	TCountsGenome readCountsUpstream; // stats about positions where reads on the minus strand overlap with the upstream ends of reads on the plus strand
	TCountsGenome readCountsDownstream; // stats about positions where reads on the minus strand overlap with the downstream ends of reads on the plus strand

	TNameStore bamNameStore;

	// read all BAM/SAM files
	for(TInputFiles::iterator inputFile = options.inputFiles.begin(); inputFile != options.inputFiles.end(); ++inputFile)
	{
		// open SAM/BAM file
		BamStream bamFile(toCString(*inputFile));
		if (!isGood(bamFile))
		{
			cerr << "Failed to open input file: " << *inputFile << endl;
			return 1;
		}

		// for every position in the genome, count the number of reads that start at a given position
		if (countReadsInBamFile(bamFile, readCountsUpstream, readCountsDownstream) != 0)
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
	}

	// find positions where reads on the minus strand overlap with the upstream ends of reads on the plus strand
	findOverlappingReads(readCountsUpstream, STRAND_MINUS, bamNameStore);
	// find positions where reads on the minus strand overlap with the downstream ends of reads on the plus strand
	findOverlappingReads(readCountsDownstream, STRAND_PLUS, bamNameStore);

	return 0;
}


