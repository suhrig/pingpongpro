// ==========================================================================
//                                pingpongpal
// ==========================================================================

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
//  - readsWithAOrUAtBase10: how many of these reads have an Adenine or a Uracil at base 10 (which is typical for piRNAs)
//  - readsWithNonTemplateBase: how many of these reads have terminal base which is different from the reference genome (also typical for piRNAs)
struct TPosition
{
        int reads;
        int readsWithAOrUAtBase10;
        int readsWithNonTemplateBase;
        TPosition():
                reads(0),
                readsWithAOrUAtBase10(0),
                readsWithNonTemplateBase(0)
        {}
};

// array of TPosition to store the above stats for each position of a contig/chromosome
typedef map< unsigned int, TPosition > TPositionMap;

// array of TPositionMap to store the above status for each contig/chromosome
typedef map< unsigned int, TPositionMap > TContigMap;

// ==========================================================================
// Functions
// ==========================================================================


// function to parse command-line arguments
ArgumentParser::ParseResult parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
        // Setup ArgumentParser.
        ArgumentParser parser("pingpongpal");
        // Set short description, version, and date.
        setShortDescription(parser, "Put a Short Description Here");
        setVersion(parser, "0.1");
        setDate(parser, "July 2012");

        // Define usage line and long description.
        addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
        addDescription(parser, "This is the application skelleton and you should modify this string.");

        addOption(parser, ArgParseOption("i", "input", "Path to the input file(s)", ArgParseArgument::INPUTFILE, "IN", true));
        vector< string > acceptedInputFormats;
        acceptedInputFormats.push_back(".sam");
        acceptedInputFormats.push_back(".bam");
        setValidValues(parser, "input", acceptedInputFormats);

        // Add Examples Section.
        addTextSection(parser, "Examples");
        addListItem(parser, "\\fBpingpongpal\\fP \\fB-v\\fP \\fItext\\fP", "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

        // Parse command line.
        ArgumentParser::ParseResult parserResult = parse(parser, argc, argv);

        if (parserResult != ArgumentParser::PARSE_OK)
                return parserResult;

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

// program entry point
int main(int argc, char const ** argv)
{
        // Parse the command line.
        AppOptions options;
        if (parseCommandLine(options, argc, argv) != seqan::ArgumentParser::PARSE_OK)
                return 1;


        TContigMap readCountsByPosition[2];

        if (options.inputFiles.size() == 0)
                options.inputFiles.push_back("/dev/stdin");

        TPosition *pposition;
        size_t seqLength;

        for(TInputFiles::iterator inputFile = options.inputFiles.begin(); inputFile != options.inputFiles.end(); ++inputFile)
        {
                // Open input stream, BamStream can read SAM and BAM files.
                BamStream bamFile(toCString(*inputFile));

                BamAlignmentRecord record;
                while (!atEnd(bamFile))
                {
                        readRecord(record, bamFile);
                        //todo: handle clipped alignments (maybe with getClippedPos?)
                        //todo: count skipped reads
                        if (record.beginPos != BamAlignmentRecord::INVALID_POS)
                        {
                                seqLength = length(record.seq);
                                if (hasFlagRC(record)) {
                                        pposition = &(readCountsByPosition[STRAND_MINUS][record.rID][record.beginPos+seqLength]);
                                } else {
                                        pposition = &(readCountsByPosition[STRAND_PLUS][record.rID][record.beginPos]);
                                }
                                pposition->reads++;
                                if (seqLength >= 10)
                                        //todo: comparison with lowercase sequences
                                        if ((record.seq[9] == 'A') || (record.seq[9] == 'T') || (record.seq[9] == 'U'))
                                                pposition->readsWithAOrUAtBase10++;
                                pposition->readsWithNonTemplateBase++;
                        }
                }
        }

        /*for (int strand = 0; strand <= 1; strand++)
        {
                for(TContigMap::iterator contig = readCountsByPosition[strand].begin(); contig != readCountsByPosition[strand].end(); ++contig)
                {
                        for(TPositionMap::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
                        {
                                cout << strand << "." << contig->first << "." << position->first << ": " << position->second << endl;
                        }
                }
        }*/

        int density[1000];
        for (int i = 0; i < 1000; i++)
                density[i] = 0;
        int minCount = 0;
        TContigMap::iterator contigOnOppositeStrand;
        TPositionMap::iterator positionOnOppositeStrand;
        for(TContigMap::iterator contig = readCountsByPosition[STRAND_PLUS].begin(); contig != readCountsByPosition[STRAND_PLUS].end(); ++contig)
        {
                contigOnOppositeStrand = readCountsByPosition[STRAND_MINUS].find(contig->first);
                if (contigOnOppositeStrand != readCountsByPosition[STRAND_MINUS].end())
                {
                        for(TPositionMap::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
                        {
                                positionOnOppositeStrand = contigOnOppositeStrand->second.find(position->first + 10);
                                if (positionOnOppositeStrand != contigOnOppositeStrand->second.end())
                                {
//                                      minCount = (position->second.reads < positionOnOppositeStrand->second.reads) ? position->second.reads : positionOnOppositeStrand->second.reads;
//                                      if (minCount > 999)
//                                              minCount = 999;
//                                      density[minCount]++;

                                if ((position->second.reads >= 100) && (positionOnOppositeStrand->second.reads >= 100))
                                                cout << contig->first << "." << position->first << ": reads+ " << position->second.reads << ", a/u+ " << position->second.readsWithAOrUAtBase10 << ", non-template bases+ " << position->second.readsWithNonTemplateBase << ", reads- " << positionOnOppositeStrand->second.reads << ", a/u- " << positionOnOppositeStrand->second.readsWithAOrUAtBase10 << ", non-template bases- " << positionOnOppositeStrand->second.readsWithNonTemplateBase << endl;

                                }
                        }
                }
        }


/*      for (int strand = 0; strand <= 1; strand++)
                for (TContigMap::iterator contig = readCountsByPosition[strand].begin(); contig != readCountsByPosition[strand].end(); ++contig)
                        for (TPositionMap::iterator position = contig->second.begin(); position != contig->second.end(); ++position)
                                cout << strand << "." << contig->first << "." << position->first << ": reads " << position->second.reads << ", a/u " << position->second.readsWithAOrUAtBase10 << ", non-template bases " << position->second.readsWithNonTemplateBase << endl;*/


//      for (int i = 0; i < 1000; i++)
//              cout << i << "\t" << density[i] << endl;

        return 0;
}


