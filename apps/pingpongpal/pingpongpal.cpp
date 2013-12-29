// ==========================================================================
//                                pingpongpal
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
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
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

typedef vector<CharString> TInputFiles;

struct AppOptions
{
        TInputFiles inputFiles;

        AppOptions()
        {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
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

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.


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
typedef map< unsigned int, TPosition > TPositionMap;
typedef map< unsigned int, TPositionMap > TContigMap;

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
                                        pposition = &(readCountsByPosition[0][record.rID][record.beginPos+seqLength]);
                                } else {
                                        pposition = &(readCountsByPosition[1][record.rID][record.beginPos]);
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
        for(TContigMap::iterator contig = readCountsByPosition[1].begin(); contig != readCountsByPosition[1].end(); ++contig)
        {
                contigOnOppositeStrand = readCountsByPosition[0].find(contig->first);
                if (contigOnOppositeStrand != readCountsByPosition[0].end())
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

