#############################################################
###           Sequence Align
#############################################################

# This python script is to align CRISPR sequening reads to sgRNA library to get the raw count of genes

import os
import sys
import gzip
import multiprocessing
import fnmatch
import glob
import argparse


### Core Functions



### Auxiliary Functions
def parseSeqFileName(fileNameList):
    '''
    Function to parse and check input file name.

    :param fileNameList: List of file names to be parsed (wildcard supported)
    :return outfileBaseList: Output file name list found
    '''

    outfileBaseList = []

    for eachFileName in fileNameList:       # loop through all input file name
        for eachFileMatched in glob.glob(eachFileName):        # for each file name, loop through each matched files
            for eachFileType in zip(*acceptedFileType)[0]:        # for each file, check if file type accepted
                if fnmatch.fnmatch(eachFileMatched, eachFileType):
                    outfileBaseList.append(os.path.split(eachFileMatched)[-1].split('.')[0])        # select the file name from file path, and select file base name from file name
                    break


    return outfileBaseList


def parseLibraryFasta(libraryFasta):
    '''
    Function to parse library sgRNA reference
    Reference format:
    >sequence id
    ATCGATCG

    :param libraryFasta: library fasta file to parsed
    :return SeqID_Dict: dict stores seq and seq id;
            IDCount_Dict: dict stores seq id and count (set 0);
            readLength: dict stores seq read length (should have same length)
    '''

    SeqID_Dict, IDCount_Dict, readLength = dict(), dict(), []

    curSeqID = ''
    curSeq = ''

    with open(libraryFasta) as infile:
        for line in infile:     # loop through each line
            if line[0] == '>':      # check if is ID or sequence
                # if line is ID
                # the following is processing the last seqID and seq read
                if curSeqID != '' and curSeq != '':
                    if curSeq not in SeqID_Dict:
                        SeqID_Dict[curSeq] = []     # if seq ID not found in dict, create a new key and value

                    SeqID_Dict[curSeq].append(curSeqID)     # add seqID value to dict
                    IDCount_Dict[curSeqID] = 0      # set the count to 0 for later alignment
                    readLength.append(len(curSeq))

                curSeqID = line.strip()
                curSeq = ''

            else:
                # if line is sequence
                curSeq += line.strip()[1:].upper()

    if len(SeqID_Dict) == 0 or len(IDCount_Dict) == 0 or readLength[0] == 0:
        raise ValueError("Couldn't parse library fasta file or no sequence detected.")
    elif max(readLength) != min(readLength):
        raise ValueError('Library reference contains sequence with different length, please make sure all reference sequences are of the same length.')

    return SeqID_Dict, IDCount_Dict, readLength[0]



### Utility Functions
def immed_print(content):
    '''
    Function to force print content on the screen immediately

    :param content:
    :return:
    '''

    print content
    sys.stdout.flush()


### Set global variables
acceptedFileType = [('*.fastq.gz', 'fqgz'),
                    ('*.fastq', 'fq'),
                    ('*.fq', 'fq'),
                    ('*.fa', 'fa'),
                    ('*.fasta', 'fa'),
                    ('*.fna', 'fa')]        # List of tuples with accepted file extension and file type

TestLines = 10000

curDir = os.path.dirname(os)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Align CRISPR sequening reads to sgRNA library to get raw gene counts')
    parser.add_argument('AnalysisName', help='Name of your experiment, use for output folder name.')
    parser.add_argument('SeqFileName', help='Name of sequence read file. Use wildcards to select multiple files. ex: "*.fasta.gz"')
    parser.add_argument('RefLibrary', help='Library fasta file of expected reads.')
    parser.add_argument('-p', '--processor', type=int, default=1, help='Number of processors to run multiple files parallelly.')
    parser.add_argument('--trim_start', type=int, help='Trim off the head sequence before the trim_start of every reads.')
    parser.add_argument('--trim_end', type=int, help='Trim off the tail sequence after the trim_end of every reads.')
    parser.add_argument('--test', action='store_true', default=False, help='Only analyze the head line of all input files for testing the pipeline.')

    args = parser.parse_args()

    ### check window #1
    print 'Start analyzing sequencing reads ......'


    numProcessor = max(args.processor, 1)

    # check input sequence files
    outfileBaseList = parseSeqFileName(args.SeqFileName)
    if len(outfileBaseList) == 0:
        sys.exit('Input error: sequence files not found.')

    # check library reference file
    try:
        SeqID_dict, IDCount_dict, expectedReadLength = parseLibraryFasta(args.RefLibrary)
        immed_print('Library reference loaded successfully:\n\t%.2E Total sequence (%.2E unique sequence)\t%dbp reads expected' % (len(IDCount_dict), len(SeqID_dict), expectedReadLength))

    except IOError:
        sys.exit('Input error: library fasta file not found')

    except ValueError as err:
        sys.exit('Input error: ' + err.args[0])

