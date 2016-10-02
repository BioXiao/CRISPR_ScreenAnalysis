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

#############################################
### Core Functions
#############################################
def alignSeq2count(infilePathList, unalignedFilePathList, countFilePathList, processPool, RefLibrary, startIdx = None, endIdx = None, test = False):
    '''
    Functions for parallel process multiple sequence files

    :param infilePathList: Input sequence files to be aligned
    :param unalignedFilePathList: Output unaligned sequence read files
    :param countFilePathList: Output sequence count files
    :param processPool:
    :param RefLibrary: Library reference file
    :param startIdx: Start index of input sequence read to be aligned
    :param endIdx: End index of input sequence read to be aligned
    :param test: If just align the first N read for testing
    :return:
    '''

    if len(infilePathList) != len(unalignedFilePathList):
        raise ValueError('Input and ouput file number do not match!')

    # assemble parameters of each file for alignment function
    fileNum = len(infilePathList)
    align_argList = zip(infilePathList, unalignedFilePathList, countFilePathList, [RefLibrary]*fileNum, [startIdx]*fileNum, [stopIdx]*fileNum, [test]*fileNum)

    # multiprocessing
    readsPerFile = processPool.map(alignArgsWrapper, align_argList)

    return zip(countFilePathList, readsPerFile)


def alignArgsWrapper(argList):
    '''
    Function to pass the argument list as seprate argument to the real aligning function

    :param argList:
    :return:
    '''

    return seq2count(*argList)


def seq2count(infilePath, outUnalignedFilePath, countFilePath, RefLibrary, startIdx = None, stopIdx = None, test = False):
    '''

    :param infilePath:
    :param outUnalignedFilePath:
    :param countFilePath:
    :param RefLibrary:
    :param startIdx:
    :param stopIdx:
    :param test:
    :return:
    '''

    immed_print('Analyzing %s' % infilePath)

    # Recognize supported file type of input file
    fileType = None
    for fileExt_Type in acceptedFileType:
        if fnmatch.fnmatch(infilePath, fileExt_Type[0]):
            fileType = fileExt_Type[1]
            break

    if fileType == 'fqgz':
        linesPerRead = 4
        infile = gzip.open(infilePath)
    elif fileType == 'fq':
        linesPerRead = 4
        infile = open(infilePath)
    elif fileType == 'fa':
        linesPerRead = 2
        infile = open(infilePath)
    else:
        raise ValueError('Sequence file type not recognized!')

    SeqID_dict, IDCount_dict, expectedReadLength = parseLibraryFasta(RefLibrary)

    curRead = 0
    numAligning = 0

    with open(outUnalignedFilePath, 'w') as unaligned_f:
        # loop through each line of input fastq file
        for i, fastqLine in enumerate(infile):
            if i % 4 != 1:      # only select the sequence line
                continue

            else:
                seq = fastqLine.strip()[startIdx:stopIdx]

                if i == 1 and len(seq) != expectedReadLength:
                    raise ValueError('Trimmed read length does not match expected reference read length!')

                if seq in SeqID_dict:
                    # seq found in reference
                    for eachSeqID in SeqID_dict[seq]:
                        IDCount_dict[eachSeqID] += 1

                    numAligning += 1

                else:
                    # seq not found in reference
                    unaligned_f.write('>%d\n%s\n' % (i, seq))

                curRead += 1

                # if test is true, than only analyze the first N reads just for testing the pipeline
                if test and curRead >= TestLines:
                    break

    with open(countFilePath, 'w') as count_f:
        for eachID_count in sorted(zip(IDCount_dict.keys(), IDCount_dict.values())):
            count_f.write('%s\t%d\n' % eachID_count)

    immed_print('Done processing %s' % infilePath)

    return curRead, numAligning, numAligning * 100.0 / curRead





############################################
### Auxiliary Functions
############################################
def parseSeqFileName(fileNameList):
    '''
    Function to parse and check input file name.

    :param fileNameList: List of file names to be parsed (wildcard supported)
    :return infilePathList: Matched input file path found
            outfileBaseList: Prepare output file name list
    '''

    infilePathList = []
    outfileBaseList = []

    for eachFileName in fileNameList:       # loop through all input file name
        for eachFileMatched in glob.glob(eachFileName):        # for each file name, loop through each matched files
            for eachFileType in zip(*acceptedFileType)[0]:        # for each file, check if file type accepted
                if fnmatch.fnmatch(eachFileMatched, eachFileType):
                    infilePathList.append(eachFileMatched)
                    outfileBaseList.append(os.path.split(eachFileMatched)[-1].split('.')[0])        # select the file name from file path, and select file base name from file name
                    break


    return infilePathList, outfileBaseList


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
        raise ValueError('Library reference contains sequence with different length, make sure all reference sequences are of the same length.')

    return SeqID_Dict, IDCount_Dict, readLength[0]


##################################################
### Utility Functions
##################################################
def immed_print(content):
    '''
    Function to force print content on the screen immediately

    :param content: content to print
    :return:
    '''

    print content
    sys.stdout.flush()


def check_makeDir(dir_path):
    '''
    Function to check and make directory if not exists.

    :param dir_path: path to make directory
    :return:
    '''

    try:
        os.makedirs(dir_path)
    except OSError:
        immed_print(dir_path + ' already exists!')
        pass


### Set global variables
acceptedFileType = [('*.fastq.gz', 'fqgz'),
                    ('*.fastq', 'fq'),
                    ('*.fq', 'fq'),
                    ('*.fa', 'fa'),
                    ('*.fasta', 'fa'),
                    ('*.fna', 'fa')]        # List of tuples with accepted file extension and file type

TestLines = 10000

curDir = os.path.dirname(os.path.realpath(__file__))        # get the dir of the current script

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
    ### check inout
    print 'Start analyzing sequencing reads ......'
    print args

    numProcessor = max(args.processor, 1)

    # check input sequence files
    infilePathList, outfileBaseList = parseSeqFileName(args.SeqFileName)
    if len(infilePathList) == 0:
        sys.exit('Input error: sequence files not found.')
    else:
        print '%d sequence read files found to be analyzed.' % len(outfileBaseList)

    # check library reference file
    try:
        SeqID_dict, IDCount_dict, expectedReadLength = parseLibraryFasta(args.RefLibrary)
        immed_print('Library reference loaded successfully:\n\t%.2E Total sequence (%.2E unique sequence)\t%dbp reads expected' % (len(IDCount_dict), len(SeqID_dict), expectedReadLength))

    except IOError:
        sys.exit('Input error: library fasta file not found')

    except ValueError as err:
        sys.exit('Input error: ' + err.args[0])

    # check output path
    unalignedDir = os.path.join(curDir, 'analysisResult', args.AnalysisName, 'unaligned_reads')
    check_makeDir(unalignedDir)     # unaligned reads file directory
    countDir = os.path.join(curDir, 'analysisResult', args.AnalysisName, 'seq_count')
    check_makeDir(countDir)     # guide sequence count file directory

    unalignedFileNameList = [eachOutFile + '_unaligned.fa' for eachOutFile in outfileBaseList]
    unalignedFilePathList = [os.path.join(unalignedDir, eachUnalignedFile) for eachUnalignedFile in unalignedFileNameList]     # all unaligned read files path
    countFileNameList = [eachOutFile + '.count' for eachOutFile in outfileBaseList]
    countFilePathList = [os.path.join(countDir, eachCountFile) for eachCountFile in countFileNameList]      # all count files path

    # pool for multiprocessing
    pool = multiprocessing.Pool(min(len(infilePathList), numProcessor))

    # start to align sequence to reference
    try:
        alignResult = anlignSeq2count(infilePathList, unalignedFilePathList, countFilePathList, pool, args.RefLibrary, args.trim_start, args.trim_end, args.test)
    except ValueError as err:
        sys.exit('Error while processing sequencing files: ' + ''.join(err.args))

    # print alignment result
    print 'Alignment result:\n'
    for eachCountFile, eachReadsPerFile in alignResult:
        print eachCountFile + ':\n\t%.2E reads\t%.2E aligned (%.2f%%)' % eachReadsPerFile

    pool.close()
    pool.join()

    immed_print('Done processing all input sequence files!')

