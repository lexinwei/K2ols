#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#################################################################
# Xin Wei, xwei97@zju.edu.cn
# Updated: 07/06/2021
debug = False
# %% import modules
import argparse, os, re, sys, subprocess, logging
from argparse import ArgumentDefaultsHelpFormatter
import numpy as np
# %% pass arguments
parser = argparse.ArgumentParser(prog = 'kraken2M', 
                                 description = 'Species classification for multiple (paired/single-end) fastq files one time with Kraken2, in order not to load the super index repeatly. Most of the following arguments are originated from Kraken2, not supporting the modification of other Kraken2 arguments, just using their default value setted by Kraken2. Please clone KrakenTools by jenniferlu717 from https://github.com/lexinwei/KrakenTools.git before running.', 
                                 formatter_class = ArgumentDefaultsHelpFormatter)

Rflag = parser.add_argument_group('REQUIRED PARAMETERS')
Rflag.add_argument('-i', '--input', type = str, required = True,
                    help = 'A directory of multiple fastq files, only the files with names ending in specific suffix (given by -s arguments) will be loaded.')
Rflag.add_argument('-s', '--suffix', action="store", type=str, required = True, default = None, 
                    help="Specify the suffix of reads files, e.g. 'R1.fastq,R2.fastq' for paired-end files, '.fq' for single-end files.")
Rflag.add_argument('-d', '--db', action="store", type=str, required = True, 
                    help= "Name for Kraken2 DB.")
Rflag.add_argument('-k', '--kraken', action="store", type=str, required = True, default = None, 
                    help= "Path of Kraken2 program.")
Rflag.add_argument('-kt', '--kraken-tools', action="store", type=str, required = True, default = None,
                    help= "Path of KrakenTools by jenniferlu717.") 
Oflag = parser.add_argument_group('OPTIONAL PARAMETERS')
Oflag.add_argument('-o', '--output', default = './kraken2_output', type = str, action = "store", required = False,
                    help = 'A directory for saving output files.')
Oflag.add_argument('-c', '--confidence', action="store", type=str, default = 0, required = False,
                    help = "Confidence score threshold, must be in [0, 1].")
Oflag.add_argument('-t', '--threads', action="store", type=str, default = 1, required = False,
                    help = "Number of threads to running Kraken2.")
Oflag.add_argument('--gzip-compressed', action="store_true", required = False,
                    help = "Input files are compressed with gzip.")

if not debug:
    opt = parser.parse_args()
    args = vars(opt)
if debug:
    parser.print_help()
    opt = parser.parse_args(['--input', '/share/home/jianglab/weixin/workspace/classified_by_kraken2/4-500-soil/cleandata', 
                            '--suffix', 'R1.fq,R2.fq',
                            '--db', '/share/home/jianglab/weixin/data/kraken2/cdb/viral',
                            '--kraken', '/share/home/jianglab/weixin/bin/kraken2/kraken2',
                            '--kraken-tools', '/share/home/jianglab/weixin/workspace/mytools/ktools/KrakenTools',
                            '--output', '/share/home/jianglab/weixin/workspace/classified_by_kraken2/4-500-soil/kraken2_output',
                            '--threads', '8'])
    args = vars(opt)
    print(args)

args['input'] = os.path.abspath(args['input'])
args['output'] = os.path.abspath(args['output'])
suffix = args['suffix'].split(',')

if not os.path.isdir(args['output']):
    os.mkdir(args['output'])
tmpDir = args['output'] + '/tmp'
if not os.path.isdir(tmpDir):
    os.mkdir(tmpDir)
if len(suffix) not in [1, 2]:
    'Bad value for --suffix argument.'
    sys.exit()
if not os.path.isdir(args['db']):
    logging.error('Unable to find the Kraken2 DB in ' + args['db'] + '.')
if args['kraken_tools'] == None or (not os.path.isdir(args['kraken_tools'])):
    logging.error('Unable to find KrakenTools.')
# %% logging DB information
log_path = os.path.join(args['output'], "log.log")
logging.basicConfig(filename = log_path, level = logging.DEBUG, force = True, 
                    filemode = 'w', format = '%(asctime)s %(levelname)s > %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S')
logging.info('Running ' + __file__)
logging.info(args)

# %% get reads header list of each file
allFileNames = []
for af in os.listdir(args['input']):
    if suffix[0] in af:
        allFileNames.append(af)
    if len(suffix) == 2:
        if suffix[1] in af:
            allFileNames.append(af)
fileNameList = [f.replace(suffix[0], '') for f in allFileNames]
mode = 'single-end'
if len(suffix) == 2:
    fileNameList = [f.replace(suffix[1], '') for f in fileNameList]
    fileNameList = list(set(fileNameList))
    fileNameList.sort()
    mode = 'paired-end'
logging.info(str(len(fileNameList)) + ' ' + mode + ' samples')

# %% concatenate reads
logging.info('*' * 15 + ' concatenate reads ' + '*' * 15)
if not os.path.isfile(tmpDir + '/' + suffix[0]):
    logging.info('start')
    for i,f in enumerate(fileNameList):
        if i == 0: 
            append_mode = '>'
        else:
            append_mode = '>>'
        command = ' '.join(['cat', args['input'] + '/' + f + suffix[0], append_mode, tmpDir + '/' + suffix[0]])
        logging.info(str(i+1) + '/' + str(len(fileNameList)) + ': ' + f + suffix[0] + ' -> ' + suffix[0])
        os.system(command)
        if len(suffix) == 2:
            command = ' '.join(['cat', args['input'] + '/' + f + suffix[1], append_mode, tmpDir + '/' + suffix[1]])
            logging.info(str(i+1) + '/' + str(len(fileNameList)) + ': ' + f + suffix[1] + ' -> ' + suffix[1])
            os.system(command)
    logging.info('end')
else:
    logging.info('The concatenated file already exit.')

# %% running kraken2
logging.info('*' * 15 + ' running kraken2 ' + '*' * 15)
if os.path.isfile(tmpDir + '/classified_seqs.fastq') or os.path.isfile(tmpDir + '/unclassified_seqs.fastq') or os.path.isfile(tmpDir + '/classified_seqs_1.fastq') or os.path.isfile(tmpDir + '/unclassified_seqs_1.fastq'):
    logging.info('Kraken2 classification looks like already done, skip and continue next part.')
else:
    if mode == 'single-end':
        command = [args['kraken2'], '--threads', args['threads'], '--db', args['db'], 
                   '--confidence', str(args['confidence']),
                   '--classified-out', tmpDir + '/classified_seqs.fastq', 
                   '--unclassified-out', tmpDir + '/unclassified_seqs.fastq',
                   '--report', tmpDir + '/' + 'report.txt',
                   '--output', tmpDir + '/' + 'output.txt',
                   '--use-names']
        if args['gzip_compressed']:
            command.append('--gzip-compressed')
        command = command.append(tmpDir + '/' + suffix[0])
    else:
        command = [args['kraken'], '--threads', args['threads'], '--db', args['db'], 
                   '--paired', tmpDir + '/' + suffix[0], tmpDir + '/' + suffix[1],
                   '--confidence', str(args['confidence']),
                   '--classified-out', tmpDir + '/classified_seqs#.fastq', 
                   '--unclassified-out', tmpDir + '/unclassified_seqs#.fastq',
                   '--report', tmpDir + '/' + 'report.txt',
                   '--output', tmpDir + '/' + 'output.txt',
                   '--use-names']
        if args['gzip_compressed']:
            command.append('--gzip-compressed')
    command = ' '.join(command)
    logging.debug('command: ' + command)
    logging.info('start')
    kstat = os.system(command)
    if kstat == 0:logging.info('end')

# %% count reads for each sample
logging.info('*' * 15 + ' count reads ' + '*' * 15)
readCounts = []
for i,f in enumerate(fileNameList):
    catMode = 'cat'
    if args['gzip_compressed']:
        catMode = 'zcat'
    command = [catMode, args['input'] + '/' + fileNameList[i] + suffix[0], '| echo $((`wc -l`/4))']
    rn = os.popen(' '.join(command)).read().strip()
    logging.info(str(i+1) + '/' + str(len(fileNameList)) + ': ' + f + suffix[0] + ' : ' + rn)
    readCounts.append(int(rn))

# %% split output.txt sample by sample
if os.path.isfile(args['output'] + '/' + fileNameList[-1] + 'out.txt'):
    logging.info('Looks like it had been splited, skip and continue to next step.')
else:
    readCountsAcc = np.cumsum(readCounts)
    logging.info('*' * 15 + ' split output.txt ' + '*' * 15)
    resF = open(tmpDir + '/output.txt', 'r')
    ind = 0
    fh = fileNameList[ind]
    logging.info(str(ind+1) + '/' + str(len(fileNameList)) + ': output.txt -> ' + fh + 'out.txt')
    spOUT = open(args['output'] + '/' + fh + 'out.txt', 'w')   
    for i,rows in enumerate(resF):
        spOUT.write(rows)
        if i+1 in readCountsAcc:
            ind += 1
            if ind == len(fileNameList): break
            fh = fileNameList[ind]
            logging.info(str(ind+1) + '/' + str(len(fileNameList)) + ': output.txt -> ' + fh + 'out.txt')
            spOUT.close()
            spOUT = open(args['output'] + '/' + fh + 'out.txt', 'w')
    spOUT.close()

# %% convert results to report
# need KrakenTools by jenniferlu717 https://github.com/lexinwei/KrakenTools.git
# make ktaxonomy
logging.info('*' * 15 + ' convert results to report ' + '*' * 15)
if not os.path.isfile(args['db'] + '/mydb_taxonomy.txt'):
    logging.info('Making ktaxonomy ...')
    command = ['python', args['kraken_tools'] + '/make_ktaxonomy.py', 
               '--node', args['db'] + '/taxonomy' + '/nodes.dmp',
               '--names', args['db'] + '/taxonomy' + '/names.dmp',
               '--seqid2taxid', args['db'] + '/seqid2taxid.map',
               '-o', args['db'] + '/mydb_taxonomy.txt']
    command = ' '.join(command)
    logging.debug('command: ', command)
    logging.info('start')
    process = subprocess.Popen(command, shell = True,
                           stdout = subprocess.PIPE, 
                           stderr = subprocess.PIPE)
    out, err = process.communicate()
    outF = open(args['db'] + '/make_ktaxonomy_out.txt', 'wb')
    outF.write(out)
    outF.close()
    errF = open(args['db'] + '/make_ktaxonomy_err.txt', 'wb')
    errF.write(err)
    errF.close()
    returncode = process.returncode
    returncode = process.returncode
    if returncode == 0: logging.info('end')
else:
    logging.info('No need to make ktaxonomy again, already exist in this DB.')

# make kreport
logging.info('start converting')
for i,s in enumerate(fileNameList):
    command = ['python', args['kraken_tools'] + '/make_kreport.py', 
            '-i', args['output'] + '/' + s + 'out.txt', 
            '-t', args['db'] + '/mydb_taxonomy.txt',
            '-o', args['output'] + '/' + s + 'kreport.txt']
    command = ' '.join(command)
    logging.info(str(i+1) + '/' + str(len(fileNameList)) + ': ' + s + 'out.txt -> ' + s + 'kreport.txt')
    stat = os.system(command)
    if not stat == 0:
        logging.warning('fail converting')
logging.info('all done')