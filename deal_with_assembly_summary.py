#!/usr/bin/env python3
# -*- coding: utf-8 -*-
debug = False
# %% import modules
import argparse
from argparse import ArgumentDefaultsHelpFormatter
from operator import index
import os
import logging
from typing import Counter, ItemsView, Tuple
from typing_extensions import Concatenate
import urllib.request
import ssl
import pandas as pd
import re
from matplotlib_venn import venn2_circles, venn2_unweighted
from matplotlib_venn import venn3_unweighted, venn3_circles
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ssl._create_default_https_context = ssl._create_unverified_context

# %% pass arguments
parser = argparse.ArgumentParser(prog = 'deal_with_assembly_summary', 
                                 description = 'A script for processing the assembly summary table from NCBI.', 
                                 formatter_class = ArgumentDefaultsHelpFormatter)
parser.add_argument('-o', '--output', default = './rsgb', type = str, action = "store", required = True,
                    help = 'A directory for saving downloaded assembly summary table. [default = %(default)s]')
parser.add_argument('-d', '--domain', action = "store", type = str, nargs = '*', 
                    default = ['archaea','bacteria','fungi','protozoa','viral','plant','invertebrate','vertebrate_mammalian','vertebrate_other'],
                    help = "Specify the domain(s) that you want to include, please use space to separate multiple domains. " \
                           "Note: all of the following domains will be included by default.")
parser.add_argument('-s', '--set', action="store", type=str, default = 'rsgb', choices={'rsgb', 'refseq', 'genbank'},
                    help="Specify refseq/genbank/rsgb, rsgb means the union of refseq and genbank.")
parser.add_argument('-r', '--rep', action="store", type=str,
                    nargs='*', default=['plant','invertebrate','vertebrate_mammalian','vertebrate_other'], 
                    help="Specify the domain(s) that you may want to use their representative genome.")
if debug:
    parser.print_help()
    # preset args for debugging
    opt = parser.parse_args(['--output', '/Users/xinwei/workspace/kraken2DBbuild/rsgb_20210706'])
    args = vars(opt)
    print(args)
else:
    opt = parser.parse_args()
    args = vars(opt)
if not os.path.isdir(args['output']):
    os.mkdir(args['output'])

# %% logging DB information
log_path = os.path.join(args['output'], "log.log")
logging.basicConfig(filename = log_path, level = logging.DEBUG, force = True, 
                    filemode = 'w', format = '%(asctime)s %(levelname)s > %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S')
logging.info('Running ' + __file__)
logging.info('*' * 10 + ' database version ' + '*' * 10)
logging.info('DB home: ' + args['output'])
logging.info('Domain: ' + ' '.join(args['domain']))
logging.info('Domain(rep): ' + ' '.join(args['rep']))
logging.info('Set type: ' + args['set'])

# %% download assembly summary
logging.info('*' * 10 + ' download assembly summary ' + '*' * 10)
pathDict = dict()
asDir = args['output'] + '/assembly_summary'
if not os.path.isdir(asDir):
    os.mkdir(asDir)
for i, d in enumerate(args['domain']):
    logging.info('Domian [' + str(i+1) + ']: ' + d)
    if args['set'] in ['rsgb', 'refseq']:
        rsURL = 'http://ftp.ncbi.nlm.nih.gov/genomes/refseq/' + d + '/assembly_summary.txt'
        rsFile = asDir + '/refseq_' + d + '_assembly_summary.txt'
        if not os.path.isfile(rsFile) or os.stat(rsFile).st_size == 0:
            logging.info('From NCBI refseq to ' + rsFile)
            urllib.request.urlretrieve(rsURL, rsFile)
            logging.info('Done.')
        else:
            logging.info(rsFile + ' is already exist.')
        pathDict[d + '_refseq'] = rsFile
    if args['set'] in ['rsgb', 'genbank']:
        gbURL = 'http://ftp.ncbi.nlm.nih.gov/genomes/genbank/' + d + '/assembly_summary.txt'
        gbFile = asDir + '/genbank_' + d + '_assembly_summary.txt'
        if not os.path.isfile(gbFile) or os.stat(gbFile).st_size == 0:
            logging.info('From NCBI genbank to ' + gbFile)
            urllib.request.urlretrieve(gbURL, gbFile)
            logging.info('Done.')
        else:
            logging.info(gbFile + ' is already exist.')
        pathDict[d + '_genbank'] = gbFile

# %% load assembly summary, meanwhile check each column
statDir = asDir + '/stat'
if not os.path.isdir(statDir):
    os.mkdir(statDir)
datDict = dict()
for k,v in pathDict.items():
    df = pd.read_csv(v, delimiter = '\t',  skiprows = 1)
    datDict[k] = df
    f = open(statDir + '/' + k + '_stat.txt', 'w')
    f.write('***** col 1: accession *****' + '\n')
    f.write('counts: ' + str(df['# assembly_accession'].count()) + '\n')
    f.write('unique: ' + str(df['# assembly_accession'].nunique()) + '\n')
    df['acc9'] = [re.findall(r"GC[AF]_(\d+)\.[0-9]+", ac)[0] for ac in df['# assembly_accession']]
    f.write('unique acc9: ' + str(df['acc9'].nunique()), file = f)
    f.write('\n***** col 5: refseq_category *****', file = f)
    f.write(df['refseq_category'].value_counts(), file = f)
    f.write('\n***** col 6: taxid *****', file = f)
    f.write('unique: ' + str(df['taxid'].nunique()), file = f)
    f.write('\n***** col 7: species_taxid *****', file = f)
    f.write('unique: ' + str(df['species_taxid'].nunique()), file = f)
    f.write('\n***** col 11: version_status *****', file = f)
    f.write(df['version_status'].value_counts(), file = f)
    f.write('\n***** col 12: assembly_level *****', file = f)
    f.write(df['assembly_level'].value_counts(), file = f)
    f.write('\n***** col 13: release_type *****', file = f)
    f.write(df['release_type'].value_counts(), file = f)
    f.write('\n***** col 14: genome_rep *****', file = f)
    f.write(df['genome_rep'].value_counts(), file = f)
    f.write('\n***** col 18: gbrs_paired_asm *****', file = f)
    f.write('paried:' + str(sum(df['gbrs_paired_asm'] != 'na')), file = f)
    f.write('no paried:' + str(sum(df['gbrs_paired_asm'] == 'na')), file = f)
    f.write('\n***** col 19: paired_asm_comparsion *****', file = f)
    f.write(df['paired_asm_comp'].value_counts(), file = f)
    f.close()

# %% union refseq and genbank
logging.info('*' * 10 + ' union refseq and genbank ' + '*' * 10)
uDict = dict()
mycol = ('#80B1D3', '#FDB462', '#B3DE68')
linewidth = 0.7
linestyle = 'dashed'
figHeight = 4 * len(args['domain'])
fig, axs=plt.subplots(len(args['domain']), 6, figsize=(24, figHeight), dpi=150)
for i, d in enumerate(args['domain']):
    logging.info('Domian [' + str(i+1) + ']: ' + d)
    rsDF = datDict[d + '_refseq']
    gbDF = datDict[d + '_genbank']
    acc2add = set(gbDF['acc9']) - set(rsDF['acc9'])
    uDF = pd.concat([rsDF, gbDF[gbDF['acc9'].isin(acc2add)]])
    uDict[d] = uDF
    if d in args['rep']:
        logging.info('Also output the table only contain representative genome.')
        uDict[d + '_rep'] = uDF[uDF['refseq_category'].isin(['representative genome', 'reference genome'])]
        uDict[d + '_rep'].to_csv(args['output'] + '/assembly_summary/' + args['set'] + '_' + d + '_rep_assembly_summary.txt', index = None, sep = '\t')
    uDF.to_csv(args['output'] + '/assembly_summary/' + args['set'] + '_' + d + '_assembly_summary.txt', index = None, sep = '\t')
    g=venn2_unweighted([set(rsDF['acc9']), set(gbDF['acc9'])], 
            set_labels = ('RefSeq', 'GenBank'), 
            set_colors = mycol,
            alpha = 0.6,
            normalize_to = 1.0,
            ax = axs[i, 0]
        )
    g=venn2_circles(subsets = (1, 1, 1), ax = axs[i, 0], 
                    linewidth = linewidth, linestyle = linestyle)
    g=venn2_unweighted([set(rsDF['taxid']), set(gbDF['taxid'])], 
            set_labels = ('RefSeq', 'GenBank'), 
            set_colors = mycol,
            alpha = 0.6,
            normalize_to = 1.0,
            ax = axs[i, 1]
        )
    g=venn2_circles(subsets = (1, 1, 1), ax = axs[i, 1], 
                    linewidth = linewidth, linestyle = linestyle)
    g=venn2_unweighted([set(rsDF['species_taxid']), set(gbDF['species_taxid'])], 
            set_labels = ('RefSeq', 'GenBank'), 
            set_colors = mycol,
            alpha = 0.6,
            normalize_to = 1.0,
            ax = axs[i, 2]
        )
    g=venn2_circles(subsets = (1, 1, 1), ax = axs[i, 2], 
                    linewidth = linewidth, linestyle = linestyle)
    title = [d, 'Accession', 'taxid', 'species_taxid']
    axs[i, 0].title.set_text(d.capitalize() + ' accession')
    axs[i, 1].title.set_text(d.capitalize() + ' taxid')
    axs[i, 2].title.set_text(d.capitalize() + ' species_taxid')

    repAcc9 = uDF['acc9'][uDF['refseq_category'].isin(['representative genome', 'reference genome'])]
    g=venn3_unweighted([set(rsDF['acc9']), set(gbDF['acc9']), set(repAcc9)], 
            set_labels = ('RefSeq', 'GenBank', 'Representative'), 
            set_colors = mycol,
            alpha = 0.6,
            normalize_to = 1.0,
            ax = axs[i, 3]
        )
    g=venn3_circles(subsets = (1, 1, 1, 1, 1, 1, 1), ax = axs[i, 3], 
                    linewidth = linewidth, linestyle = linestyle)
    
    repAcc9 = uDF['taxid'][uDF['refseq_category'].isin(['representative genome', 'reference genome'])]
    g=venn3_unweighted([set(rsDF['taxid']), set(gbDF['taxid']), set(repAcc9)], 
            set_labels = ('RefSeq', 'GenBank', 'Representative'), 
            set_colors = mycol,
            alpha = 0.6,
            normalize_to = 1.0,
            ax = axs[i, 4]
        )
    g=venn3_circles(subsets = (1, 1, 1, 1, 1, 1, 1), ax = axs[i, 4], 
                    linewidth = linewidth, linestyle = linestyle)
    repAcc9 = uDF['species_taxid'][uDF['refseq_category'].isin(['representative genome', 'reference genome'])]
    g=venn3_unweighted([set(rsDF['species_taxid']), set(gbDF['species_taxid']), set(repAcc9)], 
            set_labels = ('RefSeq', 'GenBank', 'Representative'), 
            set_colors = mycol,
            alpha = 0.6,
            normalize_to = 1.0,
            ax = axs[i, 5]
        )
    g=venn3_circles(subsets = (1, 1, 1, 1, 1, 1, 1), ax = axs[i, 5], 
                    linewidth = linewidth, linestyle = linestyle)
    title = [d, 'Accession', 'taxid', 'species_taxid']
    axs[i, 3].title.set_text(d.capitalize() + ' accession')
    axs[i, 4].title.set_text(d.capitalize() + ' taxid')
    axs[i, 5].title.set_text(d.capitalize() + ' species_taxid')
plotPath = args['output'] + '/assembly_summary/plot'
if not os.path.isdir(plotPath):
    os.mkdir(plotPath)
plt.savefig(plotPath + '/' + args['set'] + '_venn.pdf')
logging.info('Venn diagram for RefSeq and GenBank accession was saved to ' + plotPath)

