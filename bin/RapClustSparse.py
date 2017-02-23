#!/usr/bin/env python
import sys
import os
import logging
import coloredlogs

import click
import argparse
import yaml
import json

#from rapclust import sparseClust
from rapclust import sparseClust_distMat as sparseClust

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--config', required=True, help='Config file describing the experimental setup')
def processQuant(config):
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("rapclust")
    coloredlogs.install(level='INFO')

    #config = args.config
    cfg = None
    with open(config, 'r') as yamCfg:
        cfg = yaml.load(yamCfg)

    # A list of all sample directories
    sampleDirs =  []
    outdir = None
    cutoff = 10.0

    if 'conditions' not in cfg:
        logging.critical("Your configuration file must contain a \"conditions\" entry!")
        sys.exit(1)
    else:
        outdir = cfg['outdir']
        conditions = cfg['conditions']
        conditionDict = {c : {} for c in conditions}
        for c in conditions:
            samples = cfg['samples'][c]
            for i, s in enumerate(samples):
                conditionDict[c][i] = s
                sampleDirs.append(s)

        if 'cutoff' in cfg:
            cutoff = cfg['cutoff']

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            logging.critical("The output directory already exists, and is not a directory!")
            sys.exit(1)
    else:
        # create it
        os.makedirs(outdir)
    # now the outdir exists
    netFile = os.path.sep.join([outdir, "mag.net"])
    # Make the network
    orphanFile = os.path.sep.join([outdir, "orphanLink"])

    with open(os.path.sep.join([sampleDirs[0], 'cmd_info.json'])) as json_file:
        cmd_info = json.load(json_file)
    
    if 'auxDir' in cmd_info:
        auxDir = cmd_info['auxDir']
    else:
        auxDir = 'aux'

    orphanFile = os.path.sep.join([outdir, "orphanLink"])

    clustFile = os.path.sep.join([outdir, "mag.clust"])
    flatClustFile = os.path.sep.join([outdir, "mag.flat.clust"])

    (numRows, numColumns, sparseDic, weightDic) = sparseClust.readFiles(sampleDirs, auxDir)
    #labels = sparseClust.doSparceCluster(numRows, numColumns, sparseMat)
    (numTrans, sparseDist) = sparseClust.buildDistMat(sampleDirs, auxDir, numRows, sparseDic, weightDic, distType='Jaccard')
    labels = sparseClust.doSparceCluster(numTrans, numTrans, sparseDist)
    sparseClust.convert2clusterFormat(labels, clustFile, flatClustFile)


    import numpy as np
    summaryFile = os.path.sep.join([outdir, "stats.json"])
    sizes = []
    with open(clustFile) as ifile:
        for l in ifile:
            toks = l.rstrip().split()
            sizes.append(len(toks))
    sizes = np.array(sizes)
    stats = {}
    stats['min clust size'] = sizes.min()
    stats['max clust size'] = sizes.max()
    stats['mean clust size'] = sizes.mean()
    stats['num clusts'] = len(sizes)

    with open(summaryFile, 'w') as ofile:
        json.dump(stats, ofile)

if __name__ == "__main__":
    #parser = argparse.ArgumentParser(description='Cluster de novo assemblies.')
    #parser.add_argument('config', type=str, help='the configuration file')
    #args = parser.parse_args()
    processQuant()
