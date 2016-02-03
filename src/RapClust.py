import sys
import os
import logging

import click
import yaml

import eqnet

@click.command()
@click.option('--config', help='Config file describing the experimental setup')
def processQuant(config):
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

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

    logging.info("Building multiple alignment graph")
    eqnet.buildNetFile(sampleDirs, netFile, cutoff)
    # Filter the network
    filtNetFile = os.path.sep.join([outdir, "mag.filt.net"])

    logging.info("Filtering multiple alignment graph")
    eqnet.filterGraph(conditionDict, netFile, filtNetFile)
    # cluster the graph
    from subprocess import call
    clustFile = os.path.sep.join([outdir, "mag.clust"])

    logging.info("Clustering multiple alignment graph")
    call(["mcl", filtNetFile, "--abc", "-o", clustFile])

if __name__ == "__main__":
    processQuant()
