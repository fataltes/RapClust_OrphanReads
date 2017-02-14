from __future__ import print_function
from __future__ import division
import time
import logging

logger = logging.getLogger("rapclust")

import itertools
import pandas as pd
import numpy as np
import os
sep = os.path.sep

from scipy.sparse import *
from sklearn.cluster import MiniBatchKMeans

tnames_inv = {}
tnames = []

def readFiles(sampdirs, auxDir):
    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampdirs]
    eqfiles = [sep.join([sd, auxDir, '/eq_classes.txt']) for sd in sampdirs]

    sumCounts = None
    firstSamp = True
    numSamp = 0
    totNumEq = 0
    totTpm = 0
    sparceMat = []
    for sffile, eqfile in itertools.izip(sffiles, eqfiles):
        quant = pd.read_table(sffile)
        quant.set_index('Name', inplace=True)

        with open(eqfile) as ifile:
            numSamp += 1
            numTran = int(ifile.readline().rstrip())
            numEq = int(ifile.readline().rstrip())
            logging.info("quant file: {}; eq file: {}; # tran = {}; # eq = {}".format(sffile, eqfile, numTran, numEq))
            if firstSamp:
                for i in xrange(numTran):
                    tnames.append(ifile.readline().rstrip())
                for i in range(len(tnames)):
                    tnames_inv[tnames[i]] = i
            else:
                for i in xrange(numTran):
                    ifile.readline()

            # for easy access to quantities of interest
            tpm = quant.loc[tnames, 'TPM'].values
            totTpm += tpm
            estCount = quant.loc[tnames, 'NumReads'].values
            efflens = quant.loc[tnames, 'EffectiveLength'].values
            epsilon =  np.finfo(float).eps
            sumCounts = np.maximum(sumCounts, estCount)

            for i in xrange(numEq):
                toks = map(int, ifile.readline().rstrip().split('\t'))
                nt = toks[0]
                tids = tuple(toks[1:-1])
                count = toks[-1]
                denom = sum([tpm[t] for t in tids])
                for t in tids:
                    sparceMat.append((t, totNumEq + i, count)) # Add count of reads for each transcript/EqClass pair
                    #sparceMat.append((t, totNumEq + i, count * (tpm[t]/denom)))
            firstSamp = False
        totNumEq += numEq
    for t in range(len(totTpm)):
        sparceMat.append((t, totNumEq, totTpm[t])) # Add sum of all tpm values for each transcript as a single column

    numColumns = totNumEq + 1
    numRows = numTran
    return numRows, numColumns, sparceMat

def doSparceCluster(rdim, cdim, sparceMat):
    # allocate a lil_matrix of size (rdim by cdim)
    # note: lil_matrix is used since we be modifying
    #       the matrix a lot.
    print(rdim, cdim)
    print(len(sparceMat))
    S = lil_matrix((rdim, cdim))

    # add data to S
    for (i, j, d) in sparceMat:
        S[i, j] = d
    # perform clustering
    labeler = MiniBatchKMeans(n_clusters=66000, verbose=1)
    # convert lil to csr format
    # note: Kmeans currently only works with CSR type sparse matrix
    labeler.fit(S.tocoo())
    return labeler.labels_

def convert2clusterFormat(labels, clustOutfile, flatClustOutfile):
    clusts = {}
    for (row, label) in enumerate(labels):
        if label not in clusts:
            clusts[label] = []
        clusts[label] += [tnames[row]]

    with open(flatClustOutfile, 'w') as fofile, open(clustOutfile, 'w') as ofile:
        for k, v in clusts.iteritems():
            cname = "cluster{}".format(k)
            new_cluster = ""
            for member in v:
                fofile.write("{}\t{}\n".format(cname, member))
                new_cluster += "{}\t".format(member)
            ofile.write("{}\n".format(new_cluster.strip()))
