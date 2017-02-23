from __future__ import print_function
from __future__ import division
import logging

logger = logging.getLogger("rapclust")

import itertools
import pandas as pd
import numpy as np
import os
sep = os.path.sep

from scipy.sparse import *
from sklearn.cluster import DBSCAN
import hdbscan

tnames_inv = {}
tnames = []

def readFiles(sampdirs, auxDir, cutoff=10.0, loadPickle=False):
    import pickle
    distFile = 'dict.pkl'
    if loadPickle:
        numRows, numColumns, sparceDic, weightDict = pickle.load(open(distFile, 'rb'))
        return numRows, numColumns, sparceDic, weightDict 
    else:
        sffiles = [sep.join([sd, 'quant.sf']) for sd in sampdirs]
        eqfiles = [sep.join([sd, auxDir, '/eq_classes.txt']) for sd in sampdirs]

        sumCounts = None
        firstSamp = True
        numSamp = 0
        totNumEq = 0
        totTpm = 0
        sparceDic = {}
        weightDict = {}
        for sffile, eqfile in itertools.izip(sffiles, eqfiles):
            quant = pd.read_table(sffile)
            quant.set_index('Name', inplace=True)

            with open(eqfile) as ifile:
                numSamp += 1
                numTran = int(ifile.readline().rstrip())
                numEq = int(ifile.readline().rstrip())
                logging.info("quant file: {}; eq file: {}; # tran = {}; # eq = {}".format(sffile, eqfile, numTran, numEq))
                if firstSamp:
                    for i in range(numTran):
                        tnames.append(ifile.readline().rstrip())
                    for i in range(len(tnames)):
                        tnames_inv[tnames[i]] = i
                    diagCounts = np.zeros(len(tnames))
                    sumCounts = np.zeros(len(tnames))
                    ambigCounts = np.zeros(len(tnames))
                else:
                    for i in range(numTran):
                        ifile.readline()

                # for easy access to quantities of interest
                tpm = quant.loc[tnames, 'TPM'].values
                totTpm += tpm
                estCount = quant.loc[tnames, 'NumReads'].values
                #efflens = quant.loc[tnames, 'EffectiveLength'].values
                epsilon =  np.finfo(float).eps
                sumCounts = np.maximum(sumCounts, estCount)

                for i in range(numEq):
                    toks = map(int, ifile.readline().rstrip().split('\t'))
                    nt = toks[0]
                    tids = tuple(toks[1:-1])
                    count = toks[-1]

                    denom = sum([tpm[t] for t in tids])
                    for t1, t2 in itertools.combinations(tids, 2):
                        w = count
                        key = (t1, t2)
                        if key in weightDict:
                            weightDict[key] += w
                        else:
                            weightDict[key] = w

                    for t in tids:
                        if t not in sparceDic:
                            sparceDic[t] = {}
                        sparceDic[t][(totNumEq + i)] = count
                        diagCounts[t] += count * (tpm[t] / denom)
                        ambigCounts[t] += count


                minWeight = 0.5
                maxWeight = 0.0
                prior = 0.1
                edgesToRemove = []

                ##
                #  Go through the weightMap and remove any edges that
                #  have endpoints with too few mapping reads
                ##
                for k, v in weightDict.iteritems():
                    c0, c1 = diagCounts[k[0]], diagCounts[k[1]]
                    a0, a1 = ambigCounts[k[0]], ambigCounts[k[1]]
                    if a0 + a1 > epsilon and a0 > cutoff and a1 > cutoff:  # and v > 3:
                        #w = (v + prior) / min((a0 + prior), (a1 + prior))
                        w = (v + prior) / (a0 + a1 - v + prior)
                        weightDict[k] = w
                        if w > maxWeight:
                            maxWeight = w
                            # if w < minWeight:
                            #    edgesToRemove.append(k)
                    else:
                        edgesToRemove.append(k)

                # Actually delete those edges
                for e in edgesToRemove:
                    del weightDict[e]

                tnamesFilt = []
                relabel = {}
                for i in range(len(estCount)):
                    if (diagCounts[i] > cutoff):
                        relabel[i] = len(tnamesFilt)
                        tnamesFilt.append(tnames[i])
                    weightDict[(i, i)] = 1.1

                firstSamp = False
                totNumEq += numEq
        #for t in range(len(totTpm)):
        #    sparceMat.append((t, totNumEq, totTpm[t])) # Add sum of all tpm values for each transcript as a single column

        numColumns = totNumEq + 1
        numRows = numTran
        #pickle.dump((numRows, numColumns, sparceDic, weightDict), open(distFile, 'wb'))        
        return numRows, numColumns, sparceDic, weightDict

def calcDist(key, sparceDic):
    dist = 0
    t1Vec = sparceDic[key[0]]
    t2Vec = sparceDic[key[1]]
    for col in t1Vec:
        if col not in t2Vec:
            dist += t1Vec[col] #Manhattan
            #dist += (t1Vec[col])**2 #Euclidean
    for col in t2Vec:
        if col not in t1Vec:
            dist += t2Vec[col] #Manhattan
            #dist += (t2Vec[col])**2 #Euclidean
    return dist

def buildDistMat(sampdirs, auxDir, numRows, sparseDic, weightDic, distType):
    eqfiles = [sep.join([sd, auxDir, '/eq_classes.txt']) for sd in sampdirs]

    distMat = {}
    numSamp = 0
    numRows = numRows
    if distType == 'Jaccard':
        for key in weightDic:
            distMat[key] = 1.1 - weightDic[key] + 1e-4
    else:
        for eqfile in eqfiles:
            with open(eqfile) as ifile:
                numSamp += 1
                numTran = int(ifile.readline().rstrip())
                numRows = max(numRows, numTran)
                numEq = int(ifile.readline().rstrip())
                logging.info("eq file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
                for i in range(numTran):
                    ifile.readline()
                    distMat[(i, i)] = 1e-4
                for i in range(numEq):
                    toks = map(int, ifile.readline().rstrip().split('\t'))
                    tids = tuple(toks[1:-1])
                    for t1, t2 in itertools.combinations(tids,2):
                        if t1 > t2:
                            t1, t2 = t2, t1
                        key = (t1, t2)
                        if key not in distMat:
                            distMat[key] = calcDist(key, sparseDic)
    return numRows, distMat

def doSparceCluster(rdim, cdim, distMat):
    # allocate a lil_matrix of size (rdim by cdim)
    # note: lil_matrix is used since we be modifying
    #       the matrix a lot.
    print(rdim, cdim)
    print(len(distMat))
    S = lil_matrix((rdim, cdim))

    # add data to S
    sumDist = 0
    cntDist = 0
    minDist = 10000000000
    maxDist = 0
    all_contigs = set(np.arange(0, rdim, 1))
    existing_contigs = set()
    for (i, j) in distMat:
        existing_contigs.add(i)
        existing_contigs.add(j)
        S[i, j] = distMat[(i, j)]
        S[j, i] = distMat[(i, j)]
        sumDist += distMat[(i, j)]
        cntDist += 1
        if distMat[(i, j)]>maxDist:
            maxDist = distMat[(i, j)]
        if distMat[(i, j)]<minDist:
            minDist = distMat[(i, j)]
    print("AVERAGE OF DISTANCE ::: {}, MIN: {}, MAX: {}".format(sumDist/cntDist, minDist, maxDist))
    print("# of unique contigs that aren't mentioned in distance matrix: {}".format(len(all_contigs-existing_contigs)))
    # convert lil to csr format
    # note: Kmeans and DBSCAN currently only work with CSR type sparse matrix
    input_data = S.tocsr()
    sel_eps = 1.09
    print("cases less than eps={}: {}".format(sel_eps, np.sum(input_data.data<=sel_eps)))
    # perform clustering
    labeler = DBSCAN(eps=sel_eps, min_samples=3, metric='precomputed')
    #labeler = hdbscan.HDBSCAN(min_cluster_size=8, metric='precomputed', min_samples=1)
    labeler.fit(input_data)
    return labeler.labels_

def convert2clusterFormat(labels, clustOutfile, flatClustOutfile):
    clusts = {}
    curClustNum = -1
    for (row, label) in enumerate(labels):
        if label == -1:
            label = curClustNum
            curClustNum-=1
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
