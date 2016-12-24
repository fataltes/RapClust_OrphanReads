from __future__ import division
import itertools
import time
from matplotlib import pyplot as plt
trinity_human = "./tx2GeneHuman.38.80.txt"
class Classification:
    TruePos, FalsePos, TrueNeg, FalseNeg = range(4)

def classType(true1, true2, pred1, pred2):
    if true1 == true2:
        if pred1 == pred2:
            return Classification.TruePos
        else: # truely the same, predicted different
            return Classification.FalseNeg
    else: # truly different
        if pred1 == pred2: #predicted same
            return Classification.FalsePos
        else:
            return Classification.TrueNeg


def accuracyExpressed(groundTruth_clust, tr_clust):
    #count true postive for each pair of transcripts O(N^2)
    tp, fp, tn, fn = 0, 0, 0, 0
    for tr_1, tr_2 in itertools.combinations(tr_clust.keys(), 2):
        if tr_1 not in groundTruth_clust or tr_2 not in groundTruth_clust:
            continue
        ct = classType(groundTruth_clust[tr_1], groundTruth_clust[tr_2], tr_clust[tr_1], tr_clust[tr_2])
        if ct == Classification.TruePos:
            tp += 1
        elif ct == Classification.TrueNeg:
            tn += 1
        elif ct == Classification.FalsePos:
            fp += 1
        elif ct == Classification.FalseNeg:
            fn += 1
    return tp, fp, tn, fn


def accuracyExpressedFast(groundTruth_clust, groundTruth_clust_inv, tr_clust, tr_clust_inv):
    #num = len(set(tr_clust.keys()) & set(groundTruth_clust.keys()))
    num = len(set(groundTruth_clust.keys()))
    tp, fp, tn, fn, orphan_cnt = 0, 0, 0, 0, 0
    for clustName, clustMems in tr_clust_inv.iteritems():
        for tr_1, tr_2 in itertools.combinations(clustMems,2):
            if tr_1 not in groundTruth_clust or tr_2 not in groundTruth_clust:
                continue
            if groundTruth_clust[tr_1] == groundTruth_clust[tr_2]:
                tp += 1
            else:
                fp += 1
    for clustName, clustMems in groundTruth_clust_inv.iteritems():
        for tr_1, tr_2 in itertools.combinations(clustMems,2):
            if tr_1 not in tr_clust or tr_2 not in tr_clust:
                continue
            if tr_clust[tr_1] != tr_clust[tr_2]:
                fn += 1

    nc2 = (num * (num-1)) / 2
    tn = nc2 - (fp + tp + fn)
    return tp, fp, tn, fn

def weightCorrelation (weightGraph, groundTruth_clust, groundTruth_clust_inv):
    class0 = []
    class1 = []
    corrupted = set()
    for key, weight in weightGraph.iteritems():
        if key[0] in groundTruth_clust and key[1] in groundTruth_clust:
            if groundTruth_clust[key[0]] == groundTruth_clust[key[1]]:
                class1 += [weight]
            else:
                class0 += [weight]
        else:
            if key[0] not in groundTruth_clust:
                corrupted.add(key[0])
            else:
                corrupted.add(key[1])
    print ("corrupted: {}".format(len(corrupted)))
    print(corrupted)
    from matplotlib import pyplot as plt
    plt.figure()
    plt.scatter(class0+class1, [0]*len(class0) + [1]*len(class1))
    plt.ylabel("same/different cluster (1:same cluster)")
    plt.xlabel("link weights")
    plt.show()

def readNetFile(fn):
    weightGraph = {}
    with open(fn) as fp:
        for l in fp:
            left, right, weight = l.split()
            weightGraph[(left, right)] = weight
    return weightGraph

def readCDHitClust(fn, filtDict=None):
    tr_clust = {}
    tr_clust_inv = {}
    fp = open(fn)
    cnum = None
    for l in fp:
        if l.startswith('>Cluster'):
            cnum = int(l.rstrip().split()[-1])
        else:
            e = l.split()[2].lstrip('>').rstrip('.')
            if not filtDict or e in filtDict:
                tr_clust[e] = cnum
    for k,v in tr_clust.iteritems():
        if v in tr_clust_inv:
            tr_clust_inv[v].append(k)
        else:
            tr_clust_inv[v] = [k]
    return tr_clust, tr_clust_inv

def readCorset(fn):
    ft = open(fn)
    clust_dict = {}
    clust_dict_inv = {}
    for line in ft:
        toks = line.rstrip().split()
        clust_dict[toks[0]] = toks[1]
    for k,v in clust_dict.iteritems():
        if v in clust_dict_inv:
            clust_dict_inv[v].append(k)
        else:
            clust_dict_inv[v] = [k]
    return clust_dict, clust_dict_inv

def readTrueLabels(fn):
    ft = open(fn)
    groundTruth_clust = {}
    groundTruth_clust_inv = {}
    gtClusterCount = {}
    for line in ft:
        tr_gn = line[:-1].split("\t")
        groundTruth_clust[tr_gn[0]] = tr_gn[1]
        if tr_gn[0] in gtClusterCount:
            gtClusterCount[tr_gn[0]] += 1
        else:
            gtClusterCount[tr_gn[0]] = 1
    for k,v in groundTruth_clust.iteritems():
        if v in groundTruth_clust_inv:
            groundTruth_clust_inv[v].append(k)
        else:
            groundTruth_clust_inv[v] = [k]
    return groundTruth_clust, groundTruth_clust_inv

def readMCLClust(fn):
    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    print (dir_path)
    fp = open(fn)
    tr_clust = {}
    tr_clust_inv = {}
    key = 1
    for cnum, line in enumerate(fp):
        same_cluster = line.rstrip().split('\t')
        for contig in same_cluster:
            tr_clust[contig] = cnum
    for k,v in tr_clust.iteritems():
        if v in tr_clust_inv:
            tr_clust_inv[v].append(k)
        else:
            tr_clust_inv[v] = [k]
    return tr_clust, tr_clust_inv

def readCommClust(fn):
    netdic = pickle.load(open(fn,'rb'))
    tr_clust = {}
    tr_clust_inv = {}
    for contig,cnum in netdic.iteritems():
        tr_clust[str(contig)] = cnum
    for k,v in tr_clust.iteritems():
        if v in tr_clust_inv:
            tr_clust_inv[v].append(k)
        else:
            tr_clust_inv[v] = [k]
    return tr_clust, tr_clust_inv

def measurePrecRecall(qsf, sp, ctype):
    if ctype == "mcl":
        tr_clust, tr_clust_inv = readMCLClust(qsf)
    elif ctype == "corset":
        tr_clust, tr_clust_inv = readCorset(qsf)

    if(sp == 'human'):
        ft = trinity_human
    elif (sp == 'yeast'):
        ft = trinity_yeast
    elif (sp == 'chicken'):
        ft = trinity_chicken
    elif (sp == 'mouse'):
        ft = trinity_mouse
    elif (sp == 'rice'):
        ft = trinity_rice
    groundTruth_clust, ground_truth_clust_inv = readTrueLabels(ft)
    weightGraph = readNetFile("transcripts/human_rapclust/mag.filt.net")
    tp, fp, tn, fn = accuracyExpressedFast(groundTruth_clust, ground_truth_clust_inv,
                                                       tr_clust, tr_clust_inv)
    weightCorrelation(weightGraph, groundTruth_clust, ground_truth_clust_inv)
    print("tp : {}, fp : {}, tn : {}, fn : {}".format(tp, fp, tn, fn))
    precision = tp / float(tp + fp)
    recall = tp / float(tp + fn)
    Fscore = 2 * ((precision*recall) / (precision+recall))
    print("prec: {}, recall: {}, F1: {}".format(precision, recall, Fscore))

import argparse

def main():
    parser = argparse.ArgumentParser(description="Give the SRR number")
    parser.add_argument('--clustfile',type = str, help="graph file")
    parser.add_argument('--sp',type = str, help="type human")
    parser.add_argument('--ctype',type = str, default="mcl", help="you do not need to enter anything here")
    parser.add_argument('--orphanfile', type = str, help="name of file containing orphan links")
    args = parser.parse_args()
    measurePrecRecall(args.clustfile, args.sp, args.ctype)


if __name__ == "__main__":
    main()
