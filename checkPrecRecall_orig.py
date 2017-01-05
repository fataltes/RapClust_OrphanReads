from __future__ import division
import itertools
import time
from matplotlib import pyplot as plt
trinity_human = "./tx2GeneNameHuman.38.80.txt"
#trinity_human = "./contigs2genes.disambiguous.txt"
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
    tp, fp, tn, fn, orphan_cnt, cls_cntr, gt_cls_cntr, rprec, rrec = 0, 0, 0, 0, 0, 0, 0, 0, 0
    for clustName, clustMems in tr_clust_inv.iteritems():
        cls_cntr += 1
        cls_cnt, cls_tp = 0, 0
        for tr_1, tr_2 in itertools.combinations(clustMems,2):
            if tr_1 not in groundTruth_clust or tr_2 not in groundTruth_clust:
                continue
            cls_cnt += 1
            if groundTruth_clust[tr_1] == groundTruth_clust[tr_2]:
                tp += 1
                cls_tp += 1
            else:
                fp += 1
        if cls_cnt != 0:
            rprec += (cls_tp/cls_cnt)
    for clustName, clustMems in groundTruth_clust_inv.iteritems():
        if len(clustMems)>1:
            gt_cls_cntr += 1
        cls_cnt, cls_tp = 0, 0
        for tr_1, tr_2 in itertools.combinations(clustMems,2):
            if tr_1 not in tr_clust or tr_2 not in tr_clust:
                continue
            cls_cnt += 1
            if tr_clust[tr_1] != tr_clust[tr_2]:
                fn += 1
            else:
                cls_tp += 1
        if cls_cnt != 0:
            rrec += (cls_tp/cls_cnt)

    nc2 = (num * (num-1)) / 2
    tn = nc2 - (fp + tp + fn)
    return tp, fp, tn, fn, rprec/cls_cntr, rrec/gt_cls_cntr

def weightCorrelation(weightGraph, groundTruth_clust, groundTruth_clust_inv):
    class0 = []
    class1 = []
    corrupted = set()
    for key, weight in weightGraph.iteritems():
        if key[0] in groundTruth_clust and key[1] in groundTruth_clust:
            if groundTruth_clust[key[0]] == groundTruth_clust[key[1]]:
                class1 += [float(weight)]
            else:
                class0 += [float(weight)]
        else:
            if key[0] not in groundTruth_clust:
                corrupted.add(key[0])
            else:
                corrupted.add(key[1])
    print ("corrupted: {}".format(len(corrupted)))
    #print(corrupted)
    print ("")
    from matplotlib import pyplot as plt
    plt.figure()
    plt.subplot(211)
    plt.hist(class0, bins=200)
    plt.title("Different Genes")
    plt.ylabel("frequency")
    plt.subplot(212)
    plt.hist(class1, bins=200)
    plt.title("Same Gene")
    plt.ylabel("frequency")
    plt.xlabel("link weights")

def orphanCorrelation(orphanPair, groundTruth_clust, groundTruth_clust_inv):
    class0 = [[] for x in range(3)]
    class1 = [[] for x in range(3)]
    corrupted = set()
    for key, orphan in orphanPair.iteritems():
        if key[0] in groundTruth_clust and key[1] in groundTruth_clust:
            if groundTruth_clust[key[0]] == groundTruth_clust[key[1]]:
                class1[0] += [float(orphan[0])]
                class1[1] += [float(orphan[1])]
                class1[2] += [float(orphan[2])]
            else:
                class0[0] += [float(orphan[0])]
                class0[1] += [float(orphan[1])]
                class0[2] += [float(orphan[2])]
        else:
            if key[0] not in groundTruth_clust:
                corrupted.add(key[0])
            else:
                corrupted.add(key[1])
    print("corrupted: {}".format(len(corrupted)))
    #print(corrupted)
    print("")
    print("# in same gene: {}, # in different gene:{}".format(len(class1[0]), len(class0[0])))

    from matplotlib import pyplot as plt
    f, ax = plt.subplots(2, 3)
    ax[0, 0].hist([c for c in class0[0] if c <= 50], bins=range(1, 51))
    ax[0, 0].set_title("Different Genes")
    ax[0, 0].set_ylabel("frequency")
    ax[1, 0].hist([c for c in class1[0] if c <= 50], bins=50)
    ax[1, 0].set_title("Same Gene")
    ax[1, 0].set_ylabel("frequency")
    ax[1, 0].set_xlabel("orphan count")

    ax[0, 1].hist([c for c in class0[1] if c <= 5000], bins=1000)
    ax[0, 1].set_title("Different Genes")
    ax[1, 1].hist([c for c in class1[1] if c <= 5000], bins=1000)
    ax[1, 1].set_title("Same Gene")
    ax[1, 1].set_xlabel("orphan dist")

    import numpy as np
    ax[0, 2].hist([np.log(c) if c > 0 else -40 for c in class0[2]], bins=100)
    ax[0, 2].set_title("Different Genes")
    ax[1, 2].hist([np.log(c) if c > 0 else -40 for c in class1[2]], bins=100)
    ax[1, 2].set_title("Same Gene")
    ax[1, 2].set_xlabel("orphan tpm")


def readOrphans(contig_file, orphan_file):
    orphan_reads = {}
    contigs_info = []
    with open(contig_file) as contigs:
        seq_cnt = int(contigs.readline().rstrip())
        for i in xrange(seq_cnt):
            tname, tlen, ttmp = contigs.readline().rstrip().split('\t')
            contigs_info.append(tname)
    with open(orphan_file) as orphan:
        for line in orphan:
            key, value = line.rstrip().split(';') # key : c1, c2 & value : count, distance, tpm_ratio
            l, r = map(int, key.rstrip().split(","))
            val = map(float, value.rstrip().split('\t'))
            orphan_reads[(contigs_info[l], contigs_info[r])] = val
    # with open(orphan_file) as orphan:
    #     seq_cnt = int(orphan.readline().rstrip())
    #     for i in xrange(seq_cnt):
    #         tname, tlen, ttmp = orphan.readline().rstrip().split('\t')
    #         contigs_info[tname] = (int(tlen), float(ttmp))
    #     for line in orphan:
    #         key, value = line.rstrip().split(';') # key : c1, c2 & value : count, distance, tpm_ratio
    #         l, r = key.rstrip().split("\t")
    #         left, right = value.split(':')
    #         left = map(int, left.rstrip().split('\t'))
    #         right = map(int, right.rstrip().split('\t'))
    #         orphan_reads[(l, r)] = (left, right)
    return orphan_reads #, contigs_info

def sortOrphanReads(orphan_file, orphan_reads):
    orphan_reads = sorted(orphan_reads.items(), key=lambda e: e[1][0])
    with open(orphan_file+"_sorted", "w") as of:
        for t in orphan_reads:
            of.write("{}, {}; {}, {}, {}\n".format(t[0][0], t[0][1], t[1][0], t[1][1], t[1][2]))

def contigSupport(orphan_reads):
    read_cnt_dist = {}
    for k, v in orphan_reads.iteritems():
        if v[0] not in read_cnt_dist:
            read_cnt_dist[v[0]] = 0
        read_cnt_dist[v[0]] += 1
    return read_cnt_dist

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
        tr_gn = map(str.strip, line[:-1].split("\t"))
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
    dir = "transcripts/human_rapclust"
    weightGraph = readNetFile("{}/mag.filt.net".format(dir))
    tp, fp, tn, fn, rprec, rrec = accuracyExpressedFast(groundTruth_clust, ground_truth_clust_inv,
                                                       tr_clust, tr_clust_inv)
    weightCorrelation(weightGraph, groundTruth_clust, ground_truth_clust_inv)
    print("tp : {}, fp : {}, tn : {}, fn : {}".format(tp, fp, tn, fn))
    precision = tp / float(tp + fp)
    recall = tp / float(tp + fn)
    Fscore = 2 * ((precision*recall) / (precision+recall))
    rfscore = 2 * (rprec*rrec) / (rprec + rrec)
    print("prec: {}, recall: {}, F1: {}".format(precision, recall, Fscore))
    print("rel. prec: {}, rel. recall:{}, rel. F1: {}".format(rprec, rrec, rfscore))
    orphan_reads = readOrphans("{}/orphanLink".format(dir), "{}/orphanLink_contigpair.txt".format(dir))
    orphanCorrelation(orphan_reads, groundTruth_clust, ground_truth_clust_inv)
    sortOrphanReads("{}/orphanLink".format(dir), orphan_reads)
    #read_cnt_dist = contigSupport(orphan_reads)
    #plt.figure()
    #plt.plot(read_cnt_dist.keys()[0:5], read_cnt_dist.values()[0:5])
    #plt.xlabel("# reads supporting contig pair")
    #plt.ylabel("frequency")
    plt.show()

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
