from __future__ import division
import itertools
import time
from matplotlib import pyplot as plt
trinity_human = "./contigs2genes.disambiguous.txt"
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
    tp_orphan = []
    for clustName, clustMems in tr_clust_inv.iteritems():
        for tr_1, tr_2 in itertools.combinations(clustMems,2):
            if tr_1 not in groundTruth_clust or tr_2 not in groundTruth_clust:
                continue
            if groundTruth_clust[tr_1] == groundTruth_clust[tr_2]:
                tp += 1
                # look if we have any orphan link for it.
                #if (tr_1, tr_2) in orphan_reads:
                #    tp_orphan += [(tr_1, tr_2)]
                #elif (tr_2, tr_1) in orphan_reads:
                #    tp_orphan += [(tr_2, tr_1)]
            else:
                fp += 1
    helpful_orphan = []
    for clustName, clustMems in groundTruth_clust_inv.iteritems():
        for tr_1, tr_2 in itertools.combinations(clustMems,2):
            if tr_1 not in tr_clust or tr_2 not in tr_clust:
                continue
            if tr_clust[tr_1] != tr_clust[tr_2]:
                fn += 1
                # look if we have any orphan link for it.
                #if (tr_1, tr_2) in orphan_reads:
                #    helpful_orphan += [(tr_1, tr_2)]
                #elif  (tr_2, tr_1) in orphan_reads:
                #    helpful_orphan += [(tr_2, tr_1)]


    nc2 = (num * (num-1)) / 2
    tn = nc2 - (fp + tp + fn)
    return tp, fp, tn, fn, helpful_orphan, tp_orphan

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
        if tr_gn[0] in gtClusterCount.keys():
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

# def readOrphans(orphan_file):
#     orphan_reads = {}
#     contigs_info = {}
#     with open(orphan_file) as orphan:
#         seq_cnt = int(orphan.readline().rstrip())
#         for i in xrange(seq_cnt):
#             tname, tlen, ttmp = orphan.readline().rstrip().split('\t')
#             contigs_info[tname] = (int(tlen), float(ttmp))
#         for line in orphan:
#             key, value = line.rstrip().split(';')
#             l, r = key.rstrip().split("\t")
#             left, right = value.split(':')
#             left = map(int, left.rstrip().split('\t'))
#             right = map(int, right.rstrip().split('\t'))
#             orphan_reads[(l, r)] = (left, right)
#     return orphan_reads, contigs_info

def measurePrecRecall(qsf, orphanfile, sp, ctype):
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
    #t1 = time.time()
    #orphan_reads, contigs_info = readOrphans(orphanfile)
    #t2 = time.time()
    #print("It took {} secs to read orphan file".format(round(t2-t1)))
    groundTruth_clust, ground_truth_clust_inv = readTrueLabels(ft)
    tp, fp, tn, fn, helpful_orphans, tp_orphans = accuracyExpressedFast(groundTruth_clust, ground_truth_clust_inv,
                                                       tr_clust, tr_clust_inv)
    print("tp : {}, fp : {}, tn : {}, fn : {}".format(tp, fp, tn, fn))
    precision = tp / float(tp + fp)
    recall = tp / float(tp + fn)
    Fscore = 2 * ((precision*recall) / (precision+recall))
    print("prec: {}, recall: {}, F1: {}".format(precision, recall, Fscore))

    #print("fn: {}, orphan_reads: {}".format(fn, len(helpful_orphans)))

    # import numpy as np
    # read_len = 101
    # read_dists = {'tot':[], 'neut':[], 'help':[], 'harm':[]}
    # cnt = {'tot':[], 'neut':[], 'help':[], 'harm':[]}
    # tpm_ratio = {'tot':[], 'neut':[], 'help':[], 'harm':[]}
    # score =  {'tot':[], 'neut':[], 'help':[], 'harm':[]}
    #
    # helpful_orphans = set(helpful_orphans)
    # tp_orphans = set(tp_orphans)
    # for k, v in orphan_reads.iteritems():
    #
    #     curr_cnt = len(v[0]) # read count
    #
    #     # tpm ratio
    #     ltpm = contigs_info[k[0]][1]
    #     rtpm = contigs_info[k[1]][1]
    #     if ltpm == 0 or rtpm == 0: # TODO, for simplicity!
    #         curr_tpm_ratio = 0
    #     else:
    #         curr_tpm_ratio = np.min([ltpm/rtpm , rtpm/ltpm])
    #
    #     # avg read distance (avg fragment length)
    #     llen = contigs_info[k[0]][0]
    #     curr_read_dists = []
    #     r_extra_nucs = 0
    #     l_extra_nucs = 0
    #     for l, r in itertools.izip(v[0], v[1]):
    #         if r < 0:
    #             r_extra_nucs = max(abs(r), r_extra_nucs)
    #         l_extra_nucs = max(read_len - (llen-l), l_extra_nucs)
    #         curr_read_dists += [llen - l + r + read_len] # we assume that two reads won't overlap at all
    #     curr_read_dists = np.mean([r + r_extra_nucs + l_extra_nucs for r in curr_read_dists])
    #
    #     cnt['tot'] += [curr_cnt]
    #     read_dists['tot'] += [curr_read_dists]
    #     tpm_ratio['tot'] += [curr_tpm_ratio]
    #     if k in helpful_orphans:
    #         cnt['help'] += [curr_cnt]
    #         read_dists['help'] += [curr_read_dists]
    #         tpm_ratio['help'] += [curr_tpm_ratio]
    #     elif k in tp_orphans:
    #         cnt['neut'] += [curr_cnt]
    #         read_dists['neut'] += [curr_read_dists]
    #         tpm_ratio['neut'] += [curr_tpm_ratio]
    #     else:
    #         cnt['harm'] += [curr_cnt]
    #         read_dists['harm'] += [curr_read_dists]
    #         tpm_ratio['harm'] += [curr_tpm_ratio]
    #
    # cnt_min = np.min(cnt['tot']); cnt_max = np.max(cnt['tot'])
    # dist_min = np.min(read_dists['tot']); dist_max = np.max(read_dists['tot'])
    # tpm_min = np.min(tpm_ratio['tot']); tpm_max = np.max(tpm_ratio['tot'])
    # lbls = []
    # points_cnt = len(cnt['tot'])
    #
    # stats_file = "bin/human_rapclust/stats"
    #
    # with open(stats_file + ".txt", 'w') as general_stats:
    #     for k in ["harm", "help", "neut"]:
    #         with open(stats_file + "_" + k + ".txt", 'w') as stats:
    #             for i in xrange(len(cnt[k])):
    #                 general_stats.write("{}, {}, {}, {}, {}\n".format(k, i, read_dists[k][i], cnt[k][i], tpm_ratio[k][i]))
    #                 score[k] += [np.log( (cnt[k][i] - cnt_min + 1)/(cnt_max-cnt_min + points_cnt) ) +
    #                              np.log((dist_max-dist_min - (read_dists[k][i] - dist_min) + 1)/(dist_max-dist_min + points_cnt)) +
    #                              np.log( (tpm_ratio[k][i] - tpm_min + 1)/(tpm_max - tpm_min + points_cnt) )]
    #             print k
    #             print "-> cnt_mean: {}, cnt_std: {}, cnt_min: {}, cnt_max: {}".\
    #                 format(round(np.mean(cnt[k]), 2), round(np.std(cnt[k]), 2), np.min(cnt[k]), np.max(cnt[k]))
    #             print "-> dist_mean: {} , dist_std: {}, dist_min:{}, dist_max:{}".\
    #                 format(round(np.mean(read_dists[k]), 2), round(np.std(read_dists[k]), 2), np.min(read_dists[k]), np.max(read_dists[k]))
    #             print "-> tpm_mean: {} , tpm_std: {}, tpm_min:{}, tpm_max:{}".\
    #                 format(round(np.mean(tpm_ratio[k]), 2), round(np.std(tpm_ratio[k]), 2), np.min(tpm_ratio[k]), np.max(tpm_ratio[k]))
    #             print "-> score_mean: {}, score_std: {}, score_min: {}, score_max: {}". \
    #                 format(round(np.mean(score[k]), 2), round(np.std(score[k]), 2), np.min(score[k]), np.max(score[k]))
    #             cb = np.arange(np.max(cnt[k]))
    #             rb = np.arange(np.max(read_dists[k]))
    #             tb = np.arange(0, np.max(tpm_ratio[k]), 0.01)
    #             sb = np.arange(np.min(score[k]), np.max(score[k]), 1)
    #
    #             stats.write("\n\n" + k + " score\n")
    #             s, x = np.histogram(score[k], bins=sb)
    #             sn, xn = np.histogram(score[k], bins=sb, normed=1)
    #             #x = (x[1:] + x[:-1])/2
    #             x = x[:-1]
    #             plt.subplot(411)
    #             lbl, = plt.plot(xn[:-1], sn, label=k)
    #             lbls += [lbl]
    #             for a, b in zip(x[s != 0], s[s != 0]):
    #                 stats.write("{} : {}\n".format(a, b))
    #
    #             stats.write("\n\n" + k + " fragment size\n")
    #             s, x = np.histogram(read_dists[k], bins=rb)
    #             sn, xn = np.histogram(read_dists[k], bins=rb, normed=1)
    #             #x = (x[1:] + x[:-1])/2
    #             x = x[:-1]
    #             for a, b in zip(x[s != 0], s[s != 0]):
    #                 stats.write("{} : {}\n".format(a, b))
    #             plt.subplot(412)
    #             plt.plot(xn[:-1], sn, label=k)
    #
    #             stats.write("\n\n" + k + " tpm ratio\n")
    #             s, x = np.histogram(tpm_ratio[k], bins=tb)
    #             sn, xn = np.histogram(tpm_ratio[k], bins=tb, normed=1)
    #             #x = (x[1:] + x[:-1])/2
    #             x = x[:-1]
    #             for a, b in zip(x[s != 0], s[s != 0]):
    #                 stats.write("{} : {}\n".format(a, b))
    #             plt.subplot(413)
    #             plt.plot(xn[:-1], sn, label=k)
    #
    #             stats.write("\n\n" + k + " count\n")
    #             s, x = np.histogram(cnt[k], bins=cb)
    #             sn, xn = np.histogram(cnt[k], bins=cb[0:10], normed=1)
    #             #x = (x[1:] + x[:-1])/2
    #             x = x[:-1]
    #             for a, b in zip(x[s != 0], s[s != 0]):
    #                 stats.write("{} : {}\n".format(a, b))
    #             plt.subplot(414)
    #             plt.plot(xn[:-1], sn, label=k)

    #plt.legend(handles=lbls)
    #plt.show()


import argparse

def main():
    parser = argparse.ArgumentParser(description="Give the SRR number")
    parser.add_argument('--clustfile',type = str, help="graph file")
    parser.add_argument('--sp',type = str, help="type human")
    parser.add_argument('--ctype',type = str, default="mcl", help="you do not need to enter anything here")
    parser.add_argument('--orphanfile', type = str, help="name of file containing orphan links")
    args = parser.parse_args()
    measurePrecRecall(args.clustfile,args.orphanfile, args.sp, args.ctype)


if __name__ == "__main__":
    main()
