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

tnames_inv = {}
tnames = []


def buildNetFile(sampdirs, netfile, orphanLink_out_file, cutoff, auxDir, writecomponents=False):


    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampdirs]
    eqfiles = [sep.join([sd, auxDir, '/eq_classes.txt']) for sd in sampdirs]

    weightDict = {}
    diagCounts = None
    sumCounts = None
    ambigCounts = None
    firstSamp = True
    numSamp = 0
    tot = 0
    eqClasses = {}
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
                diagCounts = np.zeros(len(tnames))
                sumCounts = np.zeros(len(tnames))
                ambigCounts = np.zeros(len(tnames))
            else:
                for i in xrange(numTran):
                    ifile.readline()

            for i in range(len(tnames)):
                tnames_inv[tnames[i]]=i
            # for easy access to quantities of interest
            tpm = quant.loc[tnames, 'TPM'].values
            estCount = quant.loc[tnames, 'NumReads'].values
            efflens = quant.loc[tnames, 'EffectiveLength'].values
            epsilon =  np.finfo(float).eps
            sumCounts = np.maximum(sumCounts, estCount)

            for i in xrange(numEq):
                toks = map(int, ifile.readline().rstrip().split('\t'))
                nt = toks[0]
                tids = tuple(toks[1:-1])
                count = toks[-1]
                if tids in eqClasses:
                    eqClasses[tids] += count
                else:
                    eqClasses[tids] = count

                # Add the contribution to the graph
                denom = sum([tpm[t] for t in tids])
                for t1, t2 in itertools.combinations(tids,2):
                    #tpm1 = tpm[t1]
                    #tpm2 = tpm[t2]
                    #w = count * ((tpm1 + tpm2) / denom)
                    w = count
                    key = (t1, t2)
                    if key in weightDict:
                        weightDict[key] += w
                    else:
                        weightDict[key] = w
                for t in tids:
                    diagCounts[t] += count * (tpm[t] / denom)
                    ambigCounts[t] += count
            firstSamp = False

    lens = quant.loc[tnames, 'Length'].values

    minWeight = 0.5
    maxWeight = 0.0
    prior = 0.1
    edgesToRemove = []

    ##
    #  Go through the weightMap and remove any edges that
    #  have endpoints with too few mapping reads
    ##
    for k,v in weightDict.iteritems():
        c0, c1 = diagCounts[k[0]], diagCounts[k[1]]
        a0, a1 = ambigCounts[k[0]], ambigCounts[k[1]]
        if a0 + a1 > epsilon and a0 > cutoff and a1 > cutoff:# and v > 3:
            w = (v+prior) / min((a0+prior), (a1+prior))
            weightDict[k] = w
            if w > maxWeight:
                maxWeight = w
            #if w < minWeight:
            #    edgesToRemove.append(k)
        else:
            edgesToRemove.append(k)

    # Actually delete those edges
    for e in edgesToRemove:
        del weightDict[e]

    # weightDict = use_orphan_reads_simple(sampdirs, auxDir, '/orphan_links.txt', lens, tpm, ambigCounts, diagCounts,
    #                                      cutoff, weightDict)

    # optimally possible
    #weightDict = optimal_prec_filter(weightDict, tnames, tpm, ambigCounts, diagCounts)

    tnamesFilt = []
    relabel = {}
    for i in xrange(len(estCount)):
        if (ambigCounts[i] > cutoff):
            relabel[i] = len(tnamesFilt)
            tnamesFilt.append(tnames[i])
            weightDict[(i, i)] = 1.1

    import networkx as nx
    G = nx.Graph() if writecomponents else None
    with open(netfile, 'w') as ofile:
        writeEdgeList(weightDict, tnames, ofile, G)

    if G is not None:
        clustFile = netfile.split('.net')[0] + '.clust'
        print("Writing connected components as clusters to {}".format(clustFile))
        with open(clustFile, 'w') as ofile:
            cc = nx.connected_component_subgraphs(G)
            for c in cc:
                ofile.write('{}\n'.format('\t'.join(c.nodes())))
    #print("#edges: {}".format(len(weightDict.keys())))
def writeEdgeList(weightDict, tnames, ofile, G):
    useGraph = G is not None
    for k,v in weightDict.iteritems():
        ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))
        if useGraph:
            G.add_edge(tnames[k[0]], tnames[k[1]])


def writePajek(weightDict, tnames, relabel, ofile):
    with open(netfile, 'w') as ofile:
        ofile.write("*Vertices\t{}\n".format(len(tnamesFilt)))
        for i, n in enumerate(tnamesFilt):
            ofile.write("{}\t\"{}\"\n".format(i, n))
        ofile.write("*Edges\n")
        print("There are {} edges\n".format(len(weightDict)))
        for k,v in weightDict.iteritems():
            ofile.write("{}\t{}\t{}\n".format(relabel[k[0]], relabel[k[1]], v))
            #ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))
            #if k[0] != k[1]:
            #    ofile.write("{}\t{}\t{}\n".format(tnames[k[1]], tnames[k[0]], v))

class EquivCollection(object):
    def __init__(self):
        self.tnames = []
        self.eqClasses = {}
        self.hasNames = False

    def setNames(self, names):
        self.tnames = names
        self.hasNames = True

    def add(self, tids, count):
        if tids in self.eqClasses:
            self.eqClasses[tids] += count
        else:
            self.eqClasses[tids] = count

def readEqClass(eqfile, eqCollection):
    with open(eqfile) as ifile:
        numTran = int(ifile.readline().rstrip())
        numEq = int(ifile.readline().rstrip())
        print("file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
        if not eqCollection.hasNames:
            tnames = []
            for i in xrange(numTran):
                tnames.append(ifile.readline().rstrip())
            eqCollection.setNames(tnames)
        else:
            for i in xrange(numTran):
                ifile.readline()

        for i in xrange(numEq):
            toks = map(int, ifile.readline().rstrip().split('\t'))
            nt = toks[0]
            tids = tuple(toks[1:-1])
            count = toks[-1]
            eqCollection.add(tids, count)

def getCountsFromEquiv(eqCollection):
    countDict = {}
    tn = eqCollection.tnames
    for tids, count in eqCollection.eqClasses.iteritems():
        for t in tids:
            if tn[t] in countDict:
                countDict[tn[t]] += count
            else:
                countDict[tn[t]] = count
    # ensure no division by 0
    for t in eqCollection.tnames:
        if t in countDict:
            countDict[t] += 1.0
        else:
            countDict[t] = 1.0
    return countDict

def flattenClusters(infile, outfile):
    with open(outfile, 'w') as ofile:
        with open(infile) as ifile:
            for i,l in enumerate(ifile):
                toks = l.rstrip().split()
                cname = "cluster{}".format(i)
                for t in toks:
                    ofile.write("{}\t{}\n".format(cname, t))

def filterGraph(expDict, netfile, ofile, auxDir):
    import os
    import pandas as pd
    import math
    import logging
    from tqdm import tqdm

    logger = logging.getLogger("rapclust")
    # Get just the set of condition names
    conditions = expDict.keys()
    logging.info("conditions = {}".format(conditions))

    #for cond in conditions:
    #    sailfish[cond] = collections.defaultdict(float)
    #    for sample in samples:
    #        fileDir = cond+sample
    #        filePath = os.path.sep.join([samplesPath, fileDir, "quant.sf"])
    #        print(filePath)
    #        with open(filePath) as f:
    #            data = pd.read_table(f, header=0)
    #            for i in range(len(data)):
    #                name = data['Name'][i]
    #                numreads = data['NumReads'][i]
    #                sailfish[cond][name] += numreads
    #    # To avoid divide by 0 error
    #    for name in sailfish[cond]:
    #        sailfish[cond][name] += 1

    eqClasses = {}
    for cond in conditions:
        print(expDict[cond])
        for sampNum, sampPath in expDict[cond].iteritems():
            if cond not in eqClasses:
                eqClasses[cond] = EquivCollection()
            eqPath = os.path.sep.join([sampPath, auxDir, "/eq_classes.txt"])
            readEqClass(eqPath, eqClasses[cond])

    ambigCounts = {cond : getCountsFromEquiv(eqClasses[cond]) for cond in conditions}

    sailfish = {}
    for cond in conditions:
        sailfish[cond] = ambigCounts[cond]

    logging.info("Done Reading")
    count = 0
    numTrimmed = 0
    infomap_ofile = ofile.split('.net')[0] + '.txt'
    mlrmcl_ofile = ofile.split('.net')[0] + '.mlrmcl.graph'
    mlrmcl = [[] for i in range(len(tnames_inv))]
    numEdges = 0
    #import pdb
    with open(netfile) as f, open(ofile, 'w') as ofile, open(infomap_ofile, 'w') as infomap_ofile:
        data = pd.read_table(f, header=None)
        for i in tqdm(range(len(data))):
            count += 1
            #print("\r{} done".format(count), end="")
            #Alternative hypo
            x = data[0][i]
            y = data[1][i]
            non_null=0
            x_all=0
            y_all=0
            for cond in conditions:
                y_g = sailfish[cond][y]
                x_g = sailfish[cond][x]
                r = y_g / x_g
                non_null += (y_g * math.log(r*x_g)) - (r*x_g)
                non_null += (x_g * math.log(x_g)) - x_g
                x_all += x_g
                y_all += y_g
            #null hypothesis
            null = 0
            r_all = y_all / x_all
            for cond in conditions:
                y_g = sailfish[cond][y]
                x_g = sailfish[cond][x]
                mean_x = (x_g + y_g) / (1+r_all)
                null += (y_g * math.log(r_all * mean_x)) - (r_all * mean_x)
                null += (x_g * math.log(mean_x)) - mean_x
            D = 2*(non_null-null)
            if D <= 20:
                ofile.write("{}\t{}\t{}\n".format(x, y, data[2][i]))
                infomap_ofile.write("{}\t{}\t{}\n".format(tnames_inv[x], tnames_inv[y], data[2][i]))
                #pdb.set_trace()
                mlrmcl[tnames_inv[x]].append((tnames_inv[y], int(data[2][i]*1000)))
                numEdges+=1
            else:
                numTrimmed += 1
    #import pdb
    #pdb.set_trace()
    with open(mlrmcl_ofile, 'w') as mlrmclf:
        mlrmclf.write("{}\t{}\t1".format(len(tnames_inv), numEdges))
        for i in range(len(mlrmcl)):
            for vertex, weight in mlrmcl[i]:
                if weight > 0:
                    mlrmclf.write("{}\t{}\t".format(vertex+1, weight))
            mlrmclf.write("{}\t1100\n".format(i+1)) #adding selfloops todo:check if you are already adding all selfloops to the graph! and then these might be excluded from the list of nodes in the graph
    logging.info("Trimmed {} edges".format(numTrimmed))


def use_orphan_reads_simple(sampdirs, auxDir, orphanFileName, lens, tpm, ambigCounts, diagCounts, cutoff, weightDict):
    # Considering Orphan reads
    def nearEnd(tup):
        txp = tup[0]
        pos = tup[1]
        moverhang = 10
        ml = 100
        if pos < -moverhang or pos > lens[txp] + moverhang:
            return False
        elif pos <= ml or pos >= lens[txp] - ml:
            return True
        else:
            return False

    orphanLinkFiles = [sep.join([sd, auxDir, orphanFileName]) for sd in sampdirs]
    haveLinkFiles = all(os.path.isfile(f) for f in orphanLinkFiles)
    if haveLinkFiles:
        numOrphanLinks = 0
        for olfile in orphanLinkFiles:
            for l in open(olfile):
                left, right = l.rstrip().split(':')
                lp = [map(int, i.split(',')) for i in left.rstrip('\t').split('\t')]
                rp = [map(int, i.split(',')) for i in right.split('\t')]
                lp = [t for t in filter(nearEnd, lp)]
                rp = [t for t in filter(nearEnd, rp)]
                #if len(lp) == 1 or len(rp) == 1:
                for a, b in itertools.product(lp, rp):
                    ltpm = tpm[a[0]] + 10 ** -11  # Laplacian Smoothing
                    rtpm = tpm[b[0]] + 10 ** -11
                    tpm_ratio = 1 - (abs(ltpm - rtpm) / (ltpm + rtpm))
                    read_dist = lens[a[0]] - a[1] + b[1]
                    if tpm_ratio >= .5: # and read_dist <= 300 and tpm[a[0]] > .5 and tpm[b[0]] > .5:
                        a = a[0]; b = b[0]
                        key = (a, b) if a < b else (b, a)
                        if ambigCounts[a] < cutoff or ambigCounts[b] < cutoff:
                            continue
                        c0, c1 = diagCounts[a], diagCounts[b]
                        a0, a1 = ambigCounts[a], ambigCounts[b]
                        if key not in weightDict:
                            numOrphanLinks += 1
                            weightDict[key] = 1.0 / min(a0, a1)
                        else:
                            weightDict[key] += 1.0 / min(a0, a1)
        logging.info("Added {} orphan link edges".format(numOrphanLinks))


    # CNT_IDX =0; READ_DIST_IDX=1; TPM_RATIO_IDX=2
    # use_orphan_reads_complicated(tnames, tpm, lens, sampdirs, orphanLink_out_file, auxDir, load_calculated_file=False)
    #
    # import sys
    # cnt_min = sys.maxsize
    # cnt_max = 0
    # dist_min = sys.maxsize
    # dist_max = 0
    # tpm_min = 1
    # tpm_max = 0
    # for k, val in orphan_pair.iteritems():
    #     cnt_min = min(val[CNT_IDX], cnt_min)
    #     cnt_max = max(val[CNT_IDX], cnt_max)
    #     dist_min = min(val[READ_DIST_IDX], dist_min)
    #     dist_max = max(val[READ_DIST_IDX], dist_max)
    #     tpm_min = min(val[TPM_RATIO_IDX], tpm_min)
    #     tpm_max = max(val[TPM_RATIO_IDX], tpm_max)
    #
    # print("count: Min = {}, Max = {}".format(cnt_min, cnt_max))
    # print("read dist: Min = {}, Max = {}".format(dist_min, dist_max))
    # print("tpm: Min = {}, Max = {}".format(tpm_min, tpm_max))
    # #points_cnt = len(cnt.keys())
    #
    # increased = 0
    # added = 0
    # min_score = 1
    # max_score = 0
    # score_cnt = 0
    # score_sum = 0
    # score_squared_sum = 0
    # scores = []
    # for k, val in orphan_pair.iteritems():
    #     vals = [(val[CNT_IDX] - cnt_min+1) / (cnt_max - cnt_min+1), \
    #             (1 - ((val[READ_DIST_IDX] - dist_min) / (dist_max - dist_min))), \
    #             (val[TPM_RATIO_IDX] - tpm_min) / (tpm_max - tpm_min)
    #             ]
    #     vals = np.sort(vals)
    #     score = vals[1]*vals[2]#*vals[0]
    #     scores += [score]
    #     score_cnt += 1
    #     score_sum += score
    #     score_squared_sum += score**2
    #     min_score = score if score < min_score else min_score
    #     max_score = score if score > max_score else max_score
    #     if score >= 0.7:
    #         if k in weightDict:
    #             weightDict[k] += score
    #             increased += 1
    #         else:
    #             weightDict[k] = score
    #             added += 1
    # score_mean = score_sum/score_cnt
    # score_std = (score_squared_sum/score_cnt - score_mean**2)**0.5
    # print ("Score: Min = {} , Max = {}, Mean = {}, STD = {}".format(min_score, max_score, score_mean, score_std))
    # print ("Links: Added = {} , Value Increased = {} ".format(added, increased))
    # #from matplotlib import pyplot as plt
    # #plt.hist(scores, bins = 200)
    # #plt.show()
    # #End of Orphan read section
    return weightDict;
def use_orphan_reads_complicated(tnames, tpm, lens, sampdirs, orphanLink_out_file, auxDir, load_calculated_file = True):
    import itertools
    import numpy as np
    import os
    sep = os.path.sep

    orphan_pair = {}

    if load_calculated_file:
        line_cntr = 0
        s1 = time.time()
        for l in open(orphanLink_out_file + "_contigpair.txt"):
            key, val = l.rstrip().split(';')
            left, right = map(int, key.split(','))
            vals = map(float, val.split('\t'))
            orphan_pair[(left, right)] = vals
            line_cntr += 1
        e1 = time.time()
        print("file: {}; # pairs = {}; time = {} secs".format(orphanLink_out_file + "_contigpair.txt", line_cntr,
                                                              round(e1 - s1)))
        return orphan_pair

    CNT_IDX =0; WEIGHTED_CNT_IDX=1; UNIPARTIAL_CNT_IDX=2; WINNER_CNT=3
    AVG_DIST_IDX=4; WEIGHTED_AVG_DIST_IDX=5; UNIPARTIAL_AVG_DIST_IDX=6; WINNER_AVG_DIST_IDX=7
    TPM_RATIO_IDX=8; TPM_TPM_RATIO_IDX=9
    R_EXTRA_NUCS_IDX=10; L_EXTRA_NUCS_IDX= 11
    read_len = 101
    orphanLinkFiles = [sep.join([sd, auxDir, '/orphan_links.txt']) for sd in sampdirs]
    haveLinkFiles = all(os.path.isfile(f) for f in orphanLinkFiles)
    if haveLinkFiles:
        orphan_reads = {}
        for olfile in orphanLinkFiles:
            s1 = time.time()
            new_cntr = 0
            orphan_pair_cntr = 0
            total_orpha_pair_cntr = 0
            for l in open(olfile):
                left, right = l.rstrip().split(':')
                lp = [map(int, i.split(',')) for i in left.rstrip('\t').split('\t')]
                rp = [map(int, i.split(',')) for i in right.split('\t')]

                read_dist = {}
                tpm_ratio = {}
                tpm_tpm_ratio = {}
                read_pos = {}

                # for each read pair, we might have multiple contig pairs
                # In the following we would choose the best based with maximum score:
                # score = (1 - standardized read distance) * standardized tpm_ratio
                # The values are just standardized among all possible contig pairs for this read pair
                for l, r in itertools.product(lp, rp): # read all possible contig pairs for each read pair
                    lidx = l[0]
                    ridx = r[0]
                    ltpm = tpm[lidx]+10**-11 # Laplacian Smoothing
                    rtpm = tpm[ridx]+10**-11
                    tpm_tpm_ratio[(lidx, ridx)] = 1 - (abs(ltpm-rtpm)/(ltpm+rtpm))
                    tpm_ratio[(lidx, ridx)] = min(ltpm/rtpm, rtpm/ltpm)
                    read_dist[(lidx, ridx)] = lens[lidx] - l[1] + r[1] + read_len
                    read_pos[(lidx, ridx)] = (l[1], r[1])
                    total_orpha_pair_cntr += 1

                pair_cnt = len(tpm_ratio)
                pair_weight_sum = sum(tpm_tpm_ratio.values())

                maxKey = (-1, -1)
                maxScore = -1
                min_read_dist = np.min(read_dist.values()) - 1
                max_read_dist = np.max(read_dist.values())
                # choose contig pair with maximum score
                for k in read_dist.keys():
                    orphan_pair_cntr += 1
                    weighted_cnt = tpm_tpm_ratio[k] / pair_weight_sum
                    unipartial_cnt = 1 / pair_cnt
                    weighted_dist = weighted_cnt * read_dist[k]
                    unipartial_dist = unipartial_cnt * read_dist[k]
                    if k not in orphan_pair:
                        new_cntr += 1
                        orphan_pair[k] = [0, 0, 0, 0,
                                          [], [], [], [],
                                          tpm_ratio[k], tpm_tpm_ratio[k], 0, 0]
                    orphan_pair[k][CNT_IDX] += 1
                    orphan_pair[k][WEIGHTED_CNT_IDX] += weighted_cnt
                    orphan_pair[k][UNIPARTIAL_CNT_IDX] += unipartial_cnt
                    orphan_pair[k][AVG_DIST_IDX] += [read_dist[k]]
                    orphan_pair[k][WEIGHTED_AVG_DIST_IDX] += [weighted_dist]
                    orphan_pair[k][UNIPARTIAL_AVG_DIST_IDX] += [unipartial_dist]
                    orphan_pair[k][R_EXTRA_NUCS_IDX] = max(-1 * read_pos[k][1],
                                                                orphan_pair[k][R_EXTRA_NUCS_IDX])
                    orphan_pair[k][L_EXTRA_NUCS_IDX] = max(read_pos[k][0] + read_len - lens[lidx],
                                                                orphan_pair[k][L_EXTRA_NUCS_IDX])
                    score = (1-(read_dist[k] - min_read_dist)/(max_read_dist-min_read_dist))*(tpm_tpm_ratio[k])
                    if score > maxScore:
                        maxKey = k
                        maxScore = score
                orphan_pair[maxKey][WINNER_CNT] += 1
                orphan_pair[maxKey][WINNER_AVG_DIST_IDX] += [read_dist[maxKey]]

            e1 = time.time()
            logger.info("file: {}; orphaned read pairs: total = {}; # chosen = {}; # new = {}; time = {} secs". \
                  format(olfile, total_orpha_pair_cntr, orphan_pair_cntr, new_cntr, round(e1 - s1)))
        # averaging over all read distances between two transcripts
        # before averaging add left extra nucleotides and right extra nucleotides to the distance for those reads
        # passing the transcript length
        for k, val in orphan_pair.iteritems():
            orphan_pair[k][AVG_DIST_IDX] = round(np.sum(
                [r + val[R_EXTRA_NUCS_IDX] + val[L_EXTRA_NUCS_IDX] for r in val[AVG_DIST_IDX]])/val[CNT_IDX])
            orphan_pair[k][WEIGHTED_AVG_DIST_IDX] = round(np.sum(
                [r + val[R_EXTRA_NUCS_IDX] + val[L_EXTRA_NUCS_IDX] for r in val[WEIGHTED_AVG_DIST_IDX]]) \
                                                    / val[WEIGHTED_CNT_IDX])
            orphan_pair[k][UNIPARTIAL_AVG_DIST_IDX] = round(np.sum(
                [r + val[R_EXTRA_NUCS_IDX] + val[L_EXTRA_NUCS_IDX] for r in val[UNIPARTIAL_AVG_DIST_IDX]]) \
                                                      / val[UNIPARTIAL_CNT_IDX])
            orphan_pair[k][WINNER_AVG_DIST_IDX] = round(np.sum(
                [r + val[R_EXTRA_NUCS_IDX] + val[L_EXTRA_NUCS_IDX] for r in val[WINNER_AVG_DIST_IDX]]) \
                                                  / val[WINNER_CNT])
            orphan_pair[k][WEIGHTED_CNT_IDX] = round(orphan_pair[k][WEIGHTED_CNT_IDX], 2)
            orphan_pair[k][UNIPARTIAL_CNT_IDX] = round(orphan_pair[k][UNIPARTIAL_CNT_IDX], 2)

            orphan_pair[k][TPM_RATIO_IDX] = round(orphan_pair[k][TPM_RATIO_IDX], 4)
            orphan_pair[k][TPM_TPM_RATIO_IDX] = round(orphan_pair[k][TPM_TPM_RATIO_IDX], 4)

        with open(orphanLink_out_file + "_contigpair.txt", 'w') as ofile:
            for k, v in orphan_pair.iteritems():
                ofile.write('{},{};{}\n'.format(tnames[k[0]], tnames[k[1]],
                                             '\t'.join(str(e) for e in v[CNT_IDX:TPM_TPM_RATIO_IDX+1])))


    return orphan_pair

def optimal_prec_filter(weightDict, tnames, tpm, ambigCounts, diagCounts):
    ##
    #  Go through the weightMap and remove any edges that
    #  are between two paralogous contigs
    ##
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
        for k, v in groundTruth_clust.iteritems():
            if v in groundTruth_clust_inv:
                groundTruth_clust_inv[v].append(k)
            else:
                groundTruth_clust_inv[v] = [k]
        return groundTruth_clust, groundTruth_clust_inv

    groundTruth_clust, ground_truth_clust_inv = readTrueLabels("../contigs2genes.disambiguous.txt")
    edgesToRemove = []

    for k,v in weightDict.iteritems():
        if (tnames[k[0]] in groundTruth_clust and tnames[k[1]] in groundTruth_clust and
                    groundTruth_clust[tnames[k[0]]] != groundTruth_clust[tnames[k[1]]]):
            edgesToRemove.append(k)
    # Actually delete those edges
    tpml = []
    countl = []
    commonReadl = []
    for e in edgesToRemove:
        ltpm = tpm[e[0]] + 10**-12
        rtpm = tpm[e[1]] + 10**-12
        tpml += [1 - (abs(ltpm - rtpm) / (ltpm + rtpm))]
        countl = [min(ambigCounts[e[0]], ambigCounts[e[1]])]
        commonReadl = [weightDict[e]*min(ambigCounts[e[0]], ambigCounts[e[1]])]
        print ("{}: {}, {}, {}".format(e, tpml[-1], countl[-1], commonReadl[-1]))
        del weightDict[e]
    from matplotlib import pyplot as plt
    f, ax = plt.subplots(3, 1)
    ax[0].hist(tpml, bins=100)
    ax[1].hist(countl, bins=100)
    ax[2].hist(commonReadl, bins=100)

    plt.show()

    return weightDict

def convert2clusterFormat(infile, clustOutfile, flatClustOutfile):
    clusts = {}
    with open(infile) as ifile:
        next(ifile); next(ifile)
        for i,l in enumerate(ifile):
            toks = l.rstrip().split()
            if toks[1] not in clusts:
                clusts[toks[1]] = []
            clusts[toks[1]] += [tnames[int(toks[0])]]

    with open(flatClustOutfile, 'w') as fofile, open(clustOutfile, 'w') as ofile:
        for k, v in clusts.iteritems():
            cname = "cluster{}".format(k)
            new_cluster = ""
            for member in v:
                fofile.write("{}\t{}\n".format(cname, member))
                new_cluster += "{}\t".format(member)
            ofile.write("{}\n".format(new_cluster.strip()))
