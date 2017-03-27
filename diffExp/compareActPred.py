def catchFalsehood(fn_genes):
    allgenes = [{}, {}]
    tp = [set(), set()]
    fp = [set(), set()]
    for i in [0, 1]:
        cntr = 0
        with open(fn_genes[i]) as fo:
            for line in fo:
                if cntr <= 10000:
                    tups = line.strip().split() # tups is cluster, tp list [, fp list]
                    if len(tups) > 1:
                        for gene in tups[1].split(';'):
                            tp[i].add(gene)
                            if gene not in allgenes[i]:
                                allgenes[i][gene] = []
                            allgenes[i][gene] += [tups[0]]
                    if len(tups) > 2:
                        for gene in tups[2].split(';'):
                            fp[i].add(gene)
                            if gene not in allgenes[i]:
                                allgenes[i][gene] = []
                            allgenes[i][gene] += [tups[0]]
                cntr += 1
    tp2fn = tp[0].difference(tp[1]) # tp in orig, but not in orphan
    tn2fp = fp[1].difference(fp[0]) # not fp in orig, but fp in orphan
    return tp2fn, tn2fp 
  
def calcStats(tp2fn, tn2fp, fn_clusts, fn_contGenMap):
    from scipy.stats import entropy
    # fill out gene to contig map
    gen2cont = {}
    cont2gen = {}
    with open(fn_contGenMap) as fp:
        for line in fp:
            contig, gene = line.strip().split()            
            if gene not in gen2cont:
                gen2cont[gene] = []
            gen2cont[gene] += [contig]
            cont2gen[contig] = gene
    # fill out contig to cluster map & cluster to contigs map
    cont2clus = [{}, {}]
    clus2cont = [{}, {}]
    clusPurity = [{}, {}]
    for i in [0, 1]:
        with open(fn_clusts[i]) as fo:
            for line in fo:
                clus, cont = line.strip().split()
                cont2clus[i][cont] = clus
                if clus not in clus2cont[i]:
                    clus2cont[i][clus] = []
                clus2cont[i][clus] += [cont]
        # Calculate purity for each cluster
        #import pdb;pdb.set_trace()
        for clus, contigs in clus2cont[i].iteritems():
            clusGenes = {}
            for cont in contigs:
                if cont in cont2gen:
                    clusGen = cont2gen[cont]
                    if clusGen not in clusGenes:
                        clusGenes[clusGen] = 0
                    clusGenes[clusGen] += 1
            clusPurity[i][clus] = entropy(clusGenes.values()) 

    cntDiff = {}
    sizeDiff = {}
    purityDiff = {}
    contigCnt = {}
    genDiff = {}
    falsehood = {'tn2fp':tn2fp, 'tp2fn':tp2fn}
    problematics = set()
    for errType in ['tn2fp', 'tp2fn']:
        cntDiff[errType] = {}; sizeDiff[errType] = {}; purityDiff[errType] = {}; contigCnt[errType] = {}; genDiff[errType] = {}
        for gen in falsehood[errType]:
            clusSet = set()
            cntDiff[errType][gen] = [0, 0]
            sizeDiff[errType][gen] = [[], []]
            purityDiff[errType][gen] = [[], []]
            genDiff[errType][gen] = [set(), set()]
            contigCnt[errType][gen] = len(gen2cont[gen])           
            for contig in gen2cont[gen]:
                if contig in cont2clus[0]:
                    for i in [0, 1]:
                        clus = cont2clus[i][contig]
                        if clus not in clusSet:
                            clusSet.add(clus)
                            cntDiff[errType][gen][i] += 1
                            sizeDiff[errType][gen][i] += [len(clus2cont[i][clus])]
                            purityDiff[errType][gen][i] += [clusPurity[i][clus]]
                            for cc in clus2cont[i][clus]:
                                if cc in cont2gen:
                                    genDiff[errType][gen][i].add(cont2gen[cc])
                                else:
                                    problematics.add(cc)
    print(problematics)
    return contigCnt, cntDiff, sizeDiff, purityDiff, genDiff

import matplotlib
matplotlib.use("Pdf")
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
def saveAndPlot(contigCnt, cntDiff, sizeDiff, purityDiff, genDiff):
    import pdb
    for errType in ['tn2fp', 'tp2fn']:
        pp = PdfPages('plots/'+errType+'.pdf')
        p1 = plt.figure(); plt.hist(contigCnt[errType].values()); plt.title('contig counts'); pp.savefig(p1)        
        # cluster cnt distribution
        orig = [x[0] for x in cntDiff[errType].values()]
        orphan = [x[1] for x in cntDiff[errType].values()]
        diff = [x[1]-x[0] for x in cntDiff[errType].values()]
        p1 = plt.figure(); plt.hist(orig, alpha=.5, label="orig"); plt.hist(orphan, alpha=.5, label="orphan"); plt.legend(loc='upper right'); plt.title('cluster cnt'); pp.savefig(p1)
        p1 = plt.figure(); plt.hist(diff); plt.title('cluster cnt difference (orphan-orig)'); pp.savefig(p1)
        # clusters avg size distribution
        orig = [np.mean(x[0]) for x in sizeDiff[errType].values()]
        orphan = [np.mean(x[1]) for x in sizeDiff[errType].values()]
        diff = [x-y for x, y in zip(orphan,orig)]
        p1 = plt.figure(); plt.hist(orig, alpha=.5, label="orig"); plt.hist(orphan, alpha=.5, label="orphan"); plt.legend(loc='upper right'); plt.title('Avg clusters size'); pp.savefig(p1)
        p1 = plt.figure(); plt.hist(diff); plt.title('Avg clusters size difference (orphan-orig)'); pp.savefig(p1)
        orig = [np.mean(x[0]) for x in purityDiff[errType].values()]
        orphan = [np.mean(x[1]) for x in purityDiff[errType].values()]
        diff = [x-y for x, y in zip(orphan,orig)]
        p1 = plt.figure(); plt.hist(orig, alpha=.5, label="orig"); plt.hist(orphan, alpha=.5, label="orphan"); plt.legend(loc='upper right'); plt.title('Avg clusters entropy'); pp.savefig(p1)
        p1 = plt.figure(); plt.hist(diff); plt.title('Avg clusters entropy difference (orphan-orig)'); pp.savefig(p1)
        pp.close()

        with open(errType + "_geneStats.txt", "wb") as outfile, open(errType + "_allgens.txt", "wb") as genfile:
            outfile.write('%15s %10s %10s %15s %15s\n' % ("Gene", "#Contigs", "#Clus", "AvgClusSize", "AvgEntropy"))
            for key in contigCnt[errType].keys():
                cnt = cntDiff[errType][key][1]-cntDiff[errType][key][0]
                size =  round(np.mean(sizeDiff[errType][key][1])-np.mean(sizeDiff[errType][key][0]), 3)
                purity =  round(np.mean(purityDiff[errType][key][1])-np.mean(purityDiff[errType][key][0]), 3)
                outfile.write('%15s %10s %10s %15s %15s\n' % (key, contigCnt[errType][key], cnt, size, purity))
                genfile.write('{}:{}\t{}\n'.format(key, ','.join(genDiff[errType][key][0]), ','.join(genDiff[errType][key][1])))
        

import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Enter the address to all your files.')
    parser.add_argument('dirAddr', type=str, help='the address to the directory containing all files')
    args = parser.parse_args()
    clust_file1 = args.dirAddr + "mag.flat.clust.orig"
    clust_file2 = args.dirAddr + "mag.flat.clust.orphan"
    clust_files = [clust_file1, clust_file2]
    #vertices, edges = compareActPred(clust_files)
    #import pdb; pdb.set_trace()
    #print(len(vertices))
    gene1 = args.dirAddr + "originalInfo.txt"
    gene2 = args.dirAddr + "orphanInfo.txt"
    cont2gen_file = args.dirAddr + "genContMap.txt"
    tp2fn, tn2fp = catchFalsehood([gene1, gene2])
    print("tp2fn:{}, tn2fp:{}".format(len(tp2fn), len(tn2fp)))
    contigCnt, cntDiff, sizeDiff, purityDiff, genDiff = calcStats(tp2fn, tn2fp, clust_files, cont2gen_file)
    saveAndPlot(contigCnt, cntDiff, sizeDiff, purityDiff, genDiff)
