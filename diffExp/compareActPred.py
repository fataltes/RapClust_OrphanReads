def catchFalsehood(fn_gene1, fn_gene2):
    tp1 = set()
    cntr = 0
    with open(fn_gene1) as fo:
        for line in fo and cntr <= 10000:
            cluster, tp, fp = line.strip().split()
            for gene in tp.split():
                tp1.add(gene)
            cntr += 1
    cntr = 0
    tp2fp = {}
    with open(fn_gene2) as fo:
        for line in fo and cntr <= 10000:
            cluster, tp, fp = line.strip().split()
            if fp is not None:
                for gene in fp.split():
                    if gene in tp1:
                        tp2fp[gene] = cluster
    return tp2fp


def compareActPred(fn_clust):
    clusts = ({}, {})
    contigs = ({}, {})
    i = 0
    # read both cluster files for two experiments (mag.flat.clust)
    for fn in fn_clust:
        with open(fn) as fp:
            for line in fp:
                tup = line.strip().split() #tup[0]:cluster_id, tup[1]:contig
                if tup[0] not in clusts[i]:
                    clusts[i][tup[0]] = []
                clusts[i][tup[0]] += [tup[1]] # cluster -> list of contigs
                contigs[i][tup[1]]=tup[0] # contig -> cluster
        i += 1
    # build the (bipartite) graph
    vertices = set()
    edges = {}
    for i in [0,1]:
        for clust, conts in clusts[i].iteritems():
            clust_name0 = str(i) + clust
            vertices.add(clust_name0)
            for contig in conts:
                clust_name1 = str(1-i) + contigs[1-i][contig]
                key = (clust_name0, clust_name1)
                if key not in edges:
                    edges[key] = 0
                edges[key] += 1
    return vertices, edges

import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Enter the address to all your files.')
    parser.add_argument('dirAddr', type=str, help='the address to the directory containing all files')
    args = parser.parse_args()
    clust_file1 = args.dirAddr + "mag.flat.clust.orig"
    clust_file2 = args.dirAddr + "mag.flat.clust.orphan"
    clust_files = [clust_file1, clust_file2]
    vertices, edges = compareActPred(clust_files)
    import pdb; pdb.set_trace()
    print(len(vertices))
