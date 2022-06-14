import math
import numpy as np
from pybiomart import Server
import tqdm
import pandas as pd

global species_genes
species_genes = {}

def load_gene_info(species="human", reload=False):
    if (species not in species_genes) or reload:
        dataset_name = 'hsapiens_gene_ensembl'
        if species == "mouse":
            dataset_name = "mmusculus_gene_ensembl"
        server = Server(host='http://www.ensembl.org')
        dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets[dataset_name])
        gene_info = dataset.query(attributes=['ensembl_gene_id', 'hgnc_symbol','chromosome_name', 'start_position', 'end_position', 'strand', 'gene_biotype'], filters={'chromosome_name': ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'], 'biotype': ['lncRNA', 'protein_coding', "miRNA"]}).dropna()
        chrs = {}
        for idx, row in gene_info.iterrows():
            start = row["Gene start (bp)"] if row["Strand"] == '1' else row["Gene end (bp)"]
            strand = 1 if row["Strand"] == '1' else -1
            if row["Chromosome/scaffold name"] in chrs:
                chrs[row["Chromosome/scaffold name"]].append((row["Gene stable ID"], row["HGNC symbol"], start, strand))
            else:
                chrs[row["Chromosome/scaffold name"]] = [(row["Gene stable ID"], row["HGNC symbol"], start, strand)]
        for k in chrs.keys():
            chrs[k].sort(key=lambda y: y[2])
        species_genes[species] = chrs
    return species_genes[species]

def find_gene(pos, arr, start, end, d=float('inf')):
    pivot = math.floor((start+end)/2)
    dist = pos-arr[pivot][2]
    md = 0
    md = dist if abs(dist) < abs(d) else d
    if dist == 0:
        return (dist, pivot)
    elif end-start < 2:
        return (md, start)
    elif dist > 0:
        return find_gene(pos, arr, pivot, end, md)
    else:
        return find_gene(pos, arr, start, pivot, md)

def map_genes(bed, upstream_distance=-2000, downstream_distance=500, species="human", reload=False):
    chrs = load_gene_info(species, reload)
    res = []
    for i in tqdm.tqdm(range(bed.shape[0])):
        if bed.iloc[i,0] in chrs:
            cc = chrs[bed.iloc[i,0]]
            rest = find_gene(bed.iloc[i,1], cc, 0, len(cc))
            res.append((*bed.iloc[i,:], rest[0]*cc[rest[1]][3], cc[rest[1]][1], cc[rest[1]][0], cc[rest[1]][2]))
    df = pd.DataFrame(res, columns =[*bed.columns, 'distance', 'symbol', 'entrez_gene_id', 'gene_start'])
    return df.iloc[np.where((df.loc[:, "distance"] > upstream_distance) & (df.loc[:, "distance"] < downstream_distance))[0], :]
    
