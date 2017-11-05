import os.path as op
import findspark
import os
findspark.init()

import pyspark
sc = pyspark.SparkContext()

import shortuuid
shortuuid.uuid()

assembly = ''

data_dir = op.expanduser("~/data")
output_dir = op.join(data_dir, assembly)   # where all of the intermediate output will be stored
base_ucsc_dir = op.join(data_dir, 'ucsc-data/{}'.format(assembly))  # where all of the files downloaded from UCSC will be stored

import shutil

# create a directory to store intermediate output files
def get_outfile(table_name):
    outfile = op.join(output_dir, 'genbank-output/{}'.format(table_name))
    if op.exists(outfile):
        shutil.rmtree(outfile)
    return outfile

#wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz
#mv chromInfo.txt.gz ~/data/ucsc-data/hg19/

# Loading the refgene data

base_dir=op.join(op.expanduser("~/data/genbank-data/"), assembly)

import time

gene2pubmed = (sc.textFile(op.join(base_dir, "gene2pubmed"))
                        .filter(lambda x: x[0] != '#')
                        .map(lambda x: x.split('\t'))
                        .map(lambda x: ((int(x[0]), int(x[1])), int(x[2])))
                     )
print(gene2pubmed.take(1))

t1 = time.time()
pubmeds_set = set([x[1] for x in gene2pubmed.collect()])
t2 = time.time()
print("time taken", t2 - t1)

from pyspark import SparkContext, SparkConf

base_dir = op.join(data_dir, 'genbank-data/')
taxid_gene_info = (sc.textFile(op.join(base_dir, 'gene_info'))
                   .filter(lambda x: x[0] != '#')
                   .map(lambda x: x.split('\t'))
                   .map(lambda x: ((int(x[0]), int(x[1])),(x[2], x[8], x[9])))
                   )
taxid_gene_info.take(1)

taxid_gene_refseq_id = (sc.textFile(op.join(base_dir, "gene2refseq"))
                        .filter(lambda x: x[0] != '#')
                        .map(lambda x: x.split('\t'))
                        .map(lambda x: ((int(x[0]), int(x[1])), (x[3].split('.')[0])))
                        )
taxid_gene_refseq_id.take(1)

import time
t1 = time.time()
taxid_gene_info_refseq = taxid_gene_info.join(taxid_gene_refseq_id)

t2 = time.time()
print("time taken", t2 - t1)

len(pubmeds_set)
year_pmid = (sc.textFile(op.join(base_dir, 'recent_pmid_year.ssv'))
                  .map(lambda x: x.split())
                  .map(lambda x: (int(x[0]), int(x[1])))
                  .filter(lambda x: x[1] in pubmeds_set))
year_pmid_collected = year_pmid.collect()

pmid_year = dict([(x[1], x[0]) for x in year_pmid_collected])

print([k for k in list(pmid_year.values())[:10]])

taxid_gene_info_pubmed = (taxid_gene_info.join(gene2pubmed)
                    .filter(lambda x: x[1][1] in pmid_year)
                                 .map(lambda x: ((x[0][0], x[0][1], x[1][1]), x[1])))
taxid_gene_info_pubmed.take(1)        

# get total citation counts over all time for each taxid, gene_id combo
gene_counts = sorted(taxid_gene_info_pubmed.map(lambda x: ((x[0][0], x[0][1]), (x[1][0], 1)))
.reduceByKey(lambda x1, x2: (x1[0], x1[1] + x2[1]))
.collect(), key=lambda x: -x[1][1])

gene_counts[:3]

taxid_gene_info_year = (taxid_gene_info_pubmed.map(lambda x: ((x[0][0], x[0][1], pmid_year[x[1][1]]), x[1]))
                        .map(lambda x: ((x[0][0], x[0][1], pmid_year[x[1][1]]), (x[1][0], 1)))
                                 .reduceByKey(lambda x1, x2: (x1[0], x1[1] + x2[1])))
print(taxid_gene_info_year.take(1))
print(taxid_gene_info_year.count())
print(taxid_gene_info_year.filter(lambda x: x[1][0][2] == 'rRNA').take(1))

genes_per_year = dict(taxid_gene_info_pubmed.map(lambda x: ((pmid_year[x[1][1]], x[0][1]), 1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .map(lambda x: (x[0][0], 1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .collect())

citations_per_year = dict((taxid_gene_info_pubmed.map(lambda x: (pmid_year[x[1][1]], 1))
                      .reduceByKey(lambda x1,x2: x1+x2)
                      .collect()))

citation_genes = sorted(taxid_gene_info_pubmed.map(lambda x: (x[1][1], 1))
 .reduceByKey(lambda x1, x2: x1 + x2)
 .collect(), key=lambda x: -x[1])

print("citation_genes:", citation_genes[:10])

gene_types_citations = (taxid_gene_info_pubmed.map(lambda x: (x[1][0][2], 1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .collect()
)
print("gene_types_citations:", gene_types_citations)

gene_types_counts = (taxid_gene_info_pubmed.map(lambda x: ((x[1][0][2], x[0][1]),1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .map(lambda x: (x[0][0], 1))
 .reduceByKey(lambda x1,x2: x1+x2)
 .collect())

for gtc in gene_types_counts:
    print(gtc[0], gtc[1])
string_values = (taxid_gene_info_year.map(lambda x: "\t".join(map(str, 
                 [x[0][0], x[0][1], x[0][2], x[1][0][0], x[1][0][1], x[1][0][2], x[1][1]
                 ])))
                 .collect())
with open('gene_info_by_year.tsv', 'w') as f:
    for sv in string_values:
        f.write(sv + "\n")
        

