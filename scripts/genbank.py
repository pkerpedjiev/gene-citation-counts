# get the gene_id -> pubmed mapping
gene2pubmed = sc.textFile("/Users/pkerp/projects/genbank/data/gene2pubmed").map(lambda x: x.split("\t")).filter(lambda x: x[0][0] != '#').map(lambda x: (int(x[1]), int(x[2]))).cache()
print gene2pubmed.take(1)
  
# the gene information
# we're only really interested in the id, the symbol and the description
gene_info = sc.textFile("/Users/pkerp/projects/genbank/data/gene_info").map(lambda x: x.split("\t")).filter(lambda x: x[0][0] != '#').map(lambda x: (int(x[1]), (x[2], x[8], x[0]))).cache()
print gene_info.take(1)
  
# combine the gene information with the reference table
gene_info_pubmed = gene_info.join(gene2pubmed).cache();
print gene_info_pubmed.take(1)

# the dates were obtained by doing entrez queries and stored as
# table of the form "date pmid"
def parse_date_pmid_line(line):
    # extract the date information and the pmid
    # and return it as tuple
    parts = line.split()
    # the year is in %Y/%m/%d format
    year = int(parts[0].split('/')[0])
    pmid = int(parts[1])
    return (pmid, year)


# extract the dates and store them as values where the key is the pmid
pmid_year = sc.textFile('data/pmid_year.ssv').map(parse_date_pmid_line).cache()
print pmid_year.take(1)
# reorder the gene/pubmed information so that the pmid is the key and
# we can join it to the date information
pmid_geneid_name_desc = gene_info_pubmed.map(lambda x: (x[1][-1], (x[0], x[1][0][0], x[1][0][1], x[1][0][2]))).cache()
print pmid_geneid_name_desc.take(1)
# join the gene information with the publciation date information
pmid_year_geneid_name_desc = pmid_year.join(pmid_geneid_name_desc).cache()
pmid_year_geneid_name_desc.take(1)
# reorder the information so that we can count by gene/year pairs
# so we calculate how many times each gene was cited each year
geneid_name_desc_year_1 = pmid_year_geneid_name_desc.map(lambda (pmid, (year, (geneid, name, desc, taxid))): ((geneid, name, desc, taxid, year),1))
geneid_name_desc_year_1.count()
geneid_name_desc_year_counts = geneid_name_desc_year_1.reduceByKey(lambda x,y: x+y)
# get rid of unpopular g    enes (cited < 10 times)
# and rearrange so that they are keyed by year
filtered_geneid_name_desc_year_counts = geneid_name_desc_year_counts.filter(lambda ((geneid, name, desc, taxid, year), count): count > 10)
year_count_geneid_name_desc = filtered_geneid_name_desc_year_counts.map(lambda ((geneid, name, desc, taxid, year), count): (year, (count, geneid, name, desc, taxid)))
year_count_geneid_name_desc.take(1)
# group according to the publication year
# so that now we have a (year, [gene1, gene2, gene3...]) dataset
grouped_year_count_geneid_name_desc = year_count_geneid_name_desc.groupByKey()
grouped_year_count_geneid_name_desc.take(1)
# sort each grouping
sorted_grouped_year_count_geneid_name_desc = grouped_year_count_geneid_name_desc.mapValues(lambda x: sorted(x, reverse=True))
sorted_grouped_year_count_geneid_name_desc.take(1)
# gimme da results!!
scounts = dict(sorted_grouped_year_count_geneid_name_desc.collect())

with open('genes_by_year_all.csv', 'w') as f:
    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("year", "pos", "count", "geneid", "symbol", "name", "taxid"))
    for key in scounts:
        for i,val in enumerate(scounts[key][:20]):
            f.write(str(key) + "\t" + str(i) + "\t" + "\t".join(map(str,val)) + "\n")
            
#######################################################################
from pyspark.sql import *

curated_gene_disease = sc.textFile('data/curated_gene_disease_associations.txt')
header = curated_gene_disease.take(1)[0]

fields = [StructField(field_name, StringType(), True) for field_name in header.split('\t')]
fields[6] = StructField('NumberOfPubmeds', IntegerType(), True)
schema = StructType(fields)
geneDisease = curated_gene_disease.filter(lambda line: line.find('geneId') != 0).map(lambda line: line.split('\t')).map(lambda x: x[:6] + [int(x[6])] + x[7:])
print geneDisease.take(1)
sqlContext = SQLContext(sc)
schemaDiseases = sqlContext.applySchema(geneDisease, schema)
schemaDiseases.registerTempTable('gene_diseases')

sqlContext.sql("SELECT * from gene_diseases where associationType = 'GeneticVariation' order by NumberOfPubmeds desc").take(1)

results = sqlContext.sql("SELECT geneId from gene_diseases where associationType = 'GeneticVariation'")
results.registerTempTable('results')

########################################################################
gene_info = sc.textFile('data/curated_gene_disease_associations.txt')


res1 = sqlContext.sql("select * from results")

# 1. count how many times each gene is cited each month
#   - create a geneid, date, citationcount table
# 2. flatmap each gene id to a set of disease ids
#   - create a geneid, diseaseid table
# 3. translate the counts to a proxy for citations relevant for a particular disease
#    each month
#   - create a diseaseid, citationcount table
# 4. sort the results by date
#   - 

print schema


