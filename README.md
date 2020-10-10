[![DOI](https://zenodo.org/badge/72544773.svg)](https://zenodo.org/badge/latestdoi/72544773)

# Calculating Gene Citation Counts using NCBI's databases

**This is no longer maintained, try at your own risk. See [this comment](https://github.com/pkerpedjiev/gene-citation-counts/issues/6#issuecomment-692344645) for a link to similar data at the NCBI**

To get an idea of how citation counts have changed over the years, we need
a list of which citations were published when. This script queries Entrez
gene for citation lists for each day.

```
python scripts/pmids_by_date.py --startdate 1990/01/01 --enddate 2017/11/05; 
```

Consolidate the publications from each day into one complete list

For easier processing, we'll collapse the per-day list of files into one
file

```
rm ~/data/genbank-data/hg19/recent_pmid_year.ssv; \
for file in $(find data/pmid_by_date -name "*.ssv"); \
    do cat $file | awk '{split($1,a,"-"); print a[1], $2}' >> ~/data/genbank-data/hg19/recent_pmid_year.ssv; \
done;
```

### Get NCBI's genes to publications table

The NCBI maintains a list of citations that mention each gene in its
database. We'll download it.

```
wget -N -P ~/data/genbank-data/ ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz; \
gunzip ~/data/genbank-data/gene2pubmed.gz; \
mkdir ~/data/genbank-data/hg19/; 
```

And then filter out the human genes:

```
cat ~/data/genbank-data/gene2pubmed | awk '{if ($1 == 9606) print;}' > ~/data/genbank-data/hg19/gene2pubmed
    cat ~/data/genbank-data/gene2pubmed  | python scripts/genes_by_species.py  - | sort -nk 2 > results/genes_by_species.tsv; \
    cp results/genes_by_species.tsv ~/Dropbox/tmp/graphic-idea-3.tsv
```

We'll get a list of chromosome sizes from the UCSC genome browser for plotting purposes:

```
wget -N -P ~/data/ucsc-data/hg19/ http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz
```

We'll also grab a list of genomic positions from UCSC to plot the locations of the genes.

```
wget -N -P ~/data/ucsc-data/hg19/ http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz
```

We need a mapping of gene IDs, which are just numbers, to more meaningful names and descriptions.

```
wget -N -P ~/data/genbank-data/ ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz; \
gunzip ~/data/genbank-data/gene_info.gz; \
mkdir ~/data/genbank-data/hg19/; 
```

Filter for human data for faster downstream processing.

```
cat ~/data/genbank-data/gene_info | awk '{if ($1 == 9606) print;}' > ~/data/genbank-data/hg19/gene_info
```

Get gene to refseq information so that we can use it to get transcript level information. NCBI's gene
database annotates only genes. Actual transcript information is stored in RefSeq.

```
wget -N -P ~/data/genbank-data/ ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz; \
mkdir ~/data/genbank-data/hg19/; \
```

Again, filter for human data.

```
zcat gene2refseq.gz | awk '{if ($1 == 9606) print;}' > ~/data/genbank-data/hg19/gene2refseq
```

All of the aggregations which count the total citations per gene as well as the citations
per per gene per year are calculated using the `summarize_gene_citations.py` script. It could
be done using command line scripts but I started with python so I'll stick with it.

    # aggregate by citation counts and by year
    python scripts/summarize_gene_citations.py; \
    python scripts/summarize_gene_citations_all_species.py; \
    cp gene_info_by_year.tsv ~/Dropbox/tmp/graphic-idea-2.tsv; \



To get a ranking of human transcripts by citation count,
I used the pipeline we created for [higlass's](http://higlass.io) gene anntations:
https://hms-dbmi.github.io/higlass-docs/gene_annotations.html#creating-the-track
This is a good file to cross-reference with the the ones above to check that the
values are correct.

```
    cp ~/data/genbank-data/hg19/geneAnnotationsExonUnions.bed ~/Dropbox/tmp/graphic-idea-1.tsv \

        awk '{if ($4 == "TP53") print}' gene_info_by_year.tsv | sort -nk 7
```

To calculate how many of the papers referencing human genes reference the top 100 genes
we first need to obtain the top 100 genes (transcripts, in this case).

```
sort -nk 5  ~/data/genbank-data/hg19/geneAnnotationsExonUnions.bed | tail -n 100 | awk '{print $8}' > ~/data/genbank-data/hg19/top100_genes
```

Calculate the total number of papers which cite this set of genes.

```
python scripts/pubs_for_genes.py ~/data/genbank-data/gene2pubmed ~/data/genbank-data/hg19/top100_genes
```
