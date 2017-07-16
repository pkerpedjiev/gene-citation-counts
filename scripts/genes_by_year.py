import pandas as pd
sorted_gene_info = pd.read_csv('gene_info_by_year.tsv', names=['taxid', 'geneid', 'year', 'symbol', 'name', 'type', 'citations'],
                                delimiter='\t')
min_year = min(sorted_gene_info['year']) 
#min_year = 2000
max_year = max(sorted_gene_info['year'])

print(min_year, max_year)
print(sorted_gene_info.head())

import collections as col

all_genes = set()
all_year_genes = col.defaultdict(lambda: col.defaultdict(int))
all_decade_genes = col.defaultdict(lambda: col.defaultdict(int))
all_genes_citations = col.defaultdict(int)

for year in range(min_year, max_year+1):
    for ix, row in  (sorted_gene_info[sorted_gene_info['year'] == year]
                    .sort_values('citations', ascending=[False]).iterrows()):
        
        all_year_genes[year][row['symbol']] = row['citations']
        all_genes.add(row['symbol'])
        all_decade_genes[year / 10][row['symbol']] += row['citations']
        all_genes_citations[row['symbol']] += row['citations']
        
        
with open('all_year_genes.tsv', 'w') as f:
    f.write("gene\t{}\n".format("\t".join(map(str, range(min_year, max_year+1)))))

    for gene in all_genes:
        str_counts = []    
        for year in range(min_year, max_year+1):
            if gene in all_year_genes[year]:
                str_counts += [str(all_year_genes[year][gene])]
            else:
                str_counts += ['']
            #print year / 10, top_half_decades[year / 10]
        out_str = "{gene}\t{counts}".format(gene=gene, counts="\t".join(str_counts))
        f.write(out_str + "\n")
        #print out_str
