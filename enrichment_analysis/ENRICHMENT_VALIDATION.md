# GO Enrichment Input Validation

## Dataset
GSE190688 — HGSOC versus FTE

## Significant gene lists

| Category | Significant Ensembl genes | Unmapped genes | Unique Entrez IDs |
|---|---:|---:|---:|
| HGSOC-upregulated | 5,573 | 1,461 | 4,107 |
| FTE-upregulated | 4,294 | 618 | 3,672 |

Unmapped Ensembl IDs were excluded from the corresponding Entrez gene lists. Duplicate Entrez identifiers were counted once.

## Background validation

- Genes in DESeq2 results: 24,746
- Genes with non-NA adjusted p-values: 24,259
- Unique mapped Entrez IDs from all DESeq2 genes: 19,338
- IDs in the enrichment background file: 19,338
- Expected IDs missing from background: 0
- Unexpected IDs in background: 0

The background file exactly matches the unique Entrez IDs mapped from all genes in the DESeq2 results, including genes with NA adjusted p-values.

## Interpretation

The difference between the significant Ensembl gene counts and the Entrez gene-list sizes is explained by unmapped identifiers and duplicate Entrez mappings.

These checks establish consistency of the enrichment input files. They do not independently validate the enrichment algorithm, parameters, or biological conclusions.
