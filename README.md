# Supporting materials for manuscript: "Overcoming circularity bias reveals a hidden linear DNA and RNA virosphere"
This repository tracks the notes about the termini-conservation approach applied to recover new linear complete viral genomes.
### Data availability
See: https://www.meta-virome.org/
### Code availability

To identify "Linear Conserved" vOTUs the user is in need of 2 tabular file:
1) a tabular file with the following UViGs metatadata [uvig, votu, lenght, uvig_topology, checkv_quality, viral_confidence, taxon_oid (sample identifier)].
2) a text file (i.e. representative.txt), contaning the uvig name of each votu representative.

To recover "New Linear Complete" vOTUs the user is in need of a third tabular file:  
3) a tabular file termini_25bp.tsv with the following data [uvig,start_25bp,end_25bp]

OPTION 1 (via seqkit):
```
seqkit fx2tab -n -s input.fasta | \
awk 'BEGIN{OFS="\t"; print "uvig","start_25bp","end_25bp"} 
     {print $1, substr($2,1,25), substr($2,length($2)-24,25)}' \
> termini_25bp.tsv
```

OPTION 2 (via Biopython):
```
from Bio import SeqIO
import pandas as pd

records = []
for record in SeqIO.parse("input.fasta", "fasta"):
    seq = str(record.seq)
    if len(seq) >= 25:
        records.append({
            "uvig": record.id,
            "start_25bp": seq[:25],
            "end_25bp": seq[-25:]
        })
    else:
        records.append({
            "uvig": record.id,
            "start_25bp": seq,
            "end_25bp": seq
        })

df = pd.DataFrame(records)
df.to_csv("termini_25bp.tsv", sep="\t", index=False)
```

For any other code/analysis inquiries, please open a github issue. Note: most of these scripts were written for Python 3.  
To retrieve metadata necessary to run the analysis, please follow Methods described in metaVR database [https://doi.org/10.1093/nar/gkaf1283]

If this code is useful, please cite: TBD
