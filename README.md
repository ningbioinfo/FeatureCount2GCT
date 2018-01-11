# FeatureCount2GCT

This is a script to convert the output from FeatureCount to GCT format expression tables

## Why we need this

[FeatureCount](http://bioinf.wehi.edu.au/featureCounts/) is quite useful and popular to analyse the gene expression level from RNA-seq datasets, and the [GTEX eQTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) requires expression tables in [GCT format](http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT).

### Dependencies:

1. Python3
2. numpy

### How to run

```
python3 FC_2_GCT.py /path/to/your/data/
```

**the /data/ folder is expected to have featureCounts outputs end with ".counts.txt"**
