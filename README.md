# FeatureCount2GCT

This is a script to convert the output from FeatureCount to GCT format expression tables

## Why we need this

[FeatureCount](http://bioinf.wehi.edu.au/featureCounts/) is quite useful and popular to analyse the gene expression level from RNA-seq datasets, and the [GTEX eQTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) requires expression tables in [GCT format](http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT).

### Dependencies:

1. Python3+
2. numpy

### Manual

```
 $ python FC_2_GCT.py -h
usage: FC_2_GCT.py [-h] [-datadir [DATADIR]] [-out OUT] [-mergedata MERGEDATA]

Script to transform the output from featurecount to .gct format for eQTL
mapping

optional arguments:
  -h, --help            show this help message and exit
  -datadir [DATADIR]    path to the directory that contains the all the output
                        files from featurecount (.counts.txt)
  -out OUT              The output file prefix, by default it is "out", i.e.
                        output will be out.gct and out.normalised.gct
  -mergedata MERGEDATA  The input dir has one count output from featureCount
                        contains all the samples instead of having one txt
                        file for each sample, specify this option to 0 to turn
                        it off

Your ideas are intriguing to me, and I wish to subscribe to your newsletter.
```

### Example
```
 $ ll featureCount/ 
total 4164520
-rw-r--r--  1 ningliu  staff    34M 17 Jan 11:44 ERR009096.counts.txt
-rw-r--r--  1 ningliu  staff   331B 17 Jan 11:45 ERR009096.counts.txt.summary
-rw-r--r--  1 ningliu  staff    34M 17 Jan 11:52 ERR009097.counts.txt
-rw-r--r--  1 ningliu  staff   331B 17 Jan 11:40 ERR009097.counts.txt.summary
-rw-r--r--  1 ningliu  staff    34M 17 Jan 11:44 ERR009099.counts.txt
-rw-r--r--  1 ningliu  staff   333B 17 Jan 11:55 ERR009099.counts.txt.summary
-rw-r--r--  1 ningliu  staff    34M 17 Jan 11:38 ERR009102.counts.txt
-rw-r--r--  1 ningliu  staff   334B 17 Jan 11:54 ERR009102.counts.txt.summary
-rw-r--r--  1 ningliu  staff    34M 17 Jan 11:42 ERR009103.counts.txt
-rw-r--r--  1 ningliu  staff   332B 17 Jan 11:48 ERR009103.counts.txt.summary
-rw-r--r--  1 ningliu  staff    34M 17 Jan 11:55 ERR009104.counts.txt
-rw-r--r--  1 ningliu  staff   334B 17 Jan 11:47 ERR009104.counts.txt.summary
-rw-r--r--  1 ningliu  staff    34M 17 Jan 11:53 ERR009105.counts.txt
....
```
```
python FC_2_GCT.py -datadir featureCount/ -mergedata 0
```

or


```
$ ll Mergedata 
total 29736
-rw-r--r--@ 1 ningliu  staff    15M  5 Apr 20:42 merge_all.counts.txt
```
```
python FC_2_GCT.py -datadir Mergedata/ -mergedata 1 
```
**the /data/ folder is expected to have featureCounts outputs end with ".counts.txt"**

### To do
Give meaningful info to the Description column.
