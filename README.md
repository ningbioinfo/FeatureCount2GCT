# FeatureCount2GCT

This is a script to convert the output from FeatureCount to GCT format expression tables

## Why we need this

[FeatureCount](http://bioinf.wehi.edu.au/featureCounts/) is quite useful and popular to analyse the gene expression level from RNA-seq datasets, and the [GTEX eQTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) requires expression tables in [GCT format](http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT).

### Dependencies:

1. Python3+
2. numpy

### Manual

```
$ python FC_2_GCT.py --help
usage: FC_2_GCT.py [-h] [--datadir [DATADIR]] [--out OUT]
                  [--mergedata MERGEDATA] [--norm NORM] [--prior PRIOR]
                  [--annotation [ANNOTATION]]

Script to transform the output from featurecount to .gct format for eQTL
mapping

optional arguments:
 -h, --help            show this help message and exit
 --datadir [DATADIR]   Path to the directory that contains the all the output
                       files from featurecount (.counts.txt)
 --out OUT             Output file prefix, by default it is "out", i.e.
                       output will be out.gct and out.normalised.gct
 --mergedata MERGEDATA
                       The input dir has one count input from featureCounts
                       containing all the samples. Set to 0 to turn it off
 --norm NORM           By default this script uses the tpm normalisation from
                       edgeR, another choice could be rpkm by speicify "rpkm"
                       or cpm by specify "cpm" in this option.
 --prior PRIOR         By default the prior.count set to 1, as it is suggest
                       by the GTEx documentation, but a large prior.count may
                       be valuable to damp down the variability of small
                       count cpm values.
 --annotation [ANNOTATION]
                       the anotation file used in your featureCount, we are
                       assuming the gene_id start with a "ENSG", aka ensembl
                       annotation.

Your ideas are intriguing to me, and I wish to subscribe to your newsletter.
```

### Input Examples

1. featurecount output with one sample a file.

What's in the directory (the summary is not necessary):

```
$ ll -h featureCount
total 512856
-rw-r--r--  1 ningliu  staff    31M  6 Apr 17:27 ERR009102.counts.txt
-rw-r--r--  1 ningliu  staff   348B  6 Apr 17:27 ERR009102.counts.txt.summary
-rw-r--r--  1 ningliu  staff    31M  6 Apr 17:28 ERR009103.counts.txt
-rw-r--r--  1 ningliu  staff   345B  6 Apr 17:28 ERR009103.counts.txt.summary
-rw-r--r--  1 ningliu  staff    31M  6 Apr 17:28 ERR009104.counts.txt
-rw-r--r--  1 ningliu  staff   347B  6 Apr 17:28 ERR009104.counts.txt.summary
-rw-r--r--  1 ningliu  staff    31M  6 Apr 17:28 ERR009105.counts.txt
-rw-r--r--  1 ningliu  staff   347B  6 Apr 17:28 ERR009105.counts.txt.summary
-rw-r--r--  1 ningliu  staff    31M  6 Apr 17:28 ERR009106.counts.txt
-rw-r--r--  1 ningliu  staff   347B  6 Apr 17:28 ERR009106.counts.txt.summary
-rw-r--r--  1 ningliu  staff    31M  6 Apr 17:28 ERR009107.counts.txt
-rw-r--r--  1 ningliu  staff   347B  6 Apr 17:28 ERR009107.counts.txt.summary
-rw-r--r--  1 ningliu  staff    31M  6 Apr 17:28 ERR009108.counts.txt
-rw-r--r--  1 ningliu  staff   347B  6 Apr 17:28 ERR009108.counts.txt.summary
-rw-r--r--  1 ningliu  staff    31M  6 Apr 17:28 ERR009109.counts.txt
-rw-r--r--  1 ningliu  staff   348B  6 Apr 17:28 ERR009109.counts.txt.summary
```

What does it look like in the file:

```
$ head -n 5 featureCount/ERR009102.counts.txt
# Program:featureCounts v1.5.2; Command:"featureCounts" "-t" "exon" "-g" "gene_id" "-a" "Ref/gencode.v27lift37.annotation.gtf" "-o" "Pipeline-output/featureCounts/ERR009102.counts.txt" "Pipeline-output/SplitNCigarReads/ERR009102.markdup.split.bam" "-T" "16"
Geneid    Chr    Start    End    Strand    Length    Pipeline-output/SplitNCigarReads/ERR009102.markdup.split.bam
ENSG00000223972.5_2    chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1    11869;12010;12179;12613;12613;12975;13221;13221;13453    12227;12057;12227;12721;12697;13052;13374;14409;13670    +;+;+;+;+;+;+;+;+    1735    4
ENSG00000227232.5_2    chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1    14404;15005;15796;16607;16858;17233;17606;17915;18268;24738;29534    14501;15038;15947;16765;17055;17368;17742;18061;18366;24891;29570    -;-;-;-;-;-;-;-;-;-;-    1351    264
ENSG00000243485.5_3    chr1;chr1;chr1;chr1;chr1    29554;30267;30564;30976;30976    30039;30667;30667;31097;31109    +;+;+;+;+    1021    0
```

1. featurecount output with all samples a file.

What's in the directory (the summary is not necessary):

```
$ ll -h Mergedata
total 29736
-rw-r--r--@ 1 ningliu  staff    15M  5 Apr 20:42 merge_all.counts.txt
```

What does it look like in the file:

```
$ head -n 5 Mergedata/merge_all.counts.txt
# Program:featureCounts v1.5.1; Command:"featureCounts" "-Q" "10" "-s" "1" "-T" "16" "-a" "/localscratch/Refs/human/hg19_GRCh37d5/Homo_sapiens.GRCh37.87.gtf" "-o" "merged_samples" "PAC006" "PAC007" "PAC008" "PAC009" "PAC010" "PAC011" "PAC012" "PAC013" "PAC014" "PAC015" "PAC016" "PAC017" "PAC018" "PAC020" "PAC021" "PAC022" "PAC023" "PAC024" "PAC025" "PAC026" "PAC027" "PAC029" "PAC030" "PAC031" "PAC032" "PAC033" "PAC034" "PAC035" "PAC036" "PAC037" "PAC038" "PAC039" "PAC040" "PAC041" "PAC042" "PAC043" "PAC044" "PAC045" "PAC046" "PAC047" "PAC048" "PAC049" "PAC050" "PAC051" "PAC052" "PAC053" "PAC054" "PAC055" "PAC056" "PAC057" "PAC058" "PAC059" "PAC060" "PAC062" "PAC063" "PAC064" "PAC065" "PAC069" "PAC070" "PAC071" "PAC072" "PAC074" "PAC075" "PAC076" "PAC077" "PAC078" "PAC083" "PAC084" "PAC086" "PAC087" "PAC088" "PAC091" "PAC093" "PAC097" "PAC098" "PAC099" "PAC100" "PAC102" "PAC103" "PAC105" "PAC107" "PAC108" "PAC109" "PAC111" "PAC114" "PAC117" "PAC118" "PAC120" "PAC122" "PAC124" "PAC127" "PAC129" "PAC131" "PAC134" "PAC139" "PAC140"
Geneid    PAC006    PAC007    PAC008    PAC009    PAC010    PAC011    PAC012    PAC013    PAC014    PAC015    PAC016    PAC017    PAC018    PAC020    PAC021    PAC022    PAC023    PAC024    PAC025    PAC026    PAC027    PAC029    PAC030    PAC031    PAC032    PAC033    PAC034    PAC035    PAC036    PAC037    PAC038    PAC039    PAC040    PAC041    PAC042    PAC043    PAC044    PAC045    PAC046    PAC047    PAC048    PAC049    PAC050    PAC051    PAC052    PAC053    PAC054    PAC055    PAC056    PAC057    PAC058    PAC059    PAC060    PAC062    PAC063    PAC064    PAC065    PAC069    PAC070    PAC071    PAC072    PAC074    PAC075    PAC076    PAC077    PAC078    PAC083    PAC084    PAC086    PAC087    PAC088    PAC091    PAC093    PAC097    PAC098    PAC099    PAC100    PAC102    PAC103    PAC105    PAC107    PAC108    PAC109    PAC111    PAC114    PAC117    PAC118    PAC120    PAC122    PAC124    PAC127    PAC129    PAC131    PAC134    PAC139    PAC140
ENSG00000223972    1    2    0    2    2    3    1    1    3    1    2    1    2    4    0    2    1    3    2    1    1    3    2    2    0
ENSG00000227232    47    55    36    31    37    39    28    57    56    72    96    85    57    65    59    106    53    53    37    69    60    55    64    40    89    32    40    46    54    46    38    51    80    39    38    56    55    46    47    43    43    57    57    49    54    55    55    58    49    60    59    36    43    100    73    79    46    65    51    70    52    67    66    84    66    44    44    73    53    40    48    58    58    69    54    78    70    59    74    59    71    62    63    69    71    60    55    87    60    81    81    79    97    55    48    32
ENSG00000243485    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
```

**the --datadir directory is expected to have featureCounts outputs end with ".counts.txt"**

### To do

- [ ] Give meaningful info to the Description column.
