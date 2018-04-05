#!/usr/bin/env python
## A python script to extract information from FeatureCount output file
## and restructre into a .GCT expression table which is required by the
## GTEX eqtl pipeline


## Author: Ning Liu

## version 1.2.0

## Libraries
import argparse
from os import listdir
from os.path import isfile, join
import os
import numpy as np
import math
import csv

## Read Files
datafiles=[]
ReadCount = []
ReadCountlist = []
Norm_ReadCount = []
Norm_ReadCountlist =[]
gene_sum = []

## Get input function
def Get_file(datadir):
    global datafiles
    for f in listdir(datadir):
        if f.endswith(".counts.txt"):
            if isfile(join(datadir,f)):
                datafiles.append(datadir +"/"+ f)
    return datafiles

## Get infomation function
## get read count
def Get_readcount(fc_file):
    global ReadCount
    with open(fc_file, "r") as fcfile:
        next(fcfile)
        next(fcfile)
        if args['mergedata'] == '1':
            for line in fcfile:
                data = line.rstrip('\n').split('\t')
                ReadCount.append([int(i) for i in data[1:]])
        elif args['mergedata'] == '0':
            for line in fcfile:
                data = line.rstrip('\n').split('\t')
                ReadCount.append(int(data[6]))

    return ReadCount

## get normalized read count
## TPM normalization
## readcount/sum*1M

def Get_Normreadcount(ReadCountlist):
    global Norm_ReadCount
    global Norm_ReadCountlist
    global gene_sum

    for sample in ReadCountlist.transpose():
        gene_sum.append(sum(sample))

    for gene in ReadCountlist:
        for col in range(len(gene)):
            norm_rc = round(math.log2(((gene[col]+0.5)/(gene_sum[col]+1.0)*1000000)),3)
            Norm_ReadCount.append(norm_rc)
        Norm_ReadCountlist.append(Norm_ReadCount)
        Norm_ReadCount = []
    Norm_ReadCountlist = np.array(Norm_ReadCountlist)
    return Norm_ReadCountlist


## Add header into the output file

def insert(originalfile,string):
    with open(originalfile,'r') as f:
        with open('newfile.txt','w') as f2:
            f2.write(string)
            f2.write(f.read())
    os.rename('newfile.txt',originalfile)

## Argument

parser = argparse.ArgumentParser(description='Script to transform a featureCounts output (cut -f1,7- removed) to .gct format for eQTL mapping',
                                 epilog='Your ideas are intriguing to me, and I wish to subscribe to your newsletter.')
parser.add_argument('--datadir', nargs='?', help='Path to the directory that contains the all the output files from featurecount (.counts.txt)')
parser.add_argument('--out', default='out', help='Output file prefix, by default it is "out", i.e. output will be out.gct and out.normalised.gct')
parser.add_argument('--mergedata', default='1', help='The input dir has one count input from featureCounts containing all the samples. Set to 0 to turn it off')


## read arguments
args = vars(parser.parse_args())


## Execution
print('START!')
Get_file(args['datadir'])

#print(datafiles)

# get geneid and description

with open(datafiles[0],"r") as fcfile:
    Gene_id = []
    Description = []
    next(fcfile)
    next(fcfile)
    for line in fcfile:
        data = line.rstrip('\n').split('\t')
        Gene_id.append(data[0])
        Description.append('gene')

#print(len(Gene_id))
#print(len(Description))

##
######################### Execution

# raw_readcount
if args['mergedata'] == '1':
    print("Extracting read counts from",datafiles[0])
    Get_readcount(datafiles[0])
    ReadCountlist = np.array(ReadCount)
elif args['mergedata'] == '0':
    for afile in datafiles:
        print("Extracting read counts from",afile)
        Get_readcount(afile)
        ReadCountlist.append(ReadCount)
        ReadCount = []
    ReadCountlist = np.array(ReadCountlist).transpose()


sample_num = len(ReadCountlist[0])
gene_num = len(Gene_id)

print('There are %d samples in your data.'%(sample_num))
print('There are %d genes in your data.' %(gene_num))



out1 = args['out'] + '.gct'
out2 = args['out'] + '.normalised.gct'

### Output file
print("Writing into gct file.")
with open(out1, 'a', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for eachgene in range(gene_num):
        data=[Gene_id[eachgene], Description[eachgene]]
        for eachsample in ReadCountlist[eachgene]:
            data.append(eachsample)
        writer.writerow(data)

# headers
header = ["Name","Description"]
if args['mergedata'] == '1':
    with open(datafiles[0]) as file:
        next(file)
        for line in file:
            name = line.strip().split('\t')[1:]
            for i in name:
                header.append(i)
            break
    header.append('\n')
elif args['mergedata'] == '0':
    for i in datafiles:
        header.append(i.lstrip(args['datadir']).rstrip(".counts.txt").lstrip('/'))
    header.append('\n')


print("The headers are:", header)

if header[int(len(header)/2)].startswith('/'):
    insert(out1,'\t'.join(i.lstrip('/') for i in header))
else:
    insert(out1,'\t'.join(i for i in header))
insert(out1,'\t'.join([str(gene_num),str(sample_num),'\n']))
insert(out1,'\t'.join(['#1.2','\n']))

# normalised_read count
print("Extracting normalised read counts.")
Get_Normreadcount(ReadCountlist)

#print(len(Norm_ReadCountlist))

### Norm_Output file
print("Writing into normalised gct file.")
with open(out2, 'a', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for eachgene in range(gene_num):
        data=[Gene_id[eachgene], Description[eachgene]]
        for eachsample in Norm_ReadCountlist[eachgene]:
            data.append(eachsample)
        writer.writerow(data)

if header[int(len(header)/2)].startswith('/'):
    insert(out2,'\t'.join(i.lstrip('/') for i in header))
else:
    insert(out2,'\t'.join(i for i in header))
insert(out2,'\t'.join([str(gene_num),str(sample_num),'\n']))
insert(out2,'\t'.join(['#1.2','\n']))
