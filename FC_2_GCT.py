#!/usr/bin/env python
## A python script to extract information from FeatureCount output file
## and restructre into a .GCT expression table which is required by the
## GTEX eqtl pipeline


## Author: Ning Liu

## version 1.2.0

## Libraries
import sys
import argparse
from os import listdir
from os.path import isfile, join
import os
import numpy as np
import math
import csv
import re

## Read Files
datafiles=[]
ReadCount = []
ReadCountlist = []
Norm_ReadCount = []
Norm_ReadCountlist =[]
gene_sum = []
gene_dict = {}

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

def Get_Normreadcount(ReadCountlist, gene_dict, gene_id):
    global Norm_ReadCount
    global Norm_ReadCountlist
    global gene_sum

    for sample in ReadCountlist.transpose():
        gene_sum.append(sum(sample))

    ave_libsize = sum(gene_sum)/len(gene_sum)
    priorcount = float(args['prior'])
    n = 0
    for gene in ReadCountlist:
        whichgene = gene_id[n]
        n = n + 1
        for col in range(len(gene)):
            if args['norm'] == 'tpm':
                RPK = gene[col]/(gene_dict[whichgene]/1000) # account for gene length
                prcountscale = float(gene_sum[col])/ave_libsize*priorcount
                libscale = (gene_sum[col]+2*prcountscale)/1000000
                norm_rc = round(math.log2((RPK+prcountscale)/libscale),4) # account for sequencing depths
                Norm_ReadCount.append(norm_rc)

            elif args['norm'] == 'rpkm':
                prcountscale = float(gene_sum[col])/ave_libsize*priorcount
                libscale = (gene_sum[col]+2*prcountscale)/1000000
                RPM = (gene[col]+prcountscale)/libscale # account for sequencing depths
                norm_rc = round(math.log2(RPM/(gene_dict[whichgene]/1000)),4) # account for gene length
                Norm_ReadCount.append(norm_rc)

            elif args['norm'] == 'cpm':
                prcountscale = float(gene_sum[col])/ave_libsize*priorcount
                libscale = (gene_sum[col]+2*prcountscale)/1000000
                norm_rc = round(math.log2((gene[col]+prcountscale)/libscale),4) # account for sequencing depths
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

## get the gene length for normalisation
def get_gene_length(gtf):
    global gene_dict
    with open(gtf, 'r') as gtf:

        next(gtf);next(gtf);next(gtf);next(gtf);next(gtf)
        for line in gtf:
            data = line.strip().split('"')
            gene_id = list(filter(re.compile('^ENSG').search, data))[0]
            gene_length = str(int(data[0].strip().split('\t')[4])-int(data[0].strip().split('\t')[3]))
            if gene_id not in gene_dict:
                gene_dict[gene_id] = gene_length
            else:
                gene_dict[gene_id] = (list(gene_dict[gene_id])+list(gene_length))

        for key,value in gene_dict.items():
            if type(value) is list:
                value = [int(x) for x in value]
                value = round(sum(value)/len(value),4)
                gene_dict[key] = value
    return gene_dict

## Argument

parser = argparse.ArgumentParser(description='Script to transform the output from featurecount to .gct format for eQTL mapping',
                                 epilog='Your ideas are intriguing to me, and I wish to subscribe to your newsletter.')
parser.add_argument('--datadir', nargs='?', help='Path to the directory that contains the all the output files from featurecount (.counts.txt)')
parser.add_argument('--out', default='out', help='Output file prefix, by default it is "out", i.e. output will be out.gct and out.normalised.gct')
parser.add_argument('--mergedata', default='1', help='The input dir has one count input from featureCounts containing all the samples. Set to 0 to turn it off')
parser.add_argument('--norm', default='tpm', help='By default this script uses the tpm normalisation from edgeR, another choice could be rpkm by speicify "rpkm" or cpm by specify "cpm" in this option.')
parser.add_argument('--prior', default='1', help='By default the prior.count set to 1, as it is suggest by the GTEx documentation, but a large prior.count may be valuable to damp down the variability of small count cpm values.')
parser.add_argument('--annotation', nargs='?', help='the anotation file used in your featureCount, we are assuming the gene_id start with a "ENSG", aka ensembl annotation.')

## read arguments
args = vars(parser.parse_args())

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit()
if args['datadir'] is None:
    print('\x1b[0;37;41m'+'ERROR: Please run the script with specifying the --datadir arugments.'+'\x1b[0m')
    parser.print_help(sys.stderr)
    sys.exit()
if args['annotation'] is None:
    print('\x1b[0;37;41m'+'ERROR: Please run the script with specifying the --annotation arguments.'+'\x1b[0m')
    parser.print_help(sys.stderr)
    sys.exit()
## Execution
print('START!')



Get_file(args['datadir'])

# check the mode
if len(datafiles)>1 and args['mergedata'] == '1':
    print('\x1b[0;37;41m'+'ERROR: Please turn off the mergedata mode when inputing multiple count.txt.'+'\x1b[0m')
    parser.print_help(sys.stderr)
    sys.exit()
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



if args['norm'] != 'cpm':
    print("Calculating gene length from the gtf file.")
    get_gene_length(args['annotation'])
    print('Finishing calculation of gene lengths.')


print("Extracting normalised read counts.")

print("You choosed %s to be your normalisation method." %(args['norm']))

Get_Normreadcount(ReadCountlist, gene_dict, Gene_id)

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
