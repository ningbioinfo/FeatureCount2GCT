#!/usr/bin/env python
## A python script to extract information from FeatureCount output file
## and restructre into a .GCT expression table which is required by the
## GTEX eqtl pipeline

## Author: Ning Liu

## Libraries
import sys
import argparse
from os import listdir
from os.path import isfile, join
import os
import numpy as np
import csv


## Read Files
datafiles=[]
ReadCount = []
ReadCountlist = []
Norm_ReadCount = []
Norm_ReadCountlist =[]



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
    with open(fc_file, "r") as file:
        global ReadCount
        next(file)
        next(file)
        for line in file:
                data = line.rstrip('\n').split('\t')
                ReadCount.append(int(data[6]))
    return ReadCount

## get normalized read count
## TPM normalization
## readcount/sum*1M

def Get_Normreadcount(rc_list):
        global Norm_ReadCount
        global Norm_ReadCountlist
        for i in range(len(rc_list)):
            s = sum(rc_list[i])
            for j in rc_list[i]:
                norm_rc = np.round((j/float(s))*1000000, decimals=3)
                Norm_ReadCount.append(norm_rc)
            Norm_ReadCountlist.append(Norm_ReadCount)
            Norm_ReadCount=[]
        return Norm_ReadCount
        return Norm_ReadCountlist


## Add header into the output file

def insert(originalfile,string):
    with open(originalfile,'r') as f:
        with open('newfile.txt','w') as f2:
            f2.write(string)
            f2.write(f.read())
    os.rename('newfile.txt',originalfile)

## Argument

parser = argparse.ArgumentParser(description='Script to transform the output from featurecount to .gct format for eQTL mapping',
                                 epilog='Your ideas are intriguing to me, and I wish to subscribe to your newsletter.')
parser.add_argument('-datadir', nargs='?', help='path to the directory that contains the all the output files from featurecount (.counts.txt)')
parser.add_argument('-out', default='out', help='The output file prefix, by default it is "out", i.e. output will be out.gct and out.normalised.gct')



## read arguments
args = vars(parser.parse_args())



## Execution


Get_file(args['datadir'])

# get geneid and description

with open(datafiles[0],"r") as file:
    Gene_id = []
    Description = []
    next(file)
    next(file)
    for line in file:
         data = line.rstrip('\n').split('\t')
         Gene_id.append(data[0])
         position = data[1] + ":" + data[2] + "-" + data[3]
         Description.append(position)


##
######################### Execution

# raw_readcount
for file in datafiles:
    print("Extracting read counts from",file)
    Get_readcount(file)
    #print(len(ReadCount))
    ReadCountlist.append(ReadCount)
    ReadCount=[] # reset the readcount

print(len(ReadCountlist[0]))
print(len(Gene_id))
print(len(Description))

out1 = args['out'] + '.gct'
out2 = args['out'] + '.normalised.gct'

### Output file
print("Writing into gct file.")
with open(out1, 'a', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for i in range(len(Gene_id)):
        data=[Gene_id[i], Description[i]]
        for j in range(len(ReadCountlist)):
            data.append(str(ReadCountlist[j][i]))
        writer.writerow(data)

# headers
header = ["Name","Description"]
for i in datafiles:
    header.append(i.lstrip(sys.argv[1]).rstrip(".counts.txt"))
header.append('\n')
print("The headers are:", header)

if header[int(len(header)/2)].startswith('/'):
    insert(out1,'\t'.join(i.lstrip('/') for i in header))
else:
    insert(out1,'\t'.join(i for i in header))
insert(out1,'\t'.join([str(len(Gene_id)),str(len(ReadCountlist)),'\n']))
insert(out1,'\t'.join(['#1.2','\n']))

# normalised_read count
print("Extracting normalised read counts.")
Get_Normreadcount(ReadCountlist)

#print(len(Norm_ReadCountlist))

### Norm_Output file
print("Writing into normalised gct file.")
with open(out2, 'a', newline='') as norm_csvfile:
    writer = csv.writer(norm_csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for i in range(len(Gene_id)):
        data=[Gene_id[i], Description[i]]
        for j in range(len(Norm_ReadCountlist)):
            data.append(str(Norm_ReadCountlist[j][i]))
        writer.writerow(data)

if header[int(len(header)/2)].startswith('/'):
    insert(out2,'\t'.join(i.lstrip('/') for i in header))
else:
    insert(out2,'\t'.join(i for i in header))
insert(out2,'\t'.join([str(len(Gene_id)),str(len(ReadCountlist)),'\n']))
insert(out2,'\t'.join(['#1.2','\n']))
