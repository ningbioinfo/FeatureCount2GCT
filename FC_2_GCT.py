## A python script to extract information from FeatureCount output file
## and restructre into a .GCT expression table which is required by the 
## GTEX eqtl pipeline

## Libraries
import sys
from os import listdir
from os.path import isfile, join
import os
import numpy as np


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

## Execution


    
Get_file(sys.argv[1])

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
    Get_readcount(file)
    #print(len(ReadCount))
    ReadCountlist.append(ReadCount)
    ReadCount=[] # reset the readcount
    
#print(len(ReadCountlist))

### Output file
for i in range(len(Gene_id)):
        print(Gene_id[i],'\t',Description[i],'\t','\t'.join(str(ReadCountlist[j][i]) for j in range(len(ReadCountlist))), file = open("output.txt","a"))

# headers
header = ["Name","Description"]
for i in datafiles:
    header.append(i.lstrip(sys.argv[1]).rstrip(".counts.txt"))
header.append('\n')
print(header)

insert("output.txt",'\t'.join(i for i in header))
insert("output.txt",'\t'.join([str(len(Gene_id)),str(len(ReadCountlist)),'\n']))
insert("output.txt",'\t'.join(['#1.2','\n']))

# normalised_read count
Get_Normreadcount(ReadCountlist)
    
#print(len(Norm_ReadCountlist))

### Norm_Output file
for i in range(len(Gene_id)):
        print(Gene_id[i],'\t',Description[i],'\t','\t'.join(str(Norm_ReadCountlist[j][i]) for j in range(len(Norm_ReadCountlist))), file = open("Normalised_output.txt","a"))

insert("Normalised_output.txt",'\t'.join(i for i in header))
insert("Normalised_output.txt",'\t'.join([str(len(Gene_id)),str(len(ReadCountlist)),'\n']))
insert("Normalised_output.txt",'\t'.join(['#1.2','\n']))