#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 18:06:28 2020

@author: Alisa
"""
#inputFile = '/Users/Alisa/Downloads/hg38.hybrid.vcf'
#varBed = '/Users/Alisa/Documents/bedLike.txt'

import pybedtools
import math
import nupack_prob
import timeit

def main():
    startTime = timeit.default_timer()
    
    inputFile = input("Please provide file name for variant list: ")
    
    bpList = filter_bps(inputFile) 
    
    varBed = listToBedLike(bpList)
    
    valid_probes = sect(varBed)
    
    #mined_probes = mine(valid_probes)
    
    print('Program complete in %f seconds' % (timeit.default_timer() - startTime))


# Extracts and filters variants from VCF file, returns list of single bp heterozygous variants
def filter_bps(inputFile):
    # Read in file contents.
    with open(inputFile, 'r') as f:
        file_read = [line.strip() for line in f]
    
    
    # Create filtered list using conditional list comprehension.
    filtered_list = [x for x in file_read if x.split('\t')[0][0:4]=="chrX" \
                     and len(x.split('\t')[3])==1 and len(x.split('\t')[4])==1]
    
    bpList = []
    
    # Filter list again for heterozygocity
    for x in range(len(filtered_list) - 1):
        if filtered_list[x].split('\t')[9] == "0|1" or filtered_list[x].split('\t')[9] == "1|0":
            bpList.append(filtered_list[x])

        
    return bpList
    
# Transfers contents of filtered list to bed-like txt file
def listToBedLike(bpList):
    varBed = input("Where would you like to write the filtered variant list: ")
    #'/Users/Alisa/Documents/bedLike.txt'
    
    # Open file where varients will be stored
    w = open(varBed, 'w')

    chrX_start = 2122646
    chrX_end = 2188636  
    
    for x in range(chrX_start, chrX_end, 1):
        s1 = bpList[x].split('\t')[0] #First column = chromosome ID
        d1 = int(bpList[x].split('\t')[1]) #Second column = start coordinate
        d2 = d1 + 1 # Third column = end coordinate
        s2 = bpList[x].split('\t')[3] # fourth column = Reference name
        s3 = bpList[x].split('\t')[4] # Fifth column = Alternate name
        w.write("%s\t%d\t%d\t%s\t%s\n" % (s1, d1, d2, s2, s3))
        
    return varBed

def sect(varBed):
    probes_all = input("Please provide a file name for the probe set to be searched: ")
    #'/Users/Alisa/Documents/chrX_sorted.txt'
    a = pybedtools.BedTool(probes_all)
    b = pybedtools.BedTool(varBed)
    
    valid_probes = input("Where would you like to write the list of probes with intersections: ")
    #'/Users/Alisa/Documents/valid_probes.txt'
    w = open(valid_probes, 'w')

    bt = a.intersect(b, wa=True, u=True)

    for line in bt:
        tabs = str(line).split('\t')
        w.write(tabs[0] + '\t' + tabs[1] + '\t' + tabs[2] + '\t' + tabs[3] + '\n')
        
    return valid_probes


def mine(inputFile):
    outputFile = '/Users/Alisa/Documents/mined.txt'
    
    #Mark start time and open appropriate files
    #startTime = timeit.default_timer()
    r = open(inputFile, 'r')
    w = open(outputFile, 'w')

    startcoords = []
    endcoords = []
    seq = []
    result = 0
    probcalls = 0
    
    
    #Load data from bed-like file into lists
    for line in r:
        t =line.split('\t')
        startcoords.append(int(t[1]))
        endcoords.append(int(t[2]))
        seq.append(t[3])
        
    
    #Establish a moving window of 10000
    for x in range(math.floor(min(startcoords)/10000)*10000, math.ceil(max(endcoords)/10000)*10000, 10000):
        #Establish a temporary dictionary to hold the probes in the current window
        temp = {'sc': [], 'ec': [], 'seq': [], 'probs': []}
        
        #Extrct all probes within the current window
        while startcoords[0] < x + 10000:
            temp['sc'].append(startcoords[0])
            temp['ec'].append(endcoords[0])
            temp['seq'].append(seq[0])
            temp['probs'].append(0)
                
            #Remove values from the front of the lists to avoid iterating over them later    
            startcoords.remove(startcoords[0])
            endcoords.remove(endcoords[0])
            seq.remove(seq[0])
                
            if len(startcoords) == 0:
                break
                
                
        #If more than two probes were extracted calculate probs for all probes in the current window        
        if len(temp['sc']) >2:
            for j in range(len(temp['sc'])):
                sq = temp['seq'][j]
                struct = ''
                
                for i in range(len(sq)):
                    struct = struct + '.' #Assign struct to be string of .'s of len(seq)
                
                prob = nupack_prob.prob(sq, struct, '74.5', '.390', 'dna1998')
                temp['probs'][j] = prob
                probcalls = probcalls + 1
    
        #Filter out probes based on prob value until only two remain
        while len(temp['probs']) > 2:
            ind = temp['probs'].index(min(temp['probs']))
            temp['sc'].remove(temp['sc'][ind])
            temp['ec'].remove(temp['ec'][ind])
            temp['seq'].remove(temp['seq'][ind])
            temp['probs'].remove(temp['probs'][ind])

        
        #Load the information associated with these probes into an output file
        for s in range(len(temp['sc'])):
            if len(temp['sc']) == 0:
                break
            w.write('chrX\t')
            w.write(str(temp['sc'][s]) + '\t')
            w.write(str(temp['ec'][s]) + '\t')
            w.write(str(temp['seq'][s]) + '\t')
            w.write(str(temp['probs'][s]) +'\n')
        
        #Optional: calculate how many probes remain after filtering
        result = result + len(temp['sc'])
    
    return outputFile