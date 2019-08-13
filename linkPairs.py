import sys

import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import time
from LinkerLib import Detection
from LinkerLib import printPercentage
from LinkerLib import Triplet
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
import LinkerLib as LL

import numpy as np
try:
   import cPickle as pickle
except:
   import pickle
import argparse
import time
import ephem

from multiprocessing import Pool, cpu_count, Manager, Queue


'''
Inputs:
detPairs is a list of tuples of the format (det, links) where det is an objid and links is a list of objids.
Size of each chunk (chunkSize)
A name for the save file (saveName)
Output:
A list of triplets (tripList)
'''
def formTriplets(detPairs, chunkSize, saveName, tripletStart=[-1,-1,-1], chunkStart = 1):
    # detPairs, queue = args
    linkDict = {}
    for pair in detPairs:
        linkDict[pair[0]] = pair[1]

    tripList = []
    time0 = time.time()
    counter = 0
    x = chunkStart
    trackCount = 0
    print('size of each chunk: ' +str(chunkSize)) 
    notFound = True
    for det in detPairs:
        counter += 1
        links = det[1] 
        if counter%1000==0:
            printPercentage(counter, len(detPairs), time.time()-time0)
        if det[0] != tripletStart[0] and tripletStart[0] != -1 and notFound:
            continue
        for link in links:
            if link != tripletStart[1] and tripletStart[1] != -1 and notFound :
                continue
            links2 = linkDict[link]
            for trip in links2:
                if trip != tripletStart[2] and tripletStart[2] != -1 and notFound:
                    continue
                elif(det[0] == tripletStart[0] and link == tripletStart[1] and trip == tripletStart[2]):
                    notFound = False
                trackid = (trackCount%chunkSize)
                triplet = (trackid, [det[0], link, trip])
               
                trackCount += 1
                tripList.append(triplet)
                if chunkSize > 0 and len(tripList) >= chunkSize:
                    pickleTriplets(tripList, 'chunk{0:06d}'.format(x) + '+' + saveName + '.pickle', False)
                    # writeTriplets(tripList, 'chunk{0:06d}'.format(x) +
                    # '+' + saveName + '.txt', False, False)
                    x += 1
                    trackCount = 0
                    tripList = []
    return tripList


def main():
    args = argparse.ArgumentParser()
    args.add_argument('linkedPairs', help='path to pickle file; ' + 'file has format detectionLinks+SNOBS_SEASON###_ML0#.pickle')
    args.add_argument('-n', '--chunkSize', help='size of chunks, leave empty for only one chunk')
    args.add_argument('-c', '--cont', help='last chunk processed')
    args = args.parse_args()
    # Three arguments, of which the last two are optional.

    tripletStart = [-1,-1,-1]
    chunkStart = 1
    if(args.cont):
        print('loading last chunk: ' + args.cont)
        fname = args.cont.split('+')[0]
        if('/' in fname):
            fname = fname.split('/')[-1]
        chunkStart = int(fname[-6:]) + 1
        with open(args.cont, 'rb') as f:
            lastChunk = pickle.load(f)
        lastTrip = lastChunk[-1]
        tripletStart = lastTrip[1]
        print('last triplet: ' + str(tripletStart))
        print('starting from chunk: ' + str(chunkStart))

    # load the pickle file with list of detections and links
    print('load pickle file')
    detPairs = pickle.load(open(args.linkedPairs, 'rb'))
    print('done loading')
    # makes a list of potential triplets

    chunkSize = -1 
    if(args.chunkSize):
        chunkSize = int(args.chunkSize)
    saveName = args.linkedPairs.split('+')[-1].split('.')[0]
    print('forming triplets')
    triplets = formTriplets(detPairs, chunkSize, saveName, tripletStart, chunkStart)

    # tripChunks = splitList(triplets, numChunks, buffered, saveName, printP)
    # triplets = multiProcessPairs(detPairs, Ncpu)
    
    # we save it earlier if argument is buffered

    x = 0
    pickleTriplets(triplets, 'chunk{0:06d}'.format(x) + '+' + saveName + '.pickle', False)
    writeTriplets(triplets, 'chunk{0:06d}'.format(x) + '+' + saveName + '.txt', False, False)
    # Writes and pickles triplets. (The second argument is the save name.)
    # The name is auto-generated and doesn't require an inputted save name.

if __name__ == '__main__':
    main()


