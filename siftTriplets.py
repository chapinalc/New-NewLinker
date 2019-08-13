# Takes triplets, performs cuts (by ChiSq), and writes the triplets to chunkName + ‘+goodtriplets+’ + saveName + ‘.txt’
# The aforementioned cuts seem quite minor, only to eliminate extreme cases.
# --- -------- - ---------------- --------- --------------
import sys
import os
LINKERDIR=os.environ["TNO_PATH"]
sys.path.insert(0, LINKERDIR)
import time

from LinkerLib import printPercentage
from LinkerLib import Triplet
from LinkerLib import Detection
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
import subprocess

from astropy.io import fits
from astropy.table import Table
'''
input: --a list of triplets with objids as dets (triplets)
       --a dictionary from objid to Detection objects (detDict)
       --chunk number of the chunk (chunkname)
       --save number (savename)
output: --name of the output file (outName)
        --a dictionary of a list of objids of detections (trackIDs)
        --a list of triplets (tripList) 
'''

def writeProcessingFile(triplets, detDict, chunkname, savename):
    trackIDs = {}
    tripList = []
    for x in range(len(triplets)):
        trip = triplets[x]
        if(isinstance(trip, Triplet)):
            trip.trackid = x
            tripdets = trip.dets
            tripList.append(trip)
        else:
            tripdets = trip[1]
 
        trackID = x
        dets = []
        for objid in tripdets:
            if isinstance(objid, int):
                dets.append(detDict[objid])
            elif isinstance(objid, Detection):
                dets.append(objid)
            else:
                print(str(objid) + ' is not an int or Detection')
        trackIDs[trackID] = dets
        
    outName = LL.writeDetToOrb(trackIDs, 
            'siftRequest+' + chunkname + '+' + savename + '.fits')
    
    return outName, trackIDs, tripList 

'''
input: --name of fits file with a list of triplets to be fit
output: --name of file to open for fits file

'''
#TODO turn this into an actual wrapper rather than just calling it with subprocess
def callBulkFit(tripletFile):
    cmd = 'BulkFit -observationFile=' + str(tripletFile) + ' -orbitFile='
    paramName = tripletFile.split('+',1)[-1].split('.')[0]+'.orbit'
    cmd = cmd + paramName
    cmd = cmd +' -bindingFactor=1'
    print('running BulkFit...')
    time0 = time.time()
    os.system(cmd)
    print('done running BulkFit after ' + str(time.time()-time0) + ' seconds')

    return paramName

'''
input: --name of orbitfile
output: writes orbital parameters in aei to orbtifile
'''

def callBulkEle(orbitfile):
    cmd = 'BulkElements -orbitFile=' + str(orbitfile)
    print('running BulkElements...')
    time0 = time.time()
    os.system(cmd)
    print('done running BulkElements after ' + str(time.time()-time0) + ' seconds')
    return


'''
input: --a fits filename with trackIDs and orbital parameters (fitsName)
       --a dictionary from trackIDs to objids (trackDict)
       --a chisq threshold (chiThresh) (optional; defaults to 50)
output: a list of triplets with their objids (tripList)
'''
# Puts orbits with ChiSq under 50 (or another threshold given) into tripList.
def siftTrips(fitsName, trackDict, chiThresh=50):
    print('only keeping orbtis with chisq under ' + str(chiThresh))
    orbits = Table.read(fitsName, format='fits')

    # Records the initial time.
    time1 = time.time()
    mask_chisq = (orbits['CHISQ']<chiThresh) & (orbits['FLAGS']==0)
    orbitIDs = orbits[mask_chisq]['ORBITID']
    # orbitIDs is the list of orbits with ChiSq under 50.
    tripList = []
    for row in orbitIDs:
        trackid = row 
        objids = trackDict[trackid]
        # trackDict is a dictionary of orbits; trackDict[trackid] is a single orbit.
        del trackDict[trackid]
        trip = Triplet(objids)
        # Finds the triplet for the orbit.
        trip.trackid = trackid
        # Creates an ID for each triplet (called trackid).
        tripList.append(trip)
        # Creates tripList, a list of all triplets.

    # Records the final time.
    time2 = time.time()
    '''
    #version 2
    tripList2 = []
    for row in orbits:
        if row['CHISQ'] < chiThresh:
            trackid = '{0:06d}'.format(row['ORBITID'])
            objids = trackDict[trackid]
            tripList2.append(Triplet(objids))
    time3 = time.time()
    '''
    print('done after ' + str(time2-time1) + ' seconds')
    
    #print('v2: ' + str(time3-time2))
   
    return tripList

def removeBadChisq(triplets):
    li = [trip for trip in triplets if trip.chiSq != -1]
    return li

def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets',
                        help='path to pickle file; file has format ' + 
                        'chunk###+SNOBS_S\EASON###_ML0#.pickle')
    args.add_argument('-d', '--detFile', help='path to csv file with detections')
    args.add_argument('-o', '--orbit', action='store_true', help='produce orbitfile only')
    args.add_argument('-x', '--suppress', action='store_true', help='do not run fitter')
    args = args.parse_args()
    
    print('\nopen pickle file ' + args.triplets) 
    with open(args.triplets, 'rb') as pickleFile:
        print('load pickle file')
        triplets = pickle.load(pickleFile)
        print('done loading')
    saveName = args.triplets.split('+')[-1].split('.')[0]
    chunkName = args.triplets.split('/')[-1].split('+')[0]

    if(args.suppress):
        lets = removeBadChisq(triplets)
    else:
        if(args.detFile is None):
            detDict = {}
        else:
            detDict = LL.objidDictionary(args.detFile)

        try:
            triplets = list(set(triplets)) 
        except TypeError:
            triplets = triplets
        outname, trackIDs, lets= writeProcessingFile(triplets, detDict, chunkName, saveName)
        fitsFile = callBulkFit(outname)
        if(args.orbit):
            callBulkEle(fitsFile)
        else:
            lets = siftTrips(fitsFile, trackIDs)
    # Saving and pickling triplets to chunkName + '+goodtriplets+' + saveName + '.txt'
    writeTriplets(lets, chunkName + '+goodtriplets+' + saveName + '.txt', args.orbit)
    pickleTriplets(lets, chunkName + '+goodtriplets+' + saveName + '.pickle', rmunbound=False)

# Runs program.
if __name__=='__main__':
    main()
