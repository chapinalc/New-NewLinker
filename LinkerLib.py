# -     - - --   --   ----
import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)

import numpy as np
import pandas as pd
from collections import namedtuple

from datetime import datetime
import time
try:
   import cPickle as pickle
except:
   import pickle

import ephem
from astropy.time import Time
from astropy.table import Table
from astropy.io import fits
from astropy.io.fits.hdu.hdulist import HDUList

from orbfitScript import degToHour, degToDMS
import GammaTPlotwStatTNOExFaster as gts
from Orbit import Orbit

Det = namedtuple('Det', 'ra dec mjd mag objid ' +
                        'expnum ccd band fakeid posErr lookAhead')
# Returns the bounds of ra and dec as a string.
class Region:
    # fields include the range for RA / Dec and list of detections
    def __init__(self, raRng, decRng):
        self.raLo = raRng[0]
        self.raHi = raRng[1]
        self.decLo = decRng[0]
        self.decHi = decRng[1]
        self.detections = []
        
    def add(self, tno):
        self.detections.append(tno)
    
    def toStr(self):
        minbounds = '(minRa: ' + str(self.raLo) + ', minDec: ' + str(self.decLo) + ')'
        maxbounds = '(maxRa: ' + str(self.raHi) + ', maxDec: ' + str(self.decHi) + ')'
        return minbounds + '  ' + maxbounds 


# Returns several variables and strings pertaining to detection information.
#ra is between -180 and 180
class Detection:
    def __init__(self, ra, dec, mjd, flux=0, objid=0, 
                expnum=0, ccd=0, band=0, lookAhead=0, fakeid=-1, magerr=0):
        self.ra = ra
        if self.ra > 180:
            self.ra = ra -360
        self.erra = 0
        self.errb = 0
        self.pa = 0
        self.dec = dec
        self.mjd = mjd
        self.mag = -2.5*np.log10(flux) + 31.4
        self.flux = flux
        self.magerr = magerr
        self.objid = objid
        self.expnum = expnum
        self.ccd = ccd
        self.band = band
        self.fakeid = fakeid
        self.posErr = 0
        self.lookAhead = lookAhead
        #describes a cone region that would contain possible pair links
        #angle in radians where 0 is along the x-axis, where change in dec is 0
        #angle1 < angle2
        # is an array of angles and radii for each nite with index corresponding
        # to the nite ahead
        self.angle1, self.angle2, self.radius = gts.calcCone(lookAhead, 
                        self.ra, self.dec, self.mjd)
        #list of linked detections (intially empty)
        self.linkedList = []

    # Computes and returns posErr (which is based upon magnitude).
    def getPosErr(self):
        err_sec = 0
        if(self.mag <= 21):
            err_sec = 0.1
        else:
            err_sec = 0.1 + 0.1*(self.mag-21)/3
        self.posErr = err_sec/3600
        return self.posErr

    # Calculates magerr
    def setMagErr(self, fluxErr):
        self.magerr = abs(-2.5/np.log(10)*fluxErr/self.flux)
    
    # Returns whether two detections could be the same object depending on the magErrs
    def magFit(self, other):
        if(not (self.band == other.band)):
            return True
        selfMin = self.mag-2*self.magerr
        selfMax = self.mag+2*self.magerr
        otherMin = other.mag-2*other.magerr
        otherMax = other.mag+2*other.magerr

        overLap = not (selfMin > otherMax or selfMax < otherMin)
        return overLap

    #so nodes can have class of Detection
    def __hash__(self):
        return hash(self.objid)
    def __eq__(self, other):
        if(isinstance(self, int) and isinstance(other, int)):
            return self == other
        #print(self)
        #print(other)
        if(self.objid == other.objid):
            return True
        if(not self.expnum != other.expnum):
            return False
        if(abs(self.ra - other.ra) < 0.000001 and abs(self.dec - other.dec) < 0.000001):
            return True
        return False
    def __str__(self):
        return self.toStr()
    # Returns a string consisting of ids, coord, exp, light, lookAhead, and err.
    def toStr(self):
        coord = '(RA: ' + str(self.ra) + ', DEC: ' + str(self.dec) + ') mjd: ' + str(self.mjd)
        ids = ('fakeid: ' + str(self.fakeid) + ', objid: ' + str(self.objid))
        exp = ('expnum: ' + str(self.expnum) + ', ccd: ' + str(self.ccd))
        light = 'mag: ' + str(self.mag) + ', band: ' + str(self.band) 
        err = 'err: ' + str(self.posErr*3600)
        return (ids + ', ' + coord + ', ' + exp + ', ' + light + ' la: ' + str(self.lookAhead) + ', ' + err)
    # Returns a string  with the Juilan Date, the time taken to reach ra in the sky, and dec as dd:mm:ss.ss
    def toDat(self):
        # formerly used str(ephem.hours(self.ra*np.pi/180))
        # str(ephem.degrees(self.dec*np.pi/180))
        s = (str(self.mjd+2400000.5) + ' ' + degToHour(self.ra) +
                ' ' + degToDMS(self.dec) + ' 0.1 807')
        return s
        
    
    #checks if another detection is inside the cone
    def withinCone(self, det):
        #time0 = time.time()
        nitesAhead = int(det.mjd - self.mjd)
        #only looking ahead a certain number of nites
        if(nitesAhead + 1 >= self.lookAhead or nitesAhead < 0):
            return False
        #time1 = time.time()
        ra = det.ra
        # make sure ra is between -180 and 180
        if ra > 180:
            ra = ra - 360
        deltaRa = ra-self.ra
        deltaDec = det.dec-self.dec
        #time2 = time.time()
        dist = np.sqrt(deltaRa**2 + deltaDec**2)
        angle = np.arctan2(deltaDec, deltaRa)
        #vector angle is between 0 and 360
        if(angle < 0):
            angle += 2*np.pi
        #checks if the cone overlaps with the 0 angle
        rev = self.angle1[nitesAhead] > self.angle2[nitesAhead]
        #time3 = time.time()
        #a one night interval
        angle1 = max(self.angle1[nitesAhead], self.angle1[nitesAhead + 1])
        angle2 = min(self.angle2[nitesAhead], self.angle2[nitesAhead + 1])
        radius = max(self.radius[nitesAhead], self.radius[nitesAhead + 1])
        #time4 = time.time()
        #print('nite Time: ' + str(time1-time0))
        #print('coord Time: ' + str(time2-time1))
        #print('math Time: ' + str(time3-time2))
        #print('angleTime: ' + str(time4-time3))
        if(rev):
            if(angle < angle1 and angle > angle2 and dist < radius):
#                print('\nself: ' + self.toStr() + '\n  them: '  + det.toStr())
#                print('dist: ' + str(dist))
                return True
        
            else:
                return False
        #the cone overlaps with the 0 angle
        else:
            if((angle < angle1 or angle > angle2) and dist < radius):
                return True
            else:
                return False

    #returns the min/max ra/dec of the cone
    def bounds(self): 
        vertex1 = [self.ra, self.dec]
        vertex2 = [self.ra + np.amax(self.radius)*np.cos(np.amax(self.angle1)),
                    self.dec + np.amax(self.radius)*np.sin(np.amax(self.angle1))]
        vertex3 = [self.ra + np.amax(self.radius)*np.cos(np.amin(self.angle2)), 
                     self.dec + np.amax(self.radius)*np.sin(np.amin(self.angle2))] 
        minRa = min(vertex1[0],vertex2[0],vertex3[0])
        maxRa = max(vertex1[0],vertex2[0],vertex3[0]) 
        minDec = min(vertex1[1],vertex2[1],vertex3[1]) 
        maxDec = max(vertex1[1],vertex2[1],vertex3[1])
        return minRa, maxRa, minDec, maxDec  

#Fake class for linkmap efficiency
class Fake(): 
    def __init__(self, fakeid, objid, mjd, exp, ccd, ra, dec, flux, band):
        self.fakeid = fakeid 
        self.objid = objid 
        self.mjd = mjd  
        self.expnum = exp   
        self.mag = -2.5*np.log10(flux)+31.4 
        self.band = band  
        self.ccd = ccd 
        self.ra = ra  
        self.dec = dec     
        self.links = []
        self.alllinks = []

    def toStr(self): 
        fakeid = "fakeid: " + str(self.fakeid) + "  " 
        objid = "objid: " + str(self.objid) + "  " 
        mjd = "mjd: " + str(self.mjd) + "  "  
        exp = "expnum: " + str(self.expnum) + " "   
        ccd = "ccd: " + str(self.ccd) + " " 
        ra = "RA: " + str(self.ra) + ' '  
        dec = "Dec: " + str(self.dec) + ' '  
        mag = "Mag: " + str(self.mag) + ' ' 
        band = "Band: " + str(self.band)  
        return (fakeid + objid + mjd + exp + ccd + ra + dec + mag + band)
    def toDat(self):
        s = (str(self.mjd+2400000.5) + ' ' + str(ephem.hours(self.ra*np.pi/180)) +
                ' ' + str(ephem.degrees(self.dec*np.pi/180)) + ' 0.1 807')
        return s
 

#class that stores triplets
class Triplet:
    #input is array of detections
    def __init__(self, dets):
        
        if(not isinstance(dets, list)):
            raise TypeError("Input 'dets' is not a list")
        if(not len(dets) == 0):
            if(not isinstance(dets[0], Detection) and not isinstance(dets[0], int)):
                raise TypeError("Input 'dets' is not a list of Detections or ints")
        #is a list
        self.dets = dets
        #reduce runtime
        self.elements = 0
        self.errs = 0
        self.orbit = 0
        self.trackid = -1
        self.cands = []
        self.chiSq = -1
        self.abg = 0
        self.cov = 0
        self.aeiCov = 0 
 #returns whether two triplets share M detections
    def shareM(self, other, M):
        count = 0
        for det1 in self.dets:
            for det2 in other.dets:
                if(det1==det2):
                    count += 1
        return count >= M
    #removes duplicate exposures from detections
    def removeDupExp(self):
        expDict = {} 
        for det in self.dets:
            expDict[det.expnum] = det
        self.dets = expDict.values()
            
    #returns whether or not this triplet is the proper subset of the other
    def isSubset(self, other, proper=False):
        if(len(self.dets)==len(other.dets) and proper):
            return False
        for x in self.dets:
            if not x in other.dets:
                return False
        return True
 

    #returns whether all detections are fakes
    def isFake(self):
        for det in self.dets:
            if(det.fakeid==0):
                return False
        return True

    # Returns true if less than half of the detections are fakes.
    def isReal(self):
        fakeCounter = 0
        for det in self.dets:
            if(det.fakeid!=0):
                fakeCounter += 1
        return fakeCounter < len(self.dets)/2


    #returns a list of detections that contains detections from both triplets
    def merge(self, other):
        temp = []
        if(isinstance(other, Triplet)):
            temp = self.dets + other.dets
        elif(isinstance(other, list)):             
            temp = self.dets + other
        elif(isinstance(other, tuple)):
            temp = self.dets[:]
            for i in other:
                temp.append(i)
        else:
            print('Error: not Triplet or list or tuple')
        return list(set(temp))
    #adds a detection if not already present
    def addDetection(self, det):
        temp = self.dets + [det]
        self.dets = list(set(temp))
        self.orbit = 0
        self.chiSq = -1
        self.elements=0
        self.errs=0
    #converts triplet to a dat format
    def toDat(self, outfile):
        s = self.dets[0].toDat()
        for x in range(len(self.dets)-1):
            s = s + '\n'
            s = s + self.dets[x+1].toDat()    
        print('writing to file: ' + outfile)
        with open(outfile, 'w+') as f:
            f.write(s)
        return s

    # Calculates an orbit if one does not already exist.
    # Returns values for elements and errs.
    def calcOrbit(self):
        if(self.orbit == 0):
            self.orbit = self.setOrbit()
        #time4 = time.time()
        #print('time1: ' + str(time1-time0))
        #print('time2: ' + str(time2-time1))
       # print('time3: ' + str(time3-time2))
        #print('time4: ' + str(time4-time3))
        #orbit.get_elements()
        self.chiSq = self.orbit.chisq
        self.elements, self.errs = self.orbit.get_elements()
        return self.elements, self.errs
    
    # determines whether there is a subset of num detections that are within daysApart of each other
    def discoverableTriplet(self, daysApart, num=3):
        self.sortByMjd()
        seqLen = 0
        prevDate = 0
        for det in self.dets:
            if det.mjd-prevDate < daysApart:
                seqLen += 1
            else:
                seqLen = 1
            prevDate = det.mjd
            if seqLen >= num:
                return True
        return False

    # Returns an object orbit.
    def setOrbit(self):
        mag = min([x.mag for x in self.dets])
        #time0 = time.time()
        errs=0.1
        if mag <= 21:
            errs = 0.15
        else:
            errs = 0.15 + (mag -21.0) / 40.0
        #time1 = time.time()
        ralist = [ephem.hours(np.deg2rad(det.ra)) for det in self.dets]
        #print('ra:' + str(time.time()-time1))
        #time2 = time.time()
        declist = [ephem.degrees(np.deg2rad(det.dec)) for det in self.dets]
        #print('dec:' + str(time.time()-time2))
        datelist = [ephem.date((Time(det.mjd,format='mjd')).datetime) for det in self.dets]
        #time3 = time.time()
        #print('date:' + str(time3-time2))
        #time2 = time.time()i
        orbit = Orbit(dates=datelist, ra=ralist, dec=declist,
            obscode=np.ones(len(self.dets), dtype=int)*807, err=errs)
        self.orbit = orbit
        self.chisq = orbit.chisq
        #time3 = time.time()
        self.elements, self.errs = orbit.get_elements()
        #time4 = time.time()
        #print('time0:' + str(time1-time0))
        #print('dateTime:' + str(time2-time1))
        #print('orbitTime:' + str(time3-time2))
        #print('getTime:' + str(time4-time3))
        return orbit

    # Sets an orbit if one doesn’t already exist.
	# Calculates and returns dist and err.
    def getDistance(self):
        if self.orbit == 0:
            self.setOrbit()
        dist, err = self.orbit.barycentric_distance()
        return dist, err

    #get the covarance matrix 
    def getCovar(self):
        if(self.orbit==0):
            self.setOrbit()
        #_, _, covar = self.orbit.get_elements_abg()
        return self.orbit.covar_abg

    #gets the chi squared of the orbit
    def getChiSq(self):
        if(len(self.dets) == 0):
            return 999999999
        if(self.chiSq == -1):
            self.setOrbit()
            self.chiSq = self.orbit.chisq 
        #return self.orbit.chisq
        return self.chiSq

    '''
    input: date in the format MJD
    output: a tuple containing (ra, dec) and a scalar error in arc seconds
    '''
    # predicts the position of the triplet at a future date
    def predictPos(self, date):
        date = ephem.date((Time(date, format='mjd')).datetime)
        if(self.orbit==0):
            self.setOrbit()
        orbit = self.orbit
        predPos = orbit.predict_pos(date)
        ra = ephem.hours(predPos['ra'])
        dec = ephem.degrees(predPos['dec'])
        erra = predPos['err']['a']
        errb = predPos['err']['b']
        pa = predPos['err']['PA']
        ra = np.degrees(ra)
        if(ra > 180):
            ra -= 360
        dec = np.degrees(dec)
        return (ra, dec), erra, errb, pa

    # return number of detections without repeats
    def realLength(self):
        self.sortByMjd()
        counter = 0
        prev = 0
        for det in self.dets:
            if(det.mjd-prev > 0.1):
                counter += 1
            prev = det.mjd
        return counter

    #checks if the objects all have the same fakeid    
    def sameFake(self, thresh=0):
        fakeids  = [x.fakeid for x in self.dets]
        fakeid = max(set(fakeids), key=fakeids.count)
        counter = 0
        if(fakeid == 0):
            return 0
        for det in self.dets: 
            if(det.fakeid != fakeid):
                counter += 1
        if(counter > thresh):
            return 0
        return fakeid 

    # Sorts detections by MJD.
    def sortByMjd(self):
        self.dets.sort(key = lambda x: x.mjd)
        self.chisq = -1

    # Returns true iff each triplet is apart by at least days.
    def apartByXDays(triplet, days):
        mjd = sorted([x.mjd for x in triplet.dets])
        for i in range(1, len(mjd)):
            if mjd[i] - mjd[i - 1] < days:
                return False
        return True
    
    #checks if a triplet is small enough to be a p9 candidate
    def p9Like(self, chiThresh = 10):
        apartDays = 7
        apartDist = 90 # in arcseconds
        if(self.getChiSq > 10):
            return False
        if(self.apartByXDays(apartDays)):
            self.sortByMjd()            
            for x in range(len(self.dets)-1):
                dist = calcDist(self.dets[x], self.dets[x+1])
                if(dist*60*60 > apartDist):
                    return False

            return True
        else:
            return False
    
    # returns whether a magnitude makes sense with this triplet
    def magFit(self, detection, thresh = 0.7):
        for det in self.dets:
            if(detection.band == det.band):
                if(abs(det.mag-detection.mag) > thresh):
                    return False
        return True

    # Returns true iff the magnitude fits.
    def magCheck(self):
        for det in self.dets:
            if(not magFit(self, det)):
                return False
        return True
    
    # returns the largest magnitude difference regardless of band
    def magDiff(self):
        mag = [det.mag for det in self.dets]
        return max(mag) - min(mag)

    # Returns true if there exists a consecutive set of detections separated by between .1 and lookAhead days at least as long as n (by default a triplet).
    def containsNlet(self, lookAhead=50, n=3):
        self.sortByMjd()
        largest = 0
        curr = 0
        for x in range(len(self.dets)-1):
            mjdDiff = self.dets[x+1].mjd - self.dets[x].mjd
            if(mjdDiff < lookAhead and mjdDiff > 0.1):
                curr += 1
                if(curr > largest):
                    largest = curr
            elif(mjdDiff < lookAhead):
                pass
            else:
                curr = 0
        return largest >= n

    # returns a copy of itself, possibly with more detections
    def makeCopy(self, additional=[]):
        dets = [Detection(x.ra, x.dec, x.mjd, x.flux, x.objid, x.expnum,
                        x.ccd, x.band, x.lookAhead, x.fakeid, x.magerr) for x in self.dets]
        dets += additional
        return Triplet(dets)


    # returns false if magnitudes for a given band doesn't make sense
    def magFilter(self, thresh=1.0):
        gMax = 0
        gMin = 30
        rMax = 0
        rMin = 30
        iMax = 0
        iMin = 30
        zMax = 0
        zMin = 30
        for det in self.dets:
            if(det.band == 'g'):
                if(det.mag < gMin):
                    gMin = det.mag
                if(det.mag > gMax):
                    gMax = det.mag
            elif(det.band == 'r'):
                if(det.mag < rMin):
                    rMin = det.mag
                if(det.mag > rMax):
                    rMax = det.mag
            elif(det.band == 'i'):
                if(det.mag < iMin):
                    iMin = det.mag
                if(det.mag > iMax):
                    iMax = det.mag
            elif(det.band == 'z'):
                if(det.mag < zMin):
                    zMin = det.mag
                if(det.mag > zMax):
                    zMax = det.mag
            else:
                print('Error: band is ' + det.band)
        if(gMax - gMin > thresh and gMin != 30):
            return False
        if(rMax - rMin > thresh and rMin != 30):
            return False
        if(iMax - iMin > thresh and iMin != 30):
            return False
        if(zMax - zMin > thresh and zMin != 30):
            return False
        return True

    # Assigns a unique number to each fake and returns the maximum of those numbers.
    def majFake(self):
        fakeDict = {}
        for det in self.dets:
            if(det.fakeid in fakeDict):
                fakeDict[det.fakeid] += 1
            else:
                fakeDict[det.fakeid] = 1
        return max(fakeDict, key=fakeDict.get)

    def __hash__(self):
        return hash(tuple(sorted([hash(x) for x in self.dets])))

    def __eq__(self, other):
        if(not isinstance(other, Triplet)):
            return False
        if other == None:
            return False
        if(isinstance(trip.dets[0], Detection)):
            selfIDs = tuple([x.objid for x in self.dets]) 
            otherIDs = tuple([x.objid for x in other.dets]) 
            return sorted(selfIDs) == sorted(otherIDs)
        elif(isinstance(trip.dets[0], int)):
            selfIDs = tuple([x for x in self.dets])
            otherIDs = tuple([x for x in other.dets])
            return sorted(selfIDs) == sorted(otherIDs)
        else:
            return False

    def __ne__(self, other):
        if other == None:
            return True
        selfIDs = tuple([x.objid for x in self.dets])
        otherIDs = tuple([x.objid for x in other.dets])
        return not (sorted(selfIDs) == sorted(otherIDs))

    def __str__(self):
        if(all(isinstance(x, Detection) for x in self.dets)):
            return self.toStr()
        else:
            s = 'TrackID: ' + str(self.trackid) + '\nDet: '
            for d in self.dets:
                s+= '\n' + str(d)
            s+= '\nCands: ' + str(self.cands)
            return s 

    def toStr(self):
        s = '\n *****Triplet*****\n'
        s += 'trackid: ' + str(self.trackid) + '\n'
        s += 'realLength: ' + str(self.realLength()) + '\n'
        for x in range(len(self.dets)):
            s += ('Det' + str(x) + ': ' + self.dets[x].toStr() + '\n')
        s+= 'elements: ' + str(self.elements) + '\n'
        s+= 'errs: ' + str(self.errs) + '\n'
        s+= 'chisq: ' + str(self.chiSq) + '\n'
       
        if(len(self.cands) > 0):
            for x in range(len(self.cands)):
                s += ('Cand' + str(x) + ': ' + str(self.cands[x]) + '\n')
        
        return s

#a class for processing detections in the csvfile, organizing them into fakeObjs
class fakeObj:
    def __init__(self, fakeid):
        self.fakeid = fakeid
        self.listobj = []
        self.campaigns = []

#saves the orbital elements from the orbitTable into the list of triplets
# tripList and orbitTable should have the same length
def saveElements(tripList, orbitTable):
    print('size of list: ' + str(len(tripList)))
    assert(len(tripList) == len(orbitTable))
    
    trackids = orbitTable['ORBITID'].tolist()
    elements = orbitTable['ELEMENTS'].tolist()
    errs = orbitTable['ELCOV'].tolist()
    chisqs = orbitTable['CHISQ'].tolist()
    flags = orbitTable['FLAGS'].tolist()
    abg = orbitTable['ABG'].tolist()
    cov = orbitTable['ABGCOV'].tolist()
    aeiCov = orbitTable['ELCOV'].tolist()

    newList = []
    trackDict = {}
    for trip in tripList:
        trackDict[trip.trackid] = trip
    for x in range(len(trackids)):
        trip = trackDict[trackids[x]]
        try:
            assert(trip.trackid == trackids[x])
        except AssertionError:
            print(trip)
            print(trackids[x])
        if(flags[x] != 0):
            trip.chiSq = -1
            trip.elements = 0
            trip.errs = 0
        else:
            trip.chiSq = chisqs[x]
            els = elements[x]
            ers = errs[x]
            trip.aeiCov = aeiCov[x]           
            trip.elements = {'a': els[0], 'e': els[1], 'i': els[2],
                            'lan': els[3], 'top': els[5], 'aop': els[4]}
            trip.abg = abg[x]
            trip.cov = cov[x]
            trip.errs = {'a': np.sqrt(ers[0]), 'e': np.sqrt(ers[7]), 'i': np.sqrt(ers[17]), 
                        'lan': np.sqrt(ers[21]), 'top': np.sqrt(ers[28]), 'aop': np.sqrt(ers[35])}
        newList.append(trip)
    return newList

#Write the triplet list to a file
def writeTriplets(tripList, outfile, writeOrbit=False, isObj=True):
    print('\nwriting triplets to: ' + outfile)
    if(writeOrbit):
        splitname = outfile.split('+')
        orbitname = splitname[0] + '+' + splitname[-1].split('.')[0] + '.orbit'
        try:
            print('reading from ' + str(orbitname))
            orbits = Table.read(orbitname, format='fits')
            tripList = saveElements(tripList, orbits)
        except IOError:
            print('orbits file does not exist: ' + str(orbitname))
        except AssertionError:
            print('orbits not updated')
    with open(outfile, 'w+') as output:
        count = 0
        time0 = time.time()
        for triplet in tripList:
            try:
                triplet.sortByMjd()
            except AttributeError:
                pass
            #printPercentage(count, len(tripList), time.time() - time0)
            '''
            if(isObj and len(triplet.dets) >= 2):
                if(triplet.elements == 0 and triplet.errs == 0 and writeOrbit == True):
                    elements, errs = triplet.calcOrbit()
                    triplet.getChiSq() 
                output.write('triplet ' + str(count) + ': ' + str(triplet))
            '''
            output.write('triplet ' + str(count) + ': ' + str(triplet))
            output.write('*****************\n\n')
            count += 1

# sets all orbit to 0 before pickling because Orbit cannot be pickled
def pickleTriplets(tripList, outfile, isObj=True, rmunbound=False):
    if(isObj):
        for trip in tripList:
            trip.orbit = 0
    if(rmunbound):
        temp = []
        for trip in tripList:
            try:
                if(not(np.isnan(trip.elements['aop']) or 
                        np.isnan(trip.elements['top']) or 
                        np.isnan(trip.elements['lan']) or
                        np.isnan(trip.elements['i']))):
                    temp.append(trip)
            except TypeError:
                pass
        tripList = temp
    with open(outfile, 'wb') as f:
        print('\npickling to ' + outfile)
        pickle.dump(tripList, f)

#Returns true iff a and b are closer than dist.
def approx(a, b, dist):
    return a-b<=dist and b-a <= dist

#Prints the percent completion.
def printPercentage(progress, finish, time):
    percent = float(progress)/float(finish) *100.0
    sys.stderr.write("\rCompleted: %1f %% after %5f sec" % (percent, time))
    sys.stderr.flush()

#calculates the distance between two detections in degreees
def calcDist(det1, det2):
    deltaRa = det1.ra - det2.ra
    deltaDec = det1.dec - det2.dec
    dist = deltaRa*deltaRa+deltaDec*deltaDec
    return np.sqrt(dist)

#takes a variable and tells whether it's a float
def isNumber(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

# Creates and returns a tripDict of detections.
# prints its progress in process.
def fakeDict(csvFile):
    detlist = wrapDets(csvFile)
    tripDict = {}
    size = len(detlist)
    time0 = time.time()
    # Records an inital time.
    counter = 0
    for det in detlist:
        i = det.fakeid
        if(i > 800000000):
            i = 0
        if(i not in tripDict):
            tripDict[i] = [det]
        else:
            tripDict[i].append(det)
        counter += 1
        printPercentage(counter, size, time.time()-time0)
    return tripDict

# Returns a dictionary of detections called expdict.
def expDictionary(dets):
    expdict = {}
    for det in dets:
        if det.expnum not in expdict:
            expdict[det.expnum] = [det]
        else:
            expdict[det.expnum].append(det)
    return expdict

# input: a pandas dataframe with the necessary headers for the Detection object
#           false: returns dictionary, true: saves to pickle file
# output: a dictionary that goes from objid to Detection object
def objidDictionary(csvFile, lookAhead=0, printP=False, toPickle=False, efficient=False):
    cand = False
    if(lookAhead == -1):
        cand = True
        lookAhead = 0
    print('start wrap')
    detList = wrapDets(csvFile, lookAhead, printP, efficient)
    detDict = {}
    size = len(detList)
    print('\nadding ' + str(size) +' objects to dictionary') 
    counter = 0
    time0 = time.time()
    for det in detList:
        if(cand):
            det.lookAhead = -1
        detDict[det.objid] = det
        if(printP):
            printPercentage(counter, size, time.time()-time0)
        counter += 1
    if(toPickle):
        savename = (csvFile.split('+')[-1] if '+' in csvFile else
                    csvFile.split('/')[-1] if '/' in csvFile else
                    csvFile)
        
        savename = 'objDict+' + savename.split('.')[0]
        with open(savename + '.pickle', 'wb') as f:
            print('\npickling to ' + savename)
            pickle.dump(detDict, f)
    else:
        return detDict

#takes a csvfile of detections and wraps them with Detection class
def wrapDets(csvFile, lookAhead=0, printP=False, efficient=False):
    df = pd.read_csv(csvFile)
    df.rename(str.lower, axis='columns', inplace=True)
    df.rename(columns={'snobjid': 'objid',
            'snfake_id': 'fakeid', 'ccdnum': 'ccd'}, inplace=True)

    size = len(df['objid'])
    raList = df['ra'].tolist()
    decList = df['dec'].tolist()
    mjdList = df['mjd'].tolist()
    if 'mag' in df.columns:
        print('using mag')
        magList = df['mag'].tolist()
        fluxList = [0]*len(magList)
    else:
        print('using flux')
        fluxList = df['flux'].tolist()
    objidList = df['objid'].tolist()
    expnumList = df['expnum'].tolist()
    ccdList = df['ccd'].tolist()
    bandList = df['band'].tolist()
    try:
        fakeidList = df['fakeid'].tolist()
    except KeyError:
        print('no fakeids')
        fakeidList = [0]*len(objidList)
    try:
        fluxErrList = df['flux_err'].tolist()
    except KeyError:
        print('no flux err')
        fluxErrList = [10]*len(objidList)
    try:
        posErrList = df['errawin_world'].tolist()
    except KeyError:
        print('no poserr')
        posErrList = [0]*len(objidList)
    detlist = []
    startT = time.time()
    print('wrapping ' + str(size) + ' objects')
    update = 60
    
    for y in range(size):
        if(time.time() - startT > update and printP):
            printPercentage(y,size, time.time()-startT)
            update += 60
        if(efficient):
            det = Det(ra=float(raList[y]), dec=float(decList[y]), mjd=float(mjdList[y]),
                        mag=float(magList[y]), objid=float(objidList[y]), 
                        expnum=float(expnumList[y]),
                        ccd=float(ccdList[y]), band=bandList[y], 
                        fakeid=float(fakeidList[y]), posErr=float(posErrList[y]), lookAhead=-1)
        else:
            det = Detection(float(raList[y]), float(decList[y]), float(mjdList[y]),
                    float(fluxList[y]), int(objidList[y]),
                    int(expnumList[y]), int(ccdList[y]), bandList[y],
                    lookAhead, int(fakeidList[y])) 
            det.setMagErr(int(fluxErrList[y]))
            det.posErr = float(posErrList[y])
            if(det.mag > 40):
                det.mag = float(magList[y])
        detlist.append(det)
    return detlist

# takes a list of detections and returns all detections in an exposure
def findDetections(exposure, detList):
    detRes = []
    for det in detList:
        if(det.expnum == exposure):
            destRes.append(det)
    return detRes


# takes a chunkname and a savename and determines the trackID header
def setTrackId(chunkname, savename):
    #assumes savename is ????+????+...+SN???_SEASON###_ML##.???
    parms = savename.split('+')[-1].split('.')[0].split('_')
    #SEASON###
    season = parms[1]
    season = season[6:9]
    #ML##
    ml = parms[2]
    ml = ml[2:4]
    #SN????
    obj = parms[0]
    obj = obj[2:]
    
    if(obj == 'FAKE'):
        obj = '3'
    elif(obj == 'OBS'):
        obj = '2'
    elif(obj == 'ALL'):
        obj = '1'
    else:
        print('Error: objtype is ' + obj)
        sys.exit()
    
    num = [s for s in chunkname if s.isdigit()]
    num = ''.join(num)
    num = '%06d' % int(num)

    TRACK_ID_HEADER = obj + season + num 
    print('trackid header: ' + TRACK_ID_HEADER)
    return TRACK_ID_HEADER

# takes in a dictionary from trackid to list of dets and writes to a file
# to feed into the orbit fitter
# trackIDs is a dictionary from a 6 digit string to a list of detections
def writeDetToOrb(trackIDs, outName='defaultTempSift.fits', verbose=True):
    ra_sum = 0
    dec_sum = 0
    mjd_sum = 0
    num_dets = 0
    if(verbose):
        print('\nnumber of triplets to be fit: ' + str(len(trackIDs)))
    for track in trackIDs.values():
        for det in track:
            ra_sum += det.ra
            dec_sum += det.dec
            mjd_sum += det.mjd
            num_dets += 1
    if(num_dets == 0):
        print('nothing to fit')
        return False
    ra_center = ra_sum / num_dets
    dec_center = dec_sum / num_dets
    mjd_center = mjd_sum / num_dets

    time0 = time.time()
    hdu = fits.PrimaryHDU()
    trackList = []
    expList = []
    raList = []
    decList = []
    errList = []
    objList = []
    for trackID in sorted(trackIDs.iterkeys()):
        dets = trackIDs[trackID]
        for det in dets:
            trackList.append(trackID)
            expList.append(det.expnum)
            raList.append(det.ra)
            decList.append(det.dec)
            if(det.posErr == 0):
                errList.append(det.getPosErr()*3600)
            else:
                errList.append(det.posErr*3600)
            objList.append(det.objid)
    #print(trackList)
    outTable = Table([trackList, expList, raList, decList, errList, objList],
                names=('ORBITID', 'EXPNUM', 'RA', 'DEC', 'SIGMA', 'OBJID'),
                dtype=('int64', 'i4', 'f8', 'f8', 'f8', 'i8'))
    
    binTable = fits.BinTableHDU(outTable)
    binTable.header['RA0'] = ra_center
    binTable.header['DEC0'] = dec_center
    binTable.header['MJD0'] = mjd_center
    hdul = HDUList([hdu, binTable])
    if(verbose):
        print('writing to: ' + outName)
    hdul.writeto(outName, overwrite=True)
    time1 = time.time()
    if(verbose):
        print('fits save time: ' + str(time1-time0))
    return outName



