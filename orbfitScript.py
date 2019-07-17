# Takes triplets from a pickled file and saves information about each detection (time, day, etc.) to a .dat file.
# --- -------- - ---------------- --------- --------------
'''
This program write a script for orbit fitter
'''
import sys
import os
tnopath = os.environ['TNO_PATH']
sys.path.insert(0, tnopath)
import argparse
import pickle
import math
from GammaTPlotwStatTNOExFaster import mjdToDate

# Takes 1, 2, etc and returns 01, 02. For two-digit numbers, returns the input.
def helpFormatTime(t):
    if len(t) == 1:
        t = '0' + t
    return t

# Formats time in hh:mm:ss.ss
def formatTime(h,m,s):
    h = helpFormatTime(h)
    m = helpFormatTime(m)
    if len(s.split('.')[0]) == 1:
        s = '0' + s
    return h + ':' + m + ':' + s

# Converts motion across the sky (degrees) into time (hh:mm:ss.sss)
def degToHour(deg):
    abval = abs(deg)
    h = int(math.floor(abval / 15))
    m = int(math.floor((abval % 15) * 4))
    s = round((abval - 15 * h - m / 4.0) * 240.0, 3)
    result = formatTime(str(h), str(m), str(s))
    if deg < 0:
        result = '-' + result
    return result

# Converts degrees into dd:mm:ss.ss
def degToDMS(degree):
    abval = abs(degree)
    deg = int(abval)
    arcmin = int((abval - deg) * 60)
    arcsec = round(3600 * (abval - deg) - 60 * arcmin, 2)
    result = formatTime(str(deg), str(arcmin), str(arcsec))
    if degree < 0:
        result = '-' + result
    return result

# Converts day and time (in hours:minutes:seconds) into day. (Returns day.)
def decimalDay(day, time):
    t = time.split(':')
    h = float(t[0])
    m = float(t[1])
    s = float(t[2])
    return str(int(day) + h / 24 + m / 1440 + s / 86400)

# Creates a file for each detection and writes information to it
def scriptWriter(triplets, saveName):
    for trip in triplets:
        # Cycles through each triplet
        saveas = saveName + '_' + str(trip.dets[0].fakeid) + '.dat'
        f = open(saveas, 'w+')
        # Creates a file named saveName (which is an input) + '_' + the fakeid + '.dat'
        for det in trip.dets:
            # Cycles through each detection
            date = mjdToDate(det.mjd).split('/')
            daytime = date[2].split('T')
            f.write(date[0] + ' ' + helpFormatTime(date[1]) + ' ' + decimalDay(daytime[0],daytime[1]) +
                    ' ' + degToHour(det.ra) + ' ' + degToHour(det.dec) + ' 0.1 807\n')
            # Writes information about each detection to the above file
        f.close()


# Finds the triplets from the pickled file and executes the scriptWriter with them.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pickle', nargs=1, help='path to pickle file')
    args = parser.parse_args()
    saveName = args.pickle[0].split('+')[-1].split('.')[0]
    triplets = pickle.load(open(args.pickle[0],'rb'))
    print(degToHour(-37.13241)) # This line may be unnecessary.
    scriptWriter(triplets, saveName)

if __name__ == '__main__':
    main()
