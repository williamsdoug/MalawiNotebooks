# coding: utf-8

import numpy as np
from scipy import signal
import math


import wfdb
from os import listdir
from os.path import isfile, join

import ctg_processing

#BASE_ctu_uhb_ctgdb = '/Volumes/Recordings/physionet/ctu-uhb-ctgdb'
BASE_ctu_uhb_ctgdb = '/Users/doug/Documents/ctg_recordings'


def downsample(sig, mask, ts, factor=4):
    # downsamle signal and timescale
    sigD = signal.decimate(sig, factor, zero_phase=True)
    tsD = ts[::factor]
    # align lengths
    #N = min(len(sigD), len(tsD))
    N = min(len(sigD), len(tsD), len(sig)//factor, len(ts)//factor)
    sigD = sigD[:N]
    tsD = tsD[:N]

    # compute mask


    maskD = np.sum(mask[:N * 4].reshape((N, 4)), axis=1) > factor / 2.0

    return sigD, maskD, tsD


def loadPhysionetRecording(src, blackoutLeft=10, blackoutRight=10, mask_method='basic',
                           recordPrefix=BASE_ctu_uhb_ctgdb):
    recordName = recordPrefix + '/' + src
    # ctgRecord = wfdb.rdsamp(fname)
    sig, fields = wfdb.srdsamp(recordName)

    # print sig.shape
    fs = fields['fs']
    signames = fields['signame']
    units = fields['units']
    comments = fields['comments']
    # pprint(fields)

    ts = np.arange(sig.shape[0]) / float(fs) / 60.0
    rawFHR = sig[:, 0]
    rawUC = sig[:, 1]

    if mask_method == 'basic':
        mask = rawFHR > 0
    else:
        mask = ctg_processing.maskBadFHR(rawFHR)
        mask = ctg_processing.dialateMaskFHR(mask, blackoutLeft=blackoutLeft, blackoutRight=blackoutRight)

    fhr = np.interp(ts, ts[mask], rawFHR[mask])

    fhrD, maskD, tsD = downsample(fhr, mask, ts, factor=4)

    ucD, _, _ = downsample(rawUC, rawUC > 0, ts, factor=4)

    filtUC = ctg_processing.filterUC(ucD, freqLow=0.025, filt_order=4, fs_Hz=1.0)

    return fhrD, maskD, tsD, ucD, filtUC


def parseToken(tok):
    try:
        return int(tok)
    except Exception:
        try:
            return float(tok)
        except Exception:
            return tok


def loadPhysionetMetadata(src, recordPrefix=BASE_ctu_uhb_ctgdb,
                          selectedFields=['pH', 'BDecf', 'pCO2', 'BE', 'Apgar1', 'Apgar5']):
    recordName = recordPrefix + '/' + src
    # ctgRecord = wfdb.rdsamp(fname)
    sig, fields = wfdb.srdsamp(recordName)

    meta = {}
    for entry in fields['comments']:
        tokens = entry.split()
        if tokens[0] in selectedFields:
            meta[tokens[0]] = parseToken(tokens[1])
    return meta


def getAllRecordingNumbers(mypath=BASE_ctu_uhb_ctgdb):
    allRecno = sorted([f.split('.')[0] for f in listdir(mypath) if f.endswith('.hea')])
    return allRecno


def getSelectedRecordings(mypath=BASE_ctu_uhb_ctgdb,
                          max_BDecf=None, min_BDecf=None,
                          max_BE=None, min_BE=None,
                          max_pH=None, min_pH=None,
                          max_Apgar1=None, min_Apgar1=None,
                          max_Apgar5=None, min_Apgar5=None,
                          maxRecordings=None):
    allRecno = getAllRecordingNumbers(mypath)

    count = 0
    for recno in allRecno:
        if maxRecordings and count == maxRecordings:
            return
        meta = loadPhysionetMetadata(recno, recordPrefix=mypath)
        if math.isnan(meta['BDecf']) and (max_BDecf or min_BDecf):
            continue
        elif max_BDecf and meta['BDecf'] > max_BDecf or min_BDecf and meta['BDecf'] < min_BDecf:
            continue

        if math.isnan(meta['BE']) and (max_BE or min_BE):
            continue
        elif max_BE and meta['BE'] > max_BE or min_BE and meta['BE'] < min_BE:
            continue

        if math.isnan(meta['pH']) and (max_pH or min_pH):
            continue
        elif max_pH and meta['pH'] > max_pH or min_pH and meta['pH'] < min_pH:
            continue

        if math.isnan(meta['Apgar1']) and (max_Apgar1 or min_Apgar1):
            continue
        elif max_Apgar1 and meta['Apgar1'] > max_Apgar1 or min_Apgar1 and meta['Apgar1'] < min_Apgar1:
            continue

        if math.isnan(meta['Apgar5']) and (max_Apgar5 or min_Apgar5):
            continue
        elif max_Apgar5 and meta['Apgar5'] > max_Apgar5 or min_Apgar5 and meta['Apgar5'] < min_Apgar5:
            continue

        count += 1
        yield recno, meta


def showMetadata(recno, meta):
    print 'Recording {} -- pH: {:0.2f}  Apgar1: {:0d}  Apgar5: {:0d}  BDecf: {:0.2f}  BE: {:0.2f}  pCO2: {:0.2f}'.format(
        recno, meta['pH'], meta['Apgar1'], meta['Apgar5'], meta['BDecf'], meta['BE'], meta['pCO2'])