#
#  Copyright Douglas Williams, 2017
#  All Rights Reserved
#

import numpy as np
from scipy import signal
import wfdb
from pprint import pprint

from scipy.interpolate import UnivariateSpline
#from ehg_ctg_uc_extraction import filterSignal
#from scipy.interpolate import InterpolatedUnivariateSpline
from os.path import expanduser
home = expanduser("~")

BASE_ctu_uhb_ctgdb = '/Volumes/Recordings/physionet/ctu-uhb-ctgdb'
#BASE_ctu_uhb_ctgdb = home + '/Documents/ctg_recordings'

def interpolateHR(sig, ts, mask, useLinear=True):
    if not useLinear:
        sp = UnivariateSpline(ts[mask], sig[mask]) #, s=smoothing)
        return sp(ts)
    else:
        return np.interp(ts, ts[mask], sig[mask])



def filterUC(sig, freqLow=0.025, filt_order=4, fs_Hz=4.0):
    """Extract Uterine Contraction from individual channel using Filter Method"""

    sigDC = filterSignal(sig, fType='lowpass', freq=freqLow,  # Filter DC
                         order=filt_order, fs_Hz=fs_Hz)

    return sigDC


def computeBaseline(estHR, mask, baselinePad=None,
                    perMin=60 * 4, minPerHist=10, pctValid=0.5, percentile=70, thresh=25):
    """Estimates baseline using 10 minute median"""

    if baselinePad is None:
        baselinePad = (5, 5)
    N = perMin * minPerHist
    minCount = N * pctValid

    # Initial estimate
    mask10 = [np.sum(mask[i:i + N]) > minCount for i in range(0, len(estHR) - N, perMin)]
    median10 = [np.percentile(estHR[i:i + N][mask[i:i + N]], percentile) if mask10[j] else 0
                for j, i in enumerate(range(0, len(estHR) - N, perMin))]
    maskMid10 = np.pad(mask10, (5, 5), 'edge')
    medianMid10 = np.pad(median10, (5, 5), 'edge')
    idx = np.arange(len(maskMid10)) * perMin

    baseline = np.interp(np.arange(len(estHR)), idx[maskMid10], medianMid10[maskMid10])

    # Now refine estimate by exclusing significant devations from baseline
    revisedMask = np.copy(mask)
    revisedMask[np.abs(estHR - baseline) > thresh] = 0

    mask10 = [np.sum(revisedMask[i:i + N]) > minCount for i in range(0, len(estHR) - N, perMin)]
    median10 = [np.percentile(estHR[i:i + N][revisedMask[i:i + N]], percentile) if mask10[j] else 0
                for j, i in enumerate(range(0, len(estHR) - N, perMin))]
    maskMid10 = np.pad(mask10, baselinePad, 'edge')
    medianMid10 = np.pad(median10, baselinePad, 'edge')
    idx = np.arange(len(maskMid10)) * perMin

    baseline2 = np.interp(np.arange(len(estHR)), idx[maskMid10], medianMid10[maskMid10])
    return baseline2


def altComputeBaseline(estHR, mask,
                       segmentSizeInMin=4, samplesPerMin=60 * 4, pctValid=0.5,
                       percentile=70, thresh=10, includePoints=False):
    """Estimates baseline using 4 minute median of prior samples after excluding regions with
       high absolute deviation from segment median

    Inputs:
       estHR - Input signal after error correction
       mask - Mask indicating valid samples
       segmentSizeInMin - Duration to be used when computing baseline
       samplesPerMin - Number of samples per 1 minute interval
       pctValid - Fraction of valid samples required in order to compute validate baseline
       percentile - Percentile absolute deviation from median when computing segment stability
       thresh - Maximum allowable deviation (rel to percentile) to be considered stable segment
       """

    # Compute medians and median absolute difference for different segments and look for stability
    perSeg = segmentSizeInMin * samplesPerMin
    perHalfSeg = perSeg // 2
    # allMedians = [(i, med, np.median(np.abs(seg - med)),
    #                np.percentile(np.abs(seg - med), percentile))
    #               #for i in range(perSeg, len(estHR), perHalfSeg // 2)
    #               for i in range(perSeg, len(estHR), samplesPerMin)
    #               for seg in [estHR[i - perSeg:i][mask[i - perSeg:i]]]
    #               for med in [np.median(seg)]
    #               if np.sum(mask[i - perSeg:i]) > pctValid * perSeg]

    allSegs = [(i, estHR[i - perSeg:i][mask[i - perSeg:i]])
                  #for i in range(perSeg, len(estHR), perHalfSeg // 2)
                  for i in range(perSeg, len(estHR), samplesPerMin)
                  for seg in [estHR[i - perSeg:i][mask[i - perSeg:i]]]
                  if np.sum(mask[i - perSeg:i]) > pctValid * perSeg]


    allMedians = [(i, med, np.median(np.abs(seg - med)),
                   np.percentile(np.abs(seg - med), percentile))
                  for i, seg in allSegs
                  for med in [np.median(seg)]]

    valid = [x for x in allMedians if x[3] < thresh]
    # print 'valid:', len(valid)
    if len(valid) > 0:
        # baseline = np.interp(np.arange(len(estHR)),
        #                      [x[0] for x in allMedians if x[3] < thresh],
        #                      [x[1] for x in allMedians if x[3] < thresh])
        baseline = np.interp(np.arange(len(estHR)),
                             [x[0] for x in valid],
                             [x[1] for x in valid])
    else:
        print 'altComputeBaseline:  Unable to compute standard baseline, using alternative'
        base = np.median([x[1] for x in allMedians])
        baseline = np.ones(len(estHR)) * base

    if includePoints:
        return baseline, [{'idx':x[0], 'val':x[1]} for x in allMedians if x[3] < thresh]
    else:
        return baseline



def processRecordingCTG(src,
                        freqLow=0.0125 / 2, # freqLow=0.0125 / 16,
                        freqHigh=0.0125, freqHigh2=0.0125 * 2, baselinePad=None,
                        blackoutLeft=10, blackoutRight=10,
                        recordPrefix=BASE_ctu_uhb_ctgdb):
    recordName = recordPrefix + '/' + src
    # ctgRecord = wfdb.rdsamp(fname)
    sig, fields = wfdb.srdsamp(recordName)

    # print sig.shape
    fs = fields['fs']
    signames = fields['signame']
    units = fields['units']
    comments = fields['comments']

    #pprint(fields)

    ts = np.arange(sig.shape[0]) / float(fs) / 60.0
    rawFHR = sig[:, 0]
    rawUC = sig[:, 1]

    # Compute Mask, applying black-out to neighboring points
    # mask = rawFHR > 0
    # mask2 = np.copy(mask)
    # for i in range(len(mask)):
    #     if mask2[i] == 0:
    #         mask[max(i - blackoutLeft, 0): min(i + blackoutRight, len(mask))] = 0

    mask = maskBadFHR(rawFHR)
    mask = dialateMaskFHR(mask, blackoutLeft=blackoutLeft, blackoutRight=blackoutRight)
    estHR = interpolateHR(rawFHR, ts, mask)  # iterpolate missing datapoints


    #baseline = computeBaseline(estHR, mask, baselinePad=baselinePad)
    baseline = altComputeBaseline(estHR, mask)


    minVal = np.min(rawFHR[mask])  # clip extreme values
    estHR[estHR < minVal] = minVal
    estHR[estHR > 200] = 200
    smoothLowHR = filterUC(estHR, freqLow=freqLow, filt_order=4, fs_Hz=4.0)
    smoothHighHR = filterUC(estHR, freqLow=freqHigh, filt_order=4, fs_Hz=4.0)
    smoothHighHR2 = filterUC(estHR, freqLow=freqHigh2, filt_order=4, fs_Hz=4.0)

    filtUC = filterUC(rawUC, freqLow=0.025, filt_order=4, fs_Hz=4.0)

    selected_params = ['pH', 'BDecf', 'pCO2', 'BE', 'Apgar1', 'Apgar5',
                       'Gest. weeks', 'Weight(g)', 'Sex',
                       'Age', 'Gravidity', 'Parity', 'Diabetes ','Hypertension', 'Preeclampsia',
                       'Liq.', 'Pyrexia', 'Meconium',
                       'Presentation', 'Induced', 'I.stage', 'NoProgress', 'CK/KP', 'II.stage', 'Deliv']

    history = "\n".join([x for x in comments if startsWithAny(x, selected_params)])
    artifacts = "\n".join(['Placeholder'])

    return (ts, rawFHR, estHR, mask, baseline, smoothLowHR, rawUC,
            smoothHighHR, smoothHighHR2, filtUC, history, artifacts)

def startsWithAny(val, all_prefix):
     for pre in all_prefix:
         if val.startswith(pre):
             return True
     return False

def oldMaskBadFHR(rawFHR, blackoutLeft=5, blackoutRight=15):
    # Compute Mask, applying black-out to neighboring points
    mask = rawFHR > 0
    mask2 = np.copy(mask)
    for i in range(len(mask)):
        if mask2[i] == 0:
            mask[max(i-blackoutLeft, 0): min(i+blackoutRight, len(mask))] = 0
    return mask

def dialateMaskFHR(mask, blackoutLeft=10, blackoutRight=10):
    # Compute Mask, applying black-out to neighboring points
    mask2 = np.copy(mask)
    for i in range(len(mask)):
        if mask2[i] == 0:
            mask[max(i-blackoutLeft, 0): min(i+blackoutRight, len(mask))] = 0
    return mask


def maskBadFHR(rawFHR, pctLow=0.65, pctHigh=1.75, samplesPerSecond=4):
    """Mask badFHR using algorithm in Dawes Redman Criteria"""
    # Compute Mask, applying black-out to neighboring points
    inputMask = rawFHR > 0
    if sum(inputMask) == 0:
        return inputMask   # no valid data

    outputMask = np.copy(inputMask)

    if not inputMask[0] or not inputMask[1]:
        # use first minute as reference:
        segLen = 60 * samplesPerSecond   # 60 seconds @ 4 samplesPerSecond
        i = 0
        while sum(inputMask[i:i+segLen]) == 0:
            if i < len(inputMask) - segLen:
                i += segLen
                print 'maskBadFHR skipping to next minute', i
            else:
                print 'maskBadFHR -- no valid data'
                return inputMask
        seg = rawFHR[i:i+segLen][inputMask[i:i+segLen]]

        medVal = np.median(seg)  # use median as starting point
        if not inputMask[0]:
            first = medVal
        if not inputMask[1]:
            second = medVal
    else:
        first = rawFHR[0]
        second = rawFHR[1]

    for i in range(2, len(rawFHR)):
        if not inputMask[i]:
            continue
        thresh = (first + second) / 2.0
        if rawFHR[i] < pctLow * thresh or rawFHR[i] > pctHigh * thresh:
            outputMask[i:min(i + 3, len(outputMask))] = 0
        else:
            first = second
            second = rawFHR[i]

    return outputMask


def filterSignal(sig, fType='lowpass', useFiltFilt=True, fs_Hz=1000.0, freq=100.0, order=1):
    """Applies lowpass or highpass filter to signal"""

    # CreateFilter
    if fType not in ['lowpass', 'highpass']:
        raise Exception('Invalid filter type: {}'.format(fType))
    fNyquist = fs_Hz / 2.0
    order = min(order, max(int(fs_Hz / freq) - 1, 1))  # limit order to avoid oscillations
    b, a = signal.butter(order, freq / fNyquist, fType)
    if not useFiltFilt:
        zi = signal.lfilter_zi(b, a)

    if len(sig.shape) > 1:
        if useFiltFilt:
            newSig = np.vstack([signal.filtfilt(b, a, sig[i])
                                for i in range(sig.shape[0])])
        else:
            newSig = np.vstack([signal.lfilter(b, a, sig[i], zi=zi * sig[i, 0])[0]
                                for i in range(sig.shape[0])])
    else:

        if useFiltFilt:
            newSig = signal.filtfilt(b, a, sig)
        else:
            newSig, _ = signal.lfilter(b, a, sig, zi=zi * sig[0])

    return newSig
