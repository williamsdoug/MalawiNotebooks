# coding: utf-8

import pickle
import os
from pprint import pprint
import datetime
from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
import traceback
import sys
#from analyze_recordings import computeBaseline

from libUltrasound import combineExtractionResults
from libUC import findUC

def selectRecordings(path, selectedDate=None):
    with open(path + '/patient_data.p', 'r') as f:
        data = pickle.load(f)

    if selectedDate:
        tMin = int(datetime.date(*selectedDate).strftime("%s"))
    else:
        tMin = 0

    subset_recordings = []
    for patient_uuid in data.keys():
        thisDir = path + '/' + patient_uuid
        stats = os.stat(thisDir)
        if stats.st_mtime < tMin:
            continue

        subset_recordings.append([stats.st_ctime, patient_uuid, ])

    subset_recordings = [x[1] for x in sorted(subset_recordings, reverse=True)]
    return subset_recordings, data



def detectUC(data):
    try:
        if 'uc_source' not in data or data['uc_source'] is None:
            # use manual annotations
            if 'uc_ann' in data and data['uc_ann'] is not None:
                allUC = [ann['time'] / 60.0 for ann in data['uc_ann'] if ann['annotation'] == 'acme']
            else:
                allUC = []
        else:
            posMin = data['uc']['posMin']
            _, allUC = findUC(data['uc']['uc'], posMin)
    except Exception, e:
        allUC = []
        print '-' * 60
        print 'Exception during detectUC', e
        traceback.print_exc(file=sys.stdout)
        print '-' * 60

    return allUC



def getRecordingsLowCostCTG(subset_recordings, data, path, min_duration_plot=10, includeSep=False, skip=0):
    patient_no = 1
    total_skipped = 0
    for patient_uuid in subset_recordings:
        patient = data[patient_uuid]

        thisDir = path + '/' + patient_uuid
        with open(thisDir + '/recordings_data.p', 'r') as f:
            recordings = pickle.load(f)

        for r in recordings.values():
            if r['duration'] < min_duration_plot:
                continue

            if skip > 0 and total_skipped < skip:
                total_skipped += 1
                patient_no += 1
                continue

            try:
                with open(thisDir + '/' + r['uuid'] + '.p', 'r') as f:
                    recording = pickle.load(f)

            except Exception, e:
                print 'Exception:', e
                continue

            if includeSep:
                print
                print '-' * 40

            print
            print 'Patient:', patient_no
            patient_no += 1
            print 'Patient:  {}, {} -- {}'.format(
                patient['last_name'], patient['first_name'], patient['uuid'])
            print 'Comment:', patient['comment']
            print
            print 'Recording - Duration: {:0.0f}m   Date: {}  {}'.format(r['duration'], r['date'], r['uuid'])
            print

            ts = recording['pitch']['pos'] / 60.0
            sigFinal, maskFinal = combineExtractionResults(recording['envelope']['hr'],
                                                           recording['pitch']['hr'], showStats=False)
            # ALign all array lengths
            N = min(len(sigFinal), len(ts), len(maskFinal))
            sigFinal = sigFinal[:N]
            ts = ts[:N]
            maskFinal = maskFinal[:N]

            indices = np.arange(len(sigFinal))
            estHR = np.interp(indices, indices[maskFinal], sigFinal[maskFinal])

            allUC = detectUC(recording)

            if 'uc_source' not in recording or recording['uc_source'] is None:
                yield estHR, maskFinal, ts, allUC, None, None, r['uuid']
            else:
                yield estHR, maskFinal, ts, allUC, recording['uc']['uc'], recording['uc']['posMin'], r['uuid']

            if includeSep:
                print
                print '*' * 40
                print



