import pickle
import os
from pprint import pprint
import datetime
import copy
import numpy as np
import scipy.signal
from matplotlib import pyplot as plt
import sys


from libUC import compressSig, findUC, normalizeUC
from libTocopatchSignal import ProcessUC, isolateUC, UC_DEFAULT_PARAMS
from libTocopatchSignal import butter_bandpass_filter, butter_lowpass_filter
from libTocopatchSignal import clipExtremeAC, clipMinDC
from libFilterUC import filterArtfactsUC, PARAMS_FILTER_UC



def showAnnotatedPlot(posMin, uc, filtered, allUC, 
                      thresh_base, thresh, sustainedUC=[],
                      color='b', title='After'):

    ymax = np.max(filtered)
    scale = np.max(uc)/ymax
    
    plt.figure(figsize=(15, 2))
    plt.title(title)
    plt.plot(posMin, uc, c=color, lw=1.5)
    plt.plot(posMin, scale*filtered, c='g', alpha=0.125)
        
    plt.plot([0, posMin[-1]], [thresh_base, thresh_base], color+'--')
    plt.plot([0, posMin[-1]], [thresh, thresh], 'r--')

    for j in range(int(posMin[-1])):
        if j % 5 == 0:
            plt.plot([j,j], [-ymax,ymax], 'k', alpha=0.25)
        else:
            plt.plot([j,j], [-ymax,ymax], 'k', alpha=0.1)

    for j in allUC:
        plt.plot([j,j], [-ymax,ymax], color, lw=2)

    for entry in sustainedUC:
        j = entry['tAcme']
        plt.plot([j,j], [-ymax,ymax], 'm', lw=4)

    plt.xlabel('time in min')
    plt.xlim(posMin[0], posMin[-1])
    #plt.ylim(0,ymax+1)
    plt.ylim(0,np.max(uc))
    plt.show()


def showAllContractions(selectedRecordings, path,
                        squelchFactor=3, squelchPercentile=50,
                        minThresh=5, **kwargs):

    print 'kwargs:', kwargs
    
    with open(path+'/patient_data.p', 'r') as f:
        data = pickle.load(f)
 
    for key, v in data.items():
        name = '{} -- {}'.format(v['last_name'].lower(), v['first_name'].lower())
        
        
        if key not in selectedRecordings or 'deleteme' in name:
            continue
        
        with open(path+'/'+key+'/recordings_data.p', 'r') as f:
            detail = pickle.load(f)
        if not detail:
            continue

        print '-'*40
        print name, ':'
        for k, v in detail.items():
            if v['uuid'] not in selectedRecordings[key]:
                continue
                
            fname = '{}/{}/{}.p'.format(path, key, v['uuid'])
            print '{}/{}.p'.format(key, v['uuid'])
            #print 'Comments:', data[key]['comment']
            with open(fname, 'r') as f:
                recording = pickle.load(f)

            if 'uc' not in recording or 'pos' not in recording['uc']:
                print
                print '** Skipped'+'*'*40
                print
                continue

            print recording['uc'].keys()
            posSec = recording['uc']['pos']
            posMin = recording['uc']['posMin']

            if len(posSec) == 0:
                print
                print '** Skipped'+'*'*40
                print
                continue

            raw = recording['uc']['raw']
            fs = 1.0/(posSec[1] - posSec[0])
            
            #UC_DEFAULT_PARAMS['pctClipMin'] = None

            filtered, uc, alt_uc = isolateUC(raw, fs,
                                             **UC_DEFAULT_PARAMS)

            thresh = np.percentile(np.abs(filtered), 95)*2
            filtered = np.abs(filtered)
            filtered[filtered > thresh] = thresh

            ymax = thresh

            plt.figure(figsize=(15, 2))
            plt.title('{} -- Filtered Raw Signal'.format(name))
            plt.plot(posMin, filtered)
            for j in range(int(posMin[-1])):
                if j % 5 == 0:
                    plt.plot([j,j], [-ymax,ymax], 'k', alpha=0.5)
                else:
                    plt.plot([j,j], [-ymax,ymax], 'k', alpha=0.1)
            plt.xlabel('time in min')
            plt.xlim(posMin[0], posMin[-1])
            plt.ylim(0,ymax+1)
            plt.show()

            allIdx, allUC = findUC(uc, posMin)
            base = np.percentile(uc, squelchPercentile)
            squelchMag = max(minThresh, base*squelchFactor)

            if minThresh > thresh*squelchFactor:
                print '*** Threshold set to min abs value -- abs: {:0.0f} vs rel: {:0.0f}'.format(
                    minThresh, base*squelchFactor)

            maxVal = np.max(uc)
            if  maxVal < squelchMag:
                print '*** All values below  threshold -- max: {:0.0f} vs thresh: {:0.0f}'.format(
                    maxVal, squelchMag)
            if  maxVal < minThresh:
                print '*** All values below minimum threshold -- max: {:0.0f} vs thresh: {:0.0f}'.format(
                    maxVal, minThresh)

            showAnnotatedPlot(posMin, uc, filtered, allUC, 
                              base, squelchMag, title='Before')

            print 'Filtering on spacing and magnitude'
            _allUC, sustainedUC = filterArtfactsUC(uc, posMin, allIdx, minMag=squelchMag, **kwargs)
            showAnnotatedPlot(posMin, uc, filtered, _allUC, 
                              base, squelchMag, sustainedUC=sustainedUC,
                              title='After')

            print
            print '*'*40
            print
            
    return
