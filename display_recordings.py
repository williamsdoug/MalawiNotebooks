# coding: utf-8

import pickle
from matplotlib import pyplot as plt
import numpy as np

from libDecel import extractDecels, extractUniformDecels, summarizeDecels, computeVariabilityRMS, filterAssociationsUC

from plot_helpers import plotAddTimescaleGuides, plotAddDecelAnnotations, showMask, showPlotWithBaseline
from plot_helpers import showPlot, showDualPlot

combinedDecelColorMap = {
    'prolonged':'m',                  # Magenta
    'variable':'#ffff00',             # Yellow
    'variable_periodic': '#FFA500',   # Orange
    'late_decel':'#ff0000',           # Red
    'mild_late_decel':'#FFA07A',      # Lifght Red/salmon
    'early_decel':'#808080',          # Gray
    'shallow_early_decel':'#D3D3D3',  # Light Gray
    'acceleration':'#008000',            # Green
    'borderline_acceleration':'#7CFC00', # Light Green
}

prolongedDecelColorMap={'prolonged':'m', 'borderline_prolonged':'y'}
variableDecelColorMap={'variable':'m', 'borderline_variable':'y'}

# currently Unused
# uniformDecelColorMap={'valid':'m', 'borderline':'y'}
# accelColorMap={'acceleration':'b', 'borderline_acceleration':'g'}


def summarizeProlongedDecels(allDecels, tStart=None, tEnd=None):
    found = False
    for decel in allDecels:
        if tStart and decel['tEnd'] < tStart:
            continue
        if tEnd and decel['tStart'] > tEnd:
            continue

        print '{:.1f}-{:.1f} min : {:5.1f} bpm  dur: {:3.1f} min  mag: {:5.1f} bpm  valid: {:0.0f}%  [{}]'.format(
            decel['tStart'], decel['tEnd'], -decel['relMag'], decel['duration'] / 60.0,
            decel['mag'], decel['pctValid'] * 100.0, decel['classification'])
        print '              variability: {:0.0f} bpm  time below -15bpm: {:.1f} min  time below 50% drop: {:.1f} min'.format(
            decel['var'], decel['time_under_15'], decel['time_under_50pct'])
        found = True

    return found


def summarizeUC(annotatedUC, minValid=0.5,
                earlyDecel=['early_decel', 'shallow_early_decel'],
                lateDecel=['late_decel', 'mild_late_decel', 'relaxed_late_decel'],
                lateLag=20,
                verbose=True, showDetail=False):
    """Also annotates decel is isEarly and isLate"""
    totalUC = 0
    totalEarly = 0
    totalLate = 0

    for entry in annotatedUC:

        if 'pctValidFHR' in entry and entry['pctValidFHR'] < minValid:
            if verbose and showDetail:
                print 'skipping {:0.2f}: pctValidFHR: {:0.1f}%'.format(
                    entry['tAcme'], entry['pctValidFHR'] * 100)
            continue

        totalUC += 1

        if entry['decel']:
            decel = entry['decel']
            # clean-up any prior characterizations -- not properties only used for plotting purposes
            if 'isLate' in decel:
                del decel['isLate']
            if 'isEarly' in decel:
                del decel['isEarly']

            if decel['classification'] in lateDecel:
                decel['isLate'] = True
                totalLate += 1
            elif decel['classification'] in earlyDecel:
                decel['isEarly'] = True
                totalEarly += 1

            if verbose and showDetail:
                print 'uc @ {:5.2f}m  decel @ {:5.2f}m  lag: {:4.1f}s  type: {}'.format(
                    entry['tAcme'], decel['tNadir'], entry['lag'], decel['classification'])

        else:
            if verbose and showDetail:
                print 'uc @ {:5.2f}m  no decel'.format(entry['tAcme'])

    if totalUC > 0:
        if verbose:
            print
            print 'Totals: Early: {:0.1f}% ({})  Late: {:0.1f}%  ({})  of total {}'.format(
                100.0 * totalEarly / totalUC, totalEarly, 100.0 * totalLate / totalUC, totalLate, totalUC)

    return totalEarly, totalLate, totalUC




def displayRecordingUniform(fhr, mask, ts, uc=np.array([]), filtUC=np.array([]),
                            normUC=np.array([]), allUC=[],
                            classifiers=[], reclassifiers=[], filterParams={},
                            colorMap={}, # currently unused
                            name='',
                            showIndividual=True, verbose=True, showDetail=False,
                            #fDetail=None, fSmooth=None, fBaselineInitial=None, fBaselineFinal=None,
                            plotIncr=25, plotOverlap=5):
    ret = extractUniformDecels(fhr, mask, ts, allUC=allUC,
                               classifiers=classifiers, reclassifiers=reclassifiers,
                               filterParams=filterParams,
                               # fDetail=fDetail, fSmooth=fSmooth, fBaselineInitial=fBaselineInitial,
                               # fBaselineFinal=fBaselineFinal
                               )

    allDecels = ret['decels']
    annotatedUC = ret['annotatedUC']
    baseline = ret['baseline']
    smoothHR = ret['smoothHR']
    detailHR = ret['detailHR']

    print
    print 'Summary of UC and Associated Periodic Decelerations'
    summarizeUC(annotatedUC, showDetail=True)

    #     return
    print

    for tStart in range(0, int(ts[-1]), plotIncr):
        tEnd = tStart + plotIncr + plotOverlap

        base = np.percentile(uc, 20)
        base = 0
        if uc is not None and len(uc) > 0:
            maxUC = np.max(uc)
        else:
            maxUC = np.max(filtUC)

        if showDetail:
            fig = plt.figure(figsize=(15, 6))
            ax1 = fig.add_subplot(411)
            plt.title("Recording: {}  Zoom {:0.0f}-{:0.0f} min".format(name, tStart, tEnd))
            plt.plot(ts, detailHR, 'r', alpha=0.25)
            plt.plot(ts, smoothHR, 'b')
            plt.plot(ts, baseline, 'm')
            plt.plot(ts, baseline - 15, 'm--')
            for decel in allDecels:
                if decel['classification'] in ['unspecified', 'other']:
                    plt.plot([decel['tNadir'], decel['tNadir']], [70, 200], 'k--')
                    plt.plot([decel['tStart'], decel['tStart']], [70, 200], 'g--')
                    plt.plot([decel['tEnd'], decel['tEnd']], [70, 200], 'r--')
                else:
                    plt.plot([decel['tNadir'], decel['tNadir']], [70, 200], 'k', lw=3)
                    plt.plot([decel['tStart'], decel['tStart']], [70, 200], 'g', lw=1.5)
                    plt.plot([decel['tEnd'], decel['tEnd']], [70, 200], 'r', lw=1.5)
            for t in allUC:
                plt.plot([t, t], [70, 200], 'm', lw=3)
            plotAddTimescaleGuides(70, 200, ts[-1])
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.ylabel('Raw FHR')
            plt.ylim(70, 200)

            ax2 = fig.add_subplot(412, sharex=ax1)
            plt.setp(ax2.get_xticklabels(), visible=False)
            plt.plot(ts, detailHR - baseline, 'r', alpha=0.25)
            plt.plot(ts, smoothHR - baseline, 'b')
            plt.plot([0, ts[-1]], [0, 0], 'k--')
            plt.plot([0, ts[-1]], [-15, -15], 'm--')
            for decel in allDecels:
                if decel['classification'] in ['unspecified', 'other']:
                    plt.plot([decel['tNadir'], decel['tNadir']], [-60, 20], 'k--')
                else:
                    plt.plot([decel['tNadir'], decel['tNadir']], [-60, 20], 'k', lw=3)
            for t in allUC:
                plt.plot([t, t], [-60, 20], 'm', lw=3)
            plotAddTimescaleGuides(-60, 20, ts[-1])
            plt.ylim(-60, 20)
            plt.ylabel('delta hr')
            plt.setp(ax2.get_xticklabels(), visible=False)

            ax3 = fig.add_subplot(413, sharex=ax1)
            plt.plot(ts, uc, 'r', alpha=0.25)
            plt.plot(ts, filtUC, 'b')
            for decel in allDecels:
                if decel['classification'] in ['unspecified', 'other']:
                    plt.plot([decel['tNadir'], decel['tNadir']], [0, maxUC], 'k--')
                else:
                    plt.plot([decel['tNadir'], decel['tNadir']], [0, maxUC], 'k', lw=3)
            for t in allUC:
                plt.plot([t, t], [0, maxUC], 'm', lw=3)
            plotAddTimescaleGuides(0, maxUC, ts[-1], alphaMajor=0.5, alphaMinor=0.25)
            plt.ylim(base, )
            plt.ylabel('Raw UC')
            plt.setp(ax3.get_xticklabels(), visible=False)

            ax4 = fig.add_subplot(414, sharex=ax1)
            plt.scatter(ts, mask, c='r', s=2)
            plt.ylim(-0.1, 0.1)
            plt.ylabel('Invalid')

            plt.xlabel('time (min)')
            plt.xlim(tStart, tEnd)
            plt.show()

        fig = plt.figure(figsize=(15, 4))
        ax1 = fig.add_subplot(211)
        plt.title("Early/Late Decel: {}  Zoom {:0.0f}-{:0.0f} min".format(name, tStart, tEnd))
        plt.plot(ts, detailHR, 'r', alpha=0.25)
        plt.plot(ts, smoothHR, 'b')
        plt.plot(ts, baseline, 'm')
        plt.plot(ts, baseline - 15, 'm--')
        for entry in annotatedUC:
            tAcme = entry['tAcme']
            decel = entry['decel']
            if decel:
                plt.plot([tAcme, tAcme], [50, 200], 'm', lw=3)
            else:
                plt.plot([tAcme, tAcme], [50, 200], 'm--')
            if not decel:
                pass
            elif decel['classification'] in ['unspecified', 'other']:
                plt.plot([decel['tNadir'], decel['tNadir']], [50, 200], 'k--')
            else:
                if 'isLate' in decel:
                    plt.plot([decel['tNadir'], decel['tNadir']], [50, 200], 'r', lw=3)
                elif 'isEarly' in decel:
                    plt.plot([decel['tNadir'], decel['tNadir']], [50, 200], 'y', lw=3)
                else:
                    plt.plot([decel['tNadir'], decel['tNadir']], [50, 200], 'k', lw=3)
        plotAddTimescaleGuides(50, 200, ts[-1])
        plt.ylim(50, 200)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.ylabel('Raw FHR')

        ax2 = fig.add_subplot(212, sharex=ax1)
        plt.plot(ts, uc, 'r', alpha=0.25)
        plt.plot(ts, filtUC, 'b')
        for entry in annotatedUC:
            tAcme = entry['tAcme']
            decel = entry['decel']
            if decel:
                plt.plot([tAcme, tAcme], [0, maxUC], 'm', lw=3)
            else:
                plt.plot([tAcme, tAcme], [0, maxUC], 'm--')
            if not decel:
                pass
            elif decel['classification'] in ['unspecified', 'other']:
                plt.plot([decel['tNadir'], decel['tNadir']], [0, maxUC], 'k--')
            else:
                if 'isLate' in decel:
                    plt.plot([decel['tNadir'], decel['tNadir']], [0, maxUC], 'r', lw=3)
                elif 'isEarly' in decel:
                    plt.plot([decel['tNadir'], decel['tNadir']], [0, maxUC], 'y', lw=3)
                else:
                    plt.plot([decel['tNadir'], decel['tNadir']], [0, maxUC], 'k', lw=3)
        plotAddTimescaleGuides(0, maxUC, ts[-1], alphaMajor=0.5, alphaMinor=0.25)
        plt.ylim(base, )
        plt.ylabel('Raw UC')

        plt.xlabel('time (min)')
        plt.xlim(tStart, tEnd)
        plt.show()

        if not showDetail:
            showMask(ts, mask, xMin=tStart, xMax=tEnd)

        print 'Decelerations with related UC:'
        # summarizeDecels(allDecels, tStart=tStart, tEnd=tEnd, exclude=['unspecified'])
        summarizeDecels(allDecels, tStart=tStart, tEnd=tEnd, filterHasUC=True)

        print
        print '-' * 10
        print




def displayRecordingVariable(fhr, mask, ts, uc=np.array([]), filtUC=np.array([]),
                             classifiers=[], filterParams={}, colorMap=variableDecelColorMap, name='',
                             #fDetail=None, fSmooth=None, fBaselineInitial=None, fBaselineFinal=None,
                             plotIncr = 25, plotOverlap = 5,
                             showFull=False, showDelta=True, showVariability=True, showIndividualDecels=False):

    extractorResults = extractDecels(fhr, mask, ts, computeProlonged=False, classifiers=classifiers,
                                     filterParams=filterParams,
                                     # fDetail=fDetail, fSmooth=fSmooth, fBaselineInitial=fBaselineInitial,
                                     # fBaselineFinal=fBaselineFinal
                                     )

    displayCommon(extractorResults, mask, ts, uc=uc, filtUC=filtUC,
                  colorMap=colorMap, name=name, plotIncr=plotIncr, plotOverlap=plotOverlap,
                  showFull=showFull, showDelta=showDelta,
                  showVariability=showVariability, showIndividualDecels=showIndividualDecels)
    return


def displayRecordingProlonged(fhr, mask, ts, uc=np.array([]), filtUC=np.array([]),
                              classifiers=[], filterParams={}, colorMap=prolongedDecelColorMap, name='',
                              #fDetail=None, fSmooth=None, fBaselineInitial=None, fBaselineFinal=None,
                              plotIncr=25, plotOverlap=5,
                              showFull=False, showDelta=False, showVariability=False, showIndividualDecels=True):

    extractorResults = extractDecels(fhr, mask, ts, computeProlonged=True, classifiers=classifiers,
                                     filterParams=filterParams,
                                     # fDetail=fDetail, fSmooth=fSmooth, fBaselineInitial=fBaselineInitial,
                                     # fBaselineFinal=fBaselineFinal
                                     )

    displayCommon(extractorResults, mask, ts, uc=uc, filtUC=filtUC,
                  colorMap=colorMap, name=name, plotIncr=plotIncr, plotOverlap=plotOverlap,
                  showFull=showFull, showDelta=showDelta,
                  showVariability=showVariability, showIndividualDecels=showIndividualDecels)
    return



def hasDecelInRegion(allDecels, tStart, tEnd, exclude=[]):
    for decel in allDecels:
        if tStart and decel['tEnd'] < tStart:
            continue
        if tEnd and decel['tStart'] > tEnd:
            continue
        if decel['classification'] in exclude:
            continue
        return True
    return False


def displayCommon(extractorResults, mask, ts, uc=np.array([]), filtUC=np.array([]), tsUC = None,
                  colorMap={}, name='', plotIncr=25, plotOverlap=5, showAltBaselines=False,
                  showFull=False, showDelta=False, showVariability=False,
                  showIndividualDecels=False,
                  showBorderline=True,
                  ignoreArtifacts=[],
                  ):

    if tsUC is None:    # legacy support
        tsUC = ts
    allDecels = extractorResults['decels']
    var = extractorResults['variability']
    varRMS = extractorResults['rmsVar']['varRMS']
    tRMS = extractorResults['rmsVar']['tRMS']
    baseline = extractorResults['baseline']
    smoothHR = extractorResults['smoothHR']
    detailHR = extractorResults['detailHR']

    fastBaseline = extractorResults['fastBaseline']
    medBaseline = extractorResults['medBaseline']
    slowBaseline = extractorResults['slowBaseline']

    if 'annotatedUC' in extractorResults:
        annotatedUC = extractorResults['annotatedUC']
    else:
        annotatedUC = []

    if ignoreArtifacts:
        filterAssociationsUC(annotatedUC, ignoreArtifacts)

    if 'sustainedUC' in extractorResults:
        sustainedUC = extractorResults['sustainedUC']
    else:
        sustainedUC = []


    if not showFull and len(allDecels) == 0:
        # print 'NO DECELS FOUND'
        return

    if annotatedUC:
        print
        print 'Summary of UC and Associated Periodic Decelerations'
        summarizeUC(annotatedUC, showDetail=True)

    if sustainedUC:
        for entry in sustainedUC:
            print 'Sustained UC @ {:4.1f} min  Width: {:4.0f}'.format(entry['tAcme'], entry['width'])


    includesEnd = False
    for tStart in range(0, int(ts[-1]), plotIncr):
        tEnd = tStart + plotIncr + plotOverlap

        if tEnd > ts[-1]:
            if includesEnd:
                continue
            includesEnd = True

            tStart = max(ts[-1] - (plotIncr + plotOverlap), 0)
            tEnd = tStart + plotIncr + plotOverlap

        if not (showFull or hasDecelInRegion(allDecels, tStart, tEnd)):
            continue

        showMask(ts, mask, xMin=tStart, xMax=tEnd)

        plt.figure(figsize=(15, 2))
        plt.title("FHR - {}  {:0.0f}-{:0.0f} min".format(name, tStart, tEnd))
        plt.plot(ts, detailHR, 'r', alpha=0.35)
        plt.plot(ts, smoothHR, 'b')
        if showAltBaselines:
            if len(fastBaseline) > 0:
                plt.plot(ts, fastBaseline, 'k--')  # 'm'
            if len(fastBaseline) > 0:
                plt.plot(ts, medBaseline, 'r--')   # 'm'
            if len(fastBaseline) > 0:
                plt.plot(ts, slowBaseline, 'm--')
        else:
            plt.plot(ts, baseline, 'k--')  # 'm'
            plt.plot(ts, baseline - 15, 'r--')
        plotAddDecelAnnotations(allDecels, 70, colorMap=colorMap, ignore=ignoreArtifacts)
        plotAddTimescaleGuides(70, 200, ts[-1])
        plt.xlim(tStart, tEnd)
        plt.show()

        if showDelta:
			plt.figure(figsize=(15, 2))
			plt.title("delta FHR - {}  {:0.0f}-{:0.0f} min".format(name, tStart, tEnd))
			plt.plot(ts, detailHR - baseline, 'r', alpha=0.25)
			plt.plot(ts, smoothHR - baseline, 'b')
			plt.plot([ts[0], ts[-1]], [0, 0], 'k--')
			plt.plot([ts[0], ts[-1]], [-15, -15], 'r--')
			plotAddDecelAnnotations(allDecels, 20, colorMap=colorMap, ignore=ignoreArtifacts)
			plotAddTimescaleGuides(-60, 20, ts[-1])
			plt.xlim(tStart, tEnd)
			plt.ylim(-60, 20)
			plt.show()

        if showVariability:
            plt.figure(figsize=(15, 1))
            plt.title('variability')
            plt.plot(ts, var)
            plt.plot([ts[0], ts[-1]], [12.5, 12.5], 'r--')
            plt.plot([ts[0], ts[-1]], [-12.5, -12.5], 'r--')
            plt.plot([ts[0], ts[-1]], [-2.5, -2.5], 'g--')
            plt.plot([ts[0], ts[-1]], [2.5, 2.5], 'g--')
            plt.plot([ts[0], ts[-1]], [0, 0], 'k--')
            plotAddDecelAnnotations(allDecels, 15, colorMap=colorMap, ignore=ignoreArtifacts)
            plt.xlim(tStart, tEnd)
            plt.ylim(-20, 20)
            plt.show()

            plt.figure(figsize=(15, 1))
            plt.title('Variability Intensity (RMS x 2sqrt(2))')
            plt.plot(tRMS, varRMS)
            plt.plot([ts[0], ts[-1]], [25, 25], 'r--')
            plt.plot([ts[0], ts[-1]], [5, 5], 'g--')
            plt.plot([ts[0], ts[-1]], [0, 0], 'k--')
            plotAddDecelAnnotations(allDecels, 35, colorMap=colorMap, ignore=ignoreArtifacts)
            plt.xlim(tStart, tEnd)
            plt.ylim(0, 40)
            plt.show()


        if filtUC is not None and len(filtUC) > 0:
            #base = np.percentile(uc, 20)
            base = 0
            if uc is not None and len(uc) > 0:
                maxUC = np.max(uc)
            else:
                maxUC = np.max(filtUC)

            plt.figure(figsize=(15, 1))
            plt.title("UC")
            if uc is not None and len(uc) > 0:
                plt.plot(tsUC, uc, 'r', alpha=0.25)
            plt.plot(tsUC, filtUC, 'b')
            # plt.plot([ts[0], ts[-1]], [base, base], 'r--')
            for auc in annotatedUC:
                plt.plot([auc['tAcme'], auc['tAcme']], [0, maxUC], 'r--')
                if 'decel' in auc and auc['decel']:
                    plt.plot([auc['decel']['tNadir'], auc['decel']['tNadir']], [0, maxUC], 'k--')

            for entry in sustainedUC:
                plt.plot([entry['tAcme'], entry['tAcme']], [0, maxUC], 'm', lw=4)
            plotAddTimescaleGuides(0, maxUC, tsUC[-1])
            plt.xlim(tStart, tEnd)
            plt.ylim(base, )
            plt.show()

        if showIndividualDecels:
            for decel in allDecels:
                if tStart and decel['tEnd'] < tStart:
                    continue
                if tEnd and decel['tStart'] > tEnd:
                    continue

                if decel['classification'] in ignoreArtifacts:
                    continue

                print
                # TODO: Remove after debug
                print decel.keys()

                plt.figure(figsize=(8, 2))
                plt.title('Zoom {:.1f}-{:.1f} min - FHR [{}]'.format(
                    decel['tStart'], decel['tEnd'], decel['classification']))
                plt.plot(ts, detailHR, 'r', alpha=0.35)
                plt.plot(ts, smoothHR, 'b')
                if 'base' in decel:
                    plt.plot([ts[0], ts[-1]], [decel['base'], decel['base']], 'k--')
                    if 'accel' in decel['classification']:
                        plt.plot([ts[0], ts[-1]], [decel['base'] + 15, decel['base'] + 15], 'r--')
                    else:
                        plt.plot([ts[0], ts[-1]], [decel['base'] - 15, decel['base'] - 15], 'r--')
                plotAddDecelAnnotations([decel], 55, colorMap=colorMap)
                plotAddTimescaleGuides(50, 200, ts[-1])
                plt.xlim(decel['tStart'] - 1, decel['tEnd'] + 1)
                plt.ylim(50, 200)
                plt.show()

                if filtUC is not None and len(filtUC) > 0:
                    if uc is not None and len(uc) > 0:
                        maxUC = np.max(uc)
                    else:
                        maxUC = np.max(filtUC)

                    plt.figure(figsize=(8, 1))
                    plt.title("UC")
                    if uc is not None and len(uc) > 0:
                        plt.plot(tsUC, uc, 'r', alpha=0.25)
                    plt.plot(tsUC, filtUC, 'b')

                    for auc in annotatedUC:
                        plt.plot([auc['tAcme'], auc['tAcme']], [0, maxUC], 'r--')
                        if 'decel' in auc and auc['decel']:
                            plt.plot([auc['decel']['tNadir'], auc['decel']['tNadir']], [0, maxUC], 'k--')

                    plotAddTimescaleGuides(0, maxUC, ts[-1])
                    plt.xlim(decel['tStart'] - 1, decel['tEnd'] + 1)
                    plt.ylim(0, )
                    plt.show()



                if 'var' in decel:
                    plt.figure(figsize=(8, 1.5))
                    plt.title('Zoom {:.1f}-{:.1f} min - Variability During Decel: {:0.0f} bpm'.format(
                        decel['tStart'], decel['tEnd'], decel['var']))
                    plt.plot(ts, var)
                    plt.plot([ts[0], ts[-1]], [12.5, 12.5], 'r--')
                    plt.plot([ts[0], ts[-1]], [-12.5, -12.5], 'r--')
                    plt.plot([ts[0], ts[-1]], [-2.5, -2.5], 'g--')
                    plt.plot([ts[0], ts[-1]], [2.5, 2.5], 'g--')
                    plt.plot([ts[0], ts[-1]], [0, 0], 'k--')
                    plotAddDecelAnnotations([decel], 18, colorMap=colorMap)
                    plt.xlim(decel['tStart'] - 1, decel['tEnd'] + 1)
                    plt.ylim(-20, 20)
                    plt.show()

                summarizeDecels(allDecels, tStart=decel['tStart'], tEnd=decel['tEnd'])
        else:
            summarizeDecels(allDecels, tStart=tStart, tEnd=tEnd, exclude=ignoreArtifacts)

        print
        print '-' * 10
        print

    print
    print '*' * 40
    print
