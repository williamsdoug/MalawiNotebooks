from matplotlib import pyplot as plt


def plotAddTimescaleGuides(yMin, yMax, xEnd, alphaMajor=0.2, alphaMinor=0.1):
    for i in range(1, int(xEnd)):
        if i % 5 == 0:
            plt.plot([i,i], [yMin, yMax], 'k', alpha=alphaMajor)
        else:
            plt.plot([i,i], [yMin, yMax], 'k--', alpha=alphaMinor)


def plotAddDecelAnnotations(allDecels, yMarker, colorMap={}, lw=8, ignore=[]):
    for decel in  allDecels:
        if decel['classification'] in colorMap and decel['classification'] not in ignore:
            color = colorMap[decel['classification']]
            plt.plot([decel['tStart'], decel['tEnd']], [yMarker, yMarker], color, lw=lw)


def showPlot(sig, ts, mask, title='Envelope', showRaw=True):
    plt.figure(figsize=(15, 2))
    plt.title(title)
    plt.plot(ts[mask], sig[mask], 'r')
    if showRaw:
        plt.plot(ts, sig, 'b', alpha=0.25)
    for x in range(0, int(ts[-1])):
        if x % 5 == 0:
            plt.plot([x, x], [50, 200], 'k', alpha=0.25)
        else:
            plt.plot([x, x], [50, 200], 'k--', alpha=0.25)
    for y in range(50, 200, 10):
        plt.plot([0, ts[-1]], [y, y], 'k--', alpha=0.25)
    plt.xlim(ts[0], ts[-1])
    plt.ylim(50, 200)
    plt.show()


def showPlotWithBaseline(sig, ts, sigBase, title='Combined With Baseline'):
    plt.figure(figsize=(15, 4))
    plt.title(title)
    plt.plot(ts, sigBase, 'k', lw=3, alpha=0.5)
    plt.plot(ts, sigBase+15, 'k--')
    plt.plot(ts, sigBase-15, 'k--')
    for x in range(0, int(ts[-1])):
        if x % 5 == 0:
            plt.plot([x, x], [50, 200], 'k', alpha=0.25)
        else:
            plt.plot([x, x], [50, 200], 'k--', alpha=0.25)
    for y in range(50, 200, 10):
        plt.plot([0, ts[-1]], [y, y], 'k--', alpha=0.25)
    plt.plot(ts, sig, 'r')
    plt.xlim(ts[0], ts[-1])
    plt.ylim(50, 200)
    plt.show()


def showDualPlot(sig, ts, mask, sig2, ts2, mask2,  title='Both'):
    plt.figure(figsize=(15, 2))
    plt.title(title)
    plt.plot(ts[mask], sig[mask], 'r')
    plt.plot(ts2[mask2], sig2[mask2], 'b', alpha=0.5)

    for x in range(0, int(ts[-1])):
        if x % 5 == 0:
            plt.plot([x, x], [50, 200], 'k', alpha=0.25)
        else:
            plt.plot([x, x], [50, 200], 'k--', alpha=0.25)
    for y in range(50, 200, 10):
        plt.plot([0, ts[-1]], [y, y], 'k--', alpha=0.25)
    plt.xlim(ts[0], ts[-1])
    plt.ylim(50, 200)
    plt.show()


def showMask(ts, mask, s=2, xMin=None, xMax=None):
    # Show Invalid Mask

    plt.figure(figsize=(15, 0.5))
    plt.title('Invalid')
    plt.scatter(ts, mask, c='r', s=s)
    for x in range(0, int(ts[-1])):
        if x % 5 == 0:
            plt.plot([x, x], [-0.1, 1.1], 'k', alpha=0.25)
        else:
            plt.plot([x, x], [-0.1, 1.1], 'k--', alpha=0.25)
    plt.xlim(xMin if xMin else ts[0], xMax if xMax else ts[-1])
    plt.ylim(-0.1, 0.1)
    plt.show()

