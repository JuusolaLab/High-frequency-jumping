'''Script to analyse narrowing grating data
Copyright 2025 Joni Kemppainen
License: GPLv3 (only)

HOW TO
------
To analyse data, evoke this script from the command line
as follows:

    python plotkeiveen.py {FOLDER_1} {FOLDER_2} ... {FOLDER_N}

    where FOLDER_i contains Biosystfiles

The script writes files in folder "results in the
current working directory.
    
    analysis.pdf
    analysis.csv
    *.csv files

MORE INFO
---------

See the repository at
    github.com/juusolalab/High-frequency-jumping
for further readme instructions.


'''

import sys
import os
import csv


import numpy as np
import matplotlib.pyplot as plt
import scipy.signal 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from biosystfiles import extract


def dostim_freqs(X, D, l0=5, l1=0.33):
    '''Creates stimulus frequences.

    Creates tohse that match the experimental

    X       Time points
    D       Stimulus parameter: Time scale
    l0      Stimulus parameter: Start wavelength
    l1      Stimulus parameter: Stop wavelength
    '''
    dt = X[1]-X[0]

    stimulus = []
    for t in X:
        l = l0 * ((l1/l0) ** (t/D))
        stimulus.append(l)
    return np.array(stimulus)



def timetransform(freqs):
    '''Frequency normalise the frequency scale

    The idea is that freqs are the X coordinates,
    and when plotted with the Y coordinates, you get
    the narrowing grating stimuli.

    But his function normalised the freqs (ie. X) so
    that when you plot X,Y you get ideally a sinusoidal
    response.

    Arguments
    ---------
    freqs : list of floats
        The spatial frequencies to be normalised.
    '''
    newx = [0]
    for i,f in enumerate(freqs[1:]):
        nx = 1/freqs[i+1] + newx[-1]
        newx.append(nx)
    return np.array(newx)


def normalize(wave, win=50):
    '''Amplitude normalise a response

    Forces the response between 0 and 1 really well 
    but does not exactly retain the waveform.

    Arguments
    ---------
    wave : list
        The resposne to amplitude normalise
    win : int
        Default 50 produces good results. Test with 1
        or with 1000 to see bad ones.
    '''
    nwave = []
    for i in range(len(wave)-win):
        point = wave[i] - np.min(wave[i:i+win])
        point = point / (np.max(wave[i:i+win]) - np.min(wave[i:i+win]))
        nwave.append(point)

    for i in range(win):
        nwave.append(point)
    return np.array(nwave)



def calc_score(data):
    '''Calculate score for normalised data

    The scores are essentially signal to noise ratio
    of the response. As the grating gets narrower and narrower,
    the structured shape of the response disappears and gives
    space to random noisiness. This calc score function
    captures exactly that.
    
    Arguments
    ---------
    data : list
        Normalised data
    '''
    scores = []
    
    # This window parameter affects the base line of the
    # analysis - Too small window, lots of noise, underestimates
    # the resolvability. Too big window, "SNR" plot becomes
    # too smooth and easily over-estimates the resolvability.
    #
    # See the paper supplement for some of the analysis on how
    # we selected the window size

    #win = 7
    #win = 14
    #win = 30
    #win = 62
    #win = 126
    #win = 252
    win = 504
    #win = 1008
    #win = 2016
    #win = 4032
    #win = 8064
    N = len(data)
    
    for i in range(N):
        score = 0

        a = i-win//2
        b = i+win//2

        if a < 0:
            #b = abs(a) + b
            a = 0
        
        if b >= N:
            #a = a-b
            #b = N-1
            b = N-1

        sect = data[a:b]

        for j in range((b-a)-1):
            score += abs(sect[j]-sect[j+1]) #* (win//2-abs(j-win//2))/(win/2)
        score /= abs(b-a)
        scores.append(-score)

    return scores



def calc_lastresolved(freqs, scores):
    '''Calculates what was the last resolved angle

    Fits a horizontal line to the end of the scores, and
    checks the first score that goes under that line, and
    reports the corresponding frequency, and the fitted hline too.

    Returns
    -------
    hline : float
        Horinzontal noise floor
    last : float
        The final resolved frequency
    '''

    if len(scores)!=len(freqs):
        raise ValueError('freqs and scores differ in len')
    
    N=len(freqs)
    hline = np.mean(scores[3*N//4:-5000])

    last = None
    for i in range(N//5,N):
        if scores[i] < hline:
            last = freqs[i]
            break

    return hline, last



def process_file(fn, pdf):
    '''Analyse a Biosyst data file and write PDF out of results

    Arguments
    fn : string
        Path to the Biosyst data file
    pdf : object
        A PdfPages object from matplotlib
    '''
    data, fs = extract(fn, 0)
    data = data.T[0]
    
    origdata = data.copy()

    X = np.linspace(0, len(data)/fs, len(data))

    freqs = dostim_freqs(X, 40)
    
    nX = timetransform( freqs )
    nX /= np.max(nX)
    nX *= X[-1]
    data = data
    
    # Amplitude normalise the data
    data = normalize(data) 
    
    chunk = 5000
    nchunks = len(data)//chunk
    hratios = [1 for i in range(nchunks+2)]
    hratios[0] = 3
    hratios[1] = 3
    fig, axes = plt.subplots(
            nchunks+2, 1,
            sharex=False, sharey=False,
            figsize=(8.3,11.7),
            height_ratios=hratios,
            )

    # Plot original data (untouched)
    axes[0].plot(X, origdata, lw=0.5, color='black')
    axes[0].axis("off")
    axes[0].set_title(fn.split('data_projector/')[-1],
                      fontsize=8)
    
    # Calculate "scores" and draw noise floor
    scores = calc_score(data)
    axes[1].plot(scores, color='black', lw=0.5)
    hline, lastresolved = calc_lastresolved(freqs, scores)
    axes[1].plot([0,len(scores)],[hline,hline], color='gray',
                  linestyle='--')
    
    # Plot last resolved anlge
    i_lastresolved = freqs.tolist().index(lastresolved)
    axes[1].plot([i_lastresolved, i_lastresolved],
                  [np.min(scores),np.max(scores)],
                  color='orange', linestyle='--')

    yh = (np.max(origdata)+np.min(origdata))/2
    axes[1].axis("off")
   

    txt = ' '+str(round(lastresolved,2))+'°'
    txt += ' (' + str(round(X[i_lastresolved],2)) + ' s)'
    axes[0].annotate(txt, (X[i_lastresolved], yh))
    axes[0].plot([X[i_lastresolved], X[i_lastresolved]],
                 [np.min(origdata), np.max(origdata)],
                 color='orange', linestyle='--')

    # Each chunk corresponds to a row in the plot, to
    # spread the data more easily visible (not packed)
    for i, ax in enumerate(axes[2:]):

        # Plot the amplitude and frequency normalised data
        ax.plot(nX[i*chunk:(i+1)*chunk],
                data[i*chunk:(i+1)*chunk], lw=0.5, color='black')
        
        # Add degree lines
        for j in range(0,5):
            J = int(j*chunk/(6-1))
            xp = nX[i*chunk+J]
            yp = freqs[i*chunk+J]
            if yp < 0.33:
                yp = 0.33
                add = ' (end)'
            else:
                add = ''
            ax.annotate(' '+str(round(yp,2))+'°'+add, (xp,1.5))
            ax.plot([xp, xp], [0,2], color='gray', lw=1)

            if yp == 0.33:
                break
            
            # Adding the last resolved
            yp = freqs[i*chunk:(i+1)*chunk]
            if yp[0] > lastresolved >= yp[-1]:
                for index, fre in enumerate(yp):
                    if fre <= lastresolved:
                        xpos = nX[i*chunk+index]
                        ax.plot([xpos,xpos],[0,1.25],
                                color='orange', lw=1, linestyle='--')
                        break

        ax.set_axis_off()
        ax.set_ylim(0,2)

    pdf.savefig()
    plt.close()

    sfn = os.path.basename(fn)
    with open(f'results/{sfn}.csv', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        for i in range(len(data)):
            a = X[i]    # X-axis
            b = nX[i]   # transformed X-axis
            c = origdata[i] # original, untouched data
            d = data[i] # normalized data
            e = scores[i] # sequence similarity
            f = freqs[i] # stimulus frequency
            #e = ndata
            writer.writerow([a,b,c,d,e,f])
    return lastresolved



def process_folder(folder, pdf):
    '''Iteratively walk through Biosyst (.mat) files in a folder
    and run the process_file on each Biosyst file there,
    writing to the given PdfPages instance.
    
    Arguments
    ---------
    folder : string
        A folder that contains Biosystfiles, or alternatively,
        the path of one Biosyst file to analyse only one file. 
    pdf : object
        A PdfPages object from matplotlib

    Returns
    -------
    lasts : dict
        Keys are the name of the data file and items are
        the last resolved angle (numerical)
    '''
    
    if os.path.isfile(folder):
        fns = [folder]
    else:
        fns = [os.path.join(folder, fn) for fn in os.listdir(
            folder) if fn.endswith('.mat')]
    
    # Deterministic ordering
    fns.sort()

    lasts = {}

    for i, fn in enumerate(fns):
        print(f'{i+1}/{len(fns)}')
        lastresolved = process_file(fn, pdf)
        lasts[fn.split('data_projector/')[-1]] = lastresolved
    return lasts



def main():
    '''Process the data folders given as the input arguments
    when evoking the script
    '''

    os.makedirs('results', exist_ok=True)
    
    folders = sys.argv[1:]

    last_resolved = {}

    with PdfPages('results/analysis.pdf') as pdf:
        for ifol, folder in enumerate(folders):

            print(f'Fly folder {ifol+1}/{len(folders)}')
            lasts = process_folder(folder, pdf)
            last_resolved = {**last_resolved, **lasts}
    
    with open('results/analysis.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        for key, value in last_resolved.items():
            writer.writerow([key, value])


if __name__ == "__main__":
    main()


