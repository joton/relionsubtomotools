#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import copy
import pandas as pd
from docopt import docopt
from ._version import __version__
import emfile
from reliontomotools.fileIO import motlToRelion4, writeStarFile

__all__ = ['motl2RelionProgram']
# % Subtom motl file To Relion Converter

def motl2RelionProgram(args=None):
    doc = """motl2relion: converts subtom motl metadata to Relion.

    Usage:
      motl2relion.py (<em_file> | --evenodd <em_even_file> <em_odd_file>)\
      -o <star_file> -p=<angpix> [--bin=<binning>]\
      [-t=<label>] [--pad_zeros=<n_zeros>]

    Arguments:
                  <em_file>  EM file containing a subtom motl list.
             <em_even_file>  EM file containing half 1 of a motl list.
              <em_odd_file>  EM file containing half 2 of a motl list.
    Options:
                  --evenodd  Provides two half files instead of one single file.
       -p --angpix=<angpix>  Pixel size at bin 1 (in angstroms)
         -b --bin=<binning>  Input coordinates Bin factor in motl file [default: 1]
    -t --tomo_label=<label>  Prefix label for rlnTomoName [default: TS_]
      --pad_zeros=<n_zeros>  Number of padding zeros in tomogram number name.
                             Negative is auto [default: -1]

                 -h --help  Show this screen.
              -v --version  Show version.

    """

    arguments = docopt(doc, version=__version__)

    pixelSize = float(arguments['--angpix'])
    binning = int(arguments['--bin'])
    padZeros = int(arguments['--pad_zeros'])
    tomoLabel = arguments['--tomo_label']

    bothmotl = not arguments['<em_file>']
    if bothmotl:
        emFile1 = arguments['<em_even_file>']
        emFile2 = arguments['<em_odd_file>']
    else:
        emFile1 = arguments['<em_file>']

    fnStarOut = arguments['<star_file>']

    _, motl1 = emfile.read(emFile1)
    motl1 = motl1[0]
    dataStar1 = motlToRelion4(motl1, pixelSize,
                              binning=binning,
                              tomoLabel=tomoLabel,
                              padZeros=padZeros)

    if bothmotl:
        _, motl2 = emfile.read(emFile2)
        motl2 = motl2[0]
        dataStar2 = motlToRelion4(motl2, pixelSize,
                                  binning=binning,
                                  tomoLabel=tomoLabel,
                                  padZeros=padZeros)
        [dataStar1['data_particles']['_rlnRandomSubset']] = '1'
        [dataStar2['data_particles']['_rlnRandomSubset']] = '2'
        dataOut = copy.deepcopy(dataStar1)
        datap = pd.concat([dataStar1['data_particles'],
                           dataStar2['data_particles']])
        refIdx = np.array([int(partId) for partId in datap['_rlnTomoParticleId']])
        idxsort = np.argsort(refIdx)
        datap = datap.iloc[idxsort]

        dataOut['data_particles'] = datap
    else:
        dataOut = dataStar1

    writeStarFile(fnStarOut, dataOut)
