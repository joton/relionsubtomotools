#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import copy
import pandas as pd
from docopt import docopt
from ._version import __version__
import emfile
from reliontomotools.fileIO import relion4ToMotl, readStarFile

__all__ = ['relion2motlProgram']
# % Subtom motl file To Relion Converter

def relion2motlProgram(args=None):
    doc = """relion2motl: converts relion starfile to subtom motl metadata.

    Usage:
      star_to_motl.py -i <star_file>  -o <em_file> [--bin=<binning>]\
       [--angpix=<angpix>] [--tomo_angles]
    
    Arguments:
               <star_file>  Relion star file containing data_optics and 
                            data_particles tables.
                 <em_file>  EM file containing converted orientations to subtom
                            motl list.
    
    Options:
           --bin=<binning>  Bin factor for output motl list [default: 1]
      -p --angpix=<angpix>  Pixel size at bin 1 (in angstroms). If negative,
                            pixel size is obtained from optics table [default: -1]

          -t --tomo_angles  Instead of rlnAngle{Rot,Tilt,AnglePsi}, angles from
                            rlnTomoSubtomogram{Rot,Tilt,AnglePsi} are used.
    
                 -h --help  Show this screen.
              -v --version  Show version.

    """

    arguments = docopt(doc, version=__version__)

    starFile = arguments['<star_file>']
    emFile = arguments['<em_file>']
    binning = int(arguments['--bin'])
    useTomoAngles = arguments['--tomo_angles']
    pixelSize = float(arguments['--angpix'])

    dataStar = readStarFile(starFile)

    print(pixelSize)
    if pixelSize < 0:
        if 'data_optics' in dataStar:
            print(pixelSize)
            datao = dataStar['data_optics']
            pixelSize = float(datao['_rlnTomoTiltSeriesPixelSize'])
        else:
            raise Exception(f"File {starFile} does not contain data_optics table "
                            "or rlnTomoTiltSeriesPixelSize parameter.")

    motl = relion4ToMotl(dataStar, pixelSize, binning, useTomoAngles)
    emfile.write(emFile, motl.astype(np.float32)[np.newaxis, ...], overwrite=True)
