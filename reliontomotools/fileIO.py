#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import io
import xmltodict
from glob import glob
from transforms3d.euler import euler2mat, mat2euler
import re

from ._version import __version__


__all__ = ['readStarFile', 'writeStarFile', 'WarpXMLHandler',
           'motlToRelion4', 'relion4ToMotl'
           ]


def getTable(lines):

    nLines = len(lines)
    lId = 0
    labels = list()

    if 'loop_' in lines[lId]:
        lId += 1
        while lines[lId][0] == '_':
            labels.append(lines[lId].split()[0])
            lId += 1

        datalist = [lines[k].split() for k in range(lId, nLines)]
        data = pd.DataFrame(datalist, columns=labels)
    else:
        values = list()
        while lId < nLines and lines[lId][0] == '_':
            label, value = lines[lId].split()
            labels.append(label)
            values.append(value)
            lId += 1
        data = pd.DataFrame(1, columns=['value'], index=labels)
        data['value'] = values

    return data


def readStarFile(fName, mytable=''):

    with open(fName, 'r') as file:
        lines = file.read().splitlines()

    # Remove empty lines
    lines = [line for line in lines if line]
    # Remove single space lines
    lines = [line for line in lines if line not in ' ']
    # Remove comment lines
    lines = [line for line in lines if line[0] not in '#']
    nLines = len(lines)

    # Get tables
    tables = list()
    tLines = list()  # lines where tables start
    for k, line in enumerate(lines):
        if 'data_' in line[:10]:
            tables.append(line[5:])
            tLines.append(k)

    nTables = len(tables)
    tRanges = [None]*nTables

    for k in range(nTables-1):
        tRanges[k] = (tLines[k] + 1, tLines[k+1])
    tRanges[-1] = (tLines[-1] + 1, nLines)

    data = dict()

    for k, table in enumerate(tables):
        if np.diff(tRanges[k]) > 1:
            data[table] = getTable(lines[slice(*tRanges[k])])
        else:
            data[table] = ''

    if mytable:
        return data[mytable]
    else:
        return data


def writeStarFile(fName, data, tableName=''):

    if isinstance(data, dict):
        mydata = data
    elif isinstance(data, pd.core.frame.DataFrame):
        mydata = {}
        mydata[tableName] = data

    file = open(fName, 'w')

    for tName, table in mydata.items():

        file.write(f'\n\n# reliontomotools {__version__}\n\n')
        file.write(f'data_{tName}\n')

        ncol = len(table.columns)
        if ncol == 1 and not isinstance(table.index[0], int):

            file.write('\n')

            for label, value in table.itertuples():
                file.write(f'{label}\t{value}\n'.expandtabs(60))

            file.write('\n')

        else:  # regular table

            labels = table.columns

            file.write('\nloop_\n')

            for k, label in enumerate(labels):
                file.write(f'{label} #{k+1}\n')

            for k, line in table.iterrows():
                linesep = "\t".join(line)
                file.write(f'{linesep}\n')

    file.close()


class WarpXMLHandler():

    def __init__(self, xmlFname):

        self._xmlFname = xmlFname
        file = io.open(xmlFname, 'r', encoding='utf-16-le')
        self.data = xmltodict.parse(file.read())['TiltSeries']

    @property
    def params(self):
        return list(self.data)

    def getSimpleAttrib(self, param):

        if param in ['@AreAnglesInverted', '@UnselectFilter']:
            return self.getBoolAttrib(param)
        elif param in ['@PlaneNormal']:
            return self.getParamArray(param, ',')
        else:
            return float(self.data[param])

    def getBoolAttrib(self, param):

        return self.data[param] == 'True'

    def getParamArray(self, param, sep='\n'):

        items = self.data[param].split(sep)

        if items[0] in ['True', 'False']:
            return np.array(items) == 'True'
        elif items[0][0] in [str(x) for x in set('-0123456789.')]:
            return np.array(items, dtype=np.float)
        else:
            return items

    def getParam(self, param):

        if param in self.data:
            if param[0] == '@':
                return self.getSimpleAttrib(param)
            else:
                return self.getParamArray(param)
        else:
            if param[0] != '@':
                return self.getParam('@'+param)
            else:
                raise Exception(f'Parameter {param} not present in '
                                '{self._xmlFname}.')

    def getParamGrid(self, param):

        gridParam = self.data[param]
        w = int(gridParam['@Width'])
        h = int(gridParam['@Height'])
        d = int(gridParam['@Depth'])

        if '@Duration' in gridParam:
            t = int(gridParam['@Duration'])
            dims = 4
            shape = (w, h, d, t)
        else:
            dims = 3
            shape = (w, h, d)

        if np.prod(shape) == 1:
            gridParamList = [gridParam['Node']]
        else:
            gridParamList = gridParam['Node']

        gridData = np.empty(shape)

        for item in gridParamList:
            x = int(item['@X'])
            y = int(item['@Y'])
            z = int(item['@Z'])
            v = float(item['@Value'])
            pos = (x, y, z)

            if dims > 3:
                tw = int(item['@W'])
                pos = pos + (tw,)

            gridData[pos] = v

        return gridData

    def getParamCTF(self, param):

        paramsCTF = self.data['CTF']['Param']

        for item in paramsCTF:
            if item['@Name'] == param:
                return float(item['@Value'])

        raise Exception(f'CTF param {param} not present in '
                        'file {self._xmlFname}')

    @property
    def pixelSize(self):
        return self.getParamCTF('PixelSize')

    @property
    def Amplitude(self):
        return self.getParamCTF('Amplitude')

    @property
    def Cs(self):
        return self.getParamCTF('Cs')

    @property
    def Voltage(self):
        return self.getParamCTF('Voltage')

    @property
    def fracDose(self):

        fracDoseV = np.diff(np.sort(self.Dose))
        p = np.flatnonzero(np.isclose(np.diff(fracDoseV), 0))[0]

        return fracDoseV[p]

    @property
    def tomoName(self):
        tiltFname = self.MoviePath[0]
        p = tiltFname.rfind('_')

        return tiltFname[:p]

    def __getattr__(self, name):

        return self.__getitem__(name)

    def __getitem__(self, name):

        # try:
        if 'Grid' in name:
            return self.getParamGrid(name)
        else:
            return self.getParam(name)
        # except Exception as ins:
            # print(ins)


def cleanDir(path):

    files = glob(os.path.join(path, '*'))

    for f in files:
        try:
            os.remove(f)
        except OSError as e:
            print("Error: %s : %s" % (f, e.strerror))


def motlToRelion4(motl, pixelSizeBin1, binning=1, tomoLabel='TS_',
                  padZeros=-1):

    nPart = len(motl)

    labels = ['_rlnTomoName', '_rlnTomoParticleId', '_rlnTomoManifoldIndex',
              # '_rlnMicrographName',
              '_rlnCoordinateX', '_rlnCoordinateY', '_rlnCoordinateZ',
              '_rlnOriginXAngst', '_rlnOriginYAngst', '_rlnOriginZAngst',
              '_rlnAngleRot', '_rlnAngleTilt', '_rlnAnglePsi',
              # '_rlnGroupNumber',
              '_rlnClassNumber']
    datapart = pd.DataFrame(index=range(nPart), columns=labels)

    if padZeros < 0:
        tomoNumW = len(str(int(motl[:, 6].max())))
    else:
        tomoNumW = padZeros

    for k in range(nPart):
        dOnePart = datapart.loc[k]
        motlOne = motl[k, :]
        partIdx = int(motlOne[3])
        tomoIdx = int(motlOne[6])
        tomoName = f'{tomoLabel}{tomoIdx:0{tomoNumW}d}'
        dOnePart['_rlnTomoName'] = f'{tomoName}'
        dOnePart['_rlnTomoParticleId'] = f'{partIdx}'
        dOnePart['_rlnTomoManifoldIndex'] = f'{int(motlOne[5])}'
        # dOnePart['_rlnMicrographName'] = f'{tomoName}'
        # dOnePart['_rlnGroupNumber'] = f'{tomoIdx}'
        dOnePart['_rlnClassNumber'] = f'{int(motlOne[19])}'

        coords = motlOne[7:10].astype(float) - 1
        shifts = coords + motlOne[10:13].astype(float)
        shifts = shifts*binning + binning//2
        coords = np.round(shifts)
        shifts = -(shifts - coords)*pixelSizeBin1
        # shifts = -1.*motlOne[10:13].astype(float)*pixelSize

        # Angles applied to rotate part -> ref (Relion ref.system)
        newAngles = -np.array(mat2euler(
                    euler2mat(*(np.pi/180*motlOne[[16, 18, 17]]),
                              'szxz'), 'rzyz'))*180/np.pi

        dOnePart['_rlnCoordinateX'] = f'{int(coords[0])}'
        dOnePart['_rlnCoordinateY'] = f'{int(coords[1])}'
        dOnePart['_rlnCoordinateZ'] = f'{int(coords[2])}'
        dOnePart['_rlnOriginXAngst'] = f'{shifts[0]:.3f}'
        dOnePart['_rlnOriginYAngst'] = f'{shifts[1]:.3f}'
        dOnePart['_rlnOriginZAngst'] = f'{shifts[2]:.3f}'
        dOnePart['_rlnAngleRot'] = f'{newAngles[2]:.2f}'
        dOnePart['_rlnAngleTilt'] = f'{newAngles[1]:.2f}'
        dOnePart['_rlnAnglePsi'] = f'{newAngles[0]:.2f}'

    data = dict()
    data['data_particles'] = datapart
    return data

def getFieldIfExists(dataframe, label):

    if label in dataframe:
        return dataframe[label]
    else:
        return 0

def relion4ToMotl(dataStar, pixelSize, binfactor=1, useTomoAngles=False):

    datapart = dataStar['data_particles']

    binRatio = 1/binfactor

    nPart = len(datapart)
    motl = np.zeros((nPart, 20))

    if useTomoAngles:
        rotLabel  = '_rlnTomoSubtomogramRot'
        tiltLabel = '_rlnTomoSubtomogramTilt'
        psiLabel  = '_rlnTomoSubtomogramPsi'
    else:
        rotLabel  = '_rlnAngleRot'
        tiltLabel = '_rlnAngleTilt'
        psiLabel  = '_rlnAnglePsi'


    for k in range(nPart):
        dOnePart = datapart.iloc[k]
        motlOne = motl[k, :]

        motlOne[3] = int(getFieldIfExists(dOnePart, '_rlnTomoParticleId'))
        motlOne[5] = int(getFieldIfExists(dOnePart, '_rlnTomoManifoldIndex'))
        tomoIdx = np.array(int(re.findall('_\d+.', dOnePart['_rlnTomoName'])[0][1:]) )
        motlOne[6] = tomoIdx
        motlOne[0] = float(getFieldIfExists(dOnePart, '_rlnMaxValueProbDistribution'))
        motlOne[19] = int(getFieldIfExists(dOnePart, '_rlnClassNumber'))

        coords = dOnePart[['_rlnCoordinateX',
                           '_rlnCoordinateY',
                           '_rlnCoordinateZ']].\
                            values.astype(float)*binRatio
        # Subtom shifts are defined opposite to Relion and in pixels
        shifts = -dOnePart[['_rlnOriginXAngst',
                           '_rlnOriginYAngst',
                           '_rlnOriginZAngst']].\
                            values.astype(float)/pixelSize*binRatio

        # Angles applied to rotate part -> ref (Relion ref.system)
        angles = dOnePart[[psiLabel,
                           tiltLabel,
                           rotLabel]].\
                            values.astype(float)

        newAngles = np.array(mat2euler(
                    euler2mat(*(-np.pi/180*angles),
                              'rzyz'), 'szxz'))*180/np.pi

        # Subtom indexes start from 1
        shifts += coords + 1
        coords = np.round(shifts)
        shifts -= coords

        motlOne[7:10] = coords
        motlOne[10:13] = shifts
        motlOne[16:19] = newAngles[[0, 2, 1]]

    return motl
