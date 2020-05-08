#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 12:42:22 2019

@author: justinpringle

calculate future shoreline positions using Jara (2015)
"""

import pandas as pd
import numpy as np
import Jara_temp as jr
import time
import sys

def calcModShoreline(waves,C,params,sTart):
    '''
    calculates the shoreline position for the modelled wave sequence.
    '''
    a = params[0]
    b=params[1]
    c= params[2]
    sMax = params[3]
    sMin = params[4]
    
    S = np.arange(sMin,sMax,1)
    E = jr.eFunc(a,b,c,S)
#    print(E)
#    print(np.min(E))
    ef = jr.eqmShoreline(sMin,sMax,a,b,c)
    
    tList = []
    s0 = sTart
    for j,hs in enumerate(waves):
        
        if j == 0:
            tList.append(sTart)
            continue
        else:
            energy = 0.8**2*hs**2/4.004**2
#            
            if energy<np.min(E):
                sInf = sMax
#           
            elif energy > np.max(E):
                sInf = sMin
            else:
                sInf = ef(energy)
            
            if s0 > sInf:
                k = C[0]
            else:
                k = C[1]
                
            s1 = jr.shorelineResponse(k,s0,sInf,a=a,b=b,dT=6*3600)[0]
            tList.append(s1)
            s0=s1
            
    return tList

def initJara(paramDict):
    '''
    initializes Jara Model
    
    Parameters
    ----------
    paramDict : dictionary continaing variables
    

    Returns
    -------
    parameters a,b,c,sMax,sMin.
    '''
    sMax=paramDict['sMax']-paramDict['meanShor']
    v = paramDict['v']
    sg = 2.6
    g = 9.81
    d50 = paramDict['d50']/1000
    ht = paramDict['ht']
    B= paramDict['B']
    phi = np.deg2rad(paramDict['phi'])
    xt = paramDict['xt']
    #calc some basic info from assumptions
    da = (g*(sg-1)/(v**2))**(1/3)*d50
    w = v/d50*((10.36**2+1.049*da**3)**0.5-10.36)
    A = 2.25*((w**2/g)**(1/3))
    # print('A: ',A)
    #volume of max S
    V = jr.Vs(sMax,ht,B,xt)        
    hbM = jr.hbMax(ht, phi, A, V, B, xt)#paramDict['hbm']#12.8
    sMin = jr.sMnMx(xt,hbM,A,phi,V,B,ht)

    a,b,c = jr.parabolicParams(sMax,sMin,hbM)
    
    return [a,b,c,sMax,sMin]

if __name__ == '__main__':
    
    #read in the shoreline csv file and extract sMax 0 MSL which is CD +1
    shorDf = pd.read_csv('../data/B6.csv',delimiter=';').replace(9999,np.nan)
    sMax = 75#shorDf.groupby('y')['x'].max()[1] 
    # sMean = shorDf.groupby('y')['x'].mean()[1] 
    meanShor = shorDf.groupby('y').agg('mean')['x'][1]
    
    paramDict = {
        'sMax':shorDf.groupby('y').get_group(1).max()['x'],
        'meanShor':meanShor,
        'v':9.5*10**-7,
        'ht':15,
        'd50':0.22,
        'B':4,
        'phi':20,
        'xt':1050,
        'hbm':11}
    
    params = initJara(paramDict)
    
    # params[-1] = mslDf['x'].min()-meanShor
    
    # params = [a,b,c,sMax,sMin]
    #C should scale with time step
    #for 3hrs O(9X10^-5 = 1/(3*3600))
    #for 6hrs O(5X10^-5 = 1/(6*3600))
    
    #accretion is O(months) - 
    C =[-5e-5,-9e-6]
    #read in future wave dataframes
    rcp = 26
    hadWaves = pd.read_csv('../data/hadleyFutureWaves_rcp%s.csv'%(rcp),index_col=0)
    #get wave info and simulate shoreline positions
    #shoreline spinUp allows model to settle at "egm" shoreline distance
    shorelineSpinUp = calcModShoreline(hadWaves[hadWaves['simNo']==1]['hs'].values,C,
                                 params,-20)
    #loop over lists simulate shorelines
    #assign new column to temp df and save as csv in temp folder
    #combine into one df at end of loop
    numLists=101
    shL = []
    for i in range(numLists):
        tS = time.time()
        sStart = shorelineSpinUp[-1]
        
        shoreline = calcModShoreline(hadWaves[hadWaves['simNo']==i+1]['hs'].values,C,
                                 params,sStart)
        
#        tempDf = hadWaves[hadWaves['simNo']==i].assign(shore=shoreline)
#        tempDf.to_csv('../temp/hadSim%s_%03d.csv'%(rcp,i+1))
        shL.extend(shoreline)
        tE = time.time()
        timD = tE-tS
        minRem = (numLists-i)*(timD)/60
        secRem = (minRem-int(minRem))*60    
        
        progress = i/numLists*100.0
        sys.stdout.write("\r%d%% time remaining %d minutes %d seconds" % (progress,minRem,secRem))
        
    hadWaves = hadWaves.assign(shor=shL)
    hadWaves.to_csv('../data/hadleyFutureWavesShor_rcp%s_Ce_5e5_Ca9e6.csv'%rcp,
                    float_format='%.2f')
#    
    
    
    
    
    
    
    
    
    
    
    
    
    
