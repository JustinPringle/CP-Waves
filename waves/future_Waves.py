#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:52:22 2019

@author: justinpringle

future wave simulation

reads in future cp .csv file as dataframe
loads cpstats dicts (stats and ARDict) can look to change these to dfs
simulates N sequences of waves M timesteps long
returns dataframe with following columns:
    
    time | CP | season | hs | tp | direction | seqNumber

"""

import numpy as np
import pandas as pd
import pickle as pc
import scipy.stats as stats
import datetime as dt
import time 
from scipy import interpolate
import time, sys
from IPython.display import clear_output
import os
import glob

##### PROGRESS ############
def update_progress(progress):
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1

    block = int(round(bar_length * progress))

    clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    print(text)
    
def F(t,y,d):
    return d/t**2*(t**-1-1)**(d-1) + y
def Fp(t,d):
    return -1*(2*d/t**3*(t**-1-1)**(d-1)+d*(d-1)*(t**-1-1)**(d-2)/t**4)        
def psi(s,d,copula):
    '''
    psi is the generator function for the archimedian copulas
    '''
    if copula=='N12':
        return (1/s-1)**d
    elif copula =='gumbel':    
        
        return (-np.log(s))**d
    
def psiDash(t,d):
    '''
    psi dash is the derivative of the generator function
    for the N12 copula
    '''
    return -d/(t**2)*((1/t-1)**(d-1))
    
def psiI(t,d,copula):
    if copula =='N12':
        return 1/(1+t**(1/d))
    elif copula=='gumbel':
        return np.exp(-t**(1/d))
    

def bisection(u0,w1,theta,q):
    '''
    Bisection method to find root of dC/dv
    starts by splitting percentiles into 0.5->1 and 0.00001->0.5
    check for roots in these two bins accordingly
    stops when error is 0.1%
    '''
    
    xL = 0.5
    xU = 1
    
    error = np.abs((xU-xL)/(xU+xL))*100
    
    #check initial bin
    fU = F(xU,w1,theta)
    fL = F(xL,w1,theta)
    #change bin if fU and fL are same sign
    if fU*fL>0:
        xU=xL
        xL = 0.0001
        
        
    while error > 0.1:
        xR = (xU+xL)/2
        fU = F(xU,w1,theta)
        fL = F(xL,w1,theta)
        fR = F(xR,w1,theta)
        
        if fL*fR<0:
            xU = xR
        elif fL*fR>0 and fU*fR>0:
            #this seems to be an issue every now and again and to do with very low variate values
            xR = xL
            #print('something is wrong u0:%0.4f, xR:%0.4f, theta:%0.2f'%(u0,xR,theta))
            break
        else:
            xL=xR
        error = np.abs((xU-xL)/(xU+xL))*100
#        print(xU,xL,xR,fU,fL,fR)
        
    return xR

def createCDF(df):
    '''
    create interpolating function for cdf data
    
    groups dataframe by cp and seas
    uses cunnane ranking to generate cdf for empirical data
    interpolate function is created using this data
    '''
    
    #group by cp and seas (drop nans) convert into dict
    to_dict = df.groupby(['cps','season'])[['hs','tp','dir']].apply(lambda g: g.dropna().to_dict(orient='list')).to_dict()
#    add interpolating function for the cdf
    for key in to_dict:
        #get waves take log normlaized value
        hs = np.log(to_dict[key]['hs'])
        #sort the data and convert to percentiles
        sortHs = np.sort(hs)
        ranks = stats.rankdata(sortHs)
        probs = (ranks-0.4)/(len(sortHs)+0.2)
        emp = interpolate.interp1d(sortHs,probs)
        iemp = interpolate.interp1d(probs,sortHs)
        #bounds checking
        maxP = probs[-1]
        minP = probs[0]
        
        to_dict[key]['iemp']=(iemp,minP,maxP)
        to_dict[key]['emp']=emp
        
    return to_dict

def getU(u0,copDict,cop=None):
    '''
    gets the next random variate from the copula dependence for the cp
    '''
    assert cop!=None , 'You need to specify a copula'
    
    #kendalls tau used in all copulas
    tau = copDict['tau'][0]
    if cop == 'N12':
        theta = 2/(3*(1-tau))
        #draw random percentile from uniform distribution
        q = np.random.uniform()
        #all copulas need a generator function, it's inverse and derivative
        #the derivative will help estimate the next u
        #the inverse will give you phi given random variable u0
        w1 = psiDash(u0,theta)/q
        #solve for when F_N12 = 0
        w = bisection(u0,w1,theta,q)
        #now get v
        if psi(w,theta,cop)-psi(u0,theta,cop)<0:
            v = psiI(np.abs(psi(w,theta,cop)-psi(u0,theta,cop)),theta,cop)
        else:
            v = psiI(psi(w,theta,cop)-psi(u0,theta,cop),theta,cop)
#        print(u0,v)
        return v
    elif cop=='gumbel':
        theta = (1-tau)**-1
        q = np.random.uniform()
        if u0==0:
            v=u0
        else:
            w1 = psiI(psi(u0,theta,cop)/q,theta,cop)        
            v = psiI((1-q)*psi(w1,theta,cop),theta,cop)
        return v
    
def exponGeoff(z,ze,l,pe):
    '''
    Gyasi-Aqyei & Pegram 2014
    '''
    return 1- (1-pe)*np.exp((-z+ze)/l)
    
def invExponGeoff(z,ze,pe,l):
    '''
    Gyasi-Aqyei & Pegram 2014
    z - random var [0,1]
    pe - value of largest variable
    ze - rel freq of largest value  
    '''
    return -1*(l*np.log((1-z)/(1-ze)))+pe

def invExpon(r,l):
    return np.log(1-r)/-l    

def getVar(u,cdf,copDict,cp,var):
    '''
    returns the inverse CDF func
    i.e. for a given prob level this function outputs real value
    
    pe - largest value 
    ze - rel freq of largest value
    pm - second largest value
    zm - rel freq of second largest value
    pexp - 'l' for expon distribution 
    '''
    if var =='hs':
        #check for extremes
        pe = copDict['pe'][0]
        ze = copDict['ze'][0]
        pm = copDict['pm'][0]
        zm = copDict['zm'][0]
        pexp = copDict['pexp'][0]
        
        #if extreme
        if u>=cdf['iemp'][2]:
            #now either use geoff or exp tail
            if not np.isnan(pexp):
                #pexp has been calculated
#                print(cp)
                value = invExpon(u,pexp)
            else:
                #use geoff
                l = (-pe +pm)/(np.log((1-ze)/(1-zm)))
                value = invExponGeoff(u,ze,pe,l)
        #below limts
        elif u<cdf['iemp'][1]:
            value = float(np.exp(cdf['iemp'][1]))
        else:
            value = float(np.exp(cdf['iemp'][0](u)))
    else:
        value = np.percentile(cdf[var],u*100)
    return value
if __name__ == '__main__':
    
    #load future cps
    rcp = 26
    cpDf = pd.read_csv('../data/hadleyFuture_rcp%s.csv'%rcp,index_col=0)
    
    #load dataframes:
    #obs data
    cpObs = pd.read_csv('../data/cpWavObs.csv',index_col=0)
    #sim dataframe
    cpSimStats = pd.read_csv('../data/simDf.csv',index_col=0)
    
    #turns out that indexing info from list is faster than dataframe
    #so the approach here is to get data out as lists and loop over it
    #create dicts for stat data with {(cp,seas):{param:[value]}} format and index from these 
    
    #create a cdf dict
    cpCDF = createCDF(cpObs)
    #create dict from copula data
    cpCop = cpSimStats.groupby(['cps','season','copula'])[['tau','pe','pm','ze','zm','pexp']].apply(lambda g: g.to_dict(orient='list')).to_dict()
    numLists = 101
    
    #lists to create res df
    
   
    
    
    tr=0
    for i in range(numLists):
        listNum = []
        timeL = []
        seasL=[]
        hsL = []
        tpL=[]
        dL = []
        #get starting values
        cp0 = cpDf['cps'][0]
        u0 = np.random.uniform()
        hs = getVar(u0,cpCDF[('%s'%cp0,'%s'%cpDf['season'][0])],cpCop[('%s'%cp0,'%s'%cpDf['season'][0],'N12')],cp0,'hs')#['iemp'](u0)
        u0Tp = getU(u0,cpCop[('%s'%cp0,'%s'%cpDf['season'][0],'gumbel')],cop='gumbel')
        tp = getVar(u0Tp,cpCDF[('%s'%cp0,'%s'%cpDf['season'][0])],cpCop[('%s'%cp0,'%s'%cpDf['season'][0],'N12')],cp0,'tp')
        
        #get direction from emp distribution
        d = getVar(np.random.uniform(),cpCDF[('%s'%cp0,'%s'%cpDf['season'][0])],cpCop[('%s'%cp0,'%s'%cpDf['season'][0],'N12')],cp0,'dir')
        #starting cp that is required for the direction simulation
        
        #append initial values
        hsL.append(hs)
        tpL.append(tp)
        dL.append(d)
        listNum.append(i+1)
        
        timeL.extend(tm for tm in cpDf.times)
        seasL.extend(s for s in cpDf.season)
        #loop over CPs
        timeS = time.time()
        for j in range(1,len(cpDf['cps'])):#len(cpDf['cps'])):
        #get CP and Seas#            
            cp = cpDf['cps'][j]            
            seas = cpDf['season'][j] 
            
            timeGuess0 =time.time()
        #get u1|u0 (hs) and get u1|u0 (tp) get direction
            u1Hs = getU(u0,cpCop[('%s'%cp,'%s'%seas,'N12')],cop='N12')
            u1Tp = getU(u0,cpCop[('%s'%cp,'%s'%seas,'gumbel')],cop='gumbel')
        #convert to hs1, Tp1, direction
            hs = getVar(u1Hs,cpCDF[('%s'%cp,'%s'%seas)],cpCop[('%s'%cp,'%s'%seas,'N12')],cp,'hs')
            tp = getVar(u1Tp,cpCDF[('%s'%cp,'%s'%seas)],cpCop[('%s'%cp,'%s'%seas,'N12')],cp,'tp')
            #check steepeness
            sop = 1/7.0
            stp1 = (2*np.pi*hs)/(9.81*tp**2)
            if stp1>sop:
                tp = np.sqrt((2*np.pi*hs)/(sop*9.81))
            #keeping direction constant for CP    
            if cp == cp0:
                d = d
            else:
                d = getVar(np.random.uniform(),cpCDF[('%s'%cp,'%s'%seas)],cpCop[('%s'%cp,'%s'%seas,'N12')],cp,'dir')
            #updates
            cp0=cp
            
            hsL.append(hs)
            tpL.append(tp)
            dL.append(d)
            listNum.append(i+1)
            u0 = u1Hs
            
#            display progress
            bar_length=1
            progress = (tr)/numLists/len(cpDf['cps'])*100
            block = int(round(bar_length * progress))

            if j==len(cpDf['cps'])-1:
                timeE = time.time()
                timD=timeE-timeS
#                timeTr.append(timeE-timeS)
#            if i >0 and j<len(cpDf['cps'])-1:
#                tim
            elif i==0:
                timeE=0
                timD=0
#            if len(timeTr)>0:
#                totTime = numLists*np.mean(timeTr)
#                currentTime = np.sum(timeTr)
#                timeRem = (totTime-currentTime)/60
#            elif len(timeTr)==0 and j<2:
            tr+=1
#            timeGuess1 = time.time()
#            timeTr.append(timeGuess1-timeGuess0)
            minRem = (numLists-i)*(timD)/60#(numLists*len(cpDf['cps'])-tr)*(np.mean(timeTr))
            secRem = (minRem-int(minRem))*60    
#            text = "Progress: [{0}] {1:.1f}%".format( "." * block + "-" * (bar_length - block), progress)
            sys.stdout.write("\r%d%% time remaining %d minutes %d seconds" % (progress,minRem,secRem))
            
            
            
        reDf = pd.DataFrame({'times':timeL,'seas':seasL,'simNo':np.asarray(listNum,dtype=int),
                             'hs':np.asarray(hsL),'tp':np.asarray(tpL),'dir':np.asarray(dL)})    
        reDf.to_csv('../data/temp/hadSim%s_%03d.csv'%(rcp,i+1),float_format='%.2f')
    
    #create and save df
    print('creating DF and cleaning directory')
    path = '../data/temp'
    all_files = glob.glob(os.path.join(path, "*.csv"))
    df_from_each_file = (pd.read_csv(f,index_col=0) for f in all_files)
    concatenated_df   = pd.concat(df_from_each_file, ignore_index=True)
    concatenated_df.to_csv('../data/hadleyFutureWaves_rcp%s.csv'%rcp)
    for f in all_files:
        os.remove(f)
    #
