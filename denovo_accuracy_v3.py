# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 19:37:04 2022

@author: amyk0
"""
#%%
#from pep_matrix_v4 import alpha_only
import re
import numpy as np
import pandas as pd
#%%
file = '4589_EThcD de novo peptides.csv'
trueSeq = 'LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES'

#%%
massDict = {'G':57,
            'A':71,
            'S':87,
            'P':97,
            'V':99,
            'T':101,
            'C':103,
            'L':113,
            'I':113,
            'N':114,
            'D':115,
            'Q':128,
            'K':128,
            'E':129,
            'M':131,
            'H':137,
            'F':147,
            'R':156,
            'Y':163,
            'W':186,
            'm':147            
            }

def SeqCumMass(seq, massDict = massDict):
    # sa = alpha_only(seq)
    # slist = list(sa)
    mOx = re.compile('M\(\+15.99\)')
    seq = mOx.sub('m',seq)
    cFlag = False
    nFlag = False
    cAm = re.compile('\(-.98\)')
    m = cAm.search(seq)
    print(m)
    if m:
        cFlag = True
        seq = cAm.sub('',seq)
    nAc = re.compile('\(\+42.01\)')
    m2 = nAc.search(seq)
    if m2:
        nFlag = True
        seq = nAc.sub('',seq)
    slist = list(seq)
    mlist = []
    for aa in slist:
        mass = massDict[aa]
        mlist.append(mass)
    if cFlag:
        mlist[-1] -= -1
    if nFlag:
        mlist[0] += 42
    sArray = np.array(mlist)
    mCum = np.cumsum(sArray)
    mCum = list(mCum)
    return slist, mCum

def denovoMatch(trueSeq, testSeq):
    trueList, trueCum = SeqCumMass(trueSeq)
    testList, testCum = SeqCumMass(testSeq)
    matchAA = ['-']*len(testList)
    for i in range(len(trueCum)):
        mTest = trueCum[i]
        if mTest in testCum:
            j = testCum.index(mTest)
            aaTrue = trueList[i]
            aaTest = testList[j]
            if aaTrue == aaTest:
                matchAA[j] = aaTest
            elif aaTrue in ['I','L'] and aaTest in ['I','L']:
                matchAA[j] = aaTest
    #trueStr = ''.join(trueList)
    matchStr = '\t'.join(matchAA)
    testStr = '\t'.join(testList)    
    return matchStr, testStr

def confMatch(testSeq, conf, thres):
    lconf = conf.split(' ')
    lconf = [int(x) for x in lconf]
    #print(lconf)
    # https://numpy.org/doc/stable/reference/generated/numpy.where.html
    aconf = np.array(lconf)
    bconf = np.where(aconf >= thres, 1, 0)
    l2conf = list(bconf)
    l2conf = [str(y) for y in l2conf]
    sconf = '\t'.join(l2conf)
    return sconf

#%%
df = pd.read_csv(file)
scanList= list(df["Scan"].values)
pepList = list(df["Peptide"].values)
localConf = list(df["local confidence (%)"].values)
thres = 95

matchList = []
testList = []
confList = []

for i in range(len(pepList)):
    testSeq = pepList[i]
    conf = localConf[i]
    matchStr, testStr = denovoMatch(trueSeq, testSeq)
    matchList.append(matchStr)
    testList.append(testStr)
    sconf = confMatch(testSeq, conf, thres)
    confList.append(sconf)
    
with open(file[:-4] + '.txt', 'w') as f:
    for j in range(len(scanList)):
        f.write(str(scanList[j]) + '\t' + testList[j])
        f.write('\n')
        f.write(str(scanList[j]) + '\t' + matchList[j])
        f.write('\n')
        f.write(str(scanList[j]) + '\t' + confList[j])
        f.write('\n')
        

                


    
    
    
    
            




        
        
    
    


    
    

    
    
    


    
    
    