import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def normPower(powers): # function to calculate Dr Coggan's NP from a list of powers
    sum4 = 0
    for i in powers:
        sum4 += i**4
    return((sum4/len(powers))**0.25)


def averagePower(powers): # function to calculate AP from a list of powers
    sumpowers = 0
    for i in powers:
        sumpowers += i
    return(sumpowers/len(powers))


# now using moving averages

# function to apply a 30 second rolling average to a list of powers (for Dr Coggan NP)
def cumMovAv(powers, windowSize = 30): # can take a window size (in seconds) but defaults to 30 seconds if not given
    pdpowers = pd.Series(powers)
    
    windows = pdpowers.rolling(windowSize, 1)
    movingAverages = windows.mean()
    
    powersMovingAverage = movingAverages.tolist()
    
    return(powersMovingAverage)


# function to apply a 25 second exponential average to a list of powers (for xPower)
def exMovAv(powers):
    pdpowers = pd.Series(powers)
    
    movingAverages = pdpowers.ewm(span=25,adjust=False).mean()
    
    powersExponentialAverage = movingAverages.tolist()
    
    return(powersExponentialAverage)
