#!/usr/bin/env python
import math
import logging
from scipy.stats import binom_test

root_logger = logging.getLogger("")
debug = root_logger.debug

def binomial(n,p):
    #calculate the expected maximum number of replicated reads at a single position

    x = 1
    pvalue = 0
    while (binom_test(x,n,p) > 0.00001):
        x = x + 1
    if x >1:
        x= x - 1
    return x 


def median(list):
    '''    Will return the median of the list of numbers '''

    list.sort()
    if len(list) % 2 == 0:
        med = (list[len(list)/2]+list[len(list)/2-1])/2
    else:
        med = list[(len(list)-1)/2]

    return med

