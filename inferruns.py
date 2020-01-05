#!/usr/bin/python3

import aligndb
import config
import os

def configureRun(r, z0):
    S = 0
    while True:
        if os.path.exists(config.rawtile(r, 0, S)):
            S = S + 1
        else:
            break
    M = 0
    while True:
        if os.path.exists(config.rawtile(r, M, 0)):
            M = M + 1
        else:
            break
    print(f'Run {R} Montages {M} Slices {S} Z0 {z0}')
    return S

r = 1
z0 = 0
            
while True:
    if os.path.exists(config.rawtile(r, 0, 0)):
        z0 += configureRun(r, z0)
        r = r + 1
    else:
        break
