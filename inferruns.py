#!/usr/bin/python3

import aligndb
import config
import os

db = aligndb.DB()

def runInfo(r):
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
    return (M, S)

def createTable():
    db.exe('drop table if exists runs')
    db.exe('create table runs ( r integer, M integer, S integer, z0 integer )')

def configureRun(r, M, S, z0):
    print(f'Run {r} Montages {M} Slices {S} z0 {z0}')
    db.exe(f'insert into runs (r,M,S,z0) values ({r}, {M}, {S}, {z0})')

#################

createTable()

r = 1
z0 = 0
            
while True:
    if os.path.exists(config.rawtile(r, 0, 0)):
        M, S = runInfo(r)
        configureRun(r, M, S, z0)
        z0 += S
        r = r + 1
    else:
        break

