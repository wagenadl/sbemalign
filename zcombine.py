#!/usr/bin/python3

import config

root = f'{config.sclroot}/q1pyramid'

scales = [ (4,1), (5,2), (6,3), (7,4), (8,5) ] # (A,B) pairs

import factory
import pathlib
import os
import numpy as np
import rawimage
import aligndb

Zmax = 9604 # max Z at B=0

db = aligndb.DB()

tbl = 'zcombine_done'

def maketable():
    db.exe(f'''create table if not exists {tbl} (
    a integer,
    b integer,
    z integer
    )''')


def droptable():
    db.nofail(f'drop table {tbl}')
    
def abpath(a, b, z=None, y=None, x=None):
    pth = f'{root}/A{a}B{b}'
    if not z is None:
        zlo = z % 100
        zhi = z // 100
        pth += f'/Z{zhi}/{zlo}'
        if y is not None:
            pth += f'/Y{y}'
            if x is not None:
                pth += f'/X{x}.jpg'
    return pth

def inpath(a, z, y, x):
    zlo = z % 100
    zhi = z // 100
    pth = f'{root}/Z{zhi}/{zlo}/A{a}/Y{y}/X{x}.jpg'
    return pth

def zcombine(a, b, x, y, z):
    '''ZCOMBINE - Combines several z tiles
    ZCOMBINE(a, b, x, y, z, b) loads 2^b tiles to construct the
    averaged tile a/b/z/y/x from a/z0..z1/y/x.
    ZCOMBINE is smart enough to divide by less than 2^b if some source
    tiles are missing and doesn't save anything when there are zero source
    tiles.
    n = ZCOMBINE(...) returns the number of tiles, or -1 if previously done.'''
    
    Z = 2**b
    n = 0
    ima = None
    for dz in range(Z):
        z1 = z*Z + dz
        img = rawimage.loadimage(inpath(a, z1, y, x))
        if img is not None:
            if n==0:
                ima = img.astype(np.float32)
            else:
                ima += img
            n += 1
    if n>0:
        ima = (ima/n).astype(np.uint8)
        pathlib.Path(abpath(a, b, z, y)).mkdir(parents=True, exist_ok=True)
        rawimage.saveimage(ima, abpath(a, b, z, y, x))

    return n

def combinelevel(a, b, z):
    print(f'Working on A{a} B{b} Z{z}')
    Ymax = int(np.ceil(125000 / 512)) # max Y tiles at A=0
    Ymax = int(np.ceil(Ymax / 2**a)) # max Y at given A
    Xmax = int(np.ceil(45000 / 512)) # max X tiles at A=0
    Xmax = int(np.ceil(Xmax / 2**a)) # max X at given A
    topany = False
    for y in range(0, Ymax):
        doneany = False
        for x in range(0, Xmax):
            if zcombine(a, b, x, y, z):
                doneany = True
        if doneany:
            topany = True
    return topany

def perhapscombinelevel(a, b, z):
    cnt = db.sel(f'select count(1) from {tbl} where a={a} and b={b} and z={z}')
    if cnt[0][0]==0:
        combinelevel(a, b, z)
        db.exe(f'insert into {tbl} (a, b, z) values ({a}, {b}, {z})')
            
def allcombine(a, b):
    Zmax1 = int(np.ceil(Zmax / 2**b)) # max Z at given B
    fac = factory.Factory(12)
    for z in range(0, Zmax1):
        fac.request(perhapscombinelevel, a, b, z)
    fac.shutdown()

def infofile():
    txt = '''{
    "type": "image",
    "data_type": "uint8",
    "num_channels": 1,
    "scales": ['''

    X = 36864
    Y = 117248
    Z = 9604
    dx = 5.5
    dz = 50
    sep = ""

    B = 0
    for A in range(9):
        txt += sep
        txt += f'''
        {{
        "key": "pub/{A}-{B}",
        "size": [{X},{Y},{Z}],
        "resolution": [{dx}, {dx}, {dz}],
        "voxel_offset": [0,0,0],
        "chunk_sizes": [[512,512,1]],
        "encoding": "jpeg"
        }}'''

    sep = ","
    X = (X+511)//2
    Y = (Y+511)//2
    dx *= 2
    if dx>dz:
        dz *= 2
        B += 1
        Z = (Z+1)//2
    txt += "\n    ]\n";
    txt += "}\n";
    return txt

if __name__ == '__main__':
    maketable()
    for a,b in scales:
        allcombine(a, b)
    
    with open(f'{root}/ng-info.txt', 'w') as f:
        f.write(infofile())
    
