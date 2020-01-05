#!/usr/bin/python3

import factory
import rawimage
import pathlib
import numpy as np
import config

def abpath(a, b, z, y=None, x=None):
    root = f'{config.sclroot}/q1pyramid'
    zlo = z % 100
    zhi = z // 100
    zbit = f'Z{zhi}/{zlo}'
    pth = root
    if b==0:
        pth += f'/{zbit}/A{a}'
    else:
        pth += f'/A{a}B{b}/{zbit}'
    if y is not None:
        pth += f'/Y{y}'
        if not x is None:
            pth += f'/X{x}.jpg'
    #print(pth)
    return pth

def outpath(a, b, z, y=None, x=None):
    root = f'{config.sclroot}/q1chunks'
    pth = f'{root}/A{a}B{b}/Z{z}'
    if y is not None:
        pth += f'/Y{y}'
    if x is not None:
        pth += f'/X{x}.jpg'
    return pth

TILEIN = 512
TILEOUT = 32
NTILES = TILEIN // TILEOUT

def onechunk(a, b, z, y, x):
    print(f'Working on A{a} B{b} Z{z} Y{y} X{x}')
    yin = y // TILEIN
    xin = x // TILEIN
    imgs = []
    got = False
    for dz in range(TILEOUT):
        img = rawimage.loadimage(abpath(a, b, z+dz, yin, xin))
        imgs.append(img)
        if img is not None:
            got = True
    print(f'Working on A{a} B{b} Z{z} Y{yin} X{xin}: Got = {got}')
    if not got:
        return
    zout = z // TILEOUT
    for dy in range(NTILES):
        yout = y//TILEOUT + dy
        parentdir = outpath(a, b, zout, yout)
        pathlib.Path(parentdir).mkdir(parents=True, exist_ok=True)
        for dx in range(NTILES):
            xout = x//TILEOUT + dx
            ima = np.zeros((TILEOUT, TILEOUT, TILEOUT), dtype=np.uint8) + 128
            for dz in range(TILEOUT):
                if imgs[dz] is not None:
                    ima[dz,:,:] = imgs[dz][TILEOUT*dy:TILEOUT*(dy+1),
                                           TILEOUT*dx:TILEOUT*(dx+1)]
            ima = np.reshape(ima, (TILEOUT*TILEOUT, TILEOUT))
            rawimage.saveimage(ima, outpath(a, b, zout, yout, xout))

def onelevel(a, b, z, X, Y):
    print(f'Working on A{a} B{b} Z{z}')
    for y in range(0, Y, TILEIN):
        for x in range(0, X, TILEIN):
            onechunk(a, b, z, y, x)
    
def perhapsonelevel(a, b, z, X, Y):
    guard = pathlib.Path(outpath(a, b, z) + '/complete')
    if not guard.exists():
        onelevel(a, b, z, X, Y)
        guard.touch()
            
def chunkify(a, b, X, Y, Z):
    fac = factory.Factory(12)
    for z in range(0, Z, TILEOUT):
        fac.request(perhapsonelevel, a, b, z, X, Y)

def chunkifyall():
    X = 36864
    Y = 117248
    Z = 9604
    dx = 5.5
    dz = 50
    
    b = 0
    for a in range(9):
        chunkify(a, b, X, Y, Z)
        X = (X+511)//2
        Y = (Y+511)//2
        dx *= 2
        if dx>dz:
            dz *= 2
            b += 1
            Z = (Z+1)//2
        
if __name__ == '__main__':
    chunkifyall()
    
