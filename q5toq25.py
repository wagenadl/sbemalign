#!/usr/bin/python3

import rawimage
import factory
import os
import config

fac = factory.Factory(20)


def ifile(z):
    z1 = z//100
    z2 = z%100
    return f'{config.root}/q5elastic/Z{z1}/{z2}.jpg'


def odir(z=None):
    od = f'{config.root}/q25elastic'
    if z is None:
        return od
    else:
        return od + f'/Z{z//100}'

def ofile(z):
    return odir(z) + f'/{z%100}.jpg'

def convert(z):
    ifn = ifile(z)
    ofn = ofile(z)
    print(f'Downscaling {ifn}')
    img = rawimage.loadimage(ifn)
    img = rawimage.iscale(img, 5)
    img = rawimage.ipad(img, 4)
    if not os.path.exists(odir(z)):
        os.mkdir(odir(z))
    rawimage.saveimage(img, ofn)

if not os.path.exists(odir()):
    os.mkdir(odir())
for z0 in range(9600):
    z = z0+4
    ofn = ofile(z)
    if not os.path.exists(ofile(z)):
        fac.request(convert, z)
fac.shutdown()
