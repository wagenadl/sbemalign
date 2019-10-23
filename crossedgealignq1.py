
import aligndb
import time
import sys
import traceback

import swiftir
import pyqplot as qp
import numpy as np

import rawimage
import factory

db = aligndb.DB()
ri = db.runinfo()

ixx = [3,  5,  7,  9, 11, 13, 15,
       18, 20, 22, 24, 26, 28, 30 ]
R = 17100//5//5
MARG = (17100 - 33*512)//2

def droptable():
    db.nofail('''drop table crossalignq1''')
    db.nofail('''drop table edgealignq1''')

def maketable():
    db.exe('''create table if not exists crossalignq1 (
    r integer,
    m1 integer,
    m2 integer,
    s integer,
    ii integer,
    x1 float,
    y1 float,
    x2 float,
    y2 float,
    dx float,
    dy float,
    sx float,
    sy float,
    snr float,
    dxb float,
    dyb float,
    sxb float,
    syb float,
    snrb float
    )''')

def maketable():
    db.exe('''create table if not exists edgealignq1 (
    r integer,
    m integer,
    m2 integer,
    s integer,
    ii integer,
    x float,
    y float,
    dx float,
    dy float,
    sx float,
    sy float,
    snr float,
    dxb float,
    dyb float,
    sxb float,
    syb float,
    snrb float
    )''')

def edgecoord(r, m1, m2):
    C = ri.ncolumns(r)
    if m2==m1+C:
        # Two montages top-to-bottom
        y1,y2 = db.vsel(f'''select y1,y2 from slicealignq5pos 
          where r={r} and m1={m1} and m2={m2} order by ix1''')
        # Take average and upscale to Q1
        y1 = int(5 * (np.mean(y1) + 4*R)) - MARG
        y2 = int(5 * np.mean(y2)) - MARG
        return (y1, y2)
    elif m2==m1+1 and m2//C==m1//C:
        # Two montages left-to-right
        x1,x2 = db.vsel(f'''select x1,x2 from slicealignq5pos 
          where r={r} and m1={m1} and m2={m2} order by iy1''')
        # Take average and upscale to Q1
        x1 = int(5 * (np.mean(x1) + 4*R)) - MARG
        x2 = int(5 * np.mean(x2)) - MARG
        return (x1, x2)
    else:
        raise ValueError(f'Request for nonexistent edge: R{r} M{m1}:{m2}')
    
