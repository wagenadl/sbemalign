#!/usr/bin/python3

# Calculates global shifts to be added to solveq5mont + sovleq5rigidtile
# to make numbers comparable between all subvolumes

import aligndb
import sys
import os
import numpy as np
import rawimage
import pyqplot as qp
import factory
import swiftir
import cv2

monttbl = 'solveq5mont'
rigidtbl = 'solveq5rigidtile'
outtbl = 'q5global'
outview = 'q5bbox'

X = Y = 684*5 # Full tile size in q5 space!
nz = 200

db = aligndb.DB()
ri = db.runinfo()

def droptable():
    db.nofail(f'drop table {outtbl} cascade')
    db.nofail(f'drop view {outview}')

def createtable():
    db.exe(f'''create table if not exists {outtbl} (
    z0 integer,
    x float,
    y float )''')
    # x, y are the components of accumalated Delta.
    # See p. 1629 of E&R.

    db.exe(f'''create view {outview} as 
    (select min(gl.x + mo.x + ri.x) as x0, 
    min(gl.y + mo.y + ri.y) as y0,
    max(gl.x + mo.x + ri.x)+3420 as x1, 
    max(gl.y + mo.y + ri.y) + 3420 as y1
    from q5global as gl 
    inner join solveq5mont as mo on gl.z0=mo.z0 
    inner join solveq5rigidtile as ri 
    on mo.z0=ri.z0 and mo.r=ri.r and mo.m=ri.m)''')

droptable()
createtable()

x0 = 0
y0 = 0
db.exe(f'insert into {outtbl} (z0, x, y) values (4, 0, 0)')

for z0 in range(104, 9500, 100):
    print(f'Working on {z0}')
    z00 = z0 - 100
    sub1 = f'''select ri.r as r, ri.s as s, 
    avg(mo.x + ri.x) as x, avg(mo.y + ri.y) as y
    from {rigidtbl} as ri 
    inner join {monttbl} as mo on ri.z0=mo.z0 and ri.r=mo.r and ri.m=mo.m
    where ri.z0={z0} group by ri.r, ri.s'''
    sub0 = f'''select ri.r as r, ri.s as s,
    avg(mo.x + ri.x) as x, avg(mo.y + ri.y) as y
    from {rigidtbl} as ri 
    inner join {monttbl} as mo on ri.z0=mo.z0 and ri.r=mo.r and ri.m=mo.m
    where ri.z0={z00} group by ri.r, ri.s'''
    sql = f'''select avg(q1.x-q0.x), avg(q1.y-q0.y) 
    from ({sub1}) as q1 inner join ({sub0}) as q0 on q1.r=q0.r and q1.s=q0.s'''
    x,y = db.sel(sql)[0]
    x0 -= x
    y0 -= y
    db.exe(f'insert into {outtbl} (z0, x, y) values ({z0},{x0},{y0})')
    
