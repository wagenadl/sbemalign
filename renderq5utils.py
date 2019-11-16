#!/usr/bin/python3

import aligndb

monttbl = 'solveq5mont'
rigidtbl = 'solveq5rigidtile'
elastbl = 'solveq5softelastic'
globtbl = 'q5global'
bboxtbl = 'q5bbox'

X = Y = 684*5 # Full tile size in q5 space!
nz = 200

db = aligndb.DB()
ri = db.runinfo()

bbox = db.sel(f'select x0,y0,x1,y1 from {bboxtbl}')[0]

_glsh = {} # Our collection of global shifts
def globalshift(z0):
    '''GLOBALSHIFT - Global shift to be added to within-subvolume positions.
    x, y = GLOBALSHIFT(z0) returns the shift to be added to subvolume-local
    coordinates to get global coordinates.
    Note: These coordinates are shifted so that they are nonnegative by 
    definition, thanks to the use of the global q5bbox.
    The results are cached in memory, so database lookup is performed only
    the first time a given subvolume is requested.'''
    if z0 in _glsh:
        return _glsh[z0]
    x, y = db.sel(f'select x, y from {globtbl} where z0={z0}')[0]
    x -= bbox[0]
    y -= bbox[1]
    _glsh[z0] = (x, y)
    return x, y

_rtp = {}
def rigidtilepositions(z0, r, s):
    '''RIGIDTILEPOSITIONS - Truly global positions of tiles
    xm, ym = RIGIDTILEPOSITIONS(z0, r, s) returns the truly global
    positions of the topleft corners of each of the tiles in slice s of
    run r according to the subvolume at z0. 
    Note: These coordinates are shifted so that they are nonnegative by 
    definition, thanks to the use of the global q5bbox.
    The results are cached in memory, so database lookup is performed only
    the first time a given slice is requested.'''
    zrs = (z0,r,s)
    if zrs in _rtp:
        return _rtp[zrs]
    
    xm, ym = db.vsel(f'''select mo.x+ri.x, mo.y+ri.y
    from {monttbl} as mo
    inner join {rigidtbl} as ri on mo.z0=ri.z0 and mo.r=ri.r and mo.m=ri.m
    where ri.z0={z0} and ri.r={r} and ri.s={s} 
    order by ri.m''')
    gs = globalshift(z0)
    xm += gs[0]
    ym += gs[1]
    _rtp[zrs] = (xm, ym)
    return xm, ym

def shiftcollection(z0, r, m, s):
    '''SHIFTCOLLECTION - Known shifts according to elastic solution
    x, y, dx, dy = SHIFTCOLLECTION(z0, r, m, s) returns all the known 
    shifts for the given (r,m,s). X, Y are returned as truly global
    coordinates. Shifts DX, DY are relative to rigid placement.'''
    gx, gy = globalshift(z0)
    x, y, dx, dy = db.vsel(f'''select x, y, dx, dy from {elastbl}
    where z0={z0} and r={r} and m={m} and s={s}''')
    x += gx
    y += gy
    return x, y, dx, dy
