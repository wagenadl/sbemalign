#!/usr/bin/python3

import aligndb
import numpy as np

monttbl = 'solveq5mont'
rigidtbl = 'solveq5rigidtile'
elastbl = 'solveq5elastic'
globtbl = 'q5global'
bboxtbl = 'q5bbox'

X = Y = 684*5 # Full tile size in q5 space!

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

def globalbbox():
    '''GLOBALBBOX - Overall bounding box of imagery
    x0, y0, x1, y1 = GLOBALBBOX() returns the overall bounding box
    for all images at Q=5.
    By definition, x0 = y0 = 0.
    x1 and y1 are rounded up to integer.'''
    x0, x1, y0, y1 = db.sel(f'''select min(x), max(x), min(y), max(y)
    from (select gs.x+mo.x+ri.x as x, gs.y+mo.y+ri.y as y
    from {globtbl} as gs
    inner join {monttbl} as mo on gs.z0=mo.z0
    inner join {rigidtbl} as ri on mo.z0=ri.z0 and mo.r=ri.r and mo.m=ri.m
    ) as tbl''')[0]
    x0 -= bbox[0]
    y0 -= bbox[1]
    x1 -= bbox[0]
    y1 -= bbox[1]
    return (int(np.floor(x0)), int(np.floor(y0)),
            int(np.ceil(x1)) + X, int(np.ceil(y1)) + Y)

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

def whence(z):
    '''WHENCE - Which subvolumes to use to render a given z-plane
    zz0, ww = WHENCE(z) returns a vector of z0s of subvolumes and 
    a vector of weights for each. (The sum of all W's is one.)'''
    if z < 104:
        # In the first slices, we have no overlap
        return np.array([4]), np.array([1])
    if z >= 9504:
        return np.array([9404]), np.array([1])
    zz0 = 4
    dz = 100
    zp = z - zz0 # Z relative to first rendered z
    k = zp // dz # "Number" of subvolume
    z0b = 4 + dz*k # ID of subvolume
    b = (z-z0b)/dz
    wb = .5 - .5*np.cos(np.pi*b)
    wa = 1 - wb
    z0a = z0b - dz
    return np.array([z0a, z0b]), np.array([wa, wb])

def _touchlines(z0, r, s):
    xl, yt = rigidtilepositions(z0, r, s)
    xr = xl + X
    yb = yt + Y
    # For X, find the average of all the right edges of column c
    # and the left edges of column c+1, independently for each row.
    # Note that we never have more than 2 columns, but that's irrelevant
    R = ri.nrows(r)
    C = ri.ncolumns(r)
    mm = np.arange(ri.nmontages(r), dtype=int)
    touchx = np.zeros((R,C-1))
    for col in range(C-1):
        for row in range(R):
            m1 = row*C + col
            m2 = row*C + col+1
            touchx[row,col] = (np.mean(xr[m1]) + np.mean(xl[m2]))/2
    # For Y, find the average of all the bottom edges in row r and the
    # top edges of row r+1, irrespective of column.
    touchy = np.zeros(R-1)
    for row in range(R-1):
        istop = (mm//C) == row
        isbot = (mm//C) == row+1
        touchy[row] = (np.mean(yb[istop]) + np.mean(yt[isbot]))/2
    return touchx, touchy

def touchlines(r, s):
    '''TOUCHLINES - Where montages touch
    xx, yy = TOUCHLINES(r, s) returns the global coordinates where
    rendering of montages in slice S of run R should border. In the horizontal
    direction, XX will be a [R x C-1]-matrix representing the edges between
    columns c and c+1, independently for each row in the slice. In the
    vertical direction, YY will be a [R-1]-vector representing the edges
    between rows r and r+1 in the slice. Note that all columns share 
    common YY edge positions.'''
    z = ri.z(r, s)
    touchx = 0
    touchy = 0
    zz0, ww = whence(z)
    for n in range(len(zz0)):
        tx, ty = _touchlines(zz0[n], r, s)
        touchx = touchx + ww[n]*tx
        touchy = touchy + ww[n]*ty
    return touchx.astype(int), touchy.astype(int)

def rigidtileposition(r, m, s, collapse=True):
    '''RIGIDTILEPOSITION - Position for a single tile
    x, y = RIGIDTILEPOSITION(r, m, s) returns global top-left coordinates
    for the given tile using a weighted average of each of the subvolumes
    that it is part of.
    x, y, w = RIGIDTILEPOSITION(r, m, s, collapse=False) returns
    global top-left coordinates for the given tile according to each
    of the subvolumes that it is part of, along with weights for those
    subvolumes.

    '''
    z = ri.z(r, s)
    zz0, ww = whence(z)
    N = len(zz0)
    x = np.zeros(N)
    y = np.zeros(N)
    for n in range(N):
        lxx, tyy = rigidtilepositions(zz0[n], r, s)
        x[n] = lxx[m]
        y[n] = tyy[m]
    if collapse:
        x = np.sum(x*ww) / np.sum(ww)    
        y = np.sum(y*ww) / np.sum(ww)
        return x, y
    else:
        return x, y, ww

def renderlimits(r, m, s):
    '''RENDERLIMITS - Boundaries for rendering of tiles in q5elastic
    xl, yt, xr, yb = RENDERLIMITS(r, m, s) returns the limits [xl, xr) and
    [yt, yb) for rendering the given tile.'''

    # Our aim is to find the common area covered by all rigid positionings
    # of the tile.    
    xx, yy, ww = rigidtileposition(r, m, s, collapse=False)
    lx = int(np.ceil(np.max(xx)))
    ty = int(np.ceil(np.max(yy)))
    rx = int(np.floor(np.min(xx))) + X
    by = int(np.floor(np.min(yy))) + Y

    # Now, let's replace any inner edges with common touch lines 
    C = ri.ncolumns(r)
    R = ri.nrows(r)
    col = m % C
    row = m // C
    tox, toy = touchlines(r, s)
    if col<C-1:
        rx = tox[row, col]
    if col>0:
        lx = tox[row, col-1]
    if row<R-1:
        by = toy[row]
    if row>0:
        ty = toy[row-1]
    return lx,ty,rx,by

def rendergrid(r, m, s):
    '''RENDERGRID - Return grid points in global space for rendering a tile
    xx, yy = RENDERGRID(r, m, s) returns a vector of x-coordinates and
    a vector of y-coordinates in global space that mark a grid and bbox
    of space to be filled by the given tile.'''
    lx, ty, rx, by = renderlimits(r, m, s)
    xrig, yrig, ww = rigidtileposition(r, m, s, collapse=False)
    xrig = int(np.round(np.sum(xrig*ww))) # we know that ww sums to 1
    yrig = int(np.round(np.sum(yrig*ww)))
    N = 7
    xx = np.zeros(N, dtype=int)
    yy = np.zeros(N, dtype=int)
    xx[0] = lx
    xx[-1] = rx
    yy[0] = ty
    yy[-1] = by
    for n in range(N-2):
        xx[n+1] = xrig + (n+.5) * (X//5)
        yy[n+1] = yrig + (n+.5) * (X//5)
    return xx, yy

def _interpolatedshift(z0, r, m, s, x, y, sc):
    '''_INTERPOLATEDSHIFT - Interpolated shifts at specific points
    dx, dy = INTERPOLATEDSHIFTS(z0, r, m, s, x, y, sc) calculates the
    shift for the given tile (r, m, s) at the given point (x, y), specified
    in global coords, according to the subvolume at z0. SC must be result
    from SHIFTCOLLECTION.
    We don't simply interpolate. Rather, we assume that 
    dx = a + b (x' - x) + c (y'- y)
    near the point (x, y) and do a least-squares optimization to find a, b, c.
    For that, we set χ² = sum_k w_k (dx(x_k,y_k) - dx_k)²
    where k iterates over all the known shifts in the collection and w_k 
    is a weight set by w_k = 1 / (1 + Δ²_k / Δ²₀) where Δ_k is the distance
    between the k-th measuring point and the point (x, y), and Δ₀ is a 
    constant.'''
    # Let's first do dx
    Delta0 = 50
    Deltaxk = sc[0] - x
    Deltayk = sc[1] - y
    Deltak2 = (Deltaxk)**2 + (Deltayk)**2
    wk = 1 / (1 + Deltak2 / Delta0**2)
    # ∂χ²/χ²a = 2 sum_k w_k (a + b (x_k - x) + c (y_k - y) - dx_k)
    # ∂χ²/χ²b = 2 sum_k w_k (a + b (x_k - x) + c (y_k - y) - dx_k) (x_k - x)
    # ∂χ²/χ²c = 2 sum_k w_k (a + b (x_k - x) + c (y_k - y) - dx_k) (y_k - y)
    # Solve this as M [a b c]' = B. Drop the factors two
    M = np.zeros((3,3))
    M[0,0] = np.sum(wk)
    M[0,1] = np.sum(wk*Deltaxk)
    M[0,2] = np.sum(wk*Deltayk)
    M[1,0] = np.sum(wk*Deltaxk)
    M[1,1] = np.sum(wk*Deltaxk**2)
    M[1,2] = np.sum(wk*Deltayk*Deltaxk)
    M[2,0] = np.sum(wk*Deltayk)
    M[2,1] = np.sum(wk*Deltaxk*Deltayk)
    M[2,2] = np.sum(wk*Deltayk**2)
    B = np.zeros(3)
    B[0] = np.sum(wk*sc[2])
    B[1] = np.sum(wk*sc[2]*Deltaxk)
    B[2] = np.sum(wk*sc[2]*Deltayk)
    #print(f'Solving R{r} M{m} S{s} x={x} y={y}')
    #print('M = ', M)
    #print('B = ', B)
    abc = np.linalg.solve(M, B)
    #print(abc)
    dx = abc[0]
    # The equation for dy is exactly the same, except that we need dy_k.
    B[0] = np.sum(wk*sc[3])
    B[1] = np.sum(wk*sc[3]*Deltaxk)
    B[2] = np.sum(wk*sc[3]*Deltayk)
    abc = np.linalg.solve(M, B)
    #print(abc)
    dy = abc[0]
    return dx, dy
    
def interpolatedshifts(r, m, s, xx, yy):
    '''INTERPOLATEDSHIFTS - Interpolated shifts at specific points
    dx, dy = INTERPOLATEDSHIFTS(r, m, s, xx, yy) calculates shifts for 
    the given tile at given points, which are specified in global coords.
    This uses a weighted average over subvolumes.'''
    zz0, ww = whence(ri.z(r, s))
    N = len(zz0)
    scc = []
    for n in range(N):
        scc.append(shiftcollection(zz0[n], r, m, s))
    dx = np.zeros(xx.shape)
    dy = np.zeros(xx.shape)
    for k in range(xx.size):
        x = xx.flat[k]
        y = yy.flat[k]
        ddx = 0
        ddy = 0
        for n in range(N):
            dx1, dy1 = _interpolatedshift(zz0[n], r, m, s, x, y, scc[n])
            #print(dy1, ww[n])
            ddx += ww[n] * dx1
            ddy += ww[n] * dy1
        dx.flat[k] = ddx
        dy.flat[k] = ddy
    return dx, dy

if __name__=='__main__':
    import pyqplot as qp
    r = 3
    m = 2
    s = 50
    
    zz0, ww = whence(ri.z(r, s))
    z0 = zz0[0]
    x, y, dx, dy = shiftcollection(z0, r, m, s)
    qp.figure('/tmp/s1', 10, 10)
    N = len(x)
    qp.marker('o')
    qp.mark(x, -y)
    F = 500
    qp.pen('k', 1)
    for n in range(N):
        qp.plot([x[n], x[n]+F*dx[n]], -np.array([y[n], y[n]+F*dy[n]]))

    if len(zz0)>1:
        z0 = zz0[1]
        x, y, dx, dy = shiftcollection(z0, r, m, s)
        N = len(x)
        qp.pen('777')
        qp.marker('o')
        qp.mark(x, -y)
        F = 500
        qp.pen(width=1)
        for n in range(N):
            qp.plot([x[n], x[n]+F*dx[n]], -np.array([y[n], y[n]+F*dy[n]]))
        
    xx, yy = rendergrid(r, m, s)
    N = len(xx)
    xxx = np.repeat(np.reshape(xx,(1,N)), N, 0) 
    yyy = np.repeat(np.reshape(yy,(N,1)), N, 1) 
    dx, dy = interpolatedshifts(r, m, s, xxx, yyy)
    qp.pen('900', 1)
    qp.mark(xxx, -yyy)
    for n in range(xxx.size):
            qp.plot([xxx.flat[n], xxx.flat[n]+F*dx.flat[n]],
                    -np.array([yyy.flat[n], yyy.flat[n]+F*dy.flat[n]]))


    '''
    qp.pen('k', 0)
    qp.xaxis(y=np.min(-yyy)-50, lim=[np.min(xxx), np.max(xxx)])
    qp.yaxis(x=np.min(xxx)-50, lim=[np.min(-yyy), np.max(-yyy)])
    '''
    qp.shrink(1,1)
