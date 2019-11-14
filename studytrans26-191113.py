#!/usr/bin/python3

import aligndb
import swiftir
import rawimage
import pyqplot as qp
import numpy as np

roughtbl = 'solveq5slice'
intbl = 'interrunq25'
inq = 25
outq = 5

X = Y = 684 # Size of tile
IX = IY = 5

r = 26
m = 6
ix = 3
iy = 4

db = aligndb.DB()
ri = db.runinfo()

img = rawimage.partialq5img(r, m, 0, ix, iy)
x0,y0 = db.sel(f'''select x,y from {roughtbl}
where r={r} and m={m} and s=0''')[0]
s1 = ri.nslices(r-1) - 1
mm1,xx1,yy1 = db.vsel(f'''select m,x,y from {roughtbl} 
where r={r-1} and s={s1}
order by m''')
dx,dy = db.sel(f'select dx+dxb, dy+dyb from {intbl} where r2={r}')[0]
dx *= inq/outq
dy *= inq/outq

# Center of tile in run space
xc = x0 + X*(ix+.5)
yc = y0 + Y*(iy+.5)

# Center of other tiles
xx1c = xx1 + (X*2.5)
yy1c = yy1 + (X*2.5)

# Presumed best montage to match
m1 = mm1[np.argmin((xc+dx - xx1c)**2 + (yc-dy - yy1c)**2)]
# Position in that montage that should match
x1c = xc-dx - xx1[m1]
y1c = yc-dy - yy1[m1]
ix1 = int(x1c/X+.5)
iy1 = int(y1c/Y+.5)
if ix1<=0:
    ix1=1
elif ix1>=5:
    ix1=4
if iy1<=0:
    iy1=1
elif iy1>=5:
    iy1=4
s = ri.nslices(r-1)-1
if r-1 == 35:
  s -= 1

img1 = rawimage.q5subimg2x2(r-1, m1, s, ix1, iy1)

# Position in second image that should match center of first
x1croi = x1c - (ix1-1)*X
y1croi = y1c - (iy1-1)*Y

qp.figure('/tmp/s1', 8, 4)
qp.subplot(1,2,1)
qp.imsc(img, xx=np.arange(0,X), yy=-np.arange(0,Y))
qp.pen('090', 2)
qp.marker('+')
qp.mark(X/2, -Y/2)
qp.xlim(-X/2,3*X/2)
qp.shrink(1,1)
qp.subplot(1,2,2)
qp.imsc(img1, xx=np.arange(0,X*2), yy=-np.arange(0,Y*2))
qp.pen('090', 2)
qp.marker('+')
qp.mark(x1croi, -y1croi)
qp.shrink(1,1)

SIZ = (512, 512)

win1 = swiftir.extractStraightWindow(img, (X/2, Y/2), SIZ)
apo1 = swiftir.apodize(win1)

dx0 = 0
dy0 = 0
cont = True
f = 2

print(f'ix1 {ix1} iy1 {iy1}: x {x1croi} y {y1croi}')


while cont:
    win2 = swiftir.extractStraightWindow(img1, (x1croi+dx0, y1croi+dy0), SIZ)
    apo2 = swiftir.apodize(win2)
    (dx,dy,sx,sy,snr) = swiftir.swim(apo1, apo2)
    print(f'x: {dx0} + {dx}  y: {dy0} + {dy}  snr: {snr}')
    qp.figure(f'/tmp/s{f}', 8, 4)
    f += 1
    qp.subplot(1,2,1)
    qp.imsc(apo1, xx=np.arange(SIZ[0]), yy=-np.arange(SIZ[1]))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(SIZ[0]/2, -(SIZ[1]/2))
    qp.shrink(1,1)
    qp.subplot(1,2,2)
    qp.imsc(apo2, xx=np.arange(SIZ[0]), yy=-np.arange(SIZ[1]))
    qp.marker('x')
    qp.pen('090', 1)
    qp.mark(SIZ[0]/2, -SIZ[1]/2)
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(SIZ[0]/2 + dx, -(SIZ[1]/2 + dy))
    qp.shrink(1,1)
    dx0 += dx
    dy0 += dy

    if dx**2 + dy**2 < 1:
        cont = False

while f<20:
    try:
        qp.close(f'/tmp/s{f}.qpt')
    except:
        break
    f += 1
    

##

win2 = swiftir.extractStraightWindow(img1, (x1croi, y1croi), SIZ)
apo2 = swiftir.apodize(win2)
shf = swiftir.alignmentImage(apo1, apo2, -.65)
qp.figure('/tmp/s99', 5, 5)
qp.imsc(shf, xx=np.arange(512), yy=-np.arange(512))
xy = swiftir.findPeak(shf)
qp.marker('o', 15, fill='brush')
qp.pen('r',2)
x = xy[0]
if x<0:
    x+=512
y = xy[1]
if y<0:
    y+=512
qp.mark(x, -y)
qp.shrink(1,1)
win = swiftir.extractFoldedAroundPoint(shf, xy, 10)
qp.figure('/tmp/s98', 5, 5)
qp.imsc(win)

#

f=1
def recursiveswim(dx0, dy0, snr0):
    global f
    win2 = swiftir.extractStraightWindow(img1, (x1croi+dx0, y1croi+dy0), SIZ)
    apo2 = swiftir.apodize(win2)
    qp.figure(f'/tmp/t{f}', 4, 4)
    qp.imsc(apo2, xx=np.arange(512), yy=-np.arange(512))
    (dx,dy,sx,sy,snr) = swiftir.swim(apo1, apo2)
    print(f'f{f}  x {dx0}+{dx} y {dy0}+{dy}  {snr}')
    if snr<snr0:
        return (dx0,dy0,0,0,snr0)
    qp.pen('r',1)
    qp.marker('+')
    qp.mark(dx+256, -(256+dy))
    f = f+1
    if dx**2+dy**2 < 1:
        return (dx0+dx,dy0+dy,sx,sy,snr)
    res = []
    res.append(recursiveswim(dx0+dx, dy0+dy, snr))
    THR = 175
    if dx<-THR:
        res.append(recursiveswim(dx0+dx+512, dy0+dy, snr))
        if dy<-THR:
            res.append(recursiveswim(dx0+dx+512, dy0+dy+512, snr))
        elif dy>THR:
            res.append(recursiveswim(dx0+dx+512, dy0+dy-512, snr))
    elif dx>THR:
        res.append(recursiveswim(dx0+dx-512, dy0+dy, snr))
        if dy<-THR:
            res.append(recursiveswim(dx0+dx-512, dy0+dy+512, snr))
        elif dy>THR:
            res.append(recursiveswim(dx0+dx-512, dy0+dy-512, snr))
    else:
        if dy<-THR:
            res.append(recursiveswim(dx0+dx, dy0+dy+512, snr))
        elif dy>THR:
            res.append(recursiveswim(dx0+dx, dy0+dy-512, snr))
    best = None
    bests = 0
    for r in res:
        if r[-1]>bests:
            best = r
            bests = r[-1]
    return (dx0+best[0],dy0+best[1], best[2], best[3], best[4])

dx,dy,sx,sy,snr = recursiveswim(0,0,0)
