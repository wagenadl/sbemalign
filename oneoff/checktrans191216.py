#!/usr/bin/python3

# Finds subvolumes mucked up by trans pug in matchpointsq5
# We need to find all subvolumes in which _any_ trans involved m1!=m2

import aligndb
import numpy as np
import sys
db = aligndb.DB()
ri = db.runinfo()

for r in range(1,60):
    r2 = r + 1
    z = ri.z0(r2)
    m,m2,x,y,dx,dy,snr = db.vsel(f'''select m,m2,avg(x),avg(y),
    avg(dx),avg(dy),max(snr)
    from transrunmontq5 where r={r} and r2={r2} group by m,m2 order by m,m2''')
    use = snr>0
    minsnr = np.min(snr[use])
    maxsnr = np.max(snr[use])
    if np.any(snr==0):
        zero = '0'
    else:
        zero = ' '
    if minsnr < .2*maxsnr:
        bad = '!!'
    else:
        bad = '  '
    print(f'Z{z-1}:{z} R{r}:{r2}  {bad} {zero}: {minsnr:.2f} — {maxsnr:.2f}')

        

