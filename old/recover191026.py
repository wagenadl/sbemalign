#!/usr/bin/python3

import aligndb
db = aligndb.DB()

r,m1,m2,s,ix1,ix2,iy1,iy2 = db.vsel('''select r,m1,m2,s,ix1,ix2,iy1,iy2
                                  from slice5bk''');

N = len(r)

for n in range(N):
    db.exe(f'''delete from slicealignq5 
    where r={r[n]} and m1={m1[n]} and m2={m2[n]} and s={s[n]}
    and ix1={ix1[n]} and ix2={ix2[n]} and iy1={iy1[n]} and iy2={iy2[n]}''')
    
######################################################################

r,m,m2,s,ix,iy = db.vsel('''select r,m,m2,s,ix,iy from relmonttouch5bk''')

N = len(r)
for n in range(N):
    print(n)
    db.exe(f'''delete from relmontattouchq5 
    where r={r[n]} and m={m[n]} and m2={m2[n]} and s={s[n]}
    and ix={ix[n]} and iy={iy[n]}''')
