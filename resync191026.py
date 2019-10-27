#!/usr/bin/python3

import aligndb
db = aligndb.DB()

r_,m_,m2_,s_,ii_ = db.vsel('''select r,m,m2,s,ii from edgealignq1bk''')

r,m1,m2,s,ii = db.vsel('''select r,m1,m2,s,ii from crossalignq1bk''')
N = len(r)

for n in range(N):
    print(n)
    db.exe(f'''insert into edgealignq1bk (select * from edgealignq1 
    where r={r[n]} and m={m1[n]} and m2={m2[n]} and s={s[n]}
    and ii={ii[n]})''')
    db.exe(f'''insert into edgealignq1bk (select * from edgealignq1 
    where r={r[n]} and m2={m1[n]} and m={m2[n]} and s={s[n]}
    and ii={ii[n]})''')
    db.exe(f'''delete from edgealignq1
    where r={r[n]} and m={m1[n]} and m2={m2[n]} and s={s[n]}
    and ii={ii[n]}''')
    db.exe(f'''delete from edgealignq1
    where r={r[n]} and m2={m1[n]} and m={m2[n]} and s={s[n]}
    and ii={ii[n]}''')
    
######################################################################

N = len(r_)
for n in range(N):
    print(n)
    a = m_[n]
    b = m2_[n]
    if a>b:
        a,b = b,a
    db.exe(f'''insert into crossalignq1bk (select * from crossalignq1 where r={r_[n]} and m1={a} and m2={b} and s={s_[n]} and ii={ii_[n]})''')
    
    db.exe(f'''delete from crossalignq1 where r={r_[n]} and m1={a} and m2={b} and s={s_[n]} and ii={ii_[n]}''')
    
