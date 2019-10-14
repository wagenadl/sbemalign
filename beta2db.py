#!/usr/bin/python3

'''This loads the data from the beta alignment into the db'''

import psycopg2
import time
import sys
import em170428.runinfo
ri = em170428.runinfo.RunInfo()

db = psycopg2.connect(database='align170428')

def exe(sql, args=None):
    '''EXE - Execute a single SQL statement in its own transaction'''
    with db:
        with db.cursor() as c:
            c.execute(sql, args)

def nofail(sql, args=None):
    '''NOFAIL - Like EXE, but catches exceptions'''
    with db:
        with db.cursor() as c:
            try:
                c.execute(sql, args)
            except Exception as e:
                print(e)
                
def sel(sql, args=None):
    '''SEL - Run a SQL statement in its own transaction and return all
             fetched rows.'''
    with db:
        with db.cursor() as c:
            c.execute(sql, args)
            return c.fetchall()

nofail('drop table info')
exe('create table info (X integer, Y integer, root text )')
exe("insert into info (X, Y, root) values (17100, 17100, '/lsi2/dw/170428')")
        
nofail('drop table runs')
exe('create table runs ( r integer, M integer, S integer, z0 integer )')
with db:
    with db.cursor() as c:
        for r0 in range(ri.runCount()):
            r = r0 + 1
            x,y,z = ri.tileCenter(r,0,0)
            c.execute('insert into runs (r,M,S,z0) values (%s,%s,%s,%s)',
                      (r, ri.montageCount(r), ri.sliceCount(r),int(z)))

res = sel('select * from runs')
    
nofail('drop table betapos')
exe('''create table betapos (
       r integer, m integer, s integer,
       xc float, yc float, z integer )''')
with db:
    with db.cursor() as c:
        for r0 in range(ri.runCount()):
            r = r0 + 1
            for m in range(ri.montageCount(r)):
                print('Working on run', r, 'montage', m)
                for s in range(ri.sliceCount(r)):
                    x,y,z = ri.tileCenter(r,m,s)
                    c.execute('''insert into betapos 
                                 (r,m,s, xc,yc,z)
                                 values(%s,%s,%s,%s,%s,%s)''',
                              (int(r),int(m),int(s),float(x),float(y),int(z)))
