#!/usr/bin/python3

'''This loads the data from the beta alignment into the db'''

import psycopg2
import time
import sys

import em170428.runinfo

ri = em170428.runinfo.RunInfo()

db = psycopg2.connect(database='align170428')

with db:
    with db.cursor() as c:
        c.execute('drop table if exists runs')
        c.execute('create table runs ( r integer, M integer, S integer, z0 integer )')
        for r0 in range(ri.runCount()):
            r = r0 + 1
            x,y,z = ri.tileCenter(r,0,0)
            c.execute('insert into runs (r,M,S,z0) values (%s,%s,%s,%s)',
                      (r, ri.montageCount(r), ri.sliceCount(r),int(z)))
