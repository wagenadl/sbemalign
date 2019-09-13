#!/usr/bin/python3

dbfn = '/home/wagenaar/realign.db'

## Open the database 
import sqlite3
db = sqlite3.connect(dbfn)

import time
import sys
n = int(sys.argv[1])
spc = ''
for k in range(n):
    spc += '   '

def insert(q):
    try:
        with db:
            db.execute('insert into spans (r) values (?)', (q,))
    except:
        print('Failed to insert', q)

t0 = time.time()
for k in range(1000):
    q = k + 1000*n
    # print(spc, q)
    insert(q)
    time.sleep(.001)
    
print(time.time() - t0)
