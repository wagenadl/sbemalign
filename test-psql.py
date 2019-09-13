#!/usr/bin/python3

import psycopg2
import time
import sys

n = int(sys.argv[1])

db = psycopg2.connect(database='test')

def insert(q):
    with db:
        with db.cursor() as c:
            c.execute('insert into bar (r) values (%s)', (q,))

t0 = time.time()

for k in range(1000):
    q = k + 1000*n + 10000
    insert(q)
    #time.sleep(.001)

print(time.time() - t0)
db.close()
