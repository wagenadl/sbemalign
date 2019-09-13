#!/usr/bin/python3

dbfn = '/home/wagenaar/realign.db'

## Open the database 
import sqlite3
db = sqlite3.connect(dbfn)
c = db.cursor()

## Create a test transform
import numpy as np
xf = np.array([[.9,.1,0], [-.1,.9,0], [0,0,1]], dtype='float')
# Note that a numpy float is a 64 bit entity, so takes up 72 bytes

spanid = c.execute('insert into spans (r,s0,s1) values (1,0,10)').lastrowid
ssid = c.execute('insert into substacks (spanid,m) values (?,4)', (spanid,)).lastrowid

c.execute('update substacks set xf=? where id==?', (xf, ssid))
# xf is implicitly converted to bytes

c.execute('select id, spanid, m, xf, stage from substacks where id==?', (ssid,))
row = c.fetchone()
np.frombuffer(row[3], dtype='float').reshape((3,3))
# That fetches the transform and converts it back to a 3x3 matrix
