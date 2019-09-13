#!/usr/bin/python3

dbfn = '/home/wagenaar/realign.db'

## Open the database 
import sqlite3
db = sqlite3.connect(dbfn)
c = db.cursor()

## Create the tables

c.execute('''create table spans (
  id integer primary key,
  r integer,
  s0 integer,
  s1 integer )
''')

c.execute('''create table substacks (
  id integer primary key,
  spanid integer,
  m integer,
  xf,
  stage integer,
  foreign key(spanid) references spans(id)
    on delete cascade
    on update cascade )
''')

c.execute('''create table tiles (
  id integer primary key,
  ssid integer,
  ds integer,
  xf,
  stage integer,
  foreign key(ssid) references substacks(id)
    on delete cascade
    on update cascade )
''')

c.execute('''create table runs (
  r integer primary key,
  S integer,
  M integer,
  xf )
''')

c.execute('''create table info (
  X integer,
  Y integer,
  root text )
''')

c.execute('''create table scales (
  scale integer )
''')

## Insert some information
c.execute('''insert into info (X, Y, root)
  values (17100, 17100, "/lsi2/dw/170428")
''')

c.execute('insert into scales values (4)')
c.execute('insert into scales values (16)')

