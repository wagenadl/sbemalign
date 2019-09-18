#!/usr/bin/python3

# This calculates the alignment wrt uct

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

exe('''create table if not exists uctalign (
         r integer, m integer, s integer,
