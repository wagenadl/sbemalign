#!/usr/bin/python3

import aligndb
db = aligndb.DB()

db.exe('create table fixr25slicealignbk as (select * from slicealignq5 where r=25 and s>140)')
db.exe('create table fixr25relmontalignbk as (select * from relmontalignq5 where r=25 and s>140)')
db.exe('create table fixr25relmontattouchbk as (select * from relmontattouchq5 where r=25 and s>140)')
db.exe('create table fixr25optimizebk as (select * from optimizeq5rel where r=25 and s>140)')
db.exe('create table fixr25interrunbk as (select * from interrunq5 where r2=26)')
db.exe('create table fixr25transrunmontbk as (select * from transrunmontq5 where r=26)')


db.exe('delete from slicealignq5 where r=25 and s>140')
db.exe('delete from relmontalignq5 where r=25 and s>140')
db.exe('delete from relmontattouchq5 where r=25 and s>140')
db.exe('delete from optimizeq5rel where r=25')
db.exe('delete from warpq5rundone where r=25')
db.exe('delete from runextentq5 where r=25')
db.exe('delete from interrunq5 where r2=26')
db.exe('delete from transrunmontq5 where r=26')
