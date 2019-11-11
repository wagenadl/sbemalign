alter table slicealignq5 add column ii integer;
update slicealignq5 set ii = (ix1-ix2)/4*iy1 + (iy1-iy2)/4*ix1;

alter table slicealignq5 add column x1 integer;
update slicealignq5 set x1 = 684/2 + dx0/2;
alter table slicealignq5 add column y1 integer;
update slicealignq5 set y1 = 684/2 + dy0/2;
alter table slicealignq5 add column x2 integer;
update slicealignq5 set x2 = 684/2 - dx0/2;
alter table slicealignq5 add column y2 integer;
update slicealignq5 set y2 = 684/2 - dx0/2;

alter table slicealignq5 drop column ix1;
alter table slicealignq5 drop column iy1;
alter table slicealignq5 drop column ix2;
alter table slicealignq5 drop column iy2;
alter table slicealignq5 drop column dx0;
alter table slicealignq5 drop column dy0;
