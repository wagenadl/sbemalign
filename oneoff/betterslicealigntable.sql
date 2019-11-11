alter table slicealignq5 add column ii integer;
update slicealignq5 set ii = (ix1-ix2)/4*iy1 + (iy1-iy2)/4*ix1;

alter table slicealignq5 add column x1 integer;
update slicealignq5 set x1 = 684/2 + dx0/2;
alter table slicealignq5 add column y1 integer;
update slicealignq5 set y1 = 684/2 + dy0/2;
alter table slicealignq5 add column x2 integer;
update slicealignq5 set x2 = 684/2 - dx0/2;
alter table slicealignq5 add column y2 integer;
update slicealignq5 set y2 = 684/2 - dy0/2;

alter table slicealignq5 drop column ix1;
alter table slicealignq5 drop column iy1;
alter table slicealignq5 drop column ix2;
alter table slicealignq5 drop column iy2;
alter table slicealignq5 drop column dx0;
alter table slicealignq5 drop column dy0;

--- I should have included subtile count into x1
--- So here goes:

--- Side by side:
update slicealignq5 set x1 = x1 + 4*684 where r>=3 and r<=54 and m2=m1+1;
update slicealignq5 set y1 = y1 + ii*684 where r>=3 and r<=54 and m2=m1+1;
update slicealignq5 set y2 = y2 + ii*684 where r>=3 and r<=54 and m2=m1+1;
--- Vertical
update slicealignq5 set x1 = x1 + ii*684 where r<=2 or r>=55 or m2>m1+1;
update slicealignq5 set x2 = x2 + ii*684 where r<=2 or r>=55 or m2>m1+1;
update slicealignq5 set y1 = y1 + 4*684 where r<=2 or r>=55 or m2>m1+1;

--- BUGFIX ---
--- I have to reconstruct y2 from ii because I had incorrectly written
--- update slicealignq5 set y2 = 684/2 - dx0/2;
--- Fortunately, I know that _always_ dy0=0 for side by side and dy0=342
--- for vertical in the present state of the db.
--- Side by side
update slicealignq5 set y2 = ii*684 + 342  where r>=3 and r<=54 and m2=m1+1;
--- Vertical
update slicealignq5 set y2 = 684/4 where r<=2 or r>=55 or m2>m1+1;
