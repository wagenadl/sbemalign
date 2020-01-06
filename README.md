The SBEMALign Alignment Process
=====================

Introduction
------------

Alignment refers to the task of taking individual images from an SBEM
data acquisition volume, determining their relative positions in
space, and creating a unified 3D volume combining all the images.
SBEMAlign includes the ability to (gently) warp images to obtain a better fit.

Terminology
-----------

Tile
: a source image at original or reduced scale

Source tile
: a source image at original scale

Scaled tile
: a source image at reduced scale

Run
: a single run of the microscope

Montage
: all tiles within a run taken at nominally the same x-y-position

Slice
: all tiles within a run taken at nominally the same z-position

Subvolume
: Any number of consecutive slices, not necessarily all within a run

Submontage
: The intersection of a montage with a subvolume

Model space
: coordinate space for global reconstruction

Run space
: coordinate space for a run; defined relative to the model space by
way of a translation

Montage space
: coordinate space for a montage; defined relative to the run space by
way of a translation

Tile space
: coordinate space for a tile; defined relative to the enclosing
montage by way of a grid of matching points that define subtiles


Creating the database
---------------------

SBEMAlign stores information in a PostgreSQL database. If you are using Ubuntu Linux, installing PostgreSQL is as easy as

    sudo apt install postgresql

You then need to create a database and give yourself access to it. Typically, this involves something like:

    psql
    => CREATE TABLE align170428;
    => GRANT ALL PRIVILEGES ON DATABASE align170428 TO wagenaar;

where I assumed that you call your database “align170428” and that your username 
is “wagenaar”. You may first have to create that user in PostgreSQL. There are 
lots of tutorials on the web to help you with these steps.

Importing basic information about your project
------------------------------------

To use the tools, you first need to edit the “config.py” file to suit
your needs. In particular, you need to define where files are stored and what 
database to use.

Most likely, you will need to change all the pathnames near the top of the file.
You will
also need to change the name of the database and database user.

Next, you need to edit the function “rawtile” to match the file naming convention of your microscope. If you look at my version, you will see that our convention is rather curious and not always consistent. Hence the flexibility of a function.

Finally, you need to specify the number of columns in each run of the microscope. SBEMAlign assumes that each run involved a simple grid of images organized into columns and rows. If you specify the number of columns, it will figure out the number of rows automatically.

When you are done with the “config.py” file, you can run “inferruns.py” to let the software figure out how many montages and slices there are in each of your runs. It is
well worth checking that the output of this program matches your expectations. If it
gets things wrong, nothing else is going to work. Most likely, you can fix problems
by adjusting the “rawtile” function in “config.py” and rerunning “inferruns.py.”

Currently meaningful programs and tables
----------------------------------------

- slicealignq5
- relmontalignq5
- relmontattouchq5
- solveq5slice
- (renderq25slice)
- interrunq25
- transrunmontq5
- solveq5mont
- solveq5rigidtile
- solveq5elastic

Details of database tables
--------------------------

These are presented here in order of construction hierarchy, toplevel
tables first.

### solveq5elastic

Solved elastic shifts inside tiles

#### Members

z0
: subvolume ID

r, m, s
: run-montage-slice for a tile

x, y
: coordinates of test point relative to subvolume

dx, dy
: shift for that test point

### solveq5rigidtile

Solved rigid positions of tiles in a subvolume.

#### Members
z0
: subvolume ID

r, m, s
: run-montage-slice for a tile

x, y
: position of that tile relative to (r, m) position from solveq5mont

### solveq5mont

Solved positions of montages in a subvolume.

#### Members

z0
: subvolume ID

r, m
: run and montage

x, y
: position of topleft corner of the montage relative to the subvolume


### transrunmontq5

Relative alignment of subtile centers in one run relative to other run.
In the end, pixel
  (X + (ix+.5)*684, Y + (iy+.5)*684) in (r, m, s)
is matched to
  (X + x + dx + dxb, Y + y + dy + dyb) in (r2, m2, s2).

#### Members:

r, m, s
: Primary tile

ix, iy
: Subtile position in primary tile

r2, m2, s2
: Secondary tile

x, y
: Base position in secondary tile

dx, dy, sx, sy, snr
: Results of first SWIM

dxb, dyb, sxb, syb, snrb
: Results of second SWIM

#### Dependencies:

- solveq5slice
- interrunq25

### interrunq25

Relative alignment of two runs at scale 1:25. We align the last slice of r1
and the first slice of r2, which must have been rendered using renderq25slice
according to the results of solveq5slice. In the end, pixel location
  (X - dx/2 - dxb/2, Y - dy/2 - dyb/2) in (r1, sfinal)
matches to
  (X + dx/2 + dxb/2, Y + dy/2 + dyb/2) in (r2, sfirst).
Coordinates are stored at 1:25 scale, not at 1:5 scale!

#### Members:

r1, r2
: The runs involved. Invariant: r1 + 1 = r2

dx, dy, sx, sy, snr
: Results of first SWIM

dxb, dyb, sxb, syb, snrb
: Results of second SWIM

#### Dependencies:

- solveq5slice
- images rendered by renderq25slice


### solveq5slice

This optimizes the relative positions of all the tiles in a slice.
It does not optimize the positions of individual points.
See E&R p. 1611ff.
This is used for transrunmontq5 and renderq25slice.

#### Members:
r, m, s
: run, montage, slice number
x, y
: position of tile relative to slice

#### Constructed by:

solveq5slice.py

#### Dependencies:

- slicealignq5

### slicealignq5

Alignment between tiles within a slice

### Members:

r, s
: run and slice ID

m1, m2
: montage IDs within the slice

ix1, iy1
: subtile IDs within the first montage (counting 0..4, inclusive)

ix2, iy2
: subtile IDs within the second montage (counting 0..4, inclusive)

dx, dy, sx, sy, snr
: results from primary “swim” run

dxb, dyb, sxb, syb, snrb
: results from secondary “swim” run

dxc, dyc, sxc, syc, snrc
: results from tertiary “swim” run

dx0, dy0
: topleft position of primary ROI in first subtile (see below)

### Notes:

- For a top-to-bottom overlap (horizontal strip of overlap), y1=4, y2=0, and
ix1=ix2.
- For a left-to-right overlap (vertical strip of overlap), x1=4, x2=0, and
iy1=iy2.
- For top-to-bottom overlap, the primary alignment is done in a 684x342
rectangle centered at (342, 342+171) in the (ix1,iy1)-th subtile of the
first montage and at (342, 342-171) in the (ix2,iy2)-th subtile of the
second montage.
- For left-to-right overlap, the primary alignment is done in a 342x684
rectangle centered at (342+171, 342) in the (ix1,iy1)-th subtile of the
first montage and at (342-171, 342) in the (ix2,iy2)-th subtile of the
second montage.
- To make this math slightly less painful, the center can be calculated
as p1 ≡ (X/2 + dx0/2, Y/2 + dy0/2) in the first montage and
p2 ≡ (X/2 - dx0/2, Y/2 - dy0/2) in the second, where X = Y = 684, and the size
of the region can be calculated as (X-dx0) x (Y-dy0).
- The secondary overlap is done in a half-sized rectangle centered at
p1 - (dx,dy)/2 in the first montage and at p2 + (dx,dy)/2 in the second.
- The tertiary overlap is done in a half-sized rectangle centered at
p1 - (dx+dxb, dy+dyb)/2 in the first montage and at p2 + (dx+dxb, dy+dyb)/2 in
the second.
- The final implied match is between the point p1 - (dx+dxb+dxc, dy+dyb+dyc)/2
and p2 + (dx+dxb+dxc, dy+dyb+dyc)/2.
- The various SNR values cannot be directly compared because they are not
across the same rectangle.

#### Constructed by:

slicealignq5.py

#### Dependencies:

None
