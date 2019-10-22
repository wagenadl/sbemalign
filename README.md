The realignment process
=======================

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

Database tables
---------------

These are presented here in order of construction hierarchy, toplevel
tables first.

### warpq5rundone

Simple bookkeeping of which slices have been rendered at Q5 so far.

#### Members:
r, s
: run and slice number

#### Associated files:
/lsi2/dw/170428/runalignq5/R{r}/S{s}.jpg
: Entire slice at 1:5 scale

#### Constructed by:
warpq5run.py

#### Dependencies:

- roughq5pos
- optimizeq5

#### See also:

- runextentq5

### runextentq5

Coordinates in run space that are included in the rendered images at Q5

#### Members:
r
: run number

x0, y0, x1, y1
: coordinates of image in run space  
    x0, y0 are inclusive; x1, y1 are exclusive.

#### Constructed by:
warpq5run.py

#### Dependencies:

- roughq5pos
- optimizeq5

### roughq5pos

Slightly misnamed, this table contains positions of montages within a run.

#### Members:
r, m
: run and montage numbers

x, y
: coordinates of given montage

#### Constructed by:
optimizeq5.py

#### Dependencies:
- slicealignq5
- montagealignq5relhp
- montagealignattouchq5relhp

#### See also:

- optimizeq5

#### Notes

By construction, the average of all x-values for montages within a run is zero.
Same for y-values.

### optimizeq5

Results of optimization within runs at Q5

#### Members:

r, m, s
: run, montage, slice number

nx, ny
: grid index within a tile.
  nx, ny run from 0 to 6, inclusive.

x, y
: montage coordinates of a grid point, see _Transformations_ below.

dx, dy
: displacement at that grid point, see _Transformations_ below.

supported
: indicates whether this grid point is supported by data. (If not, it
  is the result of interpolation.)

#### Coordinate transformations

- To transform from montage coordinates to run coordinates, add the
  x-, y-values stored in roughposq5.
- To transform from run coordinates to pixels in the "runalignq5" images,
  subtract the x0-, y0-values stored in runextentq5.
- To transform between run coordinates and pixel coordinates in a given
  scaled source tile, use the runToQuad() and quadToRun() function of
  the Transformer class in warpq5run. Not sure which quad a pixel belongs to?
  Use the findInRun() or findInMontage() functions, or the rawToRun() function.

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

#### See also:

- slicealignq5pos

### slicealignq5pos

View onto slicealignq5 that precalculates the final matching points

#### Members

r,s
: as in slicealignq5

m1, m2
: as in slicealignq5

ix1, iy1, ix2, iy2
: as in slicealignq5

dx0, dy
: as in slicealignq5

x1, y1
: calculated as p1 - (dx+dxb+dxc, dy+dyb+dyc)/2 from slicealignq5

x2, y2
: calculated as p2 + (dx+dxb+dxc, dy+dyb+dyc)/2 from slicealignq5

#### Constructed by:

slicealignq5.py

