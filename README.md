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

These are organized in order of construction hierarchy, toplevel tables first.

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