1300-neuroglancer
=================

This folder takes the pyramid and combines adjacent z-slices to maintain
approximate isotropic resolution.

The x-y resolution of the pyramid is A0: 5.5 nm, 1: 11, 2: 22, 3: 44, 4: 88,
5: 176, 6: 352, 7: 704, 8: 1408. The z resolution is always 50 nm. 

I will introduce a new "B" parameter that gives z resolution: B0: 50 nm,
1: 100 nm, 2: 200 nm, 3: 400 nm, 4: 800 nm, 5: 1600 nm. I will then
create imagery at A4/B1, A5/B2, A6/B3, A7/B4, A8/B5. 

Note that it would be possible to use z resolutions of 350, 700, 1400, but
I don't think I will do that, even though it would be closer to isotropy.

Inside the /md0/dw/170428-SBEM/pyramid folder, I will create new subfolders
A4B1 A5B2 A6B3 A7B4 and A8B5 to contain the new tifs. Even though at B5, Z
only runs up to 300, I will still maintain the Zzzz/zz structure.

Also in here, I will create chunked volumes that are more useful to 
neuroglancer.
At A4B1 through A8B5, I will create 64x64x64 voxel chunks. At A3B0, I will
create 64x64x32 voxel chunks. At A2B0, 128x128x16. If disk space permits,
I may also create 128x128x8 at A1B0 and perhaps 256x256x4 at A0B0, but I
doubt I have space for that.

