#!/usr/bin/python

rawroot = '/lsi1/push/170428-SBEM'
sclroot = '/lsi2/dw/170428'
tmproot = '/home/wagenaar/tmp'

def rawtile(r, m, s):
    '''RAWTILE - Filename for raw image
    fn = RAWTILE(r, m, s) returns the path of the raw image for given
    run/montage/slice.'''
    root = rawroot
    if r==40 and s<143:
        return f'{root}/Run{r}/Montage_{m:03}/Prefix_OnPoint_{s:04}.tif'
    elif r==11:
        return f'{root}/Run{r}/Montage_{m:03}/Run{r}_OnPoint_OnPoint_{s:04}.tif'
    else:
        return f'{root}/Run{r}/Montage_{m:03}/Run{r}_OnPoint_{s:04}.tif'

