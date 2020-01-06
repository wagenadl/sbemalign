#!/usr/bin/python

# Directory for all output and most intermediate files
root = '/lsi2/dw/170428'

# Directory for intermediate files for movie creation
tmproot = '/home/wagenaar/tmp'

# Information about postgres database
database = 'align170428'
dbhost = 'localhost'
dbuser = 'wagenaar'

def rawtile(r, m, s):
    '''RAWTILE - Filename for raw image
    fn = RAWTILE(r, m, s) returns the path of the raw image for given
    run/montage/slice.'''
    root = '/lsi1/push/170428-SBEM'
    if r==40 and s<143:
        return f'{root}/Run{r}/Montage_{m:03}/Prefix_OnPoint_{s:04}.tif'
    elif r==11:
        return f'{root}/Run{r}/Montage_{m:03}/Run{r}_OnPoint_OnPoint_{s:04}.tif'
    else:
        return f'{root}/Run{r}/Montage_{m:03}/Run{r}_OnPoint_{s:04}.tif'

def ncolumns(r):
    if r>=51 or r<=2:
        return 1
    else:
        return 2
