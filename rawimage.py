#!/usr/bin/python3

import numpy as np
import cv2
import config
import re
import os

def rawtile(r, m, s):
    '''RAWTILE - Filename for raw image
    fn = RAWTILE(r, m, s) returns the path of the raw image for given
    run/montage/slice.
    img = LOADIMAGE(RAWTILE(r, m, s)) returns the actual image as a
    16-bit 17100x17100-pixel image.'''
    return config.rawtile(r, m, s)

def scaledtile(r, m, s, q):
    '''SCALEDTILE - Filename for scaled raw image
    fn = SCALEDTILE(r, m, s, q) returns the path of the scaled raw image
    run/montage/slice at scale Q.
    If Q is 1, returns RAWTILE filename.'''
    if q==1:
        return rawtile(r, m, s)
    return f'{config.sclroot}/scaled/Q{q}/R{r}/M{m}/S{s}.tif'

def partialq5tile(r, m, s, ix, iy):
    '''PARTIALQ5TILE - Filename for scaled raw subtile
    fn = PARTIALQ5TILE(r, m, s, ix, iy) returns the path of the scaled raw 
    subtile image (ix, iy) of run/montage/slice. IX and IY run from 0
    to 4 inclusive.'''
    return f'{config.sclroot}/scaled/Q5/R{r}/M{m}/S{s}.{ix}{iy}.tif'

def loadimage(url):
    '''LOADIMAGE - Download an image from anywhere
    img = LOADIMAGE(url) downloads an image from anywhere.
    Result is grayscale. Bit depth is as from source.
    URL may also be a local file name (even without file:// prefix).
    For unknown reasons, opencv's imdecode is slower than imread, so
    file:// is much faster than http://, even if the file system is a
    remote sshfs mount.'''
    if url.startswith('file://'):
        url = url[7:]
    #return skimage.io.imread(url)
    if re.compile('^[a-z]+://').match(url):
        resp = urllib.request.Request(url=url)
        with urllib.request.urlopen(resp) as f:
            dat = f.read()
            dat = np.fromfile(f, dtype=np.uint8)
        return cv2.imdecode(dat, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)
    else:
        return cv2.imread(url, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)

def loadlocaltif(fn):
    return cv2.imread(fn, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)

def ihist(img):
    '''IHIST - Calculate image histogram
    hst = IHIST(img) where IMG is either an 8-bit or a 16-bit image, calculates
    the image histogram. The result is a vector with 256 or 65536 entries.
    Histogram calculation is surprisingly slow for large images, despite
    the use of opencv.'''
    if img.dtype==np.uint8:
        hst = cv2.calcHist([img], [0], None, [256], [0, 256])
        return hst[:,0]
    elif img.dtype==np.uint16:
        hst = cv2.calcHist([img], [0], None, [65536], [0, 65536])
        return hst[:,0]
    else:
        raise ValueError('Only 8 or 16 bit images are supported')

def to8bit(img, stretch=.1, ignorethr=None, subst=None):
    '''TO8BIT - Convert 16-bit image to 8 bits
    im8 = TO8BIT(im16) converts a 16-bit image to 8-bit format, clipping
    extreme values.
    Optional argument STRETCH determines what percentage of pixels gets
    clipped. Default is 0.1%.
    Alternatively, STRETCH can be a dict mapping source percentiles to
    8-bit values. For instance, STRETCH = { 50: 130, 75: 172 } means
    that the 50th percentile of the 16-bit image gets mapped to 8-bit
    value 130 and the 75th percentile to 172. Precisely two percentiles
    must be specified in this case.
    If optional argument IGNORETHR is given, pixels blacker than that
    value are ignored in the histogram calculation. Furthermore, if SUBST
    is given, those ignored pixels are assigned output value SUBST.
    Calling TO8BIT on an 8-bit image has no effect; the contrast is not
    stretched in that case and substitutions are not applied.'''
    if img.dtype==np.uint8:
        return img
    elif img.dtype==np.uint16:
        if ignorethr is None:
            hst = cv2.calcHist([img], [0], None, [65536], [0, 65536])
        else:
            hst = cv2.calcHist([img[img>=ignorethr]], [0], None, [65536],
                               [0, 65536])
        hst = np.cumsum(hst)
        scl = hst[-1]
        if type(stretch)==dict:
            if len(stretch)!=2:
                raise ValueError('Must have exactly two percentiles')
            kk = list(stretch.keys())
            kk.sort()
            p1 = np.searchsorted(hst, .01*kk[0]*scl)
            p2 = np.searchsorted(hst, .01*kk[1]*scl)
            v1 = stretch[kk[0]]
            v2 = stretch[kk[1]]
        else:
            p1 = np.searchsorted(hst, .01*stretch*scl)
            p2 = np.searchsorted(hst, (1-.01*stretch)*scl)
            v1 = 0
            v2 = 256
        lut = np.arange(65536)
        lut = v1 + (v2-v1) * (lut - p1) / (1 + p2 - p1)
        lut[lut<0] = 0
        lut[lut>255] = 255
        lut = lut.astype(np.uint8)
        res = lut[img]
        if subst is not None:
            res[img<ignorethr] = subst
        return res
    else:
        raise ValueError('Only 8 or 16 bit images are supported')

def iscale(img, n):
    '''ISCALE - Downscale an image
    im = ISCALE(img, n) downscales the image IMG by a factor N in both
    x- and y-directions using "area" interpolation.'''
    y,x = img.shape
    Ky = y//n
    Kx = x//n
    if n*Ky<y or n*Kx<x:
        img = img[0:n*Ky,0:n*Kx]
    return cv2.resize(img, (Kx, Ky), interpolation=cv2.INTER_AREA)

def ipad(img, pad=512, padc=0):
    '''IPAD - Pad an image
    im = IPAD(img, pad) zero pads an image so that its size is an integral
    multiple of PAD in either direction. Optional argument PADC specifies
    an alternate padding color.'''
    Y, X = img.shape
    KY = pad * ((Y+pad-1)//pad)
    KX = pad * ((X+pad-1)//pad)
    if KY==Y and KX==X:
        return img
    res =  np.empty((KY,KX), dtype=img.dtype)
    res[Y:,:] = padc
    res[:Y,X:] = padc
    res[:Y,:X] = img
    return res

def saveimage(img, ofn, qual=None):
    '''SAVEIMAGE - Save an image to disk
    SAVEIMAGE(img, fn) saves the image to disk.
    Optional argument QUAL (0..100) specifies quality for
    jpeg output.
    Unlike LOADIMAGE, SAVEIMAGE can only save to local file, not to URL.'''
    if qual is None:
        cv2.imwrite(ofn, img)
    else:
        cv2.imwrite(ofn, img, (cv2.IMWRITE_JPEG_QUALITY, qual))

def scaledraw(r, m, s, q=25):
    url = f'http://leechem.caltech.edu:9090/scaledraw/r{r}/m{m}/s{s}/q{q}.jpg'
    return loadimage(url)

def betaimg(z, a=6):
    x0um = 70
    y0um = 70
    wum = 170
    hum = 610
    x0 = int(x0um/.0055) // (2**a)
    y0 = int(y0um/.0055) // (2**a)
    w = int(wum/.0055) // (2**a)
    h = int(hum/.0055) // (2**a)
    url = f'http://leechem.caltech.edu:9090/roi_pix/z{z}/x{x0}/y{y0}/w{w}/h{h}/a{a}.jpg'
    print(x0, y0, w, h, url)
    return loadimage(url)

def partialq5img(r, m, s, ix, iy):
    return loadimage(partialq5tile(r, m, s, ix, iy))

def q5subimg2x2(r, m, s, ix, iy):
    X = 684
    img = np.zeros((X*2, X*2), dtype=np.uint8)
    img[:X,:X] = partialq5img(r,m,s,ix-1,iy-1)
    img[:X,X:] = partialq5img(r,m,s,ix,iy-1)
    img[X:,:X] = partialq5img(r,m,s,ix-1,iy)
    img[X:,X:] = partialq5img(r,m,s,ix,iy)
    return img

def fullq5img(r, m, s):
    X = 684
    img = np.zeros((5*X,5*X), dtype=np.uint8)
    for ix in range(5):
        for iy in range(5):
            im1 = partialq5img(r,m,s,ix,iy)
            if im1 is None:
                return None
            img[iy*X:(iy+1)*X,ix*X:(ix+1)*X] = im1
    return img

def q25img(r, m, s):
    #url = f'http://leechem.caltech.edu:9090/scaledraw/Q25/R{r}/M{m}/S{m}.ppm'
    return loadimage(scaledtile(r, m, s, 25))

def fullq1img(r, m, s, stretch=.1):
    ifn = rawtile(r, m, s)
    if not os.path.exists(ifn):
        return None
    img = loadimage(ifn)
    return to8bit(img, stretch)

    
