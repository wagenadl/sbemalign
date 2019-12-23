#!/usr/bin/python3

import numpy as np
import cv2
import config
import re
import os
#import skimage.io

def rawtile(r, m, s):
    '''RAWTILE - Filename for raw image
    fn = RAWTILE(r, m, s) returns the path of the raw image for given
    run/montage/slice.
    img = LOADIMAGE(RAWTILE(r, m, s)) returns the actual image as a
    16-bit 17100x17100-pixel image.'''
    root = config.rawroot
    if r==40 and s<143:
        pat = '%s/Run%i/Montage_%03i/Prefix_OnPoint_%04i.tif'
        return pat % (root, r, m, s)
    if r==11:
        pat = '%s/Run%i/Montage_%03i/Run%i_OnPoint_OnPoint_%04i.tif'
    else:
        pat = '%s/Run%i/Montage_%03i/Run%i_OnPoint_%04i.tif'
    return pat % (root, r, m, r, s)

def rawsubtile(r, m, s, ix, iy):
    '''RAWSUBTILE - Filename for subtile out of raw image
    fn = RAWSUBTILE(r, m, s, ix, iy) returns the path of a tiff subtile
    out of the raw data (at full scale). There are 33x33 subtiles for
    each tile, numbered 0..32 (inclusive). Each is 512x512 pixels wide.
    This means that at each of the four edges of a raw image, 102
    pixels are not incorporated in the subtiles.'''
    dr = config.oldroot + f'/unaligned/R{r:03d}/M{m:03d}/S{s:04d}'
    fn = f'X{ix:02d}Y{iy:02d}.tif'
    return f'{dr}/{fn}'

def scaledtile(r, m, s, q):
    '''SCALEDTILE - Filename for scaled raw image
    fn = SCALEDTILE(r, m, s, q) returns the path of the scaled raw image
    run/montage/slice at scale Q.
    If Q is 1, returns RAWTILE filename.'''
    if q==1:
        return rawtile(r, m, s)
    root = config.sclroot
    pat = '%s/scaled/Q%i/R%i/M%i/S%i.tif'
    return pat % (root, q, r, m, s)

def partialq5tile(r, m, s, ix, iy):
    '''PARTIALQ5TILE - Filename for scaled raw subtile
    fn = PARTIALQ5TILE(r, m, s, ix, iy) returns the path of the scaled raw 
    subtile image (ix, iy) of run/montage/slice. IX and IY run from 0
    to 4 inclusive.'''
    return f'/lsi2/dw/170428/scaled/Q5/R{r}/M{m}/S{s}.{ix}{iy}.tif'

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
    Calling TO8BIT on an 8-bit image has no effect; the contrast is not
    stretched in that case.
    If optional argument IGNORETHR is given, pixels blacker than that
    value are ignored in the histogram calculation. Furthermore, if SUBST
    is given, those ignored pixels are assigned output value SUBST.'''
    if img.dtype==np.uint8:
        return img
    elif img.dtype==np.uint16:
        if ignorethr is None:
            hst = cv2.calcHist([img], [0], None, [65536], [0, 65536])
        else:
            hst = cv2.calcHist([img[img>=ignorethr]], [0], None, [65536],
                               [0, 65536])
        hst = np.cumsum(hst)
        hst = hst / hst[-1]
        lowb = np.argmax(hst>=.01*stretch)
        upb = np.argmax(hst>=1-.01*stretch)

        def mklut(lowb, upb):
            lut = np.arange(65536)
            lut.clip(lowb, upb, out=lut)
            lut -= lowb
            np.floor_divide(lut, (upb-lowb+1) / 256, out=lut, casting='unsafe')
            return lut.astype(np.uint8)
        res = mklut(lowb, upb)[img]
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

def q1subimg(r, m, s, ix, iy):
    return loadlocaltif(rawsubtile(r, m, s, ix, iy))

def q1subimg2x2(r, m, s, ix, iy):
    img = np.zeros((1024, 1024), dtype=np.uint8)
    img[:512,:512] = q1subimg(r,m,s,ix-1,iy-1)
    img[:512,512:] = q1subimg(r,m,s,ix,iy-1)
    img[512:,:512] = q1subimg(r,m,s,ix-1,iy)
    img[512:,512:] = q1subimg(r,m,s,ix,iy)
    return img

def q1roi(r, m, s, x, y, w, h):
    '''Q1ROI - Load an arbitrary rectangular patch from a raw tile
    img = Q1ROI(r, m, s, x, y, w, h) loads an WxH-sized image patch
    from the "unaligned" tiles with top-left at (X,Y). Note that this
    corresponds to (X+102, Y+102) in the 16-bit raw tiffs.'''
    print('q1roi', r, m, s, x, y, w, h)
    R = 512
    x0 = x // R # Leftmost column to read, inclusive
    x1 = (x+w+R-1) // R # Rightmost column to read, exclusive
    y0 = y // R
    y1 = (y+h+R-1) // R
    dx = x - R*x0
    dy = y - R*y0
    img = np.zeros((h,w), dtype=np.uint8)
    Y0 = 0
    for iy in range(y0, y1):
        yidx = slice(0, R)
        if iy==y0:
            yidx = slice(dy, R)
        elif iy==y1-1 and dy>0:
            yidx = slice(0, dy)
        X0 = 0
        for ix in range(x0, x1):
            xidx = slice(0, R)
            if ix==x0:
                xidx = slice(dx, R)
            elif ix==x1-1 and dx>0:
                xidx = slice(0, dx)
            sub = q1subimg(r, m, s, ix, iy)
            sh,sw = sub.shape
            if sh!=R or sw!=R:
                 raise Exception(f'Wrong size subimg R{r} M{m} S{s} {ix},{iy}')
            sub = sub[yidx, xidx]
            H,W = sub.shape
            if X0+W > w:
                 raise Exception(f'Bug X0+W>w at R{r} M{m} S{s} {ix},{iy}')
            if Y0+H > h:
                 raise Exception(f'Bug Y0+H>h at R{r} M{m} S{s} {ix},{iy}')
            img[Y0:Y0+H, X0:X0+W] = sub
            X0 += W
        Y0 += H
    print('q1roi done', r, m, s, x, y, w, h)
    return img

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

    
