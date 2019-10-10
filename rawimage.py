#!/usr/bin/python3

import numpy as np
import cv2
import config
import re
import skimage.io

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

def loadimage(url):
    '''LOADIMAGE - Download an image from anywhere
    img = LOADIMAGE(url) downloads an image from anywhere.
    Result is grayscale. Bit depth is as from source.
    URL may also be a local file name (even without file:// prefix).
    For unknown reasons, opencv's imdecode is slower than imread, so
    file:// is much faster than http://, even if the file system is a
    remote sshfs mount.'''
    #if url.startswith('file://'):
    #    url = url[7:]
    return skimage.io.imread(url)
    #if re.compile('^[a-z]+://').match(url):
    #    return resp = urllib.request.Request(url=url)
    #    with urllib.request.urlopen(resp) as f:
    #        dat = f.read()
    #        dat = np.fromfile(f, dtype=np.uint8)
    #    return cv2.imdecode(dat, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)
    #else:
    #    return cv2.imread(url, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)

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

def to8bit(img, stretch=.1):
    '''TO8BIT - Convert 16-bit image to 8 bits
    im8 = TO8BIT(im16) converts a 16-bit image to 8-bit format, clipping
    extreme values.
    Optional argument STRETCH determines what percentage of pixels gets
    clipped. Default is 0.1%.
    Calling TO8BIT on an 8-bit image has no effect; the contrast is not
    stretched in that case.'''
    if img.dtype==np.uint8:
        return img
    elif img.dtype==np.uint16:
        hst = cv2.calcHist([img], [0], None, [65536], [0, 65536])
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
        return mklut(lowb, upb)[img]
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
    return cv2.resize(img, (Ky, Kx), interpolation=cv2.INTER_AREA)

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
    ifn = f'/lsi2/dw/170428/scaled/Q5/R{r}/M{m}/S{s}.{ix}{iy}.tif'
    return loadimage(ifn)
