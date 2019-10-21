#!/usr/bin/python3

import numpy as np
import cv2

def getPerspective(xmodel, ymodel, ximage, yimage):
    '''GETPERSPECTIVE - Obtain transformation matrix from matching points
    T_mdl2img = GETPERSPECTIVE(xmodel, ymodel, ximage, yimage) where
    (xmodel, ymodel) are four points in model space and (ximage, yimage) are
    four matching points in the space of a source image, returns a perspective
    matrix that transforms model coordinates to image coordinates.'''
    xymodel = np.stack((xmodel, ymodel), 1).astype('float32')
    xyimage = np.stack((ximage, yimage), 1).astype('float32')
    mdl2img = cv2.getPerspectiveTransform(xymodel, xyimage)
    return mdl2img

def warpPerspective(img, mdl2img):
    '''WARPPERSPECTIVE - Copy an image to model space
    mdl, x0, y0 = WARPPERSPECTIVE(img, mdl2img) fills a portion of model
    space that fully covers the source IMG based on the perspective
    transform MDL2IMG. It returns the resulting image (MDL) as well as
    the coordinates (X0, Y0) where the top-left corner of that image
    should live in model space. (X0, Y0) are guaranteed to be integers.
    Pixels outside of the defined region are set to black (0).'''
    Y,X = img.shape
    imgcorners = np.array([[[0,0], [X,0], [0,Y], [X,Y]]], dtype='float32')
    img2mdl= np.linalg.inv(mdl2img)
    mdlcorners = cv2.perspectiveTransform(imgcorners, img2mdl)
    xmdlcorners = mdlcorners[0,:,0]
    ymdlcorners = mdlcorners[0,:,1]
    x0 = int(np.min(xmdlcorners))
    y0 = int(np.min(ymdlcorners))
    x1 = np.max(xmdlcorners)
    y1 = np.max(ymdlcorners)
    W = int(np.ceil(x1 - x0))
    H = int(np.ceil(y1 - y0))
    translate = np.array([[1,0,x0], [0,1,y0], [0,0,1]])
    xform = np.matmul(mdl2img, translate)
    mdl = cv2.warpPerspective(img, xform, (W,H), flags=cv2.WARP_INVERSE_MAP)
    return (mdl, x0, y0)

def quadToImageBox(xmdlbox, ymdlbox, mdl2img, shp):
    '''QUADTOIMAGEBOX - Find rectangle in image needed to cover model quad
    x0,y0,x1,y1 = QUADTOIMAGEBOX(xmodel, ymodel, mdl2img) finds the (integer)
    edges of a bounding rectangle that fully covers the quadrilateral 
    specified by XMODEL and YMODEL given the transformation matrix MDL2IMG.'''
    mdlcorners = np.reshape(np.stack((xmdlbox, ymdlbox), 1), (1, 4, 2))
    imgcorners = cv2.perspectiveTransform(mdlcorners.astype(np.float32), mdl2img)
    x0 = int(np.min(imgcorners[:,:,0])-1)
    x1 = int(np.max(imgcorners[:,:,0])+1+1)
    y0 = int(np.min(imgcorners[:,:,1])-1)
    y1 = int(np.max(imgcorners[:,:,1])+1+1)
    if x0<0:
        x0 = 0
    if y0<0:
        y0 = 0
    if y1>=shp[0]:
        y1 = shp[0]
    if x1>=shp[1]:
        x1 = shp[1]
    return (x0, y0, x1, y1)

def roundedquad(xx, yy):
    # xx,yy must be presented in order tl, tr, bl, br
    # results are presented in order tl, tr, br, bl
    return (np.array([int(xx[0]), int(xx[1]+1), int(xx[3]+1), int(xx[2])]),
            np.array([int(yy[0]), int(yy[1]), int(yy[3]+1), int(yy[2])+1]))

def createClipMask(xmodel, ymodel, x0,y0, w,h):
    # Box must be tl, tr, bl, br
    xxx,yyy = roundedquad(xmodel, ymodel)
    msk = np.zeros((h,w), dtype=np.uint8)
    pts = np.stack((xxx - x0, yyy - y0), 1)
    cv2.fillPoly(msk, np.array([pts]), color=255)
    return msk

def warpPerspectiveBoxed(img, xmdlbox, ymdlbox, mdl2img):
    (x0b,y0b,x1b,y1b) = quadToImageBox(xmdlbox, ymdlbox, mdl2img, img.shape)
    subimg = img[y0b:y1b,x0b:x1b] # This works by reference
    translate = np.array([[1,0,-x0b], [0,1,-y0b], [0,0,1]])
    xform = np.matmul(translate, mdl2img)
    mdlimg, mx0, my0 = warpPerspective(subimg, xform)
    h, w = mdlimg.shape
    msk = createClipMask(xmdlbox, ymdlbox, mx0, my0, w, h)
    return mdlimg, msk, mx0, my0

def copyWithMask(mdl, img, msk, x0, y0):
    if mdl.dtype != np.uint8:
        raise Exception('model must be uint8')
    if img.dtype != np.uint8:
        img = img.astype(np.uint8)
    h,w = img.shape
    H,W = mdl.shape
    if x0+w > W or y0+h > H or x0<0 or y0<0:
        # doesn't quite fit
        if x0<0:
            w += x0
            x0 = 0
        if y0<0:
            h += y0
            y0 = 0
        if x0+w>W:
            w = W - x0
        if y0+h>H:
            h = H - y0
        if w<=0 or h<=0:
            return
        inx = slice(0, w)
        iny = slice(0, h)
        outxx = slice(x0, x0+w)
        outyy = slice(y0, y0+h)
        ms = msk[iny,inx]
        np.bitwise_or(np.bitwise_and(img[iny,inx], ms),
                      np.bitwise_and(mdl[outyy,outxx],
                                     np.bitwise_not(ms)),
                      mdl[outyy,outxx])
    else:
        outxx = slice(x0, x0+w)
        outyy = slice(y0, y0+h)
        np.bitwise_or(np.bitwise_and(img, msk),
                      np.bitwise_and(mdl[outyy,outxx],
                                     np.bitwise_not(msk)),
                      mdl[outyy,outxx])

if __name__=='__main__':
    import swiftir
    import pyqplot as qp
    ifn = 'test-top-10.jpg'
    img = swiftir.loadImage(ifn)
    # Following points define the perspective transform
    xmodel = np.array([380, 390, 680, 650])
    ymodel = np.array([1150, 1380, 1390, 1140])
    ximage = np.array([140, 160, 380, 330])
    yimage = np.array([130, 320, 380, 130])
    # Following points define the area to copy
    xmdlbox = np.array([350, 360, 700, 700])
    ymdlbox = np.array([1100, 1400, 1400, 1100])

    qp.figure('/tmp/s1');
    Y,X = img.shape
    qp.imsc(img, xx=np.arange(X), yy=np.arange(Y))
    qp.pen('r', 1)
    qp.marker('+', 4)
    qp.mark(ximage, yimage)

    mdl2img = getPerspective(xmodel, ymodel, ximage, yimage)
    
    mdl, x0, y0 = warpPerspective(img, mdl2img)

    qp.figure('/tmp/s2')
    Y,X = mdl.shape
    qp.imsc(mdl, xx=np.arange(X), yy=np.arange(Y))
    qp.pen('r', 1)
    qp.marker('+', 4)
    qp.mark(xmodel - x0, ymodel - y0)
    qp.marker('x', 4)
    qp.mark(xmdlbox - x0, ymdlbox - y0)
    
    
    qp.figure('/tmp/s3')
    msk = createClipMask(xmodel, ymodel, x0,y0, X,Y)
    qp.imsc(msk, xx=np.arange(X), yy=np.arange(Y))
    
    dst = np.zeros((Y,X), dtype='uint8') + 128
    dst = 255-mdl
    qp.figure('/tmp/s4')
    copyWithMask(dst, mdl, msk, 0, 0)
    qp.imsc(dst, xx=np.arange(X), yy=np.arange(Y))

    qp.figure('/tmp/s5')
    mdl, msk, x0, y0 = warpPerspectiveBoxed(img, xmdlbox, ymdlbox, mdl2img)
    Y,X = mdl.shape
    qp.imsc(mdl, xx=np.arange(X), yy=np.arange(Y))
    qp.figure('/tmp/s6')
    mdl, msk, x0, y0 = warpPerspectiveBoxed(img, xmdlbox, ymdlbox, mdl2img)
    qp.imsc(msk, xx=np.arange(X), yy=np.arange(Y))

    ds2 = np.zeros((1500,1500), dtype=np.uint8) + 128
    Y,X = ds2.shape
    qp.figure('/tmp/s7')
    copyWithMask(ds2, mdl, msk, x0, y0)
    qp.imsc(ds2, xx=np.arange(X), yy=np.arange(Y))
    
