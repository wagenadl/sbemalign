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

def createClipMask(xmodel, ymodel, x0,y0, w,h):
    msk = np.zeros((h,w), dtype=np.uint8)
    pts = np.stack((xmodel - x0, ymodel - y0), 1)
    cv2.fillPoly(msk, np.array([pts]), color=255)
    return msk

def copyWithMask(mdl, img, msk):
    if mdl.dtype != np.uint8:
        mdl = mdl.astype(np.uint8)
    if img.dtype != np.uint8:
        img = img.astype(np.uint8)
    return np.bitwise_or(np.bitwise_and(img, msk),
                          np.bitwise_and(mdl, np.bitwise_not(msk)))

if __name__=='__main__':
    import swiftir
    import pyqplot as qp
    ifn = '/home/wagenaar/Desktop/stuff/vsdbuildup/top-10.jpg'
    img = swiftir.loadImage(ifn)
    xmodel = np.array([380, 390, 680, 650])
    ymodel = np.array([1150, 1380, 1390, 1140])
    ximage = np.array([140, 160, 380, 330])
    yimage = np.array([130, 320, 380, 130])

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
    
    qp.figure('/tmp/s3')
    msk = createClipMask(xmodel, ymodel, x0,y0, X,Y)
    qp.imsc(msk)
    
    dst = np.zeros((Y,X), dtype='uint8') + 128
    qp.figure('/tmp/s4')
    qp.imsc(copyWithMask(dst, mdl, msk))
