#!/usr/bin/python3

r = 21
zoffset_21 = 3712 # Slice number
q = 25
qroot = '/home/wagenaar/tmp/R21'
uctre = '/home/wagenaar/tmp/uctre'
uctz_scale = 10
uctxy_scale = 100
rel_scale = uctxy_scale / q
s = 150
m = 0
z = zoffset_21 + s
uctz = z // uctz_scale

import swiftir
import em170428.runinfo
import pyqplot as qp
import numpy as np

ri = em170428.runinfo.RunInfo()

em = swiftir.loadImage('%s/M%i/S%i.tif' % (qroot, m, s))
uct = swiftir.loadImage('%s/uct-%i.ppm' % (uctre, uctz))

center_of_em_x, center_of_em_y, center_of_em_z = ri.tileCenter(r, m, s)
# In unscaled pixels
center_in_uct_x = center_of_em_x / uctxy_scale
center_in_uct_y = center_of_em_y / uctxy_scale

Y,X = em.shape
Rem = 512
Ruct = int(Rem/rel_scale)
roi_em = swiftir.extractROI(em, ((X-Rem)//2, (Y-Rem)//2, Rem, Rem))
roi_uct = swiftir.extractROI(uct, (int(center_in_uct_x - Ruct/2),
                                   int(center_in_uct_y - Ruct/2),
                                   Ruct, Ruct))

qp.figure('/tmp/s_em')
qp.imsc(roi_em)
qp.figure('/tmp/s_uct')
qp.imsc(roi_uct)

# OK. Now see what SWIM thinks.
apo_em = swiftir.apodize(roi_em)
apo_uct = swiftir.apodize(roi_uct)

qp.figure('/tmp/s_em1')
qp.imsc(apo_em)
qp.figure('/tmp/s_uct1')
qp.imsc(apo_uct)

fft_em = swiftir.fft(apo_em)
fft_uct = swiftir.fft(apo_uct)

qp.figure('/tmp/s_em1')
qp.imsc(np.log(fft_em[:,:,0]**2 + fft_em[:,:,1]**2+1))
qp.figure('/tmp/s_uct1')
qp.imsc(np.log(fft_uct[:,:,0]**2 + fft_uct[:,:,1]**2+1))

keep = np.hstack((np.arange(Ruct//2), np.arange(Ruct//2)+Rem-Ruct//2))
fft_em2 = fft_em[keep,:,:]
fft_em2 = fft_em2[:,keep,:]

qp.figure('/tmp/s_em1')
qp.imsc(np.log(fft_em2[:,:,0]**2 + fft_em2[:,:,1]**2+1))

(dx, dy, sx, sy, snr) = swiftir.swim(fft_uct, fft_em2)

print(dx, dy, sx, sy, snr)
