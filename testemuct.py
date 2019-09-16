#!/usr/bin/python3

r = 21
zoffset_21 = 3712 # Slice number
q = 25
qroot = '/home/wagenaar/tmp/R21'
uctre = '/home/wagenaar/tmp/uctre'
uctz_scale = 10
uctxy_scale = 100
rel_scale = uctxy_scale / q
s = 50
m = 7
z = zoffset_21 + s
uctz = z // uctz_scale

import swiftir
import em170428.runinfo

em = swiftir.loadImage('%s/M%i/S%i.tif' % (qroot, m, s))
uct = swiftir.loadImage('%s/uct-%i.jpg' % (uctre, uctz))

center_of_em_x, center_of_em_y, center_of_em_z = ri.tileCenter(r, m, s)
# In unscaled pixels
center_in_uct_x = center_of_em_x / uctxy_scale
center_in_uct_y = center_of_em_y / uctxy_scale

Y,X = em.shape
R = 512
roi_em = swiftir.extractROI(em, ((X-R)//2, (Y-R)//2, R, R))
roi_uct = swiftir.extractROI(uct, (int(center_in_uct_x - R/8),
                                   int(center_in_uct_y - R/8),
                                   R//4, R//4))

qp.figure('/tmp/s_em')
qp.imsc(roi_em)
qp.figure('/tmp/s_uct')
qp.imsc(200-roi_uct)

# OK. Now see what SWIM thinks.
