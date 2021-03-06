#!/usr/bin/python3

# Core facilities for optimizign the relative positions of all the
# montages in a subvolume at q5.
# See E&R p. 1611ff.

crosstbl = 'slicealignq5'
intratbl = 'relmontalignq5'
edgetbl  = 'relmontattouchq5'
intertbl = 'interrunq5'
transtbl = 'transrunmontq5'

IX = IY = 5
X = Y = 684

import aligndb
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

db = aligndb.DB()
ri = db.runinfo()

def dynamicthreshold(snrs):
    return 0.5 * np.max(snrs)

def keepsome(idx, *args):
    res = []
    for a in args:
        res.append(a[idx])
    return res

class MatchPoints:
    def __init__(self):
        self.r1 = None
        self.m1 = None
        self.s1 = None
        self.r2 = None
        self.m2 = None
        self.s2 = None
        self.xx1 = None
        self.yy1 = None
        self.xx2 = None
        self.yy2 = None

    def x(self, ax):
        # ax=0 for x or 1 for y
        # This is "x" on E&R p. 1611
        if ax:
            return self.yy1
        else:
            return self.xx1

    def xp(self, ax):
        # ax=0 for x or 1 for y
        # This is "x'" on E&R p. 1611
        if ax:
            return self.yy2
        else:
            return self.xx2

    def combine(self, other):
        if self.r1 != other.r1 or self.r2 != other.r2:
            raise Exception('Cannot combine: run mismatch')
        if self.m1 != other.m1 or self.m2 != other.m2:
            raise Exception('Cannot combine: montage mismatch')
        if self.s1 != other.s1 or self.s2 != other.s2:
            raise Exception('Cannot combine: slice mismatch')
        self.xx1 = np.hstack((self.xx1, other.xx1))
        self.yy1 = np.hstack((self.yy1, other.yy1))
        self.xx2 = np.hstack((self.xx2, other.xx2))
        self.yy2 = np.hstack((self.yy2, other.yy2))

    def __repr__(self):
        n = f'-- {len(self.xx1)} pts'
        if self.r1==self.r2:
            if self.m1==self.m2:
                return f'R{self.r1} M{self.m1} S{self.s1}:{self.s2} {n}'
            else:
                if self.s1 is None:
                    return f'R{self.r1} M{self.m1}:{self.m2} {n}'
        if self.s1 is None:
            return f'R{self.r1}.M{self.m1} : R{self.r2}.M{self.m2} {n}'
        s = f'R{self.r1}.M{self.m1}.S{self.s1}'
        return s + f' : R{self.r2}.M{self.m2}.S{self.s2} {n}'

    def cross(r, m1, m2, s0=0, s1=None, thr=None, perslice=False):
        # Returns a single MatchPoints with combined data for all slices,
        # unless PERSLICE is True, in which case a list of MatchPoints
        # for individual slices is returned.
        swhere = f' and s>={s0}'
        if s1 is not None:
            swhere += f' and s<{s1}'
        #print(r, m1, m2, s0, s1, thr, perslice, swhere)
        (s, x1,y1, x2,y2, snr) = db.vsel(f'''select
        s,
        x1-dx/2-dxb/2-dxc/2,
        y1-dy/2-dyb/2-dyc/2,
        x2+dx/2+dxb/2+dxc/2,
        y2+dy/2+dyb/2+dyc/2,
        snrc from {crosstbl}
        where r={r} and m1={m1} and m2={m2} {swhere}
        order by ii''')
        if perslice:
            mpp = []
            for s1 in np.unique(s):
                X1,Y1,X2,Y2,SNR = keepsome(s==s1, x1,y1,x2,y2,snr)
                if thr is None:
                    dthr = dynamicthreshold(SNR)
                else:
                    dthr = thr
                mp = MatchPoints()
                mp.r1 = mp.r2 = r
                mp.m1 = m1
                mp.m2 = m2
                mp.s1 = mp.s2 = s1
                mp.xx1, mp.yy1, mp.xx2, mp.yy2 = keepsome(SNR>dthr, X1,Y1,X2,Y2)
                mpp.append(mp)
            return mpp
        else:
            if thr is None:
                dthr = dynamicthreshold(snr)
            else:
                dthr = thr
            mp = MatchPoints()
            mp.r1 = mp.r2 = r
            mp.m1 = m1
            mp.m2 = m2
            mp.s1 = mp.s2 = None
            mp.xx1, mp.yy1, mp.xx2, mp.yy2 = keepsome(snr>dthr, x1,y1,x2,y2)
            return mp
    def allcross(r, s0=0, s1=None, thr=None, perslice=False):
        mpp = []
        C = ri.ncolumns(r)
        R = ri.nrows(r)
        for col in range(C):
            for row in range(R-1):
                mp = MatchPoints.cross(r, row*C+col, (row+1)*C+col,
                                       s0, s1, thr, perslice)
                if perslice:
                    mpp += mp
                else:
                    mpp.append(mp)
        for col in range(C-1):
            for row in range(R):
                mp = MatchPoints.cross(r, row*C+col, row*C+col+1,
                                       s0, s1, thr, perslice)
                if perslice:
                    mpp += mp
                else:
                    mpp.append(mp)
        return mpp
        
    def anytrans(r1, m1, r2, thr=None, perslice=False):
        # We return a list of MatchPoints, one for each existing m2.
        # If perslice is True, s1 and s2 are stored in the MatchPoints,
        # otherwise, None.
        # r2 must be r1 ± 1.
        # "Forward" is when r1>r2.
        # Note that we are changing our nomenclature here to match the table.
        # If THR is None, we use a straightforward dynamic threshold.
        # If THR < 0, we use a dynamic threshold that can not be less
        # than |THR|.
        (s1,s2,m2, x1,y1, x2,y2, snr) = db.vsel(f'''select
        s,s2,m2,
        (ix+0.5)*{X}-dx/2-dxb/2,
        (iy+0.5)*{Y}-dy/2-dyb/2,
        x+dx/2+dxb/2,
        y+dy/2+dyb/2,
        snrb
        from {transtbl}
        where r={r1} and m={m1} and r2={r2}''')
        if thr is None:
            thr = dynamicthreshold(snr)
        elif thr<0:
            thr = -thr
            thr1 = dynamicthreshold(snr)
            if -thr > thr1:
                thr = thr1
        keep = snr>thr
        s1 = s1[keep]
        s2 = s2[keep]
        m2 = m2[keep]
        x1 = x1[keep]
        y1 = y1[keep]
        x2 = x2[keep]
        y2 = y2[keep]
        snr = snr[keep]
        if m2.size==0:
            msg = f'Trans failed >= {thr} at R{r1}:R{r2} M{m1}:~'
            raise Exception(msg)
        mm2 = np.unique(m2)
        mpp = []
        for m in mm2:
            mp = MatchPoints()
            mp.r1 = r1
            mp.m1 = m1
            mp.r2 = r2
            mp.m2 = m
            if perslice:
                mp.s1 = s1[0]
                mp.s2 = s2[0]
            else:
                mp.s1 = mp.s2 = None
            mp.xx1 = x1[m2==m]
            mp.yy1 = y1[m2==m]
            mp.xx2 = x2[m2==m]
            mp.yy2 = y2[m2==m]
            mpp.append(mp)
        return mpp
    def forwardtrans(r1, m1, thr=None, perslice=False):
        # Implicitly, r2=r1+1, s2=0, s1=S(r1)-1.
        return MatchPoints.anytrans(r1, m1, r1+1, thr, perslice)
    def backtrans(r2, m2, thr=None, perslice=False):
        return MatchPoints.anytrans(r2, m2, r2-1, thr, perslice)
    def alltrans(r1, r2, thr=None, perslice=False):
        if r2 != r1+1:
            raise Exception('r2 must be r1+1')
        fwd = {} # keys are (m1,m2)
        bck = {} # keys are (m1,m2) [not the other way around]
        minthr = thr
        if minthr is None:
            snr = db.vsel(f'''select snrb from {transtbl}
            where r={r1} and r2={r2}''')
            minthr = np.max((10, 0.2 * dynamicthreshold(snr)))
            thr = -minthr
        for m1 in range(ri.nmontages(r1)):
            try:
                for mp in MatchPoints.forwardtrans(r1, m1, thr, perslice):
                    fwd[(m1, mp.m2)] = mp
            except Exception as e:
                print(f'WARNING: Failed to find any points>{thr} in trans R{r1}M{m1} : R{r2}: ', e)
        for m2 in range(ri.nmontages(r2)):
            try:
                for mp in MatchPoints.backtrans(r2, m2, thr, perslice):
                    bck[(mp.m2, m2)] = mp
            except Exception as e:
                print(f'WARNING: Failed to find any points>{thr} in trans R{r2}M{m2} : R{r1}: ', e)
        mpp = []
        for mm, mp in fwd.items():
            if mm in bck:
                mp1 = bck[mm]
                mp.xx1 = np.concatenate((mp.xx1, mp1.xx2))
                mp.xx2 = np.concatenate((mp.xx2, mp1.xx1))
                mp.yy1 = np.concatenate((mp.yy1, mp1.yy2))
                mp.yy2 = np.concatenate((mp.yy2, mp1.yy1))
            mpp.append(mp)
        for mm, mp in bck.items():
            if mm not in fwd:
                mp.r1, mp.r2 = mp.r2, mp.r1
                mp.m1, mp.m2 = mp.m2, mp.m1
                mp.s1, mp.s2 = mp.s2, mp.s1
                mp.xx1, mp.xx2 = mp.xx2, mp.xx1
                mp.yy1, mp.yy2 = mp.yy2, mp.yy1
                mpp.append(mp)
        return mpp

    def intra(r, m, s0, s1, thr=None):
        # Returns a list of MatchPoints with data from the intra table
        # for each of the slice pairs in [s0, s1).
        if s1==s0+1:
            return []
        (s, x2,y2, x1,y1, snr) = db.vsel(f'''select
        s,
        (ix+0.5)*{X}-dx/2-dxb/2,
        (iy+0.5)*{Y}-dy/2-dyb/2,
        (ix+0.5)*{X}+dx/2+dxb/2,
        (iy+0.5)*{Y}+dy/2+dyb/2,
        snrb
        from {intratbl}
        where r={r} and m={m} and s>{s0} and s<{s1}''')
        mpp = []
        for s2 in range(s0+1, s1):
            X1,Y1,X2,Y2,SNR = keepsome(s2==s, x1,y1,x2,y2,snr)
            mp = MatchPoints()
            mp.r1 = mp.r2 = r
            mp.m1 = mp.m2 = m
            mp.s1 = s2 - 1
            mp.s2 = s2
            if thr is None:
                dthr = dynamicthreshold(SNR)
            else:
                dthr = thr
            mp.xx1, mp.yy1, mp.xx2, mp.yy2 = keepsome(SNR>dthr, X1, Y1, X2, Y2)
            mpp.append(mp)
        return mpp

    def edge(r, m, s0, s1, thr=None):
        # Returns a list of MatchPoints with data from the edge table
        # for each of the slice pairs in [s0, s1).
        if s1==s0+1:
            return []
        (s, x2,y2, x1,y1, snr) = db.vsel(f'''select
        s,
        ix*{X}+x-dx/2-dxb/2,
        iy*{Y}+y-dy/2-dyb/2,
        ix*{X}+x+dx/2+dxb/2,
        iy*{Y}+y+dy/2+dyb/2,
        snrb
        from {edgetbl}
        where r={r} and m={m} and s>{s0} and s<{s1}''')
        mpp = []
        for s2 in range(s0+1, s1):
            X1,Y1,X2,Y2,SNR = keepsome(s2==s, x1,y1,x2,y2,snr)
            mp = MatchPoints()
            mp.r1 = mp.r2 = r
            mp.m1 = mp.m2 = m
            mp.s1 = s2 - 1
            mp.s2 = s2
            if thr is None:
                dthr = dynamicthreshold(SNR)
            else:
                dthr = thr
            mp.xx1, mp.yy1, mp.xx2, mp.yy2 = keepsome(SNR>dthr, X1, Y1, X2, Y2)
            mpp.append(mp)
        return mpp

def combine(mpi, mpe):
    mpr = []
    if len(mpi) != len(mpe):
        raise Exception('Combine needs equal length lists')
    for k in range(len(mpi)):
        mp = mpi[k]
        mp.combine(mpe[k])
        mpr.append(mp)
    return mpr

def find(mpp, r1=None, m1=None, s1=None, r2=None, m2=None, s2=None):
    res = []
    for mp in mpp:
        if r1 is None or mp.r1==r1:
            if m1 is None or mp.m1==m1:
                if s1 is None or mp.s1==s1:
                    if r2 is None or mp.r2==r2:
                        if m2 is None or mp.m2==m2:
                            if s2 is None or mp.s2==s2:
                                res.append(mp)
    return res
    
def index(mpp):
    # Given a list of MatchPoints objects, construct an index for matrixing
    k = 0
    idx = {} # Map from (r,m,s) to k
    for mp in mpp:
        k1 = (mp.r1, mp.m1, mp.s1)
        k2 = (mp.r2, mp.m2, mp.s2)
        if k1 not in idx:
            idx[k1] = k
            k += 1
        if k2 not in idx:
            idx[k2] = k
            k += 1
    return idx

def assignK(mpp):
    '''Given a list of MatchPoint objects, modifies each to assign a unique
    "k" value to each point within a tile.'''
    kkk = {} # Map from (r,m,s) to last-used k value.
    for mp in mpp:
        rms = mp.r1, mp.m1, mp.s1
        if rms not in kkk:
            kkk[rms] = 0
        kk = np.zeros(mp.xx1.shape, dtype=int)
        for n in range(len(mp.xx1)):
            kk[n] = kkk[rms]
            kkk[rms] += 1
        mp.kk1 = kk

        rms = mp.r2, mp.m2, mp.s2
        if rms not in kkk:
            kkk[rms] = 0
        kk = np.zeros(mp.xx2.shape, dtype=int)
        for n in range(len(mp.xx2)):
            kk[n] = kkk[rms]
            kkk[rms] += 1
        mp.kk2 = kk

def allpoints(mpp):
    '''Returns a map from (r,m,s) to a pair of x,y vectors organized by k.
    You must call assignK first.'''
    res = {} # map from rms to pairs of vectors
    for mp in mpp:
        rms = mp.r1, mp.m1, mp.s1
        if rms not in res:
            res[rms] = [np.zeros(0), np.zeros(0)]
        if len(mp.kk1)==0:
            print(f'Warning: no points in R{mp.r1} M{mp.m1} S{mp.s1} : R{mp.r2} M{mp.m2} S{mp.s2}  #1')
            K = 0
        else:
            K = np.max(mp.kk1)+1
        if res[rms][0].size < K:
            res[rms][0].resize(K)
            res[rms][1].resize(K)
        for n in range(len(mp.xx1)):
            res[rms][0][mp.kk1[n]] = mp.xx1[n]
            res[rms][1][mp.kk1[n]] = mp.yy1[n]

        rms = mp.r2, mp.m2, mp.s2
        if rms not in res:
            res[rms] = [np.zeros(0), np.zeros(0)]
        if len(mp.kk2)==0:
            print(f'Warning: no points in R{mp.r1} M{mp.m1} S{mp.s1} : R{mp.r2} M{mp.m2} S{mp.s2}  #2')
            K = 0
        else:
            K = np.max(mp.kk2)+1
        if res[rms][0].size < K:
            res[rms][0].resize(K)
            res[rms][1].resize(K)
        for n in range(len(mp.xx2)):
            res[rms][0][mp.kk2[n]] = mp.xx2[n]
            res[rms][1][mp.kk2[n]] = mp.yy2[n]
    return res

def elasticindex(mpp):
    # Given a list of MatchPoint objects, construct an index for matrixing
    # individual points. The index is from (r,m,s,k) to a matrix index.
    # See E&R p. 1615 for k. (On that page, (r,m,s) is summarized as "t".)
    p = 0 # Confusingly, the equivalent variable is called "k" in INDEX,
    # ... but I am reserving k here for the localized variable.
    idx = {} # Map from (r,m,s,k) to p, where p is a matrix entry
    for mp in mpp:
        for n in range(len(mp.xx1)):
            p1 = mp.r1, mp.m1, mp.s1, mp.kk1[n]
            p2 = mp.r2, mp.m2, mp.s2, mp.kk2[n]
            idx[p1] = p
            p += 1
            idx[p2] = p
            p += 1
    return idx

def elasticdeindex(idx, xx):
    '''idx must be a map of (r,m,s,k) to p, and xx must be a p-vector.
    Result is a map of (r,m,s) to k-vectors.'''
    res = {}
    for rmsk,p in idx.items():
        rms = rmsk[0], rmsk[1], rmsk[2]
        k = rmsk[3]
        if rms not in res:
            res[rms] = np.zeros(k+1)
        if res[rms].size < k + 1:
            res[rms].resize(k+1)
        res[rms][k] = xx[p]
    return res

def deindex(idx, xx):
    res = {}
    for k,v in idx.items():
        res[k] = xx[v]
    return res

def weightdx(mp, ax):
    w = len(mp.xx1) # Could be changed, of course
    if w==0:
        loc = f'R{mp.r1}M{mp.m1}S{mp.s1}:R{mp.r2}M{mp.m2}S{mp.s2}'
        if (mp.r1==35 and mp.s1==130 and mp.m1>=6) \
           or (mp.r2==35 and mp.s2==130 and mp.m2>=6):
            print(f'Ignoring lack of matchpoints {j} for {loc}')
            Dx = 0 
        else:
            raise Exception(f'matchpoints {j} empty for {loc}')
    else:
        Dx = np.mean(mp.xp(ax) - mp.x(ax))
    return w, Dx

def matrix(mpp, idx, ax):
    EPSILON = 1e-6
    K = len(idx)
    A = np.eye(K) * EPSILON
    b = np.zeros(K)
    j = 0
    for mp in mpp:
        k = idx[(mp.r1, mp.m1, mp.s1)]
        kp = idx[(mp.r2, mp.m2, mp.s2)]
        w, Dx = weightdx(mp, ax)
        A[k,k] += w
        A[kp,kp] += w
        A[k,kp] -= w
        A[kp,k] -= w
        b[k] += w*Dx
        b[kp] -= w*Dx
        j += 1
    return A, b

def tension(mpp, idx, xm, ym):
    res = {}
    for mp in mpp:
        rms1 = (mp.r1, mp.m1, mp.s1)
        rms2 = (mp.r2, mp.m2, mp.s2)
        k = idx[rms1]
        kp = idx[rms2]
        w, Dx = weightdx(mp, 0)
        w, Dy = weightdx(mp, 1)
        dx = xm[k] - xm[kp]
        dy = ym[k] - ym[kp]
        sx = dx - Dx
        sy = dy - Dy
        res[(rms1, rms2)] = (np.max(sx**2), np.max(sy**2))
    return res

def elastictension(mpp, idx, xm, ym):
    res = {}
    for mp in mpp:
        rms1 = (mp.r1, mp.m1, mp.s1)
        rms2 = (mp.r2, mp.m2, mp.s2)
        sxx = []
        syy = []
        for n in range(len(mp.xx1)):
            p = idx[(mp.r1, mp.m1, mp.s1, mp.kk1[n])]
            pp = idx[(mp.r2, mp.m2, mp.s2, mp.kk2[n])]
            print(rms1, rms2, n)
            Dx = mp.xp(0)[n] - mp.x(0)[n]
            Dy = mp.xp(1)[n] - mp.x(1)[n]
            dx = xm[p] - xm[pp]
            dy = ym[p] - ym[pp]
            sx = dx - Dx
            sy = dy - Dy
            sxx.append(sx)
            syy.append(sy)
        res[(rms1, rms2)] = (np.max(np.array(sxx)**2), np.max(np.array(syy)**2))
    return res

def elasticmatrix(mpp, idx, ap, ax):
    EPSILON = 1e-9
    K = len(idx)
    # Fill A with E_stable
    A = scipy.sparse.diags([np.ones(K) * EPSILON], [0], (K,K), "dok")
    b = np.zeros(K)
    j = 0
    # Next, add E_point
    for mp in mpp:
        if mp.s1==mp.s2 and mp.r1==mp.r2:
            w = 50
        else:
            w = 5
        for n in range(len(mp.xx1)):
            p = idx[(mp.r1, mp.m1, mp.s1, mp.kk1[n])]
            pp = idx[(mp.r2, mp.m2, mp.s2, mp.kk2[n])]
            Dx = mp.xp(ax)[n] - mp.x(ax)[n]
            #print(f'point {mp.r1},{mp.m1},{mp.s1}:{mp.r2},{mp.m2},{mp.s2} {mp.xx1[n]},{mp.yy1[n]} [{ax}]: {Dx}') 
            A[p,p] += w
            A[pp,pp] += w
            A[p,pp] -= w
            A[pp,p] -= w
            b[p] += w*Dx
            b[pp] -= w*Dx
    # Finally, add E_elast
    Q = 4
    def distfoo(dist2):
        D0 = 684 # Anything within D0 px should be taken very seriously
        return 1/(dist2/(D0*D0) + .1)
    for mp in mpp:
        rms = (mp.r1, mp.m1, mp.s1)
        xx = ap[rms][0]
        yy = ap[rms][1]
        for n in range(len(mp.xx1)):
            k = mp.kk1[n]
            #print('TEST', mp, n, k, xx[k], mp.xx1[n])
            dx = xx[k] - xx
            dy = yy[k] - yy
            dst = dx**2 + dy**2
            # Now, we need to find the Q points with least dst, not
            # counting k itself.
            if len(dst)>Q+1:
                kk = np.argpartition(dst, Q+1)
                kk = kk[:Q+1]
            else:
                kk = np.arange(len(dst))
            kk = kk[kk!=k]
            p = idx[(mp.r1, mp.m1, mp.s1, k)]
            for kstar in kk:
                pstar = idx[(mp.r1, mp.m1, mp.s1, kstar)]
                w = distfoo(dst[kstar])
                A[p,p] += w
                A[pstar,pstar] += w
                A[p,pstar] -= w
                A[pstar,p] -= w
                #print(f'elast {mp.r1},{mp.m1},{mp.s1}:{mp.r2},{mp.m2},{mp.s2} {mp.xx1[n]},{mp.yy1[n]}+{dx[kstar]},{dy[kstar]} [{ax}]: {w}') 

        rms = (mp.r2, mp.m2, mp.s2)
        xx = ap[rms][0]
        yy = ap[rms][1]
        for n in range(len(mp.xx2)):
            k = mp.kk2[n]
            dx = xx[k] - xx
            dy = yy[k] - yy
            dst = dx**2 + dy**2
            # Now, we need to find the Q points with least dst, not
            # counting k itself.
            kk = np.argsort(dst)
            kk = kk[1:]
            if len(kk)>Q:
                kk = kk[:Q]
            p = idx[(mp.r2, mp.m2, mp.s2, k)]
            for kstar in kk:
                pstar = idx[(mp.r2, mp.m2, mp.s2, kstar)]
                w = distfoo(dst[kstar])
                A[p,p] += w
                A[pstar,pstar] += w
                A[p,pstar] -= w
                A[pstar,p] -= w
    return A, b
    
