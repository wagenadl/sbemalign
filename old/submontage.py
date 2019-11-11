#!/usr/bin/python3

import aligndb

db = aligndb.DB()
ri = db.runinfo()

class SubMontage:
    def __init__(self, r, m, s0=0, s1=None):
        self.r = r
        self.m = m
        self.s0 = s0
        if s1 is None:
            self.s1 = ri.nslices(r)
        else:
            self.s1 = s1
    def __repr__(self):
        return f'SubMontage(R{self.r} M{self.m} [{self.s0}–{self.s1})'
    def fromSubVolume(r, m, z0, nz):
        z0r = ri.z0(r)
        S = ri.nslices(r)
        s0 = z0 - z0r
        if s0<0:
            s0 = 0
        if s0>=S:
            raise ValueError(f'Empty intersection R{r} Z{z}+{nz}')
        s1 = z0+nz - z0r
        if s1>S:
            s1 = S
        if s1<=0:
            raise ValueError(f'Empty intersection R{r} Z{z}+{nz}')
        return SubMontage(r, m, s0, s1)
    def allSubMontages(z0, nz):
        res = []
        for r0 in range(ri.nruns()):
            r = r0 + 1
            z0r = ri.z0(r)
            S = ri.nslices(r)
            if z0 < z0r+S and z0+nz > z0r:
                for m in range(ri.nmontages(r)):
                    res.append(SubMontage.fromSubVolume(r, m, z0, nz))
        return res
