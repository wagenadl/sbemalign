#!/usr/bin/python3

import produceq25
import copro
import em170428.runinfo
import queue
import rawimage
from pathlib import Path

ri = em170428.runinfo.RunInfo()

nthr = 4
threads = []
for n in range(nthr):
    t = threading.Thread(target=worker)
    t.start()
    threads.append(t)
q = queue.Queue()

def worker():
    c = produceq25.ProduceQ25()
    while True:
        rms = q.get()
        c.ensure(rms)
        q.task_done()

def producesomeq25(r=None, m=None, s=None):
    if r is None:
        base = Path(rawimage.scaledtile(1,1,1,25)).parent.parent.parent
        if base.joinpath('complete').exists():
            return
        for r in range(ri.runCount()):
            print('Working on run %i', r+1)
            producesomeq25(r+1, m, s)
        q.join()
        base.joinpath('complete').touch()
    elif m is None:
        base = Path(rawimage.scaledtile(r,1,1,25)).parent.parent
        if base.joinpath('complete').exists():
            return
        for m in range(ri.montageCount(r)):
            print('Working on run %i montage %i', r+1, m)
            producesomeq25(r, m, s)
        q.join()
        base.joinpath('complete').touch()
    elif s is None:
        base = Path(rawimage.scaledtile(r,1,1,25)).parent
        if base.joinpath('complete').exists():
            return
        for m in range(ri.sliceCount(r)):
            print('Working on run %i montage %i slice %i', r+1, m, s)
            producesomeq25(r, m, s)
        q.join()
        base.joinpath('complete').touch()
    else:
        q.put((r,m,s))
