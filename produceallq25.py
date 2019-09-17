#!/usr/bin/python3

import produceq25
import copro
import em170428.runinfo
import queue
import rawimage
from pathlib import Path
import threading

ri = em170428.runinfo.RunInfo()

q = queue.Queue()
c0 = produceq25.ProduceQ25()

def worker():
    c = produceq25.ProduceQ25()
    while True:
        rms = q.get()
        try:
          c.ensure(rms)
        except Exception as e:
          print(e)
          print('FAILED TO PRODUCE', rms)
        q.task_done()

nthr = 4
threads = []
for n in range(nthr):
    t = threading.Thread(target=worker)
    t.start()
    threads.append(t)

def producesomeq25(r=None, m=None):
    if r is None:
        base = Path(rawimage.scaledtile(1,1,1,25)).parent.parent.parent
        if base.joinpath('complete').exists():
            return True
        ok = True
        for r in range(ri.runCount()):
            if not producesomeq25(r+1, m):
                ok = False
        q.join()
        if ok:
            base.joinpath('complete').touch()
        return ok
    elif m is None:
        print('Q25 - Working on run %i' % r)
        base = Path(rawimage.scaledtile(r,1,1,25)).parent.parent
        if base.joinpath('complete').exists():
            return True
        ok = True
        for m in range(ri.montageCount(r)):
            if not producesomeq25(r, m): 
                ok = False
        q.join()
        if ok:
            base.joinpath('complete').touch()
        return ok
    else:
        print('Q25 - Working on run %i montage %i' % (r, m))
        base = Path(rawimage.scaledtile(r,m,1,25)).parent
        if base.joinpath('complete').exists():
            return True
        mustcheck = False
        for s in range(ri.sliceCount(r)):
            if not c0.check((r,m,s)):
                print('Q25 - Queuing run %i montage %i slice %i' % (r, m, s))
                q.put((r,m,s))
                mustcheck = True
        q.join()
        if mustcheck:
            ok = True
            for s in range(ri.sliceCount(r)):
                if not c0.check((r,m,s)):
                    ok = False
                    break
        else:
           ok = True
        if ok:
            base.joinpath('complete').touch()
        return ok

producesomeq25()

