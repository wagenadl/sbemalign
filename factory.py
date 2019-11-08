#!/usr/bin/python3

'''FACTORY - Wrapper around queue and threading'''

import threading
import queue
import sys
import traceback

class Factory:
    def __init__(self, nthr=4):
        self.nthr = nthr
        self.fails = {}
        def wrkr():
            try:
                while True:
                    tsk = self.qu.get()
                    if tsk is None:
                        self.qu.task_done()
                        #print('Worker done')
                        return
                    try:
                        #print('Retrieved task', tsk)
                        self.produce(*tsk)
                    except Exception  as e:
                        ei = sys.exc_info()
                        tb = ei[2]
                        self.fails[tsk] = (e, tb)
                        print("Traceback:")
                        traceback.print_tb(tb)
                        print("Exception:\n ", e)
                        print('FAILED TO PRODUCE', tsk)
                    self.qu.task_done()
            except KeyboardInterrupt:
                sys.exit(1)
        self.thr = []
        self.qu = queue.Queue(nthr)
        for n in range(nthr):
            t = threading.Thread(target=wrkr)
            t.start()
            self.thr.append(t)

    def __enter__(self):
        while len(self.thr) < self.nthr:
            t = threading.Thread(target=wrkr)
            t.start()
            self.thr.append(t)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        #print('Exit factory')
        if exc_type is not None:
            print('Factory shutting down with exception', exc_type, exc_val)
        self.shutdown()
        
    def __del__(self):
        #print('Delete factory')
        self.shutdown()
        
    def shutdown(self):
        #print('Sending None tasks')
        for n in range(len(self.thr)):
            self.qu.put(None)
        #print('Waiting for queue to clear')
        self.qu.join()
        for t in self.thr:
            #print('Joining thread')
            t.join()
        self.thr = []

    def produce(self, *tsk):
        '''Subclasses should override this to do the work.
        Signal failure by raising an exception.'''
        tsk = list(tsk)
        foo = tsk.pop(0)
        foo(*tsk)

    def request(self, *tsk):
        '''Clients call this to add a task to the queue.'''
        if self.nthr==0:
            self.produce(*tsk)
        else:
            try:
                self.qu.put(tsk)
            except KeyboardInterrupt:
                sys.exit(1)

    def failures(self):
        '''At the end, return dict of failed requests.'''
        return self.fails

if __name__ == '__main__':
    import time

    ################## Function based approach ##########################
    def myfunc(a, b):
        print(f'Working on adding {a} and {b}.')
        time.sleep(.1)
        c = a + b
        print(f'Function says {a} plus {b} makes {c}.')
        
    # Scoping using a "with" statement ensures that all tasks are
    # completed at the end of the block:
    with Factory() as f:
        for k in range(5):
            for l in range(5):
                f.request(myfunc, k, l)

    ################### Class based approach ######################
    class MyFactory(Factory):
        def produce(self, a, b):
            print(f'Working on adding {a} and {b}.')
            time.sleep(.1)
            c = a + b
            print(f'Class says {a} plus {b} makes {c}.')

    with MyFactory() as f:
        for k in range(5):
            for l in range(5):
                f.request(k, l)

    #################### Capturing failures ######################
    def failfunc(a, b):
        print(f'Working on adding {a} and {b}.')
        time.sleep(.1)
        if a==1 and b==1:
            raise Exception(f'I can never remember how to add {a} and {b}.')
        c = a + b
        print(f'Error-prone function says {a} plus {b} makes {c}.')
    f1 = Factory()
    with f1:
        for k in range(5):
            for l in range(5):
                f1.request(failfunc, k, l)
    print('Failures:')
    for k,v in f1.failures().items():
        print(f'  {k}:', v[0])
            
    # Of course, failures can be captured the same way in the class-based
    # approach.
