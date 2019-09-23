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
            while True:
                tsk = self.qu.get()
                if tsk is None:
                    self.qu.task_done()
                    print('Worker done')
                    return
                try:
                    print('Retrieved task', tsk)
                    self.produce(tsk)
                except Exception  as e:
                    ei = sys.exc_info()
                    tb = ei[2]
                    self.fails[tsk] = (e, tb)
                    print("Traceback:")
                    traceback.print_tb(tb)
                    print("Exception:\n ", e)
                    print('FAILED TO PRODUCE', tsk)
                self.qu.task_done()
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

    def produce(self, tsk):
        '''Subclasses should override this to do the work.
        Signal failure by raising an exception.'''
        pass

    def request(self, tsk):
        '''Clients call this to add a task to the queue.'''
        print('Scheduling', tsk)
        self.qu.put(tsk)

    def failures(self):
        '''At the end, return dict of failed requests.'''
        return self.fails

if __name__ == '__main__':
    import time

    # Define a subclass with an actual task
    class MyFactory(Factory):
        def produce(self, tsk):
            print('Dummy producing', tsk)
            time.sleep(.2*tsk) # This is a really easy task
            if tsk==5:
                # Just to show that failure is handled gracefully
                raise Exception('Insomnia')
            print('Done dummy producing', tsk)

    f = MyFactory()
    with f:
        # Scoping using a "with" statement ensures that all tasks are
        # completed at the end of the block.
        for k in range(10):
            f.request(k)
    print('Failures:')
    for k,v in f.failures().items():
        print(f'  {k}:', v[0])
    # If you don't anticipate any failures, you could also simply do:
    #   with MyFactory() as f:
    #     for k in range(10):
    #       f.request(k)

