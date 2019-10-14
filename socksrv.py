#!/usr/bin/python3

import socket
import queue
import os
import threading
import time

t00 = time.time()

fn = "/tmp/socket.sock"
try:
    os.unlink(fn)
except OSError:
    if os.path.exists(fn):
        raise
    
sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
sock.bind(fn)
sock.listen(4)

def worker():
    while True:
        item = q.get()
        con, adr, k = item
        print('Got item', k, time.time() - t00)
        txt = con.recv(256)
        print('Got: ', str(txt, 'utf-8'))
        print('sleeping')
        time.sleep(10)
        print('woke up')
        con.send(bytes('Reply\n', 'utf-8'))
        print('closed')
        con.close()
        print('done', k, time.time() - t00)
        q.task_done()

q = queue.Queue()
nthr = 4
threads = []
for n in range(nthr):
    t = threading.Thread(target=worker)
    t.start()
    threads.append(t)

k = 0    
while True:
    print('Waiting for connection')
    con, adr = sock.accept()
    k += 1
    print('Got connection', k, time.time() - t00)
    item = (con, adr, k)
    q.put(item)
    print('Queued')
    
