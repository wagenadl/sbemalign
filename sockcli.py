import socket

fn = "/tmp/socket.sock"
print('prepping socket')
sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
print('done')
sock.connect(fn)
print('connected')
sock.sendall(bytes('hello\n', 'utf-8'))
print('sent')
rv = sock.recv(256)
print('recvd', str(rv, 'utf-8'))
sock.close()
