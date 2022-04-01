# -*- coding: utf-8 -*-
import socket

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect(("127.0.0.1",6666))
print s.getsockname()

s.send('shi')
"""
# 接收数据:
buffer = []

while True:
    # 每次最多接收1k字节:
    d = s.recv(1024)
    if d:
        buffer.append(d)
        print "yes"
    else:
        break
data = ''.join(buffer)

s.close()

print data
"""
print s.recv(1024)
print s.recv(1024)
s.send('exit')
s.close()