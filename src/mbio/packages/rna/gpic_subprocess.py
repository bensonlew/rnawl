#!/mnt/ilustre/users/sanger/app/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "liubinxu"

import gevent
import time
import gipc
import gevent.subprocess as subprocess

spc = subprocess.Popen("ls -rtl", shell=True)
ret = spc.wait()
print ret

'''
with gipc.pipe() as (r, w):
    p = gipc.start_process(target=child_process, args=(r, ))
    p = gipc.start_process(target=random_process)
    wg = gevent.spawn(writegreenlet, w)
    try:
        p.join()
    except KeyboardInterrupt:
        wg.kill(block=True)
        p.terminate()
    p.join()
'''
