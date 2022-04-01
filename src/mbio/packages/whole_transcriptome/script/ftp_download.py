# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import logging
import os

from ftplib import FTP

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

ftp = FTP()
ftp.set_debuglevel(2)
ftp.connect('193.62.192.7', '21')
ftp.login('anonymous', '')
ftp.cwd('vol1/fastq/SRR103/000/SRR1039510/')
filenamelist = ftp.nlst()
for filename in filenamelist:
    localfilepath = os.path.join(os.getcwd(), filename)
    with open(localfilepath, 'w') as handle:
        ftp.retrbinary('RETR {}'.format(filename), handle.write, 1024)
else:
    logging.info('succeed in downloading {} to {}'.format(os.getcwd()))
