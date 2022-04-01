# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import re

RESULTS = {
    'A3SS.MATS.JC.txt',
    'A3SS.MATS.JCEC.txt',
    'A5SS.MATS.JC.txt',
    'A5SS.MATS.JCEC.txt',
    'JC.raw.input.A3SS.txt',
    'JC.raw.input.A5SS.txt',
    'JC.raw.input.MXE.txt',
    'JC.raw.input.RI.txt',
    'JC.raw.input.SE.txt',
    'JCEC.raw.input.A3SS.txt',
    'JCEC.raw.input.A5SS.txt',
    'JCEC.raw.input.MXE.txt',
    'JCEC.raw.input.RI.txt',
    'JCEC.raw.input.SE.txt',
    'MXE.MATS.JC.txt',
    'MXE.MATS.JCEC.txt',
    'RI.MATS.JC.txt',
    'RI.MATS.JCEC.txt',
    'SE.MATS.JC.txt',
    'SE.MATS.JCEC.txt',
    'fromGTF.A3SS.txt',
    'fromGTF.A5SS.txt',
    'fromGTF.MXE.txt',
    'fromGTF.RI.txt',
    'fromGTF.SE.txt',
    'fromGTF.novelEvents.A3SS.txt',
    'fromGTF.novelEvents.A5SS.txt',
    'fromGTF.novelEvents.MXE.txt',
    'fromGTF.novelEvents.RI.txt',
    'fromGTF.novelEvents.SE.txt'
}

os.system(r'grep 命令 Rmats*/*.err > command.txt')

pwd = os.getcwd()
if '/sanger/' in pwd:
    p = 'SANGER'
elif 'isanger' in pwd:
    p = 'ISANGER'

heads = [
    '#!/bin/bash\n',
    '#SBATCH -c 20\n',
    '#SBATCH -p {}\n'.format(p),
    '#SBATCH --mem=200G\n',
]

commands = set()

for line in open('command.txt'):
    m = re.search('.*命令内容为(.*)', line)
    if m:
        cmd = m.group(1)
        mm = re.search(r'--od (\S+) -t', cmd)
        if mm:
            od = mm.group(1)
            if RESULTS == RESULTS & set(os.listdir(od)):
                print od
            else:
                cmd = cmd.replace('--nthread 10 --tstat 10', '--nthread 20 --tstat 20')
                commands.add(cmd + '\n')

open('run_rmats.slurm.sh', 'w').writelines(heads + list(commands))
