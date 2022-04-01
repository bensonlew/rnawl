# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

import os
import re

os.system('grep 命令 */*.err > command.txt')
heads = [
    '#!/bin/bash\n',
    '#SBATCH -c 20\n',
    '#SBATCH -p ISANGER\n',
    '#SBATCH --mem=80G\n',
]
commands = list()
for line in open('command.txt'):
    m = re.match(r'^.*?(/mnt/i?lustre.*)$', line)
    if m:
        command = m.group(1)
        mm = re.match(r'.*-1\s(.*?)\s-2.*', command)
        if mm:
            r1_path = mm.group(1)
            r1_fq = os.path.basename(r1_path)
            r1_name = r1_fq[:-14]
            commands.append('{} {}.sam --un-conc {}\n'.format(
                command.replace('-p 8', '-p 20')[:-26], r1_name, r1_name)
            )
        else:
            raise Exception('ERROR: can not process {}'.format(line))
else:
    with open('run_hisat2.slurm.sh', 'w') as f:
        f.writelines(heads + commands)
