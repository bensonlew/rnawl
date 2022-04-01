import os
import sys

wf = sys.argv[1]
qstat = os.path.join(wf, 'HiseqReadsStat')
qstat1 = os.path.join(wf, 'HiseqReadsStat1')

def rename(dir):
    print(dir)
    for file in os.listdir(dir):
        print(file)
        if os.path.isdir(os.path.join(dir, file)):
            rename(os.path.join(dir, file))
        else:
            if len(os.path.basename(file).split('.')) > 2:
                os.rename(os.path.join(dir,file), os.path.join(dir,os.path.basename(file).split('.')[-3]+'.' + os.path.basename(file).split('.')[-2]+'.' + os.path.basename(file).split('.')[-1]))

rename(qstat)
rename(qstat1)