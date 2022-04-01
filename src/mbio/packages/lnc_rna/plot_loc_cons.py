# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
from concurrent.futures import ThreadPoolExecutor
import os
import subprocess

parser = OptionParser(description='Draw location conservatation line plot from bwtool out file')
parser.add_option('-i', '--input', dest='input', help='input bwtool out file')
parser.add_option('-t', '--thread', dest='thread', type=int, help='number of plotting threads to launch')
parser.add_option('-n', '--cutoff', dest='cutoff', type=float, help='drop records with valid points less than cutoff rate')
parser.add_option('-r', '--interpreter', dest='interpreter', help='program interpreter for Rscript')
parser.add_option('-s', '--script', dest='script', help='script file used for drawing')
parser.add_option('-c', '--color', dest='color', help='color of line in figure')
parser.add_option('-k', '--retain', dest='retain', help='output file containing left name')
parser.add_option('-o', '--output', dest='output', help='output directory containing plots')
(opts, args) = parser.parse_args()

def main(bwtool_out, thread, cutoff, file_out, dir_out):
    pools = list()
    with ThreadPoolExecutor(max_workers=thread) as pool, open(file_out, 'w') as f:
        print 'INFO: start reading {}'.format(bwtool_out)
        for n, line in enumerate(open(bwtool_out)):
            items = line.strip().split('\t')
            scores = items[7].split(',')
            if scores.count('NA') > len(scores) * (1 - cutoff):
                print 'INFO: skip {} on {}'.format(items[3], items[0])
            else:
                pools.append(pool.submit(draw, items))
                f.write('{}\n'.format(items[3]))
    print 'INFO: succeed in exporting {} figures in {}'.format(len(pools), dir_out)
    if os.path.getsize(file_out) > 0:
        print 'INFO: succeed in exporting {}'.format(file_out)

def draw(items):
    tsv = os.path.join(opts.output, '{}.tsv'.format(items[3]))
    with open(tsv, 'w') as f:
        for n, i in enumerate(items[7].split(',')):
            f.write('{}\t{}\n'.format(n + 1, i))
    pdf = os.path.join(opts.output, '{}.pdf'.format(items[3]))
    cmd = '{} {}'.format(opts.interpreter, opts.script)
    cmd += ' -i {}'.format(tsv)
    cmd += ' -t {}'.format(items[3])
    cmd += ' -x {}'.format('"position: {} {}-{}"'.format(items[0], items[1], items[2]))
    cmd += ' -y {}'.format('"phastCons score by site"')
    cmd += ' -c {}'.format(opts.color)
    cmd += ' -o {}'.format(pdf)
    spc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = spc.wait()
    if ret:
        print 'WARN: {}'.format(cmd)
        raise Exception('ERROR: failed to export {}'.format(pdf))
    else:
        print 'INFO: succeed in exporting {}'.format(pdf)

if __name__ == '__main__':
    if opts.input and opts.thread and opts.cutoff and opts.interpreter and opts.script and opts.color and opts.retain and opts.output:
        main(opts.input, opts.thread, opts.cutoff, opts.retain, opts.output)
    else:
        parser.print_help()
