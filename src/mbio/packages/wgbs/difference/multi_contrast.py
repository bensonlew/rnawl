# -*- coding:utf-8 -*-

import ConfigParser
import logging
import os
import subprocess
import sys

import pandas as pd

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

rcp = ConfigParser.RawConfigParser()
rcp.read(os.path.join(os.environ['HOME'], 'biocluster/src/biocluster/main.conf'))
rcp.get('Command', 'software_dir')

RS_BIN = os.path.join(rcp.get('Command', 'software_dir'), 'bioinfo/rconda/bin/Rscript')

RS = r'''library(DSS)
require(bsseq)
dat.case <- read.table("{case_fp}", header = T)
dat.ctrl <- read.table("{ctrl_fp}", header = T)
BSobj <- makeBSseqData(list(dat.case, dat.ctrl), c("{case_name}", "{ctrl_name}"))
dmlResult <- DMLtest(BSobj, group1 = c("{case_name}"), group2 = c("{ctrl_name}"), smoothing = T)
dmls <- callDML(dmlResult, p.threshold = 0.05)
write.table(dmls, "{out_fp}", quote = F, sep = "\t", row.names = F)
'''


def parse_ini(filepath):
    p_dict = dict()
    for line in open(filepath):
        eles = line.strip().split('\t')
        if len(eles) == 2:
            p_dict[eles[1]] = eles[0]
    return p_dict


def make_pairs(groups):
    pairs = list()
    for i in range(0, len(groups)):
        case = groups[i]
        for j in range(i + 1, len(groups)):
            ctrl = groups[j]
            pairs.append({'case': case, 'ctrl': ctrl})
    return pairs


def do_contrast(pair, p_dict, out_dir):
    case_name = pair['case']
    ctrl_name = pair['ctrl']
    out_fp = os.path.join(out_dir, '{case}_vs_{ctrl}.txt'.format(**pair))
    kwargs = {'case_fp': p_dict[case_name], 'ctrl_fp': p_dict[ctrl_name], 'case_name': case_name,
              'ctrl_name': ctrl_name, 'out_fp': out_fp}
    if not os.path.isfile(out_fp):
        rscript_fp = os.path.join(out_dir, '{case}_vs_{ctrl}.R'.format(**pair))
        open(rscript_fp, 'w').write(RS.format(**kwargs))
        proc = subprocess.Popen([RS_BIN, rscript_fp], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        outs, errs = proc.communicate()
        returncode = proc.poll()
        logger.error(errs)
        if returncode:
            raise subprocess.CalledProcessError(returncode, ' '.join((RS_BIN, rscript_fp)), outs)
    df = pd.read_table(out_fp)
    dml_set = set(df['chr'] + '_' + df['pos'].astype(str))
    logger.debug('Get {} DML(s)'.format(len(dml_set)))
    return dml_set


def export_result(src_table, dml_set, out_fp):
    df = pd.read_table(src_table)
    df['idx'] = df['chrom'] + '_' + df['pos'].astype(str)
    df.set_index('idx', inplace=True)
    dml_df = df.reindex(dml_set)
    dml_df.to_csv(out_fp, sep='\t', index=False)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('usage: {} <list FILE> <result FILE> <output DIR>'.format(os.path.basename(sys.argv[0])))
    list_fp = sys.argv[1]
    result_fp = sys.argv[2]
    output_dir = sys.argv[3]

    logger.debug('Parse {}'.format(list_fp))
    path_dict = parse_ini(list_fp)
    logger.debug('Generate contrast pairs')
    vs_pairs = make_pairs(path_dict.keys())
    dmls = set()
    logger.debug('Start call DML for each pair')
    for i, vs_pair in enumerate(vs_pairs):
        logger.debug('Call DML of {} at {}/{}'.format(vs_pair, i + 1, len(vs_pairs)))
        dmls.update(do_contrast(vs_pair, path_dict, output_dir))
    logger.debug('Export {}'.format(os.path.join(output_dir, 'dml.result.txt')))
    export_result(result_fp, dmls, os.path.join(output_dir, 'dml.result.txt'))
