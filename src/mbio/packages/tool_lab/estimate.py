import os
from biocluster.config import Config
from multiprocessing import Pool
import subprocess
import pandas as pd


# def run_script(codes):
#     subprocess.call("Rscript {}".format(codes), shell=True)

class estimate(object):
    def __init__(self, exp_table, id, platform, output):
        self.exp_table = exp_table
        self.id = id
        self.output = output
        self.platform = platform
        software_dir = Config().SOFTWARE_DIR
        self.rscript = software_dir + "/program/R-3.3.1/bin/Rscript"
        # self.rscript = '/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/miniconda/bin/Rscript'
    # def run(self, script_list):
    #     pool = Pool(self.pool_size)
    #     pool.map(run_script, script_list)
    #     pool.close()
    #     pool.join()

    def estimate(self, output=None):
        if output is None:
            output = os.getcwd()
        exp_matrix_new = os.path.join(output, 'exp_matrix_new.txt')
        estimate_gct = os.path.join(self.output, "estimate.gct")
        estimate_score_gct = os.path.join(self.output, "estimate_score.gct")
        script_name = os.path.join(output, 'estimate.r')
        f = open(script_name, 'w')
        f.write('library(estimate)\n')
        f.write('exp_table <- read.table("{}", sep="\t",  header=TRUE, row.names=1)\n'.format(self.exp_table))
        f.write('write.table(file="{}", exp_table, sep="\t", row.names=T, col.names=T, quote=F)\n'.format(exp_matrix_new))
        f.write('filterCommonGenes(input.f="{}",output.f="{}",id="{}")\n'.format(exp_matrix_new, estimate_gct, self.id))
        f.write('estimateScore("{}","{}",platform="{}")\n'.format(estimate_gct, estimate_score_gct, self.platform))
        f.close()
        os.system('{} estimate.r'.format(self.rscript))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to identify pirna from piRbase database.')
    parser.add_argument('-e', type=str, help='exp file')
    parser.add_argument('-i', type=str, help='GeneSymbol or Entrez')
    parser.add_argument('-p', type=str, help='affymetrix or agilent or illumina')
    parser.add_argument('-o', type=str, help='output path')
    args = parser.parse_args()
    estimate = estimate(args.e, args.i, args.p, args.o)
    estimate.estimate()
