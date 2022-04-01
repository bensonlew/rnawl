import os
from biocluster.config import Config
from multiprocessing import Pool
import subprocess

# def run_script(codes):
#     subprocess.call("Rscript {}".format(codes), shell=True)

class nomogram(object):
    def __init__(self, surv_table, factor_list, time_list, output):
        self.surv_table = surv_table
        self.factor_list = factor_list.strip().split(';')
        self.time = time_list.strip().split(';')
        self.output = output
        software_dir = Config().SOFTWARE_DIR
        self.rscript = software_dir + "/program/R-3.5.1/bin/Rscript"

    # def run(self, script_list):
    #     pool = Pool(self.pool_size)
    #     pool.map(run_script, script_list)
    #     pool.close()
    #     pool.join()

    def nomogram(self, output=None):
        if output is None:
            output = os.getcwd()
        script_name = os.path.join(output, 'nomogram.r')
        f = open(script_name, 'w')
        f.write('library(survival)\n')
        f.write('library(rms)\n')
        f.write('surv_cox <- read.table("{}", sep="\t", header=TRUE)\n'.format(self.surv_table))
        f.write('dd <- datadist(surv_cox)\n')
        f.write('options(datadist="dd")\n')
        factor = '+'.join(self.factor_list)
        f.write('res.cox <- cph(Surv(time, status) ~ {}, data =  lung, surv = T, x=T, y=T)\n'.format(factor))
        f.write('survival <- Survival(res.cox)\n')
        sur = list()
        for i in self.time:
            f.write('survival{} <- function(x)survival({}, x)\n'.format(i, i))
            sur.append('survival{}'.format(i))
        sur_str = ','.join(sur)
        f.write('nom = nomogram(res.cox, fun=list({}), fun.at = c(0.05, seq(0.1,0.9,by=0.05), 0.95), funlabel = c("1 year", "2 year"))\n'.format(sur_str))
        # f.write('write.table(nom, file = "{}", sep="\t")\n'.format(self.output))
        f.write('pdf(file="nomogram.pdf")\n')
        f.write('a = plot(nom)\n')
        f.write('dev.off()\n')
        f.close()
        os.system('{} nomogram.r'.format(self.rscript))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to identify pirna from piRbase database.')
    parser.add_argument('-i', type=str, help='surv file')
    parser.add_argument('-fl', type=str, help='factor list')
    parser.add_argument('-tl', type=str, help='time list')
    parser.add_argument('-o', type=str, help='output path')
    args = parser.parse_args()
    nomogram = nomogram(args.i, args.fl, args.tl, args.o)
    nomogram.nomogram()
