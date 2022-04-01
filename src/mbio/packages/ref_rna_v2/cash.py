import os
import subprocess
import re
import argparse

class cash():
    def __init__(self, case_bam=str, control_bam=str,output=None,
                 java_path=None, cash_path=None, group_file= None,
                 gtf=None, samtools_path=None):
        self.case_bam = case_bam
        self.control_bam = control_bam
        self.gtf = gtf
        self.output = output
        self.java_path = java_path
        self.cash_path = cash_path
        self.samtool_path = samtools_path
        self.group_file = group_file

#    def get_group_dict(self):
#       group_id = ''
#        sample_id_list = []
#        group_dict = {group_id: sample_id_list}
#        with open(self.group) as f:
#            for line in f:
#                group_id = line.split()[0]
#                if group_id not in group_dict.keys():
#                    sample_id_list = line.split()[0:1]
#                    group_dict[group_id] = sample_id_list
#                else:
#                    sample_id_list.append(line.split()[0:1])
#                    group_dict[group_id] = sample_id_list
#        return group_dict
    def get_group_dict(self):
        group_dict = {}
        with open(self.group_file,'r') as f:
            for line in f:
                if line.startswith('#'):
                    pass
                else:
                    group_dict[line.split()[0]] = line.split()[1]
        if os.path.basename(self.control_bam.split(",")).split(".")[0] not in group_dict.keys():
            raise Exception("bam not match group file")
        if os.path.basename(self.case_bam.split(",")).split(".")[0] not in group_dict.keys():
            raise Exception("bam not match group file")
        return group_dict


    def no_name_gtf(self):
        id = ''
        detail = ''
        dict = {id : detail}
        g2n_all = []
        fw = open('temp.gtf','w')
        with open(self.gtf, 'r') as fr:
            for line in fr:
                new_line = re.sub('gene_name.*','',line)
                fw.write(new_line)
                for items in line.split("\t")[8].split('";'):
                    id = items.split(' "')[0]
                    detail = items.split(' "')[1]
                    dict[id] = detail
                if not dict.has_key('gene_name'):
                    g2n = {dict['gene_id']: "_"}
                else:
                    g2n = {dict['gene_id']: dict['gene_name']}
                if g2n not in g2n_all:
                    g2n_all.append(g2n)
        fw.close()
        return g2n_all

    def check_bai(self):
        for bam_A in self.case_bam.split(","):
            if os.path.exists(bam_A + ".bai"):
                pass
            else:
                cmd = '{} '.format(self.samtool_path)
                cmd += 'index {}'.format(bam_A)
                subprocess.check_call(cmd,shell=True)
        for bam_B in self.control_bam.split(","):
            if os.path.exists(bam_B + ".bai"):
                pass
            else:
                cmd = '{} '.format(self.samtool_path)
                cmd += 'index {}'.format(bam_B)
                subprocess.check_call(cmd,shell=True)

    def run_cash_cmd(self):
        group_dict = self.get_group_dict()
        control_name = group_dict[os.path.basename(self.control_bam.split(",")).split(".")[0]]
        case_name = group_dict[os.path.basename(self.case_bam.split(",")).split(".")[0]]
        cmd = '{} '.format(self.java_path)
        cmd += '-jar -Xmx10g {} '.format(self.cash_path)
        cmd += '--Case:{} '.format(case_name)
        cmd += '{} '.format(self.case_bam)
        cmd += '--Control:{} '.format(control_name)
        cmd += '{} '.format(self.control_bam)
        cmd += '--GTF {} '.format('temp.gtf')
        cmd += '--Output {}'.format(self.output)
        subprocess.check_call(cmd,shell=True)
        return control_name, case_name


    def outfile_generat(self):
        g2n_all = self.no_name_gtf()
        control, case = self.run_cash_cmd()
        as_even_list = ['Cassette', 'Cassette_multi', 'A5SS', 'A3SS',
                        'AltEnd', 'AltStart', 'MXE', 'IR']
        f1 = open(self.output + '.' + control + 'vs' +
                 case + '.alldiff.statistics.txt','r')
        with open(control + 'vs' + case
                  + 'as_stat.txt','w') as f_stat:
            f_stat.write("AS_type\tDiff_num\tAll_num\n")
            for line in f1:
                if line.split()[0] in as_even_list:
                    f_stat.write(line)
        f1.close()
        f2 = open(self.output + '.' + case + 'vs' +
                  control + '.alldiff.detail.txt','r')
        with open(control + 'vs' + case
                  + 'as_detail.txt','w') as f_detail:
            f_detail.write("AS_ID\tGene_ID\tGene_Name\tAStype\tchr\tevent_start\tevent_end\t"
                           + case + "_Junc_Inclusive\t" + case
                           + "_Junc_Exclusive\t" + control + "_Junc_Inclusive\t"
                           + control + "_Junc_Exclusive\t" + "delta_PSI"
                           + "P-Value\tFDR")
            _ = f2.readline()
            for (num,line) in enumerate(f2):
                line_list = re.split('[\t\s:]+',line)
                new_list = list(line_list[-1] + num)
                new_list += line_list[0:1]
                new_list += list(g2n_all[line_list[0]])
                new_list += line_list[-1:] +line_list[1:2]
                new_list += list(line_list[2].split("-"))
                new_list += line_list[4:15]
                f_detail.write("\t".join(new_list) + "\n")
            f2.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', type=str, metavar='group_file', required=True,
                        help="the same format as i-sanger groupfile")
    parser.add_argument('-case_bam',type=str,metavar='path_of_case_bam',required=True,
                        help="path of case bams splited by comma")
    parser.add_argument('-control_bam',type=str,metavar='path_of_control_bam',required=True,
                        help="path of control bams splited by comma")
    parser.add_argument('-gtf',type=str,metavar='gtf_file',required=True,
                        help="gtf file can be obtained by gffread or database")
    parser.add_argument('-java',type=str,metavar='jave_path',required=True,
                        default="/mnt/ilustre/users/sanger-dev/app/program/sun_jdk1.8.0/bin/java")
    parser.add_argument('-tool',type=str,metavar='cash_path',required=True,
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/cash_v2.2.1/cash.jar")
    parser.add_argument('-o',type=str,metavar='outfile_prefix',required=True,
                        default="output")
    parser.add_argument('-samtools',type=str,metavar='path_of_samtools',required=True,
                        default="/mnt/ilustre/users/sanger-dev/app/program/Python/bin/samtools")
    args = parser.parse_args()
    toolbox=cash(
        case_bam=args.case_bam,
        control_bam=args.control_bam,
        output=args.o,
        java_path=args.java,
        cash_path=args.tool,
        group_file=args.g,
        gtf=args.gtf,
        samtools_path=args.samtools
    )
    toolbox.check_bai()
    toolbox.outfile_generat()