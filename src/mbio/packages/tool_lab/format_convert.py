# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
import argparse
from Bio import SeqIO
from biocluster.config import Config
from multiprocessing import Pool
import subprocess


software_dir = Config().SOFTWARE_DIR
seqkit_path = os.path.join(software_dir, "bioinfo/meta/seqkit/seqkit")
#transeq = os.path.join(software_dir, "bioinfo/seq/EMBOSS-6.6.0/bin/transeq")

def change_data(infile, out):
    infile_path = infile
    all_files = os.listdir(infile_path)
    sample_list = []

    for file in all_files:
        if file.endswith('.fna'):
            file_name = file.split('.fna')[0]
            fna_tmp = infile_path + '/' + file
            faa_tmp = infile_path+'/'+file_name+'.faa'
            if not os.path.exists(faa_tmp):
                # cmd = '%s -sequence %s -table %s -trim -outseq %s'%(transeq, fna_tmp,11 ,faa_tmp)
                # os.system(cmd)
                r_seq = SeqIO.parse(fna_tmp, "fasta")  #parse
                with open(faa_tmp,'w') as fw:
                    for e in r_seq:
                        try:
                            faa_seq = e.seq.translate(table="Bacterial", cds=True)
                            fw.write('>'+e.description+'\n')
                            fw.write(str(faa_seq)+'\n')
                        except Exception as e:
                            print(e)


    new_all_files = os.listdir(infile_path)
    for file in new_all_files:
        file_path = os.path.join(infile_path, file)
        #old_file_name = file_path.split('/')[-1]
        if file.endswith('.faa'):
            file_name = file.split('.faa')[0]
            print file_name
            if file_name not in sample_list:
                sample_list.append(file_name)
            outfile_path = os.path.join(out, file_name + '.pep')

            #下面是为了生成新的文件.function文件，pgap软件需要检查
            function_path = os.path.join(out, file_name + '.function')
            with open(function_path, 'w') as w:
                for seq_record in SeqIO.parse(file_path, 'fasta'):
                    seq_id = seq_record.id
                    w.write('{}\t{}\n'.format(seq_id, '-'))
            #下面是为了生成新的pep文件，将文件名称加到序列id前面
            with open(outfile_path, 'w') as outf:
                for seq_record in SeqIO.parse(file_path, 'fasta'):
                    seq_id = seq_record.id
                    seq = seq_record.seq
                    outf.write('>{}|{}\n{}\n'.format(file_name,seq_id, seq))
        elif file.endswith('.fna'):
            file_name = file.split('.fna')[0]
            outfile_path = os.path.join(out, file_name + '.nuc')
            #下面是为了生成新的nuc文件，将文件名称加到序列id前面
            with open(outfile_path, 'w') as outm:
                for seq_record in SeqIO.parse(file_path, 'fasta'):
                    seq_id = seq_record.id
                    seq = seq_record.seq
                    outm.write('>{}|{}\n{}\n'.format(file_name,seq_id, seq))

    return sample_list


def seqkit_command(file1, file2, output):
    """
    seqkit软件得到相同id的核酸和蛋白文件
    :param file1: 文件1核酸
    :param file2: 文件2蛋白
    :param output:结果目录输出
    :return:
    """
    cmd1 = "{} common {} {} -o {} -w 0".format(seqkit_path, file1, file2, output + ".nuc")
    print(cmd1)
    try:
        subprocess.check_output(cmd1, shell=True)
        print("seqkit运行完成!")
    except subprocess.CalledProcessError:
        raise Exception("seqkit运行失败!")
    cmd2 = "{} common {} {} -o {} -w 0".format(seqkit_path, file2, file1, output + ".pep")
    print(cmd2)
    try:
        subprocess.check_output(cmd2, shell=True)
        print("seqkit运行完成!")
    except subprocess.CalledProcessError:
        raise Exception("seqkit运行失败!")


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[input_dir]',required=True, help='input dir')
    parse.add_argument('-out', metavar='[output_dir]', required=True, help='median out dir')
    parse.add_argument('-o', metavar='[output_dir]', required=True, help='input output dir')
    parse.add_argument('-t', metavar='[thread]', required=True, help='input thread number')
    args = parse.parse_args()
    infile = args.i
    out = args.out
    output = args.o
    thread = args.t
    sample_list = change_data(infile, out)
    p = Pool(int(thread))
    for i in sample_list:
        file1_path = os.path.join(out, i + ".nuc")
        file2_path = os.path.join(out, i + ".pep")
        output_path = os.path.join(output, i)
        p.apply_async(seqkit_command, args=(file1_path,file2_path,output_path, ))
    p.close() # 需要先关闭，才能join
    p.join()


if __name__ == '__main__':
    main()
