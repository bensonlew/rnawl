# -*- coding:utf-8 -*-
__author__ = 'guanqing.zou'

import os
import sys
import datetime

def divide_by_seq_num(infile,unit_num=1000, output_dir='divide_result_dir'):
    divide_files_list = []
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    n = 0
    out_file = output_dir+'/'+str(n)+'.fasta'
    fw = open(out_file,'w')
    divide_files_list.append(out_file)
    with open(infile) as f:
        for line in f:
            if line[0] == '>':
                n += 1
                if n > unit_num:
                    fw.close()
                    out_file = output_dir+'/'+str(n)+'.fasta'
                    fw = open(out_file,'w')
                    divide_files_list.append(out_file)
            fw.write(line)

    fw.close()
    return divide_files_list


def divide_fun_by_size(infile,unit_size,output_dir='divide_result_dir'):
    divide_files_list = []
    max_size = unit_size  #1M
    all_size = os.path.getsize(infile)
    divide_num = all_size/float(max_size)
    if divide_num < int(divide_num)+0.5 :
        divide_num = int(divide_num)
    else:
        divide_num = int(divide_num) + 1

    if divide_num > 1:
        divide_dir = output_dir
        if not os.path.exists(divide_dir):
            os.makedirs(divide_dir)

        divide_size = int(all_size/divide_num)
        with open(infile,'rb') as f:
            touch_end = False
            for i in range(divide_num):
                seq = f.read(divide_size)
                extract_len = seq.decode('utf-8').rfind('>')
                seq_len = len(seq.decode('utf-8'))
                if seq_len < divide_size:
                    f.seek(-seq_len, os.SEEK_CUR)
                    extract_seq = f.read()
                    out_dir = divide_dir +'/' + str(i)
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir)
                    out_file = out_dir +'/' + str(i) + '.fasta'
                    with open(out_file, 'w') as fw:
                        fw.write(extract_seq.decode('utf-8'))
                    divide_files_list.append(out_dir)
                    break
                else:
                    f.seek(-divide_size, os.SEEK_CUR)

                if i == (divide_num - 1):
                    extract_seq = f.read()
                else:
                    if extract_len == 0:  #如果前面是一整条序列，则往后取完
                        add_read = int(max_size/4)
                        f.seek(divide_size,os.SEEK_CUR)
                        cur_all_offsize = divide_size
                        read_len = divide_size
                        while True:
                            after = f.read(add_read)
                            len_after = len(after.decode("utf-8"))
                            cur_all_offsize += len_after
                            if len_after != add_read:
                                touch_end = True
                                break

                            cur_first_index = after.decode('utf-8').find('>')
                            if cur_first_index != -1:
                                read_len += cur_first_index
                                print('read_len %s'%read_len)
                                break
                            else:
                                read_len += len_after
                        f.seek(-cur_all_offsize,os.SEEK_CUR)
                        if touch_end:
                            extract_seq = f.read()
                        else:
                            extract_seq = f.read(read_len)
                    else:
                        extract_seq = f.read(extract_len)

                out_dir = divide_dir +'/' + str(i)
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                out_file = out_dir +'/' + str(i) + '.fasta'
                with open(out_file, 'w') as fw:
                    fw.write(extract_seq.decode('utf-8'))
                divide_files_list.append(out_dir)

                if touch_end:
                    break
    else:
        print('Not Divide')
        divide_files_list.append(infile)

    return divide_files_list


if __name__ == '__main__':
    infile = sys.argv[1]
    unit_num = int(sys.argv[2])
    type = sys.argv[3]
    s_time = datetime.datetime.now()
    if type=='num':
        divide_by_seq_num(infile,unit_num)
    else:
        divide_fun_by_size(infile,unit_num)
    end_time = datetime.datetime.now()
    with open('time_log','w') as fw:
        fw.write('start: %s\n'%str(s_time))
        fw.write('end: %s\n'%str(end_time))
        fw.write('spend: %s'%(end_time - s_time).seconds)
