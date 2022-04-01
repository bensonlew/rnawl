#!/usr/bin/python
# -*- coding: utf-8 -*-
# !/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = konghualei
# last modify: 20170119
import os
import re
from collections import defaultdict
import re
import os
import subprocess
from mbio.packages.ref_rna.express.express_distribution import distribution


def prepare(input_file, gtf_file, count_matrix, gene_length):
    # 对featureCounts生成的文件进一步处理，求fpkm.tpm表达量
    #:param input_file, infile, file produced by featureCounts
    #:param gene_length, outfile, 两列，第一列为基因名，第二列为基因外显子长度
    #:param count_matrix, outfile, 生成表达量矩阵
    #:param software, string, 计算表达量的软件
    #:param feature_type, string, 基因
    #:param gtf_file, infile, 默认是无。当gtf_type为merge_cufflinks或merge_stringtie时，才会设置
    #:param gtf_type, string, ref merge_cufflinks merge_stringtie 三种参数
    import re
    import os
    """
    if gtf_type == 'merge_stringtie':
        data=class_code_type(gtf_file=gtf_file, class_code_file=class_code_file, feature_type='Genes',gtf_type=gtf_type)
        #print data
    """
    with open(input_file, 'r+') as f1, open(gene_length, 'w+') as f2, open(count_matrix, 'w+') as f3:
        _head = f1.readline().strip().split("\t")[6:]
        f3.write("\t" + "\t".join(_head) + "\n")
        f2.write("Genes\tLength\n")
        ss = 0
        for line in f1:
            line1 = line.strip().split('\t')
            _value = [str(i) for i in line1[6:]]
            """
            z=0
            for ss in _value:
                if ss == '0':
                    z+=1
            if z==len(_value):
                break
            if gtf_type=='merge_cufflinks':   #只有cufflinks的结果
                if data:
                    if line1[0] in data.keys():
                        if "ref_gene_id" in data[line1[0]].keys():
                            f3.write(data[line1[0]]['ref_gene_id']+'\t'+'\t'.join(_value)+'\n')
                            f2.write(data[line1[0]]['ref_gene_id']+'\t'+line1[5]+'\n')
                            ss+=1
                        else:
                            f3.write(line1[0]+'\t'+'\t'.join(_value)+'\n')
                            f2.write(line1[0]+'\t'+line1[5]+'\n')
                else:
                    raise Exception("没有生成对应的class_code信息，请检查！")
            else:
            """
            f3.write(line1[0] + '\t' + '\t'.join(_value) + '\n')
            f2.write(line1[0] + '\t' + line1[5] + '\n')
    print 'end'
    print ss


def prepare_kallisto(input_file_dir, count_matrix, tpm_matrix, gene_length):
    """
    :param input_file_dir, 多个样本kallisto表达量的文件夹
    :count_matrix, 生成count 的路径和文件名
    :tpm_matrix, 生成tpm 的路径和文件名
    :gene_length, 生成gene_length 的路径和文件名
    """
    if os.path.isdir(input_file_dir):
        pass
    else:
        raise Exception("请输入多个样本kallisto表达量的文件夹！")
    data = {}
    sample_name = []
    _gene_length = {}
    for files in os.listdir(input_file_dir):
        file_path = os.path.join(input_file_dir, files)
        print file_path
        _sample_name = str(files.split(".tsv")[0])
        print _sample_name
        sample_name.append(_sample_name)
        with open(file_path, "r+") as f:
            f.readline()
            for line in f:
                line1 = line.strip().split("\t")
                if str(line1[0]) not in data.keys():
                    data[str(line1[0])] = {}
                if str(_sample_name) not in data[str(line1[0])].keys():
                    data[str(line1[0])][_sample_name] = {}
                    data[str(line1[0])][_sample_name]["count"] = line1[-2]
                    data[str(line1[0])][_sample_name]["tpm"] = line1[-1]
                if str(line1[0]) not in _gene_length.keys():
                    _gene_length[str(line1[0])] = line1[1]
    sample_num = len(sample_name)
    with open(count_matrix, 'w+') as f1, open(tpm_matrix, 'w+') as f2, open(gene_length, "w+") as f3:
        f1.write("Transcript_ID" + "\t" + "\t".join(sample_name) + "\n")
        f2.write("Transcript_ID" + "\t" + "\t".join(sample_name) + "\n")
        f3.write("Transcript_ID" + "\t" + "Length" + "\n")
        for transcript in _gene_length.keys():
            f3.write(transcript + "\t" + _gene_length[transcript] + "\n")
        for keys in data.keys():
            values = data[keys]
            count = []
            tpm = []
            for single_name in sample_name:
                if single_name not in values.keys():
                    count.append(str(0))
                    tpm.append(str(0))
                else:
                    count.append(str(values[single_name]["count"]))
                    tpm.append(str(values[single_name]["tpm"]))
            f1.write(keys + "\t" + "\t".join(count) + "\n")
            f2.write(keys + "\t" + "\t".join(tpm) + "\n")
    print "end"


def all_gene_list(file_path, all_gene_list_path):
    """
    生成全部基因的名称, 没有列名
    :param file_path infile 表达量矩阵
    :param all_gene_list_path outfile 输出全部基因的名称，提供给下游差异分析
    """
    if os.path.exists(file_path):
        pass
    else:
        raise Exception("{}表达量矩阵不存在".format(file_path))
    with open(file_path, "r+") as f1, open(all_gene_list_path, "w+") as f2:
        f1.readline()  # 除去标题
        for line in f1:
            line1 = line.strip().split("\t")[0]
            f2.write(line1 + "\n")


def event_dict_data(file_path):
    """
    根据样本名和基因名，生成字典格式的数据
    :param file_path, 表达量矩阵  count、fpkm、tpm
    :param sample_name, 列表格式——样本名
    """
    _data = {}
    with open(file_path, 'r+') as f:
        _sample_name = f.readline().strip().split('\t')[1:]
        if len(_sample_name) > 1:
            for line in f:
                line1 = line.strip().split('\t')
                _data[line1[0]] = {}
                for sample_id in range(len(_sample_name)):
                    _data[line1[0]][_sample_name[sample_id]] = line1[1:][sample_id]
    print _sample_name
    return _data, _sample_name


def class_code_type(gtf_file, class_code_file, feature_type, gtf_type):
    """
    提取基因或转录本的class_code信息，注释类型选择merged.gtf时原型展示是否为新基因或新转录本
    备注：gtf文件没有第一行标题
    :param gtf_file, infile, gtf文件
    :param class_code_file, outfile_path, class_code 文件生成路径
    :param feature_type, 输出基因或转录本class_code信息
    :param gtf_type, string 选定gtf类型，并输出相应的class_code信息 merge_cufflinks merge_stringtie 两种参数
    """
    import os
    if not os.path.exists(gtf_file):
        raise Exception("gtf文件不存在！")
    # class_code_file=os.path.join(os.path.split(gtf_file)[0], 'class_code')
    with open(gtf_file, 'r+') as f, open(class_code_file, 'w+') as w:
        class_code = {}
        t = 0
        for line in f:

            line1 = line.strip().split('\t')
            _feature1 = line1[8].split(';')
            _feature = [i.replace(' ', "") for i in _feature1]
            ref_gene_id = None
            # print _feature
            # for ss in _feature:
            #    print ss.split('\"')
            for feat in _feature:
                if re.search(r'^gene_id', feat):
                    _gene_id = feat.split('\"')[1]
                    # print _gene_id
                elif re.search(r'transcript_id', feat):
                    _transcript_id = feat.split('\"')[1]
                elif re.search(r'class_code', feat):
                    _class_code = feat.split('\"')[1]
                """
                if gtf_type == "merge_cufflinks":
                    if re.search("oId",feat):
                        #print feat
                        #print 'haha'
                        ref_gene_id=feat.split('\"')[1]
                        #print ref_gene_id
                """
            if feature_type == 'Genes':
                if _gene_id:
                    class_code[_gene_id] = {}
                if _transcript_id:
                    class_code[_gene_id]["transcript_id"] = _transcript_id
                if _class_code:
                    class_code[_gene_id]['class_code'] = _class_code
                if ref_gene_id:
                    class_code[_gene_id]["ref_gene_id"] = ref_gene_id
                    t += 1
                    # if i%100==1:
                    # print i
                    #    pass
        w.write("Genes" + "\t" + "class_code" + "\n")
        for keys in class_code.keys():
            if "ref_gene_id" in class_code[keys].keys():
                w.write(str(keys) + "\t" + str(class_code[keys]["class_code"]) + '\t' + str(
                    class_code[keys]['ref_gene_id']) + "\n")
            else:
                w.write(str(keys) + "\t" + str(class_code[keys]["class_code"]) + "\n")
        print "class_code分类信息已经完成！"
    return class_code

def check(old_gtf_path, new_gtf_path):
    import re, os
    with open(old_gtf_path, 'r+') as f1, open(new_gtf_path, 'w+') as f2:
        for lines in f1:
            line = lines.strip().split('\t')
            if re.search(r'gene_id', line[8]):
                f2.write(lines)
            else:
                break
    return True


# def class_code(infile, outfile):
#     """
#     :param infile: 输入的gtf文件
#     :param outfile: 生成的class_code文件
#     """
#     import time, re, subprocess
#     start = time.time()
#     awk = """ awk -F "\\t" '{print $9}' %s """ % (infile)
#     info = subprocess.check_output(awk, shell=True)
#     output = open(outfile, 'w+')
#     for line in info.strip().split("\n"):
#         gene_id = re.search(r'gene_id\s*"(\S+)";', line).group(1)
#         transcript_id = re.search(r'transcript_id\s*"(\S+)";', line).group(1)
#         class_code = re.search(r'class_code\s*"(\S+)";', line).group(1)
#         if 'gene_name' not in line:
#             print line
#             print "没有gene_name信息！"
#             break
#         else:
#             gene_name = re.search(r'gene_name\s*"(\S+)";', line).group(1)
#             output.write(gene_id + "\t" + transcript_id + "\t" + gene_name + "\t" + class_code + "\n")
#     output.close()
#     end = time.time()
#     time = round(end - start, 4)
#     print "本次运算共消耗 %s s" % (time)
#     infile_length = subprocess.check_output("""wc -l %s """ % (infile), shell=True)
#     print "输入文件长度为 %s" % (str(infile_length.split(" ")[0]))
#     if os.path.exists(outfile):
#         outfile_length = subprocess.check_output("""wc -l %s """ % (outfile), shell=True)
#         print "输出文件长度为 %s" % (str(outfile_length.split(" ")[0]))
#     else:
#         raise Exception("没有生成输出文件，计算失败！")
#     print 'end'


def gff3(infile, outfile, out_path, _type='merged'):
    """
    :param infile:  输入的gtf文件
    :param outfile:
    :param out_path: 生成的class_code文件路径
    :param _type: :替换ID后的merged.gtf文件类型  或参考基因组的gtf文件
    :return: gene2transcript map信息
    """
    import re, subprocess, os, shutil
    new_file = out_path + "/new.%s" % (os.path.basename(infile))
    trans_gene = """grep 'ID=transcript:' %s > %s""" % (infile, new_file)
    print trans_gene
    os.system(trans_gene)
    awk = """awk -F "\\t" '{print $9}' %s""" % (infile)
    info = subprocess.check_output(awk, shell=True)
    data = {}
    tmp = os.path.split(outfile)[0] + "/tmp"
    with open(tmp, 'w+') as f1:
        for lines in info.strip().split("\n"):
            trans_gene = re.search(r'ID=transcript:(\S+);Parent=gene:(\S+);Name=(.*?)', lines)
            if trans_gene:
                transcript_id = trans_gene.group(1)
                gene_id = trans_gene.group(2)
                f1.write(gene_id + "\t" + transcript_id + "\n")
    os.system(""" uniq %s > %s""" % (tmp, outfile))
    os.remove(tmp)
    print 'end'


def gtf(infile, outfile, _type="merged"):
    import time
    # start = time.time()
    import re, subprocess, os
    awk = """awk -F "\\t" '{print $9}' %s""" % (infile)
    info = subprocess.check_output(awk, shell=True)
    print "awk提取gtf文件成功！"
    data = {}
    tmp = os.path.split(outfile)[0] + "/tmp"
    with open(tmp, 'w+') as f1:
        i = 0
        for line in info.strip().split("\n"):
            i += 1
            if i % 10000 == 1:
                print "已提取{}行".format(str(i))
            if "gene_id" in line and "transcript_id" in line:
                gene_id = re.search(r'gene_id\s*"(\S+)";', line).group(1)
                transcript_id = re.search(r'transcript_id\s*"(\S+)";', line).group(1)
                # if transcript_id not in data.keys() and gene_id:
                # data[transcript_id]=gene_id
                f1.write(gene_id + "\t" + transcript_id + "\n")
    # time = round(end - start, 5)
    # print "本次计算时长为{} s".format(time)
    os.system(""" uniq %s > %s""" % (tmp, outfile))
    os.remove(tmp)
    print 'end'


def add_gene_name(old_express, new_express, class_code, type='gene'):
    with open(old_express, 'r+') as f1, open(class_code, 'r+') as f2, open(new_express, 'w+') as f3:
        gene_trans_info = {}
        f2.readline()
        for lines in f2:
            line = lines.strip().split("\t")
            if type == "transcript":
                gene_trans_info[line[0]] = {"gene_name": line[3]}
            if type == "gene":
                gene_trans_info[line[1]] = {"gene_name": line[3]}
        title = f1.readline()
        f3.write(title)
        for ll in f1:
            line = ll.strip().split("\t")
            seq_id = line[0]
            # if type == "transcript":
            #    seq_id = line[0]
            # if type == "gene":
            #    seq_id = line[1]
            if seq_id in gene_trans_info.keys():
                f3.write(seq_id + "," + gene_trans_info[seq_id]['gene_name'] + "\t" + "\t".join(line[1:]) + "\n")
            else:
                f3.write(ll)


def group_express(old_fpkm, new_fpkm, old_count, new_count, sample_group_info, filename, outputfile):
    """
    :params: old_fpkm, fpkm表
    :params: new_fpkm, 新生成的groupfpkm表
    :params: sample_group_info, 样本的分组信息, 字段格式
    :params: rfile, express_distribution.r文件路径
    :params: filename, string格式，文件名
    :params: outputfile, 生成文件路径
    """

    def group_assembly(old_file_path, new_file_path, sample_group_info):
        with open(old_file_path, 'r+') as f1, open(new_file_path, 'w+') as f2:
            sample = f1.readline().strip().split("\t")
            print sample
            group_name = sample_group_info.keys()
            data = {}
            for lines in f1:
                line = lines.strip().split("\t")
                seq_id = line[0]
                data[seq_id] = {}
                for i in range(len(sample)):
                    print i
                    data[seq_id][sample[i]] = float(line[i + 1])
            f2.write("\t" + "\t".join(group_name) + "\n")
            for keys in data.keys():
                grp_data = []
                for grp in group_name:
                    sample_id = sample_group_info[grp]
                    sum_value = 0
                    try:
                        for sam in sample_id:
                            for keys1 in data[keys].keys():
                                m_ = re.search(r'{}'.format(sam), keys1)
                                if m_:
                                    sum_value += data[keys][keys1]
                    except Exception:
                        print data[keys]
                        print sam
                        print sum_value
                        raise Exception("error！")
                    average_value = round(float(sum_value) / len(sample_id), 6)
                    grp_data.append(str(average_value))
                f2.write(keys + "\t" + "\t".join(grp_data) + "\n")

    group_assembly(old_fpkm, new_fpkm, sample_group_info)
    group_assembly(old_count, new_count, sample_group_info)
    # if os.path.exists(new_fpkm):
    #     input_matrix = new_fpkm
    #     distribution(rfile, input_matrix, outputfile, filename)
    # if os.path.exists(new_count):
    #     input_matrix = new_count
    #     distribution(rfile, input_matrix,outputfile, filename+"_count")
    print 'end'


def filter_ref_gene_or_transcript(input_file, class_code, output_path, query_type=None,gene_list=True):
    """
    :param input_file: 输入的表达量文件
    :param output_path: 输出的文件夹路径
    :param sequence_type: gene or transcript
    :param assembly_method: stringtie 或cufflinks拼接
    :param query_type gene or transcript
    :return: 2个文件，ref的表达量文件，及gene list
    """
    class_code_info = {}
    with open(class_code,'r+') as f1:
        f1.readline()
        for lines in f1:
            line=lines.strip().split("\t")
            if query_type == 'gene':
                if line[2] == '=':
                    if line[1] not in  class_code_info.keys():
                        class_code_info[line[1]]=line[2]
            if query_type == 'transcript':
                if line[2] == '=':  # ref transcript
                    if line[0] not in class_code_info.keys():
                        class_code_info[line[0]] = line[2]  # 把transcript 放入class_code_info里面

    filename = os.path.basename(input_file)
    if gene_list:
        list_file = open(output_path + "/{}_list".format(query_type), 'w+')
    with open(input_file, 'r+') as f1, open(output_path + "/{}".format(filename), 'w+') as f2:
        f2.write(f1.readline())
        for lines in f1:
            line = lines.strip().split("\t")
            if line[0] in class_code_info.keys():
                f2.write(lines)
                if gene_list:
                    list_file.write(line[0] + "\n")
            else:
                continue
    if gene_list:
        list_file.close()

def new_gene_location(changed_id_gtf, output_path=None,filename=None):
    if not os.path.exists(changed_id_gtf):
        raise Exception("{}文件不存在".form(changed_id_gtf))
    gene_location_path=os.path.join(output_path,filename)
    gene_info = {}
    with open(changed_id_gtf,'r+') as f1, open(gene_location_path,'w+') as f2:
        for lines in f1:
            if re.search(r'class_code "u"',lines) and  re.search(r'\t+transcript\t+', lines):
                line=lines.strip().split("\t")
                chrom = line[0]
                start = int(line[3])
                end = int(line[4])
                strand = line[6]
                txpt_id_m = re.search(r'transcript_id\s+\"(\S+)\";', lines)
                gene_id_m = re.search(r'gene_id\s+\"(\S+)\";', lines)
                if txpt_id_m and gene_id_m:
                    trans_id = txpt_id_m.group(1)
                    genes_id = gene_id_m.group(1)
                    if genes_id in gene_info.keys():
                        old_start = gene_info[genes_id]['start']
                        old_end = gene_info[genes_id]['end']
                        if int(end -start) > int(old_end - old_start):
                            gene_info[genes_id]['start'] = start
                            gene_info[genes_id]['end'] = end
                        else:
                            continue
                    else:
                        gene_info[genes_id] = {"chr":chrom,"start":start,"end":end,"strand":strand}
        f2.write("gene_id\tchr\tstrand\tstart\tend\n")
        for keys,values in gene_info.items():
            tmp = [values['chr'],values['strand'],str(values['start']),str(values['end'])]
            f2.write(keys+"\t"+"\t".join(tmp)+"\n")
    print 'end'


if __name__ == "__main__":
    # outdata=prepare(input_file="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/feature/sample4_sample1_sample2_sample3",
    #        count_matrix='feature_count.xls', gene_length='feature_gene_length.xls')
    """
    single_sample(input_file="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/feature/sample2_2_sample2_1_sample1_2_sample1_1",\
        count_file='/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/feature/feature_count.xls',\
        fpkm_file='/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/feature/out.fpkm.xls', \
        tpm_file='/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/feature/out.tpm.xls',\
        single_sample_dir='new_sample_dir', _gene=True)

    prepare_kallisto(input_file_dir="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/kallisto/Merge",
            count_matrix="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/kallisto/test_count",
            tpm_matrix="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/kallisto/test_tpm",
            gene_length="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/kallisto/test_gene_length")

    data1=class_code_type(gtf_type="merge_stringtie",gtf_file="/mnt/ilustre/users/sanger-dev/workspace/20170210/Single_assembly_module_hisat_stringtie/Assembly/output/assembly_newtranscripts/data/new_new_transcript.gtf", class_code_file="/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/ref_rna/class_code",feature_type='Genes')
    print data1
    #class_code_type(gtf_file="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/kallisto/lalalal.gtf", class_code_file="/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/tools/kallisto/class_code", feature_type='Transcripts')

    data=class_code_type(gtf_file="/mnt/ilustre/users/sanger-dev/workspace/20170210/Single_assembly_module_hisat_stringtie/Assembly/output/assembly_newtranscripts/_new_transcript.gtf",
                        class_code_file="/mnt/ilustre/users/sanger-dev/workspace/20170210/Single_express_featurecounts_34/Express/Featurecounts/output/class_code",
                        feature_type='Transcripts',gtf_type="merge_stringtie")
    i=0
    for keys,values in data.items():
        i+=1
        if i<=20:
            if 'ref_gene_id' in values.keys():
                print keys,values['class_code'],values['ref_gene_id']
            else:
                print keys,values['class_code']
    prepare(input_file="/mnt/ilustre/users/sanger-dev/workspace/20170215/Single_feature_sample_v1111/Featurecounts/output/sample2_2vssample2_1vssample1_2vssample1_1",
            gtf_file="/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/Eukaryote/tophat2/cufflinks/merge_123456/merged.gtf",
            gtf_type="merge_cufflinks", class_code_file="/mnt/ilustre/users/sanger-dev/workspace/20170215/Single_feature_sample_v1111/Featurecounts/output/new_new_class_code",
            count_matrix="/mnt/ilustre/users/sanger-dev/workspace/20170215/Single_feature_sample_v1111/Featurecounts/output/new_new_count",
            gene_length="/mnt/ilustre/users/sanger-dev/workspace/20170215/Single_feature_sample_v1111/Featurecounts/output/new_new_gene_length",
            software="featurecounts", feature_type="Genes")
    convert("/mnt/ilustre/users/sanger-dev/workspace/20170216/Single_feature_sample_v1111/Featurecounts/output/gtf/ref.gtf", "/mnt/ilustre/users/sanger-dev/workspace/20170216/Single_feature_sample_v1111/Featurecounts/output/gtf/merged.gtf",\
           "/mnt/ilustre/users/sanger-dev/workspace/20170216/Single_feature_sample_v5555/Featurecounts/output/new_new_merged.gtf")
    convert("/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/module/ref.gtf", "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/module/merged_10.gtf",
            "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/module/new_new_new_new.gtf")
    #convert("/mnt/ilustre/users/sanger-dev/workspace/20170216/Single_feature_sample_v6666/Featurecounts/output/ref.gtf","/mnt/ilustre/users/sanger-dev/workspace/20170216/Single_feature_sample_v6666/Featurecounts/output/merged.merged.gtf","/mnt/ilustre/users/sanger-dev/workspace/20170216/Single_feature_sample_v6666/Featurecounts/output/new_new_new_merged.gtf")
       infile = "/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/merged.gtf"
    outfile = "/mnt/ilustre/users/sanger-dev/workspace/20170331/Single_express_small_cufflinks_2/Express/output/new_new_class_code"
    class_code(infile, outfile)
    gtf("/mnt/ilustre/users/sanger-dev/workspace/20170330/Single_rsem_zebra_tools_1/a.gtf",\
        "/mnt/ilustre/users/sanger-dev/workspace/20170330/Single_rsem_zebra_tools_1/gene2transcript")
    gtf("/mnt/ilustre/users/sanger-dev/workspace/20170214/Single_assembly_module_tophat_cufflinks/Assembly/output/Cuffmerge/merged.gtf",\
        "/mnt/ilustre/users/sanger-dev/workspace/20170331/Single_express_small_cufflinks_4/Express/Rsem1/newnew.gene2transcript")
    add_gene_name(old_express = "/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_express_small_cufflinks_4/Express/Rsem1/output/A2_1.genes.results",\
        new_express = "/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_express_small_cufflinks_4/Express/Rsem1/output/new.A2_1.genes.results",\
        class_code = "/mnt/ilustre/users/sanger-dev/workspace/20170401/Single_express_small_cufflinks_4/Express/output/class_code",\
        type = "gene")
    des = "/mnt/ilustre/users/sanger-dev/workspace/20170405/Single_express_small_cufflinks_6/Express/output"
    add_gene_name(old_express = des + "/rsem/A2_1.isoforms.results", new_express = des + '/rsem/newA2_1.isoforms.results',\
                class_code = des+'/class_code', type ='transcript')
    gtf("/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_assembly_module_tophat_stringtie_zebra/Assembly/assembly_newtranscripts/merged.gtf",\
    "/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_rsem_stringtie_zebra_3/Express/Rsem1/newgene2transcript"
    )

    fpkm_path = "/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_rsem_stringtie_zebra_7/Express/MergeRsem/output/transcripts.TMM.fpkm.matrix"
    new_fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_rsem_stringtie_zebra_7/Express/MergeRsem/output/new_transcript.fpkm"
    sample_group_info = {"a1":["ERR1621569_sickle_l","ERR1621480_sickle_l"],"a2":["ERR1621658_sickle_l","ERR1621391_sickle_l"]}
    rfile = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/module/express_distribution.r"
    # input_matrix = "/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_rsem_stringtie_zebra_7/Express/MergeRsem/output/new_transcript.fpkm"
    filename = "transcriptGroup"
    outputfile = "/mnt/ilustre/users/sanger-dev/workspace/20170410/Single_rsem_stringtie_zebra_7/Express/MergeRsem/output"
    group_express(old_fpkm = fpkm_path, new_fpkm = new_fpkm, sample_group_info = sample_group_info, rfile = rfile, filename=filename, outputfile=outputfile)

    old_express = "/mnt/ilustre/users/sanger-dev/workspace/20170425/Single_merge_rsem_fpkm_9/MergeRsem/output/genes.counts.matrix"
    new_express = "/mnt/ilustre/users/sanger-dev/workspace/20170425/Single_merge_rsem_fpkm_9/MergeRsem/output/new.genes.counts.matrix"
    class_code = "/mnt/ilustre/users/sanger-dev/workspace/20170425/Single_merge_rsem_fpkm_9/MergeRsem/class_code"
    add_gene_name(old_express,new_express,class_code,"gene")
    print 'end'

    fpkm_path = "/mnt/ilustre/users/sanger-dev/workspace/20170504/Single_feature_sample_v2/Featurecounts/output/fpkm_tpm.fpkm.xls"
    new_fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170504/Single_feature_sample_v2/Featurecounts/output/new_fpkm.fpkm.xls"
    sample_group_info = {"B":['sample2_2', 'sample2_1'], "A":['sample1_2', 'sample1_1']}
    rfile = "/mnt/ilustre/users/sanger-dev/workspace/20170504/Single_feature_group_v3/Featurecounts/express_distribution.r"
    filename="genes"
    outputfile = "/mnt/ilustre/users/sanger-dev/workspace/20170504/Single_feature_group_v3/Featurecounts/group"
    group_express(old_fpkm = fpkm_path, new_fpkm = new_fpkm, sample_group_info = sample_group_info, rfile = rfile, filename=filename, outputfile=outputfile)

    fpkm_path = "/mnt/ilustre/users/sanger-dev/workspace/20170510/Single_feature_Truestringtie_sample_1/Express/Featurecounts/output/fpkm_tpm.fpkm.xls"
    new_fpkm = "/mnt/ilustre/users/sanger-dev/workspace/20170510/Single_feature_Truestringtie_sample_1/Express/Featurecounts/output/group_fpkm"
    # compare_column_specimen["CD|HGD"] = [ "CL1", "CL2", "CL5","HGL1", "HGL3","HGL4"]
    # compare_column_specimen["CD|HFD"] = ["CL1", "CL2", "CL5", "HFL4","HFL6","HFL3"]
    # compare_column_specimen["HGD|HFD"] = ["HGL1", "HGL3","HGL4","HFL4","HFL6","HFL3"]
    sample_group_info = {"CD":["CL1", "CL2", "CL5"],"HGD":["HGL1", "HGL3","HGL4"],"HFD":["HFL4","HFL6","HFL3"]}
    filenames="genes"
    outputfile = "data"
    count_path = "/mnt/ilustre/users/sanger-dev/workspace/20170510/Single_feature_Truestringtie_sample_1/Express/Featurecounts/output/count.xls"
    new_count = "/mnt/ilustre/users/sanger-dev/workspace/20170510/Single_feature_Truestringtie_sample_1/Express/Featurecounts/output/group_count"
    group_express(fpkm_path,new_fpkm,count_path,new_count, sample_group_info,"rdata",filenames,outputfile)
    """
    changed_id_gtf = "/mnt/ilustre/users/sanger-dev/workspace/20170616/Single_assembly_module_tophat_stringtie_true_file_1/RefrnaAssemble/NewTranscripts/output/change_id_merged.gtf"
    output_path = "/mnt/ilustre/users/sanger-dev/sg-users/konghualei/ref_rna/class_code"
    filename = "new_gene_location"
    new_gene_location(changed_id_gtf, output_path, filename)