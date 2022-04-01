# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'

from api_base import ApiBase
import datetime
import os
import re
from collections import defaultdict
import math


class Assembly(ApiBase):
    def __init__(self, bind_object):
        """
        wgs中所有数据组装中的导表
        __author__ = HONGDONG
        __last_modify__ = 20180426
        :param bind_object:
        """
        super(Assembly, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_qc_file(self, file_path, assem_id, possname=None):
        """
        :param file_path:
        :param assem_id:
        :return:
        """
        data_list = []
        assem_id = self.check_objectid(assem_id)
        if not os.path.isdir(file_path):
            raise Exception("{}不是一个路径！".format(file_path))
        file_reads = os.listdir(file_path)
        for m in file_reads:
            n = re.match(r"(.*)\.stat$", m)
            if n:
                if self.get_qc_info(os.path.join(file_path, m), assem_id, possname):
                    data_list.append(self.get_qc_info(os.path.join(file_path, m), assem_id, possname))
        if data_list:
            self.col_insert_data("sg_assembly_reads", data_list)
        else:
            self.bind_object.logger.info("请检查{}文件夹中有没有stat文件！".format(file_path))

    def get_qc_info(self, file_path, assem_id, possname=None):
        self.check_exists(file_path)
        insert_data = {}
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for line in data:
                # print line
                temp = line.strip().split('\t')
                if temp[0].split('.')[1] in ['unmap', 'unmapping']:
                    start = "--"
                    end = "--"
                else:
                    try:
                        if len(temp[0].split('.')[1].split(':')[1].split('-')) == 2:
                            start = int(temp[0].split('.')[1].split(':')[1].split('-')[0])
                            end = int(temp[0].split('.')[1].split(':')[1].split('-')[1])
                        else:
                            start = int(temp[0].split('.')[1].split(':')[1].split('-')[0])
                            end = '--'
                    except:
                        start = "1"
                        end = '--'
                insert_data = {
                    "assembly_id": assem_id,
                    "specimen_id": temp[0].split('.')[0],
                    "chr": temp[0].split('.')[1].split(':')[0],
                    "start": start,
                    "end": end,
                    "reads_num": int(temp[1])/2,
                    "base_num": int(temp[2])/2,
                    "q30_rate": float(temp[4]),
                    "gc_rate": float(temp[3])
                }
                if possname:
                    if temp[0].split('.')[1] in ['unmap', 'unmapping']:
                        mapper_reads = '--'
                        mapper_bases = '--'
                        unmapping_reads = int(temp[1]) / 2
                        unmapping_bases = int(temp[2]) / 2
                    else:
                        mapper_reads = int(temp[1]) / 2
                        mapper_bases = int(temp[2]) / 2
                        unmapping_reads = '--'
                        unmapping_bases = '--'
                    posname = temp[0].split('.')[1]
                    if posname in possname.keys():
                        region_id = possname[posname]
                    elif posname == 'unmapping':
                        region_id = "unmapping"
                    else:
                        region_id = posname
                    insert_data.update({
                        "mapper_reads": mapper_reads,
                        "mapper_bases": mapper_bases,
                        "unmapping_reads": unmapping_reads,
                        "unmapping_bases": unmapping_bases,
                        'region_id': region_id
                    })
        return insert_data if insert_data else ""

    def add_pc_file(self, file_path, assem_id, possname=None):
        """
        :param file_path:
        :param assem_id:
        :param possname:  每个区域对应的自主命名
        :return:
        文件夹中stat.xls文件导表，尚未添加assem_id参数。  AH03_unmapping.20190311_160720431_4201.denovo.stat.xls
        """
        data_list = []
        if not os.path.isdir(file_path):
            raise Exception("{}不是一个路径！".format(file_path))
        file_reads = os.listdir(file_path)
        for i in file_reads:
            n = re.match(r"(.*)\.denovo\.stat\.xls$", i)
            # print n
            if n:
                sp = i.split('.')[0]
                """
                sp是样本的名称，从excel的名字中提取
                """
                data_list.append(self.get_pc_info(os.path.join(file_path, i), assem_id, sp, possname))

        if data_list:
            self.col_insert_data("sg_assembly_stat", data_list)
        else:
            self.bind_object.logger.info("请检查{}文件夹中有没有stat文件！".format(file_path))

    def get_pc_info(self, file_path_1, assem_id, file_name_1, possname=None):
        """
        YC_bulk_chr3:1-20000
        :param file_path_1:
        :param assem_id:
        :param file_name_1:
        :param possname:
        :return:
        """
        insert_data = {}
        region_id = ''
        assem_id = self.check_objectid(assem_id)
        sample_name = file_name_1
        if possname:
            tt = file_name_1.split(':')[0].split("_")
            sample_name = "_".join(tt[:-1])
            posname = file_name_1.split('_')[-1]
            if posname in possname.keys():
                region_id = possname[posname]
            elif posname == 'unmapping':
                region_id = "unmapping"
            else:
                region_id = posname
        with open(file_path_1, 'r') as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split('\t')
                insert_data = {
                    "assembly_id": ("--" if assem_id == '' else assem_id),
                    "specimen_id": ("--" if sample_name == '' else sample_name),
                    "largest_scaf": ("--" if temp[0] == '' else temp[0]),
                    "largest_len": ("--" if temp[1] == "" else int(temp[1])),
                    "large_scaf": ("--" if temp[2] == '' else temp[2]),
                    "bases_scaf": ("--" if temp[3] == '' else int(temp[2])),
                    "n50_len": ("--" if temp[5] == '' else int(temp[5])),
                    "n50_scaf": ("--" if temp[4] == '' else int(temp[4])),
                    "n90_scaf": ("--" if temp[6] == '' else int(temp[6])),
                    "n90_len": ("--" if temp[7] == '' else int(temp[7])),
                    "gc_content": ("--" if temp[8][:-1] == '' else float(temp[8][:-1])),
                    "n_rate": ("--" if temp[9][:-1] == '' else float(temp[9][:-1])),
                }
                if possname:
                    insert_data.update({'region_id': region_id})
        return insert_data if insert_data else ""

    def get_anno_db(self):
        """
        直接从sg_anno_db数据库中查找nr，go，kegg等库文件，该条记录只有一个，永远用最新的记录
        :return:
        """
        result = self.db["sg_anno_db"].find({}).sort([("created_ts", -1)])
        if result.count() == 0:
            raise Exception("sg_anno_db中没有找到对应信息，请检查库文件！")
        else:
            try:
                nr_db = result[0]['nr_db']
                uniport_db = result[0]['uniport_db']
                go_db = result[0]['go_db']
                kegg_db = result[0]['kegg_db']
                nog_db = result[0]['nog_db']
                pfam_db = result[0]['pfam_db']
            except Exception as e:
                raise Exception("sg_anno_db中nr_db，uniport_db，go_db，kegg_db，nog_db查找不到，请检查！{}".format(e))
            else:
                self.bind_object.logger.info("从sg_anno_db中查找相应信息成功!")
        return nr_db, uniport_db, go_db, kegg_db, nog_db, pfam_db

    def get_ref_db(self, species_version_id):
        species_version_id = self.check_objectid(species_version_id)
        result = self.col_find_one("sg_species_version", {"_id": species_version_id})
        if result:
            ref_db = os.path.join(os.path.dirname(result['ref']), "makedbblast/ref")
        else:
            raise Exception("sg_species_version中没有找到{}对应信息！".format(species_version_id))
        return ref_db

    def add_anno_file(self, file_path, assem_id, possname=None):
        """
        :param file_path:
        :param assem_id:
        :return:
        """
        # data_list = []
        if not os.path.isdir(file_path):
            raise Exception("{}不是一个路径！".format(file_path))
        file_reads = os.listdir(file_path)
        for i in file_reads:
            n = re.match(r"(.*)\.final\.summary$", i)
            # print n
            if n:
                sp = i.split('.')[0]
                """
                sp是样本的名称，从excel的名字中提取
                """
                self.get_anno_info(os.path.join(file_path, i), assem_id, sp, possname)

    def get_anno_info(self, file_path_1, assem_id, file_name_1, possname=None):
        """
        chr start end结果文件中命令：chr1;chr1:1;chr1:1-2000
        前端传过来的格式为：chr1:-; chr1:1-; chr1:1-2000(或者chr1:-2000)
        :param file_path_1:
        :param assem_id:
        :param file_name_1:
        :param possname:
        :return:
        """
        insert_data = {}
        data_list1 = []
        region_id = ''
        sample_name = file_name_1
        assem_id = self.check_objectid(assem_id)
        sample_region_id = ''
        if possname:
            tt = file_name_1.split(':')[0].split("_")
            sample_name = "_".join(tt[:-1])
            posname = file_name_1.split('_')[-1]
            if posname in possname.keys():  # chr1:1-2000这里直接判断了完成的chr start end
                region_id = possname[posname]
            elif posname == 'unmapping':
                region_id = "unmapping"
            else:
                if re.match('(.*):(.*)', posname):
                    if re.match('(.*):([0-9]*)$', posname):  # chr1:1  -- chr1:1-
                        region = posname + '-'
                        if region in possname.keys():
                            region_id = possname[region]
                        else:
                            self.set_error("没有找到具体位置与region之间的对应关系")
                    elif re.match("(.*):(.*)-(.*)", posname):  # chr1:-2000 -- chr1:1-2000; chr1:1-2000--chr1:1-2000
                        region = ':'.join([posname.split(':')[0], ''.join(['-', posname.split(':')[1].split('-')[1]])])
                        if region in possname.keys():
                            region_id = possname[region]
                        else:
                            self.set_error("没有找到具体位置与region之间的对应关系")
                else:
                    region = posname + ':-'
                    if region in possname.keys():
                        region_id = possname[region]
                    else:
                        self.set_error("没有找到具体位置与region之间的对应关系")
            sample_region_id = sample_name + "_" + region_id
        with open(file_path_1, 'r') as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split('\t')
                # if temp[-3] == '--':
                #     continue
                insert_data = {
                    "scaffold_id": temp[0],
                    "assembly_id": assem_id,
                    "sample_region_id": sample_region_id,
                    "specimen_id": sample_name,
                    "chr": temp[-3],
                    "start": temp[-2],
                    "end": temp[-1],
                    "nr_id": temp[1],
                    "nr_desc": temp[2],
                    "ko_id": temp[5],
                    "ko_desc": temp[6],
                    "go_id": temp[7],
                    "go_desc": temp[8],
                    "eggnog_id": temp[9],
                    "eggnog_desc": temp[10],
                    "uniprot_id": temp[3],
                    "uniprot_desc": temp[4]
                }
                if possname:
                    # noinspection PyBroadException
                    try:
                        pfam_acc = temp[11]
                        pfam_ann = temp[12]
                        interpro_acc = temp[13]
                        interpro_ann = temp[14]
                    except:
                        pfam_acc = "--"
                        pfam_ann = "--"
                        interpro_acc = "--"
                        interpro_ann = "--"
                    insert_data.update({
                        'region_id': region_id,
                        "pfam_acc": pfam_acc,
                        "pfam_ann": pfam_ann,
                        "interpro_acc": interpro_acc,
                        "interpro_ann": interpro_ann
                    })
                data_list1.append(insert_data)
        if data_list1:
            self.col_insert_data("sg_assembly_anno", data_list1)
        else:
            self.bind_object.logger.info("文件{}为空！".format(file_path_1))
        # print "sg_specimen导表成功"
        return insert_data if insert_data else ""


if __name__ == "__main__":
    # project_sn = 'wgs_test'
    # task_id = 'wgs_test'
    # member_id = 'wgs_test'
    # data = Assembly(None)
    # file_path = "/mnt/ilustre/users/sanger-dev/workspace/20180509/Assembly_wgs_test_0509133840_3539_9226/" \
    #             "output/assembly_stat"
    # data.add_pc_file(file_path,'5aea749da4e1af445ef0b521')
    data = Assembly(None)
    # file_path = '/mnt/ilustre/users/sanger-dev/workspace/20180508/Single_summary_stat_20180508/SummarySet/output'
    # data.add_anno_file(file_path, "5aea749da4e1af445ef0b521")
    # file_path = '/mnt/ilustre/users/sanger-dev/workspace/20180517/Assembly_wgs_test_0517155225570453_1242/Assembly/SoapDenovo/output/JY102.denovo.scafSeq'
    # file_path_1 = "/mnt/ilustre/users/sanger-dev/workspace/20180517/Assembly_wgs_test_0517155225570453_1242/Assembly/AssemblyStat/output/JY102.20180517_155427333_4479.denovo.stat.xls"
    # data.add_bt_file(file_path)
    # data.add_qc_file("/mnt/ilustre/users/sanger-dev/workspace/20190311/Assembly_wgs_v2_03111600"
    #                  "47688909_9833/output/qc_stat", "5cab0ca417b2bf53a302dc2a",
    #                  {u'chr1:10000-20000': u'region01', u'chr2:10000-40000': u'region02'})
    data.add_anno_file("/mnt/ilustre/users/sanger-dev/workspace/20190430/Assembly_sanger_85433_0430122"
                       "036663892_4602/output/final_anno", "5cc7cd1417b2bf33b2700b8c", {u'chr3:-20000': u'region01'})