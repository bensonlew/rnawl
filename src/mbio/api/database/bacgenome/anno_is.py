# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20201031
from biocluster.api.database.base import Base, report_check
import os,re
import datetime
import types
import json
from Bio import SeqIO
from bson.son import SON
from collections import defaultdict
from bson.objectid import ObjectId


class AnnoIs(Base):
    def __init__(self, bind_object):
        super(AnnoIs, self).__init__(bind_object)
        self._project_type = "bacgenome"
        self.seq_convert = {
            "Chromosome" : 'Chr',
            "Chromosome1" : "Chr1",
            "Chromosome2" : "Chr2",
            "Chromosome3" : "Chr3",
            "Plasmid": 'p',
            "PlasmidA" :'pA',
            "PlasmidB" : 'pB',
            "PlasmidC" :'pC',
            "PlasmidD" : 'pD',
            "PlasmidE" :'pE',
            "PlasmidF" : 'pF',
            "PlasmidG" :'pG',
            "PlasmidH" : 'pH',
            "PlasmidI" : 'pI',
            "PlasmidJ" :'pJ',
            "PlasmidK" : 'pK',
            "PlasmidL" :'pL',
            "PlasmidM" : 'pM',
            "PlasmidN" : 'pN',
            "PlasmidO" :'pO',
            "PlasmidP" : 'pP',
        }

    @report_check
    def add_anno_is(self, task_id=None, project_sn=None,params=None, name=None, specimen_id=None):
        if task_id is None:
            task_id = self.bind_object.sheet.id
        project_sn = project_sn if project_sn else self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '插入序列分析',
            'created_ts': created_ts,
            'name': name if name else 'IS_Origin',
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            "software" : "ISEScan_v1.7.2.1",
            # 'specimen_id': specimen_id,
        }
        collection = self.db['is']
        main_id = collection.insert_one(insert_data).inserted_id
        collection.update({'_id': main_id},{'$set':{'main_id':main_id}})
        return main_id

    @report_check
    def add_anno_is_stat(self, detail_file, main_id=None,sample=None):
        """
        导入样本统计表
        :param detail_file:
        :param main_id:
        :param sample:
        :return:
        """
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        data_list = []
        collection_detail = self.db['is_stat']
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            if sample:
                for line in lines[1:]:
                    if sample in line:
                        line = re.split(r"\s+", line)
                        data = [("is_id", main_id),
                                ("sample", sample),
                                ("is_num", int(line[1]))]
            else:
                spline = lines[-1].strip().split(" ")
                data = [("is_id", main_id),
                        ("sample", spline[0].split(".fna")[0]),
                        ("is_num", int(spline[1]))]
            data_son = SON(data)
            data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入is_stat%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入is_stat%s信息出错" )
        else:
            self.bind_object.logger.info("导入is_stat%s信息成功!" % detail_file)

    @report_check
    def add_anno_is_detail(self, detail_file, seq_file, main_id=None, sample=None, gff_file=None, transposon_file=None,blast_file=None, analysis=None):
        """
        导入is的统计表，以is维度进行导表
        :param detail_file:
        :param seq_file:
        :param main_id:
        :param sample:
        :param gff_file:
        :param transposon_file:
        :param blast_file:
        :return:
        """
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        gff_dict = defaultdict(list)
        if gff_file:## 根据过滤后的is表导入最终的结果
            with open(gff_file, 'r') as f:
                f.readline()
                for line in f:
                    spline = line.strip().split("\t")
                    location = spline[1]
                    if location not in gff_dict:
                        gff_dict[location] = [spline]
                    else:
                        aa_list = gff_dict[location]
                        aa = spline
                        if aa not in aa_list:
                            aa_list.append(aa)
                        gff_dict[location] = aa_list
        transposon_dict = defaultdict(list)
        if transposon_file: ##转座酶的信息结果，判断is包含哪些基因
            with open(transposon_file, 'r') as f:
                f.readline()
                for line in f:
                    spline = line.strip().split("\t")
                    location = spline[1]
                    if location not in transposon_dict:
                        transposon_dict[location] = [spline]
                    else:
                        aa_list = transposon_dict[location]
                        aa = spline
                        if aa not in aa_list:
                            aa_list.append(aa)
                        transposon_dict[location] = aa_list
        blast_dict = defaultdict(list)
        if blast_file: ##is的blast的结果，用于挑选is_name，挑选原则是从5条比对结果中挑选与is预测得到的family和group相同的一条结果
            with open(blast_file, 'r') as f2:
                f2.readline()
                for line in f2:
                    spline = line.strip().split("\t")
                    gene = spline[0]
                    new_gene = "IS" + str(gene.split("IS")[1]).zfill(3)
                    loca = spline[1]
                    if new_gene not in blast_dict:
                        blast_dict[new_gene] = [loca]
                    else:
                        aa_list = blast_dict[new_gene]
                        if loca not in aa_list:
                            aa_list.append(loca)
                        blast_dict[new_gene] = aa_list

        data_list = []
        collection_detail = self.db['is_detail']
        seq_dict = {}
        for seq_record in SeqIO.parse(seq_file, 'fasta'):
            id = seq_record.id
            seq = str(seq_record.seq)
            if id not in seq_dict:
                seq_dict[id] = seq
        with open(detail_file, 'r') as f:
            lines = f.readlines()
            is_num2 = 0
            location_list = []
            for line in lines[2:]:
                line = re.split(r"\s+", line)
                cluster = line[2]
                if line[0].strip() in gff_dict:
                    aa_list = gff_dict[line[0]]
                    start = int(line[3])
                    end = int(line[4])
                    for aa in aa_list:## 从众多的IS结果中挑出IS预测的结果
                        new_aa_list = aa
                        origin_start = int(new_aa_list[2])
                        origin_end = int(new_aa_list[3])
                        result = self.compare_pos(origin_start,origin_end,start,end)
                        pick_is_id = new_aa_list[0].split("_")[1] + str(new_aa_list[0].split("_")[2]).zfill(3)
                        if result: ## 从is的起始终止位置挑出来相关的gene，为了计算同一个is包含多少基因
                            if line[0] not in location_list:
                                is_num2 = 1
                                location_list.append(line[0])
                            else:
                                is_num2 += 1
                            is_id = "IS"+ str(is_num2).zfill(3)
                            all_gene_list = []
                            if line[0].strip() in transposon_dict:
                                aa_gene_list = transposon_dict[line[0].strip()]
                                for gene_list in aa_gene_list:
                                    start = int(gene_list[2])
                                    end = int(gene_list[3])
                                    result = self.compare_pos(origin_start,origin_end,start,end)
                                    if result:
                                        all_gene_list.append(gene_list)
                            cds_num = len(all_gene_list)
                            is_name = '-'
                            if "Accessory_Gene" in cluster:
                                cluster = cluster.rstrip("_Accessory_Gene")
                                if len(cluster.split("_")) >= 3:
                                    if is_name in ['unknown', '-', '']:
                                        is_name = '-'
                                    else:
                                        is_name = line[1]
                                    group = cluster.split("_")[1]
                                    family = cluster.split("_")[2]
                                else:
                                    is_name = "-"
                                    family = cluster.split("_")[1]
                                    group = cluster.split("_")[0]
                            elif "Passenger_Gene" in cluster:
                                cluster = cluster.rstrip("_Passenger_Gene")
                                if len(cluster.split("_")) >= 3:
                                    if is_name in ['unknown', '-', '']:
                                        is_name = '-'
                                    else:
                                        is_name = line[1]
                                    group = cluster.split("_")[1]
                                    family = cluster.split("_")[2]
                                else:
                                    family = cluster.split("_")[1]
                                    group = cluster.split("_")[0]
                            else:
                                if len(cluster.split("_")) >= 3:
                                    if is_name in ['unknown', '-', '']:
                                        is_name = '-'
                                    else:
                                        is_name = line[1]
                                    group = cluster.split("_")[1]
                                    family = cluster.split("_")[2]
                                else:
                                    is_name = "-" ## fix by qingchen.zhang @20201224 应产品线需求改成 -
                                    family = cluster.split("_")[1]
                                    group = cluster.split("_")[0]
                            new_pick_is_id = line[0] + "_"+ pick_is_id
                            if pick_is_id in blast_dict:
                                is_name_list = blast_dict[pick_is_id]
                                is_name = self.pick_is_name(family, group, is_name_list)
                            else:
                                is_name = is_name

                            data = [("is_id", main_id),
                                    ("sample", sample),
                                    ("is", is_id),
                                    ("is_name", is_name),
                                    ("family", family),
                                    ("group", group),
                                    ("start", int(line[3])),
                                    ("end", int(line[4])),
                                    ("len", int(line[4]) - int(line[3])),
                                    ("score", int(line[11])),
                                    ("evalue", float(line[19])),
                                    ("cds_num", cds_num),
                                    ]
                            if new_pick_is_id in seq_dict:
                                data.append(("seq", seq_dict[new_pick_is_id]))
                            else:
                                data.append(("seq", "-"))
                            if line[21] in ['p']:
                                data.append(("type", "partial"))
                            else:
                                data.append(("type", "complete"))
                            if analysis:
                                data.append(("location", line[0]))
                            else:
                                data.append(("location", line[0]))
                            data_son = SON(data)
                            data_list.append(data_son)
        try:
            collection_detail.insert_many(data_list)
            s_collection = self.db['is']
            s_collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
        except Exception as e:
            self.bind_object.logger.error("导入is_detail%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入is_detail%s信息出错" )
        else:
            self.bind_object.logger.info("导入is_detail%s信息成功!" % detail_file)
        ## 下面代码更新circos表的字段，供前端判断使用
        task_id = "_".join(self.bind_object.sheet.id.split("_")[0:2])
        detail_table = self.db['circos_table_detail']
        try:
            type_main = self.db['is']
            type_main_id = type_main.find_one({"_id": main_id})['_id']
            type_main_detail = self.db['is_detail']
            main_table = self.db['circos_table']
            circos_table_id = main_table.find_one({"task_id": task_id})["_id"]
            if not isinstance(circos_table_id, ObjectId):
                circos_table_id = ObjectId(circos_table_id)
            for one in type_main_detail.find({"is_id": type_main_id}):
                if one['location'].startswith('Scaffold'):
                    detail_table.update_one(
                        {"circos_table_id": circos_table_id, "specimen_id": one['sample']},
                        {'$set': {'IS': 1}}, upsert=True)
                else:
                    detail_table.update_one(
                        {"circos_table_id": circos_table_id, "location": one['location'], "specimen_id": one['sample']},
                        {'$set': {'IS': 1}}, upsert=True)
        except Exception as e:
            self.bind_object.logger.info('没有插入序列信息')
        else:
            self.bind_object.logger.info('导入插入序列信息')

    @report_check
    def add_anno_is_summary(self, detail_file, seq_file, sample, main_id=None, transposon_file=None, gene_file=None,analysis=None):
        """
        导入is的各个元素表的情况
        :param detail_file:is.xls
        :param seq_file:_sequence.fna
        :param sample:sample
        :param main_id:main_id
        :param transposon_file:
        :param gene_file:
        :return:
        """
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        gff_dict = defaultdict(list)
        # if gff_file: 后面改掉了，原因是用不到此表了
        #     with open(gff_file, 'r') as f:
        #         for line in f:
        #             spline = line.strip().split("\t")
        #             location = spline[1]
        #             start = spline[2]
        #             end = spline[3]
        #             if location not in gff_dict:
        #                 gff_dict[location] = [location + "\t" + start + "\t" + end]
        #             else:
        #                 aa_list = gff_dict[location]
        #                 aa = location + "\t" + start + "\t" + end
        #                 if aa not in aa_list:
        #                     aa_list.append(aa)
        #                 gff_dict[location] = aa_list
        transposon_dict = defaultdict(list)
        if transposon_file: ##转座酶的信息结果，将所有的基因标记上转座酶的结果，如果没有则标记protein结果
            with open(transposon_file, 'r') as f:
                f.readline()
                for line in f:
                    spline = line.strip().split("\t")
                    location = spline[1]
                    if location not in transposon_dict:
                        transposon_dict[location] = [spline]
                    else:
                        aa_list = transposon_dict[location]
                        aa = spline
                        if aa not in aa_list:
                            aa_list.append(aa)
                        transposon_dict[location] = aa_list
        seq_dict = {}
        seq_descrip = {}
        for seq_record in SeqIO.parse(seq_file, 'fasta'): ## 主要是提取左端和右端的序列
            id = seq_record.id
            seq = str(seq_record.seq)
            description = str(seq_record.description)
            if id not in seq_dict:
                seq_dict[id] = seq
            if id not in seq_descrip:
                seq_descrip[id] = description
        gene_dict = {}
        if gene_file: ## 提取基因的序列
            for seq_record in SeqIO.parse(gene_file, 'fasta'):## 用此来提取基因的序列
                id = seq_record.id
                gene_seq = str(seq_record.seq)
                if id not in gene_dict:
                    gene_dict[id] = gene_seq
        data_list = []
        collection_detail = self.db['is_component']
        location_list = []
        with open(detail_file, 'r') as f:## 循环遍历得到的is.xls的结果
            lines = f.readlines()
            is_num = 0
            for line in lines[1:]:
                line = line.strip().split("\t")
                is_num_name = line[0].split("_")[1] + line[0].split("_")[2].zfill(3)
                location = line[1]
                origin_start = int(line[2])
                origin_end = int(line[3])
                is_gene_list = []
                if location in transposon_dict:
                    all_gene_list = transposon_dict[location]
                    for gene_list in all_gene_list:
                        start = int(gene_list[2])
                        end = int(gene_list[3])
                        result = self.compare_pos(origin_start,origin_end,start,end)
                        if result:
                            is_gene_list.append(gene_list)
                # location = self.seq_convert[location]
                if location not in location_list:
                    is_num = 1
                    location_list.append(line[1])
                else:
                    is_num += 1
                is_id = "IS" + str(is_num).zfill(3)
                is_id_left = location + "_" + is_num_name + "_left"
                is_id_right = location + "_" + is_num_name + "_right"
                component_left = ""
                component_right = ""
                if is_id_left in seq_dict:
                    component_left = "IRL"
                if is_id_right in seq_dict:
                    component_right = "IRR"
                location = location
                if len(is_gene_list) != 0: ##用于判断is中有多少基因，以及酶的类型，序列
                    for gene_list in is_gene_list:
                        gene_id = gene_list[0]
                        gene_location = gene_list[1]
                        gene_start = gene_list[2]
                        gene_end = gene_list[3]
                        gene_starnd = gene_list[4]
                        anno_type = gene_list[5]
                        if gene_location in self.seq_convert.keys():
                            if gene_location in ["Chromosome", "Chromosome1", "Chromosome2", "Chromosome3"]:
                                origin_gene_id = gene_id
                            else:
                                origin_gene_id = gene_location + "_" + gene_id
                        else:
                            origin_gene_id = gene_id
                        data = [("is_id", main_id),
                        ("sample", sample),
                        ("location", location),
                        ("is", is_id),
                        ("strand", gene_starnd),
                        ("component", gene_id)]
                        if anno_type in ["transposase"]:
                            anno = "transposase"
                            anno_type = "Transposase"
                        else:
                            anno = "protein"
                            anno_type = "tetm"
                        data.append(("anno", anno))
                        data.append(("start", int(gene_start)))
                        data.append(("end", int(gene_end)))
                        data.append(("length", abs(int(gene_start) - int(gene_end))))
                        data.append(("anno_type", anno_type))
                        if origin_gene_id in gene_dict:
                            data.append(("seq", gene_dict[origin_gene_id]))
                        else:
                            data.append(("seq", "-"))
                        data_son = SON(data)
                        data_list.append(data_son)
                for compo in [component_left, component_right]: ### 用于判断重复序列的类型和插件类型
                    data = [("is_id", main_id),
                        ("sample", sample),
                        ("location", location),
                        ("is", is_id),
                        ("component", compo)]
                    if compo in [component_left]:
                        anno = "repeat"
                        seq_id = is_id_left
                        if line[7] != "-":
                            left_start = line[7].split("..")[0]
                            lefs_end = line[7].split("..")[1]
                            if int(left_start) < int(lefs_end):
                                data.append(("strand", "+"))
                            else:
                                data.append(("strand", "-"))
                            data.append(("anno", anno))
                            data.append(("start", int(left_start)))
                            data.append(("end", int(lefs_end)))
                            data.append(("length", abs(int(left_start) - int(lefs_end))))
                            data.append(("anno_type", "DR"))
                            if seq_id in seq_dict:
                                data.append(("seq", seq_dict[seq_id]))
                            else:
                                data.append(("seq", "-"))
                            data_son = SON(data)
                            data_list.append(data_son)
                    elif compo in [component_right]:
                        anno = "repeat"
                        seq_id = is_id_right
                        if line[8] != "-":
                            right_start = line[8].split("..")[0]
                            right_end = line[8].split("..")[1]
                            if int(right_start) < int(right_end):
                                data.append(("strand", "+"))
                            else:
                                data.append(("strand", "-"))
                            data.append(("anno", anno))
                            data.append(("start", int(right_start)))
                            data.append(("end", int(right_end)))
                            data.append(("length", abs(int(right_start) - int(right_end))))
                            data.append(("anno_type", "Rev"))
                            if seq_id in seq_dict:
                                data.append(("seq", seq_dict[seq_id]))
                            else:
                                data.append(("seq", "-"))
                            data_son = SON(data)
                            data_list.append(data_son)

        try:
            if len(data_list) >1:
                collection_detail.insert_many(data_list)
        except Exception as e:
            self.bind_object.logger.error("导入is_component%s信息出错:%s" % (detail_file, e))
            self.bind_object.set_error("导入is_component%s信息出错" )
        else:
            self.bind_object.logger.info("导入is_component%s信息成功!" % detail_file)

        try:
            collection = self.db["is"]
            result = collection.find_one({"_id": main_id})
            if result:
                if "specimen_id" in result:
                    specimen_id_list = str(result['specimen_id']).split(",")
                    if sample not in specimen_id_list:
                        specimen_id_list.append(sample)
                    new_specimen_id = ",".join(specimen_id_list)
                else:
                    new_specimen_id = sample
            else:
                new_specimen_id = sample
            collection.update({"_id": main_id}, {"$set": {"specimen_id": new_specimen_id}})
        except Exception as e:
            self.bind_object.logger.error("更新主表is主表信息出错:%s" % (e))
            self.bind_object.set_error("更新主表is主表信息出错:%s" )
        else:
            self.bind_object.logger.info("更新主表is主表信息成功!")


    def compare_pos(self, origin_start, origin_end, start,end):
        """
        比较大小,用于判断is是否需要导表
        :return:
        """
        if origin_start <= start and origin_end >= end:
            return True
        else:
            return False

    def pick_is_name(self, family, group, is_name_list):
        """
        功能是从blast的结果中挑出与family和group相同的一条is_name
        :param family:
        :param group:
        :param is_name_list:
        :return:
        """
        is_name = "-"
        for all_name in is_name_list:
            all_list = all_name.strip().split("_")
            if family in [all_list[1]] and group in [all_list[2]]:
                is_name = all_list[0]
                break
            else:
                continue
        return is_name

    def add_main_id(self, main_id=None):
        try:
            s_collection = self.db['is']
            s_collection.update({'_id': main_id}, {'$set': {'main_id': main_id}})
        except Exception as e:
            self.bind_object.logger.error("导入is main_id出错:%s" % (e))
            self.bind_object.set_error("导入is main_id出错" )



