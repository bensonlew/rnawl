#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    作者：kefei.huang
    时间：20171206
    说明：用来计算个体识别。原理是检查两个样品之间有没有一样的基因型。
    1. 首先检查mongo数据库，确认是否已经有结果，如果有结果，有哪些样品是已经计算过的。
    2. 构建一个计算样品的tab文件列表，绝对路径为单位。目前存放tab文件的路径是固定的。
    3. 计算相似度的结果，与早先计算好的结果组合（如果有），按照从大到小排序，取前10个（取更多也没问题）
    4. 更新mongo数据表

    注意：目前来说，只支持在父本数据。而且必须是已经算完，在数据库中可以查到的数据。暂时不支持给出自定义的tab文件。如果需要自定义查询，暂时线下操作

"""
from __future__ import print_function
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
import re
import codecs
import datetime


class SampleSimilarityAgent(Agent):
    def __init__(self, parent):
        super(SampleSimilarityAgent, self).__init__(parent)
        options = [
            {"name": "sample_id", "type": "string"}
        ]
        self.add_option(options)
        self.library_path = self.config.SOFTWARE_DIR + "/" + "database/human/pt_ref/tab_all"
        self.logger.info("options finished")

    def check_options(self):
        if not self.option("sample_id"):
            raise OptionError("请输入样品名称")
        # ##检查两个数据，确认这个数据是否存在
        abs_sample = "{}/{}.tab".format(self.library_path, self.option("sample_id"))
        if not os.path.exists(abs_sample):
            raise OptionError("{}文件不存在,请确认".format(abs_sample))
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "1G"

    def end(self):
        super(SampleSimilarityAgent, self).end()


class SampleSimilarityTool(Tool):
    def __init__(self, config):
        super(SampleSimilarityTool, self).__init__(config)
        self.sample_id = self.option("sample_id")
        self.library_path = str(self.config.SOFTWARE_DIR + "/" + "database/human/pt_ref/tab_all")
        self.sample_tab = "{}/{}.tab".format(self.library_path, self.option("sample_id"))
        self.sample_sort_file = "{}.outfile.xls".format(self.option("sample_id"))
        self.logger.info("init finished")

    def run(self):
        super(SampleSimilarityTool, self).run()
        group_return = self.get_new_sample()    
        father_list = group_return[0]
        old_similarity = group_return[1]
        self.Identification(father_list, old_similarity)
        self.dump_similarity(old_similarity)
        self.end()

    def make_dict(self, input):
        dic = {}
        with codecs.open(input, "r", "utf-8") as FA:
            for line in FA:
                if len(line) < 8:    # 过滤掉列数小于8的位点
                    continue
                else:
                    line = line.strip()
                    tmp = line.split("\t")
                    sample_name, chrom, pos, ref, alt, dp, ref_dp, alt_dp = tmp
                    if alt == ".":
                        alt = ref
                    if int(ref_dp) + int(alt_dp) < 30:    # 过滤掉深度小于30的位点
                        continue
                    alt1 = alt[0]       # 取alt的第一个碱基作为突变碱基
                    YeChun = "%s/%s" % (ref, ref)
                    TuChun = "%s/%s" % (alt1, alt1)
                    ZaHe = "%s/%s" % (ref, alt1)
                    rf = float(ref_dp)/(float(ref_dp) + float(alt_dp))
                    if len(ref) > 1:
                        continue         # 去除ref碱基数超过1的位点
                    if rf >= 0.9:
                        tmp.append(YeChun)
                    elif rf <= 0.1:
                        tmp.append(TuChun)
                    else:
                        tmp.append(ZaHe)
                    sample_name, chrom, pos, ref, alt, dp, ref_dp, alt_dp, genetype = tmp
                    content = sample_name+'\t'+chrom+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+dp+'\t'+ref_dp+'\t'+alt_dp+'\t'+genetype+'\n'
                    dic[chrom+"_"+pos] = genetype           # 以染色体及位置作为键,基因型作为值,生成字典
        return dic

    def get_new_sample(self):
        father_list = []
        mother_list = []
        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        collection = db['sg_sample_similarity']
        dir_result = os.walk(self.library_path)
        for root, dirs, files in dir_result:
            for each in files:
                print(each)
                if re.search(r"-F", each) is not None:
                    father_list.append(each.replace('.tab', ''))
                elif re.search(r"-M", each) is not None:
                    mother_list.append(each.replace('.tab', ''))
        # ###获取新增样品
        merge_list = father_list
        merge_list.extend(mother_list)
        
        sample_id_mongo = collection.find_one({'sample_id': self.sample_id})
        if sample_id_mongo is not None:
            new_sample = list(set(merge_list).difference(set(sample_id_mongo["dec"])))  # #只要参考文件夹内，与monggo
            # 数据库中不同的样品，都算是新增样品
            old_similarity = sample_id_mongo["result"]
        else:
            new_sample = merge_list
            old_similarity = "empty"

        abs_sample = []
        for each in new_sample:
            if re.search(r'{}'.format(self.sample_id), each) is None:
                abs_name = "{}/{}.tab".format(self.library_path, each)
                abs_sample.append(abs_name)
        return [abs_sample, old_similarity]

    def Identification(self, sample_list, old_similarity):
        run_num = 1
        f1 = codecs.open("{}.tmpfile.xls".format(self.sample_id), 'w', "utf-8")
        current_dict = self.make_dict(self.sample_tab)           # 生成当前样本的字典
        # ##写入中间文件###
        if old_similarity != "empty":
            for each in old_similarity:
                content = "\t".join(each)
                print(content, file=f1)
        for sn in sample_list:
            sn = sn.strip()
            other_dict = self.make_dict(sn)                              # 生成所需要对比的字典
            same_gene = 0
            all_gene = 0
            for cur_k, cur_v in current_dict.iteritems():
                if other_dict.has_key(cur_k):
                    all_gene += 1                              # 当对比的两个样本的位置一致时,累加,作为分母
                    if cur_v == other_dict[cur_k]:
                        same_gene += 1                          # 当位置一致的时候,基因型一致时,累加,作为分子
                    else:
                        continue
                else:
                    continue
            if all_gene > 0:
                percent = '{:.4f}'.format(100*float(same_gene)/float(all_gene))+"%"
            else:
                percent = 0
            content = [str(sn), str(all_gene), str(same_gene), str(percent)]
            content = "\t".join(content)
            print(content, file=f1)
            # print("{}/{}".format(run_num,overall_num))
            run_num += 1                                       # 将计算出的每一行结果写入临时文件中,后续对匹配度进行排序工作
        f1.close()
        cmd = "less {}.tmpfile.xls|sort -nrk4 > {}.outfile.xls".format(self.sample_id, self.sample_id)
        os.system(cmd)

    def dump_similarity(self, old_similarity):
        count = 0
        similarity_data = []
        similarity_sample = []
        with codecs.open(self.sample_sort_file, "r", "utf-8") as IN:
            for line in IN:
                line = line.strip()
                linetemp = line.split("\t")
                linetemp[0] = linetemp[0].split("/")[-1].replace('.tab', '')
                similarity_sample.append(linetemp[0])
                count += 1
                if count <= 10:
                    similarity_data.append(linetemp)

        insert_table = {
            "sample_id": self.sample_id,
            "upload_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "result": similarity_data,
            "dec": similarity_sample
        }
        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        collection = db['sg_sample_similarity']
        collection.update({"sample_id": self.sample_id}, {'$set': insert_table}, upsert=True)
