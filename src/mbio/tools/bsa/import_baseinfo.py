# # !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "HONGDONG"

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class ImportBaseinfoAgent(Agent):
    """
    该tool用于bsa导入前面的所有的基础数据，这些数据是不需要计算的是产品线直接进行导入的
    version v1.0
    author: hongdongxuan
    last_modify: 2018.6.11
    最后一次更新是为了对接wgs项目结果
    """
    def __init__(self, parent):
        super(ImportBaseinfoAgent, self).__init__(parent)
        options = [
            {"name": "bsa_path", "type": "infile", "format": "bsa.bsa_path"},   # 产品线上传的数据的基础路径，该路径下还有其他的子路径，要进行下判断
            {"name": "task_id", "type": "string"},   # 任务id
            {"name": "project_sn", "type": "string"},  # 项目id
            {"name": "member_id", "type": "string"}  # 会员id
        ]
        self.add_option(options)
        self.step.add_steps("ImportBaseinfo")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.ImportBaseinfo.start()
        self.step.update()

    def stepfinish(self):
        self.step.ImportBaseinfo.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("bsa_path"):
            raise OptionError("必须输入bsa_path！", code="31500401")
        if not self.option('task_id'):
            raise OptionError('必须有任务id', code="31500402")
        if not self.option('project_sn'):
            raise OptionError('必须有项目id', code="31500403" )
        if not self.option('member_id'):
            raise OptionError('必须有会员id', code="31500404")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 3
        self._memory = '10G'

    def end(self):
        super(ImportBaseinfoAgent, self).end()


class ImportBaseinfoTool(Tool):
    """
    按照模块进行导表操作
    """
    def __init__(self, config):
        super(ImportBaseinfoTool, self).__init__(config)
        self._version = '1.0.1'
        self.base_info = self.api.api("bsa.base_info")
        self.api_base = self.api.api("bsa.api_base")

    def import_qc(self):
        """
        导入qc信息+qc质控图片
        :return:
        """
        stat_path = os.path.join(self.option("bsa_path").prop['path'], '01.fastq_qc/cleandata_qc/qc.xls')
        if os.path.exists(stat_path):
            self.logger.info("将使用wgs v1的结果去进行导表")
            self.base_info.add_sg_specimen_qc(self.option('task_id'), self.option('project_sn'), stat_path)
            self.logger.info("导入样本的信息以及样本质控信息成功！")
            atgc_path = os.path.join(self.option("bsa_path").prop['path'], '01.fastq_qc/cleandata_qc/atgc')
            for files in os.listdir(atgc_path):
                result = self.get_sample_qc_id(files, "sg_specimen_qc")
                if result:
                    origin_id = result['_id']
                    self.base_info.add_qc_atgc_curve(self.option('task_id'), origin_id, os.path.join(atgc_path, files))
                    self.logger.info("导入{}碱基含量分布图成功".format(files))
                else:
                    self.set_error("样本%s在sg_specimen_qc没有查找到！", variables =(files), code="31500401")
            qual_path = os.path.join(self.option("bsa_path").prop['path'], '01.fastq_qc/cleandata_qc/qual')
            for files in os.listdir(qual_path):
                result = self.get_sample_qc_id(files, "sg_specimen_qc")
                if result:
                    origin_id = result['_id']
                    self.base_info.add_qc_qual_curve(self.option('task_id'), origin_id, os.path.join(qual_path, files))
                    self.logger.info("导入{}碱基错误率分布图成功".format(files))
                else:
                    self.set_error("样本%s在sg_specimen_qc没有查找到！", variables =(files), code="31500402")
        else:
            json_dir = os.path.join(self.option("bsa_path").prop['path'], '01.fastq_qc/qc_stat')
            self.logger.info('开始导入sg_specimen表格')
            for m in os.listdir(json_dir):
                specimen_id = m.split(".json")[0].split("-")[0]
                self.base_info.add_sg_specimen(self.option('project_sn'), self.option('task_id'), specimen_id,
                                               specimen_id, '')
            self.logger.info('导入sg_specimen表格成功！')
            self.base_info.add_fastp_json_stat(self.option('project_sn'), self.option('task_id'),  json_dir)
            self.logger.info("导入质控导表成功！")

    def import_map(self):
        """
        导入比对信息, 这里的复杂的地方在于，每张图的数据 是要叠加的，这样X轴要取最大值，然后几个样本的X轴的值要补齐，
        并且对应value值要用NaN替代
        :return:
        """
        total_map_path = os.path.join(self.option("bsa_path").prop['path'],
                                      '03.map_stat/result.stat/Total.mapped.detail.xls')
        main_id = self.base_info.add_sg_mapping(self.option('task_id'), self.option('project_sn'),
                                                self.option("member_id"), "origin_no_params")
        self.logger.info("导入mapping主表成功！")
        self.base_info.add_sg_mapping_detail(total_map_path, main_id)
        self.logger.info("导入mapping细节表成功！")

        insert_path = os.path.join(self.option("bsa_path").prop['path'], '03.map_stat/insert')
        insert_x, insert_categories = self.get_file_xmax(insert_path, "insert")
        insert_curve_id = self.api_base.sg_curve(self.option('task_id'), main_id, "", insert_categories, 1,
                                                 "mapping_fragment")
        for files in os.listdir(insert_path):
            self.base_info.add_insert_curve(insert_curve_id, os.path.join(insert_path, files), insert_x-1)
            self.logger.info("导入{}插入片段分布图成功！".format(files))

        depth_path = os.path.join(self.option("bsa_path").prop['path'], '03.map_stat/depth')
        depth_x, dep_categories = self.get_file_xmax(depth_path, "dep")
        depth_curve_id = self.api_base.sg_curve(self.option('task_id'), main_id, "", dep_categories, 1,
                                                "mapping_deep", "", "")
        for files in os.listdir(depth_path):
            self.base_info.add_depth_curve(depth_curve_id, os.path.join(depth_path, files))
            self.logger.info("导入{}测序深度分布图成功！".format(files))

        coverage_path = os.path.join(self.option("bsa_path").prop['path'], '03.map_stat/coverage')
        ref_change = os.path.join(self.option("bsa_path").prop['path'], '02.reference/ref.chrlist')
        for files in os.listdir(coverage_path):
            name = files.split('.')[0]
            area_id = self.base_info.add_sg_area(self.option('project_sn'), self.option('task_id'), name, main_id)
            self.base_info.add_sg_area_detail(area_id, ref_change, os.path.join(coverage_path, files))

    def get_file_len(self, file_path):
        """
        获取文件的行数
        :param file_path:
        :return:
        """
        return len(open(file_path, 'rU').readlines())

    def get_file_xmax(self, dir_path, types="dep"):
        """
        获取mapping insert或者depth文件夹中所有文件的最大行数，x轴最大值
        :param dir_path:
        :param types:
        :return:
        """
        xmax = 0
        categories = []
        for files in os.listdir(dir_path):
            temp_len = self.get_file_len(os.path.join(dir_path, files))
            if temp_len > xmax:
                xmax = temp_len
        if types == "dep":
            for m in range(1, 1001):
                categories.append(str(m))
            categories.append(">1000")
        else:
            categories = range(0, xmax)
        return xmax, categories

    def get_sample_qc_id(self, files, collection):
        """
        模糊匹配查找到样本对应的质控主表的id
        :return:
        """
        sample_name = files.strip().split(".")[0].split('-')[0]
        # result = self.api.api('bsa.api_base').col_find_one(collection,  {"task_id": self.option('task_id'),
        #                                                                  "specimen_id":
        # {'$regex': "^" + sample_name}})
        result = self.api.api('bsa.api_base').col_find_one(collection, {"task_id": self.option('task_id'),
                                                                        "specimen_id": sample_name})
        # if not result:
        #     self.set_error("样本{}在{}中没有找到对应信息！".format(sample_name, collection))
        return result

    def import_variant(self):
        """
        导入突变信息
        :return:
        """
        variant_stat_path = os.path.join(self.option("bsa_path").prop['path'], '04.snp_indel/variant_stat')
        anno_stat_path = os.path.join(self.option("bsa_path").prop['path'], '04.snp_indel/anno_stat')
        call_id = self.base_info.add_sg_snp_call(self.option('task_id'), self.option('project_sn'),
                                                 self.option("member_id"), "origin_no_params")
        self.base_info.add_sg_snp_call_stat(variant_stat_path + '/snp.stat', call_id)
        snp_anno_id = self.base_info.add_sg_snp_anno(self.option('task_id'), self.option('project_sn'),
                                                     self.option("member_id"), "origin_no_params")
        self.logger.info("1")
        self.base_info.add_sg_snp_indel_anno_stat(anno_stat_path + '/snp.stat', snp_anno_id, "snp")
        indel_id = self.base_info.add_sg_indel_call(self.option('task_id'), self.option('project_sn'),
                                                    self.option("member_id"), "origin_no_params")
        self.logger.info("2")
        self.base_info.add_sg_indel_call_stat(variant_stat_path + '/indel.stat', indel_id)
        indel_anno_id = self.base_info.add_sg_indel_anno(self.option('task_id'), self.option('project_sn'),
                                                         self.option("member_id"), "origin_no_params")
        self.base_info.add_sg_snp_indel_anno_stat(anno_stat_path + '/indel.stat', indel_anno_id, "indel")
        self.logger.info("3")
        self.base_info.add_indel_length_bar(self.option('task_id'), indel_id, variant_stat_path + '/indel.len')
        self.logger.info("导入indel长度分布")

    def import_annovar(self):
        """
        导入注释信息
        :return:
        """
        self.logger.info("开始导入注释信息")
        anno_file = os.path.join(self.option("bsa_path").prop['path'], '05.annovar/anno_count/pop.summary')
        version_file = os.path.join(self.option("bsa_path").prop['path'], '02.reference/info.log')
        anno_id = self.base_info.add_sg_anno(self.option('task_id'), self.option('project_sn'),
                                             self.option("member_id"), "origin_no_params")
        self.base_info.add_sg_anno_detail(anno_file, version_file, anno_id, self.option("task_id"))
        self.logger.info("导入注释信息成功")

    def run(self):
        super(ImportBaseinfoTool, self).run()
        self.import_qc()
        self.import_map()
        self.import_variant()
        self.import_annovar()
        self.end()
