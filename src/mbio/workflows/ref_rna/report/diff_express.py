# -*- coding: utf-8 -*-
# __author__ = 'zhangpeng'
# modify by khl

"""有参转录组表达差异分析"""
import os
import re
import pandas as pd
import json
from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.objectid import ObjectId
import shutil
import glob


class DiffExpressWorkflow(Workflow):
    """
    报告中调用组间差异性分析检验时使用
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(DiffExpressWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name="express_file", type="string"),
            dict(name="count", type="outfile", format="rna.express_matrix"),
            dict(name="group_detail", type="string"),
            dict(name="fc", type="float"),
            dict(name="group_id", type="string"),
            dict(name="group_id_id", type="string"),
            dict(name="express_method", type="string"),
            dict(name="update_info", type="string"),
            dict(name="control_file", type="string"),
            dict(name="class_code", type="string"),
            dict(name="diff_express_id", type="string"),
            dict(name="diff_method", type="string", default="DESeq2"),
            dict(name="type", type="string"),
            dict(name="log", type="string"),
            dict(name="express_level", type="string"),  # 对应fpkm/tpm
            dict(name="pvalue_padjust", type="string", default="padjust"),
            dict(name="pvalue", type="float", default=0.05),
            dict(name="class_code_type", type="string", default="express_diff"),  # 传给export_class_code参数
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.diff_exp = self.add_tool("rna.diff_exp")
        self.diff_exp_ref = self.add_tool("rna.diff_exp")
        # self.output_dir = self.diff_exp.output_dir
        self.group_spname = dict()

    def get_samples(self):
        """筛选样本名称"""
        edger_group_path = self.option("group_id")
        self.logger.info(edger_group_path)
        samples = []
        with open(edger_group_path, 'r+') as f1:
            f1.readline()
            for lines in f1:
                line = lines.strip().split("\t")
                samples.append(line[0])
        return samples

    def fpkm(self, samples):
        fpkm_path = self.option("express_file").split(",")[0]
        count_path = self.option("express_file").split(",")[1]
        fpkm = pd.read_table(fpkm_path, sep="\t")
        count = pd.read_table(count_path, sep="\t")
        print "heihei1"
        print fpkm.columns[1:]
        print 'heihei2'
        print samples
        no_samp = []
        sample_total = fpkm.columns[1:]

        for sam in sample_total:
            if sam not in samples:
                no_samp.append(sam)
        self.logger.info("heihei3")
        self.logger.info(no_samp)
        if no_samp:
            new_fpkm = fpkm.drop(no_samp, axis=1)
            new_count = count.drop(no_samp, axis=1)
            print new_fpkm.columns
            print new_count.columns
            self.new_fpkm = self.diff_exp.work_dir + "/fpkm"
            self.new_count = self.diff_exp.work_dir + "/count"
            header = ['']
            header.extend(samples)
            new_fpkm.columns = header
            new_count.columns = header
            new_count.to_csv(self.new_count, sep="\t", index=False)
            new_fpkm.to_csv(self.new_fpkm, sep="\t", index=False)
            return self.new_fpkm, self.new_count
        else:
            return fpkm_path, count_path

    def filter_new(self, class_code, old_count, old_fpkm, new_count, new_fpkm):
        new_id = []
        with open(class_code, 'r+') as f1:
            f1.readline()
            for lines in f1:
                line = lines.strip().split("\t")
                # if self.option("type") == "gene":
                if line[2] == '=':
                    new_id.append(line[0])
                    # if self.option("type") == "transcript":
                    #     if line[2] != 'u':
                    #         new_id.append(line[0])
        with open(old_count, 'r+') as f1, open(old_fpkm, 'r+') as f2, open(new_count, 'w+') as f3, open(new_fpkm,
                                                                                                        'w+') as f4:
            if not new_id:
                raise Exception("生成的class_code信息为空!")
            else:
                f3.write(f1.readline())
                f4.write(f2.readline())
                for lines in f1:
                    line = lines.strip().split("\t")
                    if line[0] in new_id:
                        f3.write(lines)
                for ll in f2:
                    ll2 = ll.strip().split("\t")
                    if ll2[0] in new_id:
                        f4.write(ll)
        self.logger.info("根据class_code提取ref基因完毕！")
        return new_count, new_fpkm

    def run_diff_exp(self, _fpkm, _count, specimen):
        options = {
            "count"       : _count,
            "fpkm"        : _fpkm,
            "control_file": self.option("control_file"),
            "method"      : self.option("diff_method"),
            "fc"          : self.option("fc"),
            "pvalue_padjust" : self.option("pvalue_padjust"),
        }
        if self.option("pvalue_padjust") == "padjust":
            options["diff_fdr_ci"] = self.option("pvalue")
        if self.option("pvalue_padjust") == "pvalue":
            options["diff_ci"] = self.option("pvalue")
        self.option("count").set_path(_count)
        if self.option("group_id_id") != "all":
            options['edger_group'] = self.option("group_id")
        self.diff_exp.set_options(options)
        self.diff_exp.run()

    def run_diff_exp_ref(self, _fpkm, _count, specimen):
        options = {
            "count"       : _count,
            "fpkm"        : _fpkm,
            "control_file": self.option("control_file"),
            # "edger_group":self.option("group_id"),
            "method"      : self.option("diff_method"),
            "fc"          : self.option("fc"),
            "pvalue_padjust" : self.option("pvalue_padjust"),
        }
        if self.option("pvalue_padjust") == "padjust":
            options["diff_fdr_ci"] = self.option("pvalue")
        if self.option("pvalue_padjust") == "pvalue":
            options["diff_ci"] = self.option("pvalue")
        self.option("count").set_path(_count)
        if self.option("group_id_id") != "all":
            options['edger_group'] = self.option("group_id")
        self.diff_exp_ref.set_options(options)
        # self.diff_exp_ref.on("end", self.set_db)
        self.diff_exp_ref.run()

    def set_db(self):
        """
        保存结果表保存到mongo数据库中
        """
        api_diff_exp = self.api.refrna_express
        diff_files = os.listdir(self.diff_exp.output_dir)
        diff_files_ref = os.listdir(self.diff_exp_ref.output_dir)
        self.logger.info("打印diff_files文件信息!")
        self.logger.info(diff_files)
        if self.option("group_id_id") == "all":
            # self.samples = self.diff_exp.option('count').prop['sample']
            self.samples = self.diff_exp.option("count").prop['sample']
            # self.samples_ref = self.diff_exp_ref.option("count").prop['sample']
            self.samples.sort()
            self.logger.info(self.samples)
            self.group_spname['all'] = self.samples
        else:
            self.group_spname = self.diff_exp.option('edger_group').get_group_spname()
            # self.group_spname = self.diff_exp_ref.option('edger_group').get_group_spname()
            edger_group = self.option("group_id")
            self.samples = []

            for keys, values in self.group_spname.items():
                values.sort()
                self.samples.extend(values)
            self.logger.info("specimenname")
            self.logger.info(self.samples)
        compare_column = list()
        workflow_compare_column = {}  # workflow中真正想要运行的control table 分组信息
        workflow_group_name = []
        with open(self.option("control_file"), 'r+') as f1:
            f1.readline()
            _compare_column = []
            for lines in f1:
                _name, _compare_name = lines.strip().split("\t")
                if "|".join([_name, _compare_name]) not in _compare_column:
                    _compare_column.append("|".join([_name, _compare_name]))
        print "打印期望分析的_compare_column："
        print _compare_column

        """导all差异分析表"""
        for f in diff_files:
            if re.search(r'_edgr_stat.xls$', f):
                # 获得所有比较的差异分析文件
                con_exp = f.split('_edgr_stat.xls')[0].split('_vs_')
                c_and_t = [con_exp[1], con_exp[0]]
                if "|".join(c_and_t) in _compare_column:
                    compare_column.append('|'.join(c_and_t))
                else:
                    continue
                # print self.output_dir + '/' + f
                name = con_exp[1]
                compare_name = con_exp[0]

                """添加diff_detail表"""
                api_diff_exp.add_express_diff_detail(name=name, compare_name=compare_name, ref_all='all',
                                                     express_diff_id=self.option("diff_express_id"),
                                                     diff_stat_path=self.diff_exp.output_dir + '/' + f, workflow=False,
                                                     class_code=self.option("class_code"),
                                                     query_type=self.option("type"),
                                                     pvalue_padjust=self.option("pvalue_padjust"))
        """添加summary表"""
        if os.path.exists(self.diff_exp.output_dir + '/merge.xls'):
            api_diff_exp.add_diff_summary_detail(diff_express_id=self.option('diff_express_id'),
                                                 count_path=self.diff_exp.output_dir + '/merge.xls', ref_all='all',
                                                 query_type=self.option('type'),
                                                 class_code=self.option('class_code'), workflow=False)
        else:
            raise Exception("此次ref+new差异分析没有生成summary表！")

        """导ref差异分析表"""
        for f in diff_files_ref:
            if re.search(r'_edgr_stat.xls$', f):
                # 获得所有比较的差异分析文件
                con_exp = f.split('_edgr_stat.xls')[0].split('_vs_')
                c_and_t = [con_exp[1], con_exp[0]]
                if "|".join(c_and_t) in _compare_column:
                    compare_column.append('|'.join(c_and_t))
                else:
                    continue
                # print self.output_dir + '/' + f
                name = con_exp[1]
                compare_name = con_exp[0]

                """添加diff_detail表"""
                api_diff_exp.add_express_diff_detail(name=name, compare_name=compare_name, ref_all='ref',
                                                     express_diff_id=self.option("diff_express_id"),
                                                     diff_stat_path=self.diff_exp_ref.output_dir + "/" + f,
                                                     workflow=False,
                                                     class_code=self.option("class_code"),
                                                     query_type=self.option("type"),
                                                     pvalue_padjust=self.option("pvalue_padjust"))
        """添加summary表"""
        if os.path.exists(self.diff_exp_ref.output_dir + "/merge.xls"):
            api_diff_exp.add_diff_summary_detail(diff_express_id=self.option('diff_express_id'),
                                                 count_path=self.diff_exp_ref.output_dir + "/merge.xls", ref_all='ref',
                                                 query_type=self.option('type'), class_code=self.option('class_code'),
                                                 workflow=False)
        else:
            raise Exception("此次ref差异分析没有生成summary表！")
        """更新主表信息"""
        # if self.ref_all == 'all':
        #     #只用更新一次主表信息即可
        self.update_express_diff(table_id=self.option('diff_express_id'), compare_column=compare_column,
                                 compare_column_specimen=self.group_spname, samples=self.samples)
        self.logger.info("更新主表成功！")
        self.end()

    def update_express_diff(self, table_id, compare_column, compare_column_specimen, samples):
        db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        #client = Config().mongo_client
        #db_name = Config().MONGODB + '_ref_rna'
        self.logger.info(db)
        self.logger.info('haha')

        collection = db['sg_express_diff']
        """方便前端取数据, 生成compare_column_specimen"""

        print "compare_column_specimen"
        print compare_column_specimen

        compare_column = sorted(compare_column)
        tmp_compare_column = []
        for col in compare_column:
            if col not in tmp_compare_column:
                tmp_compare_column.append(col)

        print "compare_column"
        print compare_column

        if self.option("group_id_id") != "all":
            new_compare_column_specimen = {}
            for keys1 in compare_column:
                new_compare = keys1.split("|")
                new_sam = []
                # print new_compare
                for _group in new_compare:
                    if _group in compare_column_specimen.keys():
                        new_sam.extend(compare_column_specimen[_group])
                # print new_sam
                new_compare_column_specimen[keys1] = new_sam
            print new_compare_column_specimen

            group_detail = {}
            group_detal_dict = json.loads(self.option("group_detail"))
            for group, samples in group_detal_dict.items():
                group_detail[group] = samples

            collection.update({'_id': ObjectId(table_id)}, {
                '$set': {'group_detail'           : group_detail, 'compare_column': tmp_compare_column,
                         'specimen'               : self.samples,
                         'compare_column_specimen': new_compare_column_specimen}})
        else:

            collection.update({'_id': ObjectId(table_id)},
                              {'$set': {'compare_column': tmp_compare_column, 'specimen': self.samples}})

    def run(self):
        self.on_rely([self.diff_exp, self.diff_exp_ref], self.set_db)
        if self.option("group_id_id").lower() != 'all':
            specimen = self.get_samples()
        else:
            with open(self.option("control_file"), 'r+') as f1:
                f1.readline()
                specimen = []
                for lines in f1:
                    line = lines.strip().split("\t")
                    for ll in line:
                        if ll not in specimen:
                            specimen.append(ll)
        _fpkm, _count = self.fpkm(specimen)
        self.run_diff_exp(_fpkm=_fpkm, _count=_count, specimen=specimen)
        self.filter_new(class_code=self.option("class_code"), old_count=_count, old_fpkm=_fpkm,
                        new_count=self.work_dir + "/count.tmp", new_fpkm=self.work_dir + "/fpkm.tmp")
        if os.path.exists(self.work_dir + "/fpkm.tmp") and os.path.exists(self.work_dir + "/count.tmp"):
            self.run_diff_exp_ref(_fpkm=self.work_dir + "/fpkm.tmp", _count=self.work_dir + "/count.tmp",
                                  specimen=specimen)
        else:
            raise Exception("没有生成count.tpm和fpkm.tpm表")
        super(DiffExpressWorkflow, self).run()

    def end(self):
        self.output_dir = self.work_dir + "/output"
        output1_dir = self.diff_exp.output_dir
        output2_dir = self.diff_exp_ref.output_dir
        # os.mkdir(self.work_dir + "/diff_exp")  # self.output_dir 为tool的dir，上传文件夹经整合后放于work_dir中
        # os.system("cp -r {} {}".format(output1_dir, self.work_dir + "/diff_exp/refandnew"))
        # os.system("cp -r {} {}".format(output2_dir, self.work_dir + "/diff_exp/ref"))
        self.move2outputdir(output1_dir, "diff_exp")
        self.move2outputdir(output2_dir, "diff_exp_ref")
        self.logger.info("设置表达量结果文件成功")
        # result = self.add_upload_dir(self.work_dir + "/diff_exp")
        # add annotation
        try:
            upload_dir = self.get_workflow_output_dir()
            big_workflow_out = os.path.join(os.path.dirname(os.path.dirname(upload_dir)),"workflow_results")
            refgenedetail_path = os.path.join(big_workflow_out,
                                              "Annotation/GeneAnnotation/AnnoOverview/refgene_anno_detail.xls")
            newgenedetail_path = os.path.join(big_workflow_out,
                                              "Annotation/GeneAnnotation/AnnoOverview/newgene_anno_detail.xls")
            reftransdetail_path = os.path.join(big_workflow_out,
                                               "Annotation/TransAnnotation/AnnoOverview/reftrans_anno_detail.xls")
            newtransdetail_path = os.path.join(big_workflow_out,
                                               "Annotation/TransAnnotation/AnnoOverview/newtrans_anno_detail.xls")
            self.logger.info(upload_dir)
            self.logger.info(big_workflow_out)
            self.logger.info(refgenedetail_path)
            self.logger.info(newgenedetail_path)
            self.logger.info(reftransdetail_path)
            self.logger.info(newtransdetail_path)
            if "DiffExp_G_" in self.output_dir:
                diff_ref = glob.glob(self.output_dir + "/diff_exp_ref/" + "*_vs_*.xls")[0]
                self.paste_annotation(diff_ref,[refgenedetail_path,],diff_ref[:-3] + "annot.xls")
                diff = glob.glob(self.output_dir + "/diff_exp/" + "*_vs_*.xls")[0]
                self.paste_annotation(diff,[refgenedetail_path,newgenedetail_path],diff[:-3] + "annot.xls")
            else:
                diff_ref = glob.glob(self.output_dir + "/diff_exp_ref/" + "*_vs_*.xls")[0]
                self.paste_annotation(diff_ref,[reftransdetail_path,],diff_ref[:-3] + "annot.xls")
                diff = glob.glob(self.output_dir + "/diff_exp/" + "*_vs_*.xls")[0]
                self.paste_annotation(diff,[reftransdetail_path,newtransdetail_path],diff[:-3] + "annot.xls")
        except:
            pass
        result = self.add_upload_dir(self.output_dir)
        result.add_relpath_rules([[".", "", "表达量差异分析结果文件"], ])
        super(DiffExpressWorkflow, self).end()

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            if mode == 'link':
                shutil.copytree(olddir, newdir, symlinks=True)
            elif mode == 'copy':
                shutil.copytree(olddir, newdir)
            else:
                raise Exception('错误的移动文件方式，必须是\'copy\'或者\'link\'')
        else:
            allfiles = os.listdir(olddir)
            oldfiles = [os.path.join(olddir, i) for i in allfiles]
            newfiles = [os.path.join(newdir, i) for i in allfiles]
            self.logger.info(newfiles)
            for newfile in newfiles:
                if os.path.isfile(newfile) and os.path.exists(newfile):
                    os.remove(newfile)
                elif os.path.isdir(newfile) and os.path.exists(newfile):
                    shutil.rmtree(newfile)
            for i in range(len(allfiles)):
                if os.path.isfile(oldfiles[i]):
                    os.link(oldfiles[i], newfiles[i])
                else:
                    os.system('cp -r {} {}'.format(oldfiles[i], newdir))

    @staticmethod
    def paste_annotation(init_table, annot_table_list, out_file):
        init_pd = pd.read_table(init_table, index_col=0, header=0, sep='\t')
        annot_list = list()
        for each in annot_table_list:
            annot_list.append(pd.read_table(each, index_col=0, header=0, sep='\t'))
        annot_pd = pd.concat(annot_list)
        merged_pd = pd.concat([init_pd, annot_pd],  axis=1, join_axes=[init_pd.index])
        merged_pd.to_csv(out_file, header=True, index=True, sep='\t')

    def get_workflow_output_dir(self):
        workflow_output = self._sheet.output
        if workflow_output.startswith('tsanger:'):
            workflow_output = workflow_output.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        else:
            workflow_output = workflow_output.replace('sanger:','/mnt/ilustre/data/')
        return workflow_output