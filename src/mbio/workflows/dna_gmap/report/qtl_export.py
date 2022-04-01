# -*- coding: utf-8 -*-
# __author__ = 'qing_mei'
# modified 20180709
# workflow

import re
import os
import tarfile
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class QtlExportWorkflow(Workflow):
    """
    标记筛选的接口
    sg_qtl_export_id
    接收到的path已经在controller里转成绝对路径
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(QtlExportWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "sg_qtl_export_id", "type": "string"},     # main_id
            {"name": "task_id", "type": "string"},
            {"name": "sg_feature_file_id", "type": "string"},
            {"name": "sg_lg_id", "type": "string"},
            {"name": "popt", "type": "string"},
            {"name": "type", "type": "string"},     # 转化的类型，str
            {"name": "trit_path", "type": 'string'},
            {"name": "sexAver_loc_path", "type": 'infile', 'format': 'bsa.vcf'},  # F1
            {"name": "sexAver_map_path", "type": 'infile', 'format': 'bsa.vcf'},  # F1
            {"name": "total_loc_path", "type": 'infile', 'format': 'bsa.vcf'},    # 非F1
            {"name": "total_map_path", "type": 'infile', 'format': 'bsa.vcf'},    # 非F1
            {"name": "total_csv_path", "type": 'infile', 'format': 'bsa.vcf'},    # 非F1
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.qtl_export = self.add_tool('dna_gmap.qtl_export')
        self.target_dir_ = ''

    def check_options(self):
        if not self.option('sg_lg_id'):
            raise OptionError("请设置sg_lg_id", code="14800501")
        if not self.option('sg_feature_file_id'):
            raise OptionError("请设置sg_feature_file_id", code="14800502")
        if not self.option("sg_qtl_export_id"):
            raise OptionError("请设置sg_qtl_export_id", code="14800503")
        if not self.option("task_id"):
            raise OptionError("请设置task_id", code="14800504")
        if not self.option("trit_path"):
            raise OptionError("性状文件不存在", code="14800505")
        if self.option('popt').lower() not in ['cp', 'f1', 'f2', 'ril', 'dh', 'bc']:
            raise OptionError("群体参数错误：%s", variables=(self.option('popt')), code="14800506"
)
        # if re.match(r'ri\d+', self.option("popt"), re.I):
        #     pass
        # else:
        #     if self.option("popt").upper() not in ["BC", "DH", "CP", "F2", "F1"]:
        #         raise OptionError("群体类型{}不属于BC, DH, CP, F2, RIL, Ri\d+".format(self.option("popt")))
        if self.option('popt').lower() in ['cp', 'f1']:
            if not self.option('sexAver_loc_path') and not self.option('sexAver_map_path'):
                raise OptionError("F1群体类型sg_lg必须传递参数sexAver_loc_path，sexAver_map_path", code="14800507")
        else:
            if not self.option('total_loc_path') and not self.option('total_map_path') and not self.option('total_csv_path'):
                raise OptionError("非F1群体类型sg_lg必须传递参数total_loc_path,total_map_path,total_csv_path", code="14800508")
        return True

    def qtl_export_run(self):
        """
        """
        opt = {
            "type": self.option('type'),
            "popt": self.option('popt'),
            "trit_path": self.option('trit_path')
        }
        if self.option('popt').lower() in ['cp', 'f1']:
            opt.update({
                "sexAver_loc_path": self.option('sexAver_loc_path').prop['path'],
                "sexAver_map_path": self.option('sexAver_map_path').prop['path']
            })
        else:
            opt.update({
                "total_csv_path": self.option('total_csv_path').prop['path'],
                "total_loc_path": self.option('total_loc_path').prop['path'],
                "total_map_path": self.option('total_map_path').prop['path']
            })
        self.qtl_export.set_options(opt)
        self.qtl_export.on('end', self.set_output, 'qtl_export')
        self.qtl_export.run()

    def set_output(self, event):
        """
        链接结果
        return：tar压缩结果
        """
        obj = event['bind_object']
        if event['data'] == 'qtl_export':
            self.linkdir(obj.output_dir, self.output_dir)
        if os.path.exists(os.path.join(self.output_dir + 'format_result.tar')):
            os.remove(os.path.join(self.output_dir + 'format_result.tar'))
            self.logger.info('已删除上次结果的压缩文件')
        allfiles = os.listdir(self.output_dir)
        self.logger.info('########################开始tar')
        self.logger.info(" ".join(allfiles))
        os.system("cd {0} && tar -czvf format_result.tar.gz {1}".format(self.output_dir, " ".join(allfiles)))
        self.set_db()
        pass

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        """
        更新主表的export_path
        """
        self.get_target_dir()
        self.logger.info("设置re路径OK 开始qtl_export的导表！")
        api = self.api.api("dna_gmap.qtl_export")
        task_id = self.option('task_id')
        main_id = self.option('sg_qtl_export_id')
        sg_lg_id = self.option('sg_lg_id')
        sg_feature_file_id = self.option('sg_feature_file_id')
        # api.updata_sg_qtl_export(main_id, self.output_dir + "/format_result.tar")
        api.updata_sg_qtl_export(main_id, self.target_dir_ + "/format_result.tar.gz")
        self.logger.info("更新sg_qtl_export的export_path成功！")
        api.add_sg_qtl_export_detail(main_id, task_id, sg_lg_id, sg_feature_file_id)
        self.logger.info("sg_qtl_export_detail导表成功！")
        self.end()

    def run(self):
        self.qtl_export_run()
        super(QtlExportWorkflow, self).run()

    def end(self):
        """
        这里后面要重新定义下文件名字
        :return:
        """
        # self.set_output_file()
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(QtlExportWorkflow, self).end()

    def get_target_dir(self):
        """
        获取远程磁盘的路径
        :return:
        """
        # self.target_dir_ = self._sheet.output.strip().split('://')[1]
        self.target_dir_ = self._sheet.output.rstrip('/')

