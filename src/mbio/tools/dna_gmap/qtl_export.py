# -*- coding: utf-8 -*-
# __author__ = 'qing_mei'
# modified 20180705
# tool

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class QtlExportAgent(Agent):
    """
    qtl格式化输出
    f1群体：rqtl和
    传进去sg_task的ril群体
    """

    def __init__(self, parent):
        super(QtlExportAgent, self).__init__(parent)
        options = [
            # {"name": "sg_lg_id", "type": "string"},  # 放在params，所以不需要
            # {"name": "sg_feature_id", "type": 'string'},  # 放在params，所以不需要
            {"name": "type", "type": 'string'},     # 页面上选择的待转化类型
            {"name": "popt", "type": "string"},
            {"name": "trit_path", "type": 'string'},
            {"name": "sexAver_loc_path", "type": 'string'},  # F1
            {"name": "sexAver_map_path", "type": 'string'},  # F1
            {"name": "total_loc_path", "type": 'string'},    # 非F1
            {"name": "total_map_path", "type": 'string'},    # 非F1
            {"name": "total_csv_path", "type": 'string'},    # 非F1
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("type"):
            raise OptionError("页面传递的转化类型不能为空", code="34801101")
        if self.option('popt').lower() not in ['cp', 'f1', 'f2', 'ril', 'dh', 'bc']:
            raise OptionError("群体参数错误：%s", variables=(self.option('popt')), code="34801102")
        if not self.option("trit_path"):
            raise OptionError("性状文件不存在", code="34801103")
        # if re.match(r"ri\d+", str(self.option("popt")), re.I).group():
        #     pass
        # else:
        #     if self.option("popt").upper() not in ["BC", "DH", "CP", "F2", "F1"]:
        #         raise OptionError("群体类型{}不属于BC, DH, CP, F2, RIL, Ri\d+".format(self.option("popt")))
        if self.option('popt').lower() in ['cp', 'f1']:
            if not self.option('sexAver_loc_path') and not self.option('sexAver_map_path'):
                raise OptionError("F1群体类型sg_lg必须传递参数sexAver_loc_path，sexAver_map_path", code="34801104")
        else:
            if not self.option('total_loc_path') and not self.option('total_map_path') and not self.option('total_csv_path'):
                raise OptionError("非F1群体类型sg_lg必须传递参数total_loc_path,total_map_path,total_csv_path", code="34801105")

    def set_resource(self):
        self._cpu = 2
        self._memory = "3G"

    def end(self):
        super(QtlExportAgent, self).end()


class QtlExportTool(Tool):
    def __init__(self, config):
        super(QtlExportTool, self).__init__(config)
        self.rscript_path = 'program/R-3.3.3/bin/Rscript'
        self.convertNOCP_path = self.config.PACKAGE_DIR + "/dna_gmap/qtl-convertNOCP.R"

    def get_popt(self, popt):
        """
        转化传进来的群体类型
        """
        if popt.lower() == 'f2':
            return 'f2'
        elif popt.lower() == 'bc' or popt.lower() == 'dh':
            return 'bc'
        else:
            return 'riself'

    def run_nocp(self, format, csv_path, trit_path):
        """
        非F1群体，运行nocp的R程序
        先转化输入的群体类型
        生成结果是一个路径：内含两个文件
        """
        popt = self.get_popt(self.option('popt'))
        cmd = "{} {}".format(self.rscript_path, self.convertNOCP_path)
        cmd += " --mark {} --trt {}".format(self.option('total_csv_path'), trit_path)
        cmd += " --pop {} --out {}".format(popt, self.output_dir + "/" + format + '_nocp')
        cmd += ' --format {}'.format(format)
        self.run_cmd(cmd, "nocp_qtlcart_format", format + '_nocp')

    def run_cmd(self, cmd, cmd_name, name=None):
        """
        执行cmd
        """
        self.logger.info(cmd)
        command = self.add_command(cmd_name, cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{}运行完成".format(cmd_name))
        else:
            self.set_error("%s运行失败", variables=(cmd_name), code="34801101")
            self.set_error("%s运行失败".format(cmd_name), code="34801104")
        # self.option(cmd_name).set_path(self.output_dir + "/" + name)

    def link_path(self, filepath1, filepath2):
        if os.path.exists(filepath2):
            os.remove(filepath2)
        os.link(filepath1, filepath2)
        return True

    def run(self):
        """
        f1群体、f2群体不用转化mapqtl、rqtl，qtlcart
        f2群体只转化qtlcart
        待确定：等宣红东确定hl程序内目前非f1群体是否都有total.csv total.loc total.map文件，
        """
        super(QtlExportTool, self).run()
        trit_path = self.option('trit_path')
        type_list = self.option('type').lower().strip().split(',')   # 得到list
        # self.qtl_export_api = self.api.api("dna_gmap.qtl_export")
        if self.option('popt').lower() in ['cp', 'f1']:
            lg_path_0 = self.option('sexAver_loc_path')   # workflow传进来
            lg_path_1 = self.option('sexAver_map_path')
            self.logger.info("cp_trit文件是{}".format(trit_path))
            self.link_path(lg_path_0, self.output_dir + "/" + "cp_total.sexAver.loc")
            self.link_path(lg_path_1, self.output_dir + "/" + "cp_total.sexAver.map")
            self.link_path(trit_path, self.output_dir + "/" + "cp_trit")
        else:
            lg_path_0 = self.option('total_csv_path')   # workflow传进来
            lg_path_1 = self.option('total_loc_path')   # workflow传进来
            lg_path_2 = self.option('total_map_path')   # workflow传进来
            self.logger.info("nocp_trit文件是{}".format(trit_path))
            if 'rqtl' in type_list or 'mapqtl' in type_list:
                if not os.path.exists(self.output_dir + "/rqtl_equal_mapqtl_nocp"):
                    os.mkdir(self.output_dir + "/rqtl_equal_mapqtl_nocp")
                self.link_path(lg_path_1, self.output_dir + "/rqtl_equal_mapqtl_nocp/" + "nocp_total.loc")
                self.link_path(lg_path_2, self.output_dir + "/rqtl_equal_mapqtl_nocp/" + "nocp_total.map")
                self.link_path(trit_path, self.output_dir + "/rqtl_equal_mapqtl_nocp/" + "nocp_trit")
            if 'qtlcart' in type_list:
                format = 'qtlcart'
                self.run_nocp(format, lg_path_0, trit_path)
                os.system('rm {}/*bak'.format(self.output_dir + "/qtlcart_nocp"))
        self.end()
