# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20171211
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from biocluster.config import Config
import gevent
import os


class MsMatchModule(Module):
    """
    找出所有的母本进行计算胎儿与母本的匹配度
    lasted modified by hongdong@20180821
    包含tool：ms_match
    """
    def __init__(self, work_id):
        super(MsMatchModule, self).__init__(work_id)
        options = [
            {"name": "dad_tab", "type": "infile", "format": "paternity_test.tab"},  # 输入F/M/S的样本ID
            {"name": "mom_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "preg_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "ref_point", "type": "infile", "format": "paternity_test.rda"},
            {"name": "err_min", "type": "int", "default": 2}
        ]
        self.add_option(options)
        self.tab_data = Config().SOFTWARE_DIR + "/database/human/pt_ref/tab_all/"
        self.api_pt = self.api.api('medical.paternity_test_v3.paternity_test_v3')
        self.ms_matchs = []

    def check_options(self):
        if not self.option("dad_tab"):
            raise OptionError("必须输入父本tab文件")
        if not self.option("ref_point"):
            raise OptionError("必须输入参考位点的bed文件")
        # if not self.option("mom_tab"):
        #     raise OptionError("必须输入母本tab文件")
        if not self.option("preg_tab"):
            raise OptionError("必须输入胎儿tab文件")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_ms_match(self):
        samples = self.get_all_mom()
        for sample in samples:
            ms_match = self.add_tool("medical.paternity_test_v3.ms_match")
            ms_match.set_options({
                "dad_tab": self.option("dad_tab"),
                "mom_tab": os.path.join(self.tab_data, sample) + ".tab",
                "preg_tab": self.option("preg_tab"),
                "ref_point": self.option("ref_point"),
                "err_min": self.option("err_min")
            })
            self.ms_matchs.append(ms_match)
        for j in range(len(self.ms_matchs)):
            self.ms_matchs[j].on('end', self.set_output, 'ms_match')
        if self.ms_matchs:
            if len(self.ms_matchs) > 1:
                self.on_rely(self.ms_matchs, self.end)
            elif len(self.ms_matchs) == 1:
                self.ms_matchs[0].on('end', self.end)
        else:
            raise Exception("ms_matchs列表为空！")
        for module in self.ms_matchs:
            gevent.sleep(1)
            module.run()

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
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'ms_match':
            self.linkdir(obj.output_dir, self.output_dir)

    def get_all_mom(self):
        """
        获取所有的母本的tab文件
        :return:
        """
        sample_list = []
        mom_id = os.path.basename(self.option("mom_tab").prop['path']).split(".")[0]
        self.logger.info(mom_id)
        samples = self.api_pt.find_sample_batch(mom_id, "M")
        if samples:
            for m in samples:
                if self.api_pt.check_exists(os.path.join(self.tab_data, m) + '.tab'):
                    sample_list.append(m)
                else:
                    self.logger.info("样本在tab_all文件夹中不存在，请确认是不是异常样本！".format(m))
        else:
            sample_list.append(mom_id + ".tab")
        return sample_list

    def run(self):
        super(MsMatchModule, self).run()
        self.run_ms_match()

    def end(self):
        self.find_match_mom()
        super(MsMatchModule, self).end()

    def find_match_mom(self):
        with open(os.path.join(self.output_dir, "ms_match.txt"), "w") as w:
            w.write("#match_mom_id\n")
            for m in os.listdir(self.output_dir):
                if m != 'ms_match.txt':
                    mom_preg = self.open_one_file(os.path.join(self.output_dir, m))
                    if mom_preg:
                        w.write("{}\n".format(mom_preg))

    def open_one_file(self, file_path):
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            temp = data[0].strip().split("\t")
            if float(temp[7]) >= 95:
                return temp[5]
            else:
                return ""
