# -*- coding: utf-8 -*-
# __author__ = 'hongdong'
# last_modify:20171211
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import re
import os


class FamilyAnalysisModule(Module):
    """
    做亲子鉴定的分析。
    包含tool：family_merge、family_analysis, 后面还可以将画图的tool添加进去
    """
    def __init__(self, work_id):
        super(FamilyAnalysisModule, self).__init__(work_id)
        self.step.add_steps('family_merge', 'family_analysis')
        options = [
            {"name": "dad_tab", "type": "infile", "format": "paternity_test.tab"},  # 输入F/M/S的样本ID
            {"name": "mom_tab", "type": "infile", "format": "paternity_test.tab"},  # fastq所在路径
            {"name": "preg_tab", "type": "infile", "format": "paternity_test.tab"},
            {"name": "ref_point", "type": "infile", "format": "paternity_test.rda"},
            {"name": "err_min", "type": "int", "default": 2},
        ]
        self.add_option(options)
        self.family_merge = None
        self.father_analysis = None
        self._end_info = 0

    def check_options(self):
        if not self.option("dad_tab"):
            raise OptionError("必须输入父本tab文件")
        if not self.option("ref_point"):
            raise OptionError("必须输入参考位点的bed文件")
        if not self.option("mom_tab"):
            raise OptionError("必须输入母本tab文件")
        if not self.option("preg_tab"):
            raise OptionError("必须输入胎儿tab文件")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def merge_run(self):
        self.family_merge = self.add_tool("medical.paternity_test_v2.family_merge")
        self.family_merge.set_options({
            "dad_tab": self.option("dad_tab"),
            "mom_tab": self.option("mom_tab"),
            "preg_tab": self.option("preg_tab"),
            "ref_point": self.option("ref_point"),
            "err_min": self.option("err_min")
        })
        self.family_merge.on('end', self.set_output, 'family_merge')
        self.family_merge.on('end', self.analysis_run)
        self.family_merge.on('start', self.set_step, {'start': self.step.family_merge})
        self.family_merge.on('end', self.set_step, {'end': self.step.family_merge})
        self.family_merge.run()

    def analysis_run(self):
        self.father_analysis = self.add_tool("medical.paternity_test_v2.father_analysis")
        results = os.listdir(self.family_merge.output_dir)
        rdata = ""
        for f in results:
            if re.match(r'.*family_joined_tab\.Rdata$', f):
                rdata = f
            else:
                print "Oops!"
        self.father_analysis.set_options({
            "tab_merged": os.path.join(self.family_merge.output_dir, rdata)
        })
        self.father_analysis.on('end', self.set_output, "family_analysis")
        self.father_analysis.on('end', self.end)
        self.father_analysis.on('start', self.set_step, {'start': self.step.family_analysis})
        self.father_analysis.on('end', self.set_step, {'end': self.step.family_analysis})
        self.father_analysis.run()

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
        if event['data'] == 'family_merge':
            self.linkdir(obj.output_dir, 'family_merge')
        elif event['data'] == 'family_analysis':
            self.linkdir(obj.output_dir, 'family_analysis')
        else:
            pass

    def run(self):
        self.merge_run()
        super(FamilyAnalysisModule, self).run()

    def end(self):
        repaths = [
            [".", "", "无创亲子鉴定结果输出目录"],
        ]
        regexps = [
            [r"FamilyMerge/output/family_merge_tab.Rdata", "Rdata", "家系合并后的表格"],
            [r"FamilyMerge/output/family_merge_tab.xlsx", "xlsx", "家系合并后的表格,excel格式"],
            [r"FamilyAnalysis/output/family_analysis.Rdata", "Rdata", "父权值等信息计算表格"]
        ]

        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        # print regexps
        super(FamilyAnalysisModule, self).end()
