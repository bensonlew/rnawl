# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200416

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class ManhattanWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ManhattanWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "snp_table", "type": "infile", "format": "tool_lab.snp_table"},
            {"name": "guide", "type": "string"},  # 不绘制传"False", 否则为"True"
            {"name": "p_threshold1", "type": "string"},
            {"name": "p_threshold2", "type": "string"},
            {"name": "color1", "type": "string"},
            {"name": "color2", "type": "string"},
            {"name": "re_axis", "type": "string"},
            {"name": "one_chr", "type": "string"},
            {"name": "chr", "type": "int"},
            {'name': 'main_id', 'type': 'string'},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.manhattan = self.add_tool("tool_lab.manhattan")
        self.pdf2img = self.add_tool("tool_lab.pdf2image")

    def check_options(self):
        if not self.option("snp_table"):
            raise OptionError("必须输入snp_table文件")
        if not self.option("color1"):
            raise OptionError("必须输入散点颜色1")
        if not self.option("color2"):
            raise OptionError("必须输入散点颜色2")
        if not self.option("guide"):
            raise OptionError("必须输入是否显示显著性参考线")
        return True

    def run_manhattan(self):
        guides2 = "yes"
        if self.option("guide") != "False":
            p_threshold1 = self.option('p_threshold1')
        else:
            p_threshold1 = 0
        if not self.option("p_threshold2"):
            guides2 = "False"
            p_threshold2 = 5e-08
        else:
            p_threshold2 = self.option("p_threshold2")
        if not self.option("re_axis"):
            re_axis = "-1"
        else:
            re_axis = self.option("re_axis")
        if self.option("one_chr") == "True":
            chr_num = self.option("chr")
        else:
            chr_num = "-1"
        options = {
            "snp_table": self.option('snp_table'),
            "guides2": guides2,
            "p_threshold1": p_threshold1,
            'p_threshold2': p_threshold2,
            "color": self.option('color1') + "," + self.option('color2'),
            "re_axis": re_axis,
            "one_chr": chr_num,
        }
        self.manhattan.set_options(options)
        self.manhattan.on("end", self.set_output, "manhattan")
        self.manhattan.run()
    def run_pdf2image(self):
        options = {
            "pdf": self.output_dir + "/manhattan/qqman.pdf"
        }
        self.pdf2img.set_options(options)
        self.pdf2img.on("end", self.set_output, "pdf2img")
        self.pdf2img.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'manhattan':
            self.linkdir(obj.output_dir, 'manhattan')
            self.run_pdf2image()
        elif event['data'] == 'pdf2img':
            self.linkdir(obj.output_dir, 'pdf2img')
            png_file = os.listdir(self.output_dir + '/pdf2img/')
            new_png_file = self.output_dir + '/pdf2img/manhattan.png'
            if os.path.exists(new_png_file):
                os.remove(new_png_file)
            os.link(self.output_dir + '/pdf2img/' + png_file[0], new_png_file)
            self.set_db()

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

    def set_db(self):
        self.logger.info("保存结果到mongo")
        api_manhattan = self.api.api("tool_lab.manhattan")
        api_manhattan.add_sg_manhattan_path(self.option("main_id"),
                                            self._sheet.output + "/pdf2img/manhattan.png")
        self.end()

    def run(self):
        self.run_manhattan()
        super(ManhattanWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(ManhattanWorkflow, self).end()
