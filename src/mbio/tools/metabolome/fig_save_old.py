# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
import json
from bson import ObjectId


class FigSaveOldAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(FigSaveOldAgent, self).__init__(parent)
        options = [
            {"name": "run_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "table_id", "type": "string"},
            {"name": "table_name", "type": "string"},
            {"name": "submit_loc", "type": "string"},
            {"name": "project", "type": "string"},
            {"name": "proj_dir", "type": "string"},
            {"name": "interaction", "type": "int", "default": 1},
        ]
        self.add_option(options)

    def check_options(self):
        pass
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(FigSaveOldAgent, self).end()

class FigSaveOldTool(Tool):
    def __init__(self, config):
        super(FigSaveOldTool, self).__init__(config)
        base_dir = self.config.SOFTWARE_DIR + '/bioinfo/figsave'
        self.python = "/miniconda2/bin/python"
        self.fig_save = base_dir + '/fig_save'
        self.set_environ(PATH=self.fig_save, PYTHONPATH=base_dir + '/lib/python2.7/site-packages',
                         LD_LIBRARY_PATH=base_dir + "/x86_64-conda-linux-gnu/sysroot/usr/lib")
        self.common_api = self.api.api("metagenomic.common_api")
        self.common_api._project_type = self.option("project")
        self.run_info = self.common_api.db['sg_task'].find_one({"task_id": self.option("task_id")})
        self.meth2out = {}
        self.subloc2col = {}
        self.get_subloc2col()

    def run(self):
        """
        运行
        :return:
        """
        super(FigSaveOldTool, self).run()
        self.tmp_dir = os.path.join(self.work_dir, "tmp_out")
        cmd = "{} {}/fig2pdf_old.py -proj {} -task_id {} -tmp_out {}".format(self.python, self.fig_save,
                                                                         self.option("project"),
                                                                         self.option("task_id"), self.tmp_dir)
        if not self.option("interaction"):
            cmd += " -interaction 0 -proj_kws '{}'".format(json.dumps({"proj_type": self.run_info["type"]}))
        else:
            self.update_status(self.option("table_id"), self.option("submit_loc"), {"status": "end"})
            cmd += " -sub_loc {} -table_id {} -table_name {}".format(self.option("submit_loc"),
                                                                    self.option("table_id"),
                                                                    self.option("table_name"))
        proj_kws = {"proj_type": self.run_info["type"], "mix": 'F'}
        if "mix" in self.run_info:
            proj_kws["mix"] = self.run_info["mix"]
        cmd += " -proj_kws '{}'".format(json.dumps(proj_kws))
        if self.config.DBVersion:
            cmd += " -db_version {}".format(self.config.DBVersion)
        if not self.command_done('fig_save'):
            command = self.add_command("fig_save", cmd).run()
            self.wait(command)
        self.set_output()
        self.end()

    def command_done(self, name):
        o_file = name + '.o'
        if os.path.exists(o_file):
            filesize = os.path.getsize(o_file)
            offset = 200 if filesize > 200 else filesize
            with open(o_file, 'rb') as r:
                r.seek(-offset, 2)
                tail_info = r.read()
                if "exitcode:0" in tail_info:
                    return True
        return False

    def update_status(self, obj_id, sub_loc, update_info):
        print(obj_id)
        print(sub_loc)
        print(self.common_api.db)
        main_col = self.common_api.db[self.subloc2col[sub_loc]]
        main_col.update_one({"_id": ObjectId(obj_id)}, {"$set": update_info})

    def get_subloc2col(self):
        fig_info = os.path.join(self.fig_save, "conf/metab_fig_info.txt")
        with open(fig_info, 'r') as r:
            for line in r:
                line = line.strip().split('\t')
                if not self.option("interaction"):
                    self.meth2out[line[0]] = line[4].replace(' ', '')
                else:
                    self.meth2out[line[0]] = line[4] and line[4].split(' ')[1]
                self.subloc2col[line[1]] = line[2]

    def set_output(self):
        if not self.option("interaction"):
            with open(os.path.join(self.work_dir, "workflow_tables.txt"), 'r') as r:
                for line in r:
                    line = line.strip().split('\t')
                    table_id = str(line[1])
                    table_name = line[2]
                    self.out(table_id, table_name, line[0])
        else:
            self.out(self.option("table_id"), self.option("table_name"),
                        self.option("submit_loc"))

    def out(self, table_id, table_name, sub_loc):
        table_type = 'pos'
        ana_name = ''
        set_name = ''
        if not self.option("interaction"):
            tb = table_name.split("_Origin_")
            if len(tb) == 2:
                table_type = tb[1]
                set_name = tb[1]
            else:
                if 'mix' in self.run_info and self.run_info["mix"] == 'T':
                    table_type = 'mix'
                elif table_name == 'raw' and sub_loc == "exp":
                    table_type = 'raw'
                else:
                    table_type = 'pos'
            if set_name in ['', 'pos', 'neg', 'mix']:
                set_name = "DiffSet_mix"
        if self.option("interaction"):
            ana_name = table_name
        image_opt = os.path.join(self.tmp_dir, table_id + '.sg_image_options.txt')
        with open(image_opt, 'r') as r:
            for line in r:
                line = line.strip().split('\t')
                params = eval(line[1])
                if "table_type" not in params or params["table_type"] == 'all':
                    params["table_type"] = table_type
                if "type" in params and sub_loc == "exp":
                    params["table_type"] = params['type']
                if "metabset_name" in params:
                    params["set_name"] = params["metabset_name"]
                else:
                    params["set_name"] = set_name
                params["ana_name"] = ana_name
                out_path = self.meth2out[params["method"]].format(**params)
                out_path = out_path.replace('__', '_')
                out_path = os.path.join(self.output_dir, out_path)
                if not os.path.exists(os.path.dirname(out_path)):
                    os.makedirs(os.path.dirname(out_path))
                if os.path.exists(out_path):
                    os.remove(out_path)
                if os.path.exists(os.path.join(self.tmp_dir, line[0] + '.pdf')):
                    os.link(os.path.join(self.tmp_dir, line[0] + '.pdf'), out_path)
        self.update_status(table_id, sub_loc, {"pdf_saved": 1})