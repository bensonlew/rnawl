# -*- coding: utf-8 -*-
# __author__ = 'hongyu'
# last_modify:20171211
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import os
import gevent


class DailyDedupModule(Module):
    """
    用于全库排查守护进程
    """
    def __init__(self, work_id):
        super(DailyDedupModule, self).__init__(work_id)
        self.step.add_steps('daily_dedup')
        options = [
            {"name": "dad_id", "type": "string"},
            {"name": "mom_id", "type": "string"},
            {"name": "preg_id", "type": "string"},
            {"name": "ref_point", "type": "string"},
            {"name": "ref_data", "type": "string"},
            {"name": "ref_all", "type": "string"},
            {"name": "err_min", "type": "int"},
            {"name": "mem", "type": "int"},
            {"name": "dad_list", "type": "string"}
        ]
        self.add_option(options)
        self.tools_dedup = []
        self.api_pt = self.api.api("medical.paternity_test_v3.paternity_test_v3")
        self.err_list = [5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50]

    def check_options(self):
        for name in ["dad_id", "mom_id", "preg_id", "ref_point", "ref_data", "ref_all", "err_min"]:
            if not self.option(name):
                raise OptionError("必须设置{}参数！".format(name))
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def check_table(self):
        self.api_pt.attempt_to_export_tab_file(self.option("mom_id"), self.output_dir)
        self.api_pt.attempt_to_export_tab_file(self.option("preg_id"), self.output_dir)

    def dedup_run(self):
        n = 0
        # for p in range(2, self.option('err_min')):
        dedup = self.add_tool("medical.paternity_test_v3.dedup_multi")
        self.step.add_steps('dedup_{}'.format(n))
        dedup.set_options({
            "dad_id": self.option("dad_id"),
            "mom_id": self.option("mom_id"),
            "preg_id": self.option("preg_id"),
            "mom_tab": self.output_dir + '/' + self.option("mom_id") + '.tab',
            "preg_tab": self.output_dir + '/' + self.option("preg_id") + '.tab',
            "ref_point": self.option("ref_point"),
            "err_min": self.option('err_min'),
            "dad_list": self.option("dad_list"),
            "father_path": self.option("ref_all"),
            "mem": self.option("mem")
        })
        step = getattr(self.step, 'dedup_{}'.format(n))
        step.start()
        dedup.on('end', self.finish_update, 'dedup_{}'.format(n))
        self.tools_dedup.append(dedup)
        # n += 1
        for j in range(len(self.tools_dedup)):
            self.tools_dedup[j].on('end', self.set_output, 'dedup')
        if len(self.tools_dedup) > 1:
            self.on_rely(self.tools_dedup, self.end)
        else:
            self.tools_dedup[0].on("end", self.end)
        for tool in self.tools_dedup:
            tool.run()

    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到指定目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
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
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                next_dirname = dirname + '/' + os.path.basename(oldfiles[i])
                self.linkdir(oldfiles[i], next_dirname)

    def set_output(self, event):
        obj = event["bind_object"]
        if event['data'] == "dedup":
            self.linkdir(obj.output_dir, self.output_dir)

    def set_db(self):
        for i in self.err_list:
            dir_path = self.output_dir + '/pt_result_' + str(i)
            if not os.path.exists(dir_path):
                continue
            result_name = self.option("dad_id") + "_" + self.option("mom_id") + "_" + self.option("preg_id")
            results = os.listdir(dir_path)
            for f in results:
                if f == result_name + '.txt':
                    father_dedup_id = self.api_pt.find_dedup_err_id(i, self.option("mom_id"), self.option("preg_id"))
                    if not father_dedup_id:
                        break
                    self.api_pt.import_dedup_data(dir_path + '/' + f, father_dedup_id)
        father_ids = self.api_pt.search_father_id(self.option("mom_id"), self.option("preg_id"))
        for father_id in father_ids:
            for j in self.err_list:
                self.api_pt.update_dedup_father(j, self.option("mom_id"), self.option("preg_id"), father_id)

    def run(self):
        self.check_table()
        self.dedup_run()
        super(DailyDedupModule, self).run()
    
    def end(self):
        # self.set_db()
        super(DailyDedupModule, self).end()
