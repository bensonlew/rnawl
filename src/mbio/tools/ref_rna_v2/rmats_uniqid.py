# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from collections import OrderedDict
import itertools
import os
import subprocess


class RmatsUniqidAgent(Agent):
    def __init__(self, parent):
        super(RmatsUniqidAgent, self).__init__(parent)
        options = [
            {"name": "rmats_detail", "type": "string", "default": None},
            {"name": "significant_value", "type": "float", "default": 0.05},
            {"name": "delta_PSI", "type": "float", "default": 0.0},
            {"name": "significant_diff", "type": "string", "default": "FDR"},
            {"name": "compare_plan", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.step.add_steps("rmats_diffcomp")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.rmats_diffcomp.start()
        self.step.update()

    def stepfinish(self):
        self.step.rmats_diffcomp.finish()
        self.step.update()

    def check_options(self):
        if self.option("rmats_detail") is None:
            raise OptionError("rmats detail files is not set", code = "33707701")
        for path in self.option("rmats_detail").split(","):
            if not os.path.exists(path):
                raise OptionError("%s not exist", variables = (path), code = "33707702")

    def set_resource(self):
        self._cpu = 2
        self._memory = '4G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["diffcomp_stat.txt", "", "AS差异事件比较统计表"]
        ])
        super(RmatsUniqidAgent,self).end()

class RmatsUniqidTool(Tool):
    def __init__(self,config):
        super(RmatsUniqidTool, self).__init__(config)

    def get_uniqid(self):
        count = 0
        uniq_event = OrderedDict()
        result_list = self.option("rmats_detail").split(",")
        for (g_num,items) in enumerate(result_list):
            with open(items,'r') as f:
                _ = f.readline()
                for (num,line) in enumerate(f):
                    if line.split()[1] == "A3SS":
                        items_list = []
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                        if not uniq_event.has_key(tuple(uniq_site)):
                            uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                            items_list.append(g_num)
                            uniq_event[tuple(uniq_site)].append(items_list)
                        else:
                            #uniq_event[tuple(uniq_site)].append(g_num)
                            uniq_event[tuple(uniq_site)][1].append(g_num)
                    if line.split()[1] == "A5SS":
                        items_list = []
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                        if not uniq_event.has_key(tuple(uniq_site)):
                            uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                            items_list.append(g_num)
                            uniq_event[tuple(uniq_site)].append(items_list)
                        else:
                            uniq_event[tuple(uniq_site)][1].append(g_num)
                            #uniq_event[tuple(uniq_site)].append(g_num)
                    if line.split()[1] == "MAX":
                        items_list = []
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[22:26]
                        if not uniq_event.has_key(tuple(uniq_site)):
                            uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                            items_list.append(g_num)
                            uniq_event[tuple(uniq_site)].append(items_list)
                        else:
                            uniq_event[tuple(uniq_site)][1].append(g_num)
                            #uniq_event[tuple(uniq_site)].append(g_num)
                    if line.split()[1] == "RI":
                        items_list = []
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[20:22]
                        if not uniq_event.has_key(tuple(uniq_site)):
                            uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                            items_list.append(g_num)
                            uniq_event[tuple(uniq_site)].append(items_list)
                        else:
                            uniq_event[tuple(uniq_site)][1].append(g_num)
                            #uniq_event[tuple(uniq_site)].append(g_num)
                    if line.split()[1] == "SE":
                        items_list = []
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[8:14]
                        if not uniq_event.has_key(tuple(uniq_site)):
                            uniq_event[tuple(uniq_site)] = list(['AS_' + str(num + count)])
                            items_list.append(g_num)
                            uniq_event[tuple(uniq_site)].append(items_list)
                        else:
                            uniq_event[tuple(uniq_site)][1].append(g_num)
                            #uniq_event[tuple(uniq_site)].append(g_num)
                            #+ line.split()[39:40] + line.split()[33:35]
            count = count + num + 1
        return uniq_event

    def diffgroup_stat(self):
        #diffgroup = []
        uniq_event_new = self.get_uniqid()
        result_list = self.option("rmats_detail").split(",")
        for items in result_list:
            #diffgroup.append(items.split('/')[-2])
            with open(items,'r') as f:
                _ = f.readline()
                for line in f:
                    if self.option("significant_diff").lower() == "fdr":
                        significant_value = line.split()[34]
                        significant_value_all = line.split()[52]
                    elif self.option("significant_diff").replace(' ','').lower() == "pvalue":
                        significant_value = line.split()[33]
                        significant_value_all = line.split()[51]
                    if line.split()[1] == "A3SS":
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                        if tuple(uniq_site) in uniq_event_new.keys():
                            if line.split()[39] == "null":
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[39])) > self.option("delta_PSI") and float(significant_value) < self.option("significant_value") :
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            else:
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                    if line.split()[1] == "A5SS":
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[14:20]
                        if tuple(uniq_site) in uniq_event_new.keys():
                            if line.split()[39] == "null":
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[39])) > self.option("delta_PSI") and float(significant_value) < self.option("significant_value") :
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            else:
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                    if line.split()[1] == "MAX":
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[22:26]
                        if tuple(uniq_site) in uniq_event_new.keys():
                            if line.split()[39] == "null":
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[39])) > self.option("delta_PSI") and float(significant_value) < self.option("significant_value") :
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            else:
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                    if line.split()[1] == "RI":
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[10:14] + line.split()[20:22]
                        if tuple(uniq_site) in uniq_event_new.keys():
                            if line.split()[39] == "null":
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[39])) > self.option("delta_PSI") and float(significant_value) < self.option("significant_value") :
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            else:
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                    if line.split()[1] == "SE":
                        uniq_site = line.split()[1:2] + line.split()[4:6] + line.split()[8:14]
                        if tuple(uniq_site) in uniq_event_new.keys():
                            if line.split()[39] == "null":
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["null"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            elif abs(float(line.split()[39])) > self.option("delta_PSI") and float(significant_value) < self.option("significant_value") :
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["yes"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
                            else:
                                if line.split()[57] == "null":
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "null"] + line.split()[57:58] + line.split()[51:53])
                                elif abs(float(line.split()[57])) >self.option("delta_PSI") and float(significant_value_all) <self.option("significant_value") :
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "yes"] + line.split()[57:58] + line.split()[51:53])
                                else:
                                    uniq_event_new[tuple(uniq_site)].append(
                                        line.split()[0:1] + ["no"] + line.split()[39:40] + line.split()[33:35] + [
                                            "no"] + line.split()[57:58] + line.split()[51:53])
        return uniq_event_new

    def get_finaltable(self):
        title_detail = []
        result_list = self.option("rmats_detail").split(",")
        uniq_event_new = self.diffgroup_stat()
        detail_list = ['_AS_ID','_diff_significant_JC','_delta_PSI_JC','_Pvalue_JC','_FDR_JC','_diff_significant_JCEC','_delta_PSI_JCEC','_Pvalue_JCEC','_FDR_JCEC']
        diffgroup = self.option("compare_plan").split(",")
        header_list = list(itertools.product(diffgroup,detail_list))
        for i in header_list:
            single = ''.join(list(i))
            title_detail.append(single)
        with open('diffcomp_stat.txt','w') as f:
            f.write("New_ID\tType\tGene_ID\tGene_name\t" + '\t'.join(title_detail) + "\n")
            for v in uniq_event_new.keys():
                comp_result = []
                for (num,items) in enumerate(result_list):
                    if num in uniq_event_new[v][1]:
                        comp_result.insert(num,uniq_event_new[v].pop(2))
                    else:
                        comp_result.insert(num,['_','_','_','_','_','_','_','_','_'])
                result_all = uniq_event_new[v][0:1] + list(v)[0:3] + list(itertools.chain(*comp_result))
                result_all_str = "\t".join(result_all)
                f.write(result_all_str + "\n")

    def set_output(self):
        cmd = """sed -i \'s/\tnull/\t_/g\' %s""" % (self.work_dir + '/diffcomp_stat.txt')
        subprocess.call(cmd, shell=True)
        diffcomp_stat = self.work_dir + '/diffcomp_stat.txt'
        if os.path.exists(self.output_dir + '/diffcomp_stat.txt'):
            os.remove(self.output_dir + '/diffcomp_stat.txt')
        os.link(diffcomp_stat, os.path.join(self.output_dir, 'diffcomp_stat.txt'))

    def run(self):
        super(RmatsUniqidTool, self).run()
        self.get_finaltable()
        self.set_output()
        self.end()
