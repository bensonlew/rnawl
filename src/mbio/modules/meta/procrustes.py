# -*- coding: utf-8 -*-


import os
import glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import unittest

class ProcrustesModule(Module):

    def __init__(self,work_id):
        super(ProcrustesModule, self).__init__(work_id)
        options = [
            {'name': 'otu_table', 'type': 'infile', 'format': 'sequence.profile_table','required':True},
            {'name': 'asso_table', 'type': 'infile', 'format': 'sequence.profile_table,meta.otu.otu_table','required':True},
            {'name': 'otu_dist', 'type': 'string', 'default': 'bray_curtis'},
            {'name': 'asso_dist', 'type': 'string', 'default': 'bray_curtis'},
            {'name': 'method', 'type': 'string', 'default': 'pcoa'}
        ]
        self.add_option(options)
        self.pca_ref_tool = self.add_tool("meta.beta_diversity.pca")
        self.pca_query_tool = self.add_tool("meta.beta_diversity.pca")
        self.discalc_ref_tool = self.add_tool("meta.beta_diversity.distance_calc")
        self.discalc_query_tool = self.add_tool("meta.beta_diversity.distance_calc")
        self.pcoa_ref_tool = self.add_tool("meta.beta_diversity.pcoa")
        self.pcoa_query_tool = self.add_tool("meta.beta_diversity.pcoa")
        self.procrustes_tool = self.add_tool("meta.procrustes")
        self.run_tools = [] #并行运行的tool

    def check_options(self):
        self.method = self.option('method')
        if self.method not in ["pca", "pcoa"]:
            raise OptionError("PARAMETERS ERROR:wrong value of method (%s), pca or pcoa expected!" , variables=( self.method), code="22701301")
        return True

    def discalc_run_ref(self):
        """
        运行计算reference距离矩阵，代谢物丰度文件
        :return:
        """
        options = {
            "otutable": self.option('otu_table').prop['path'],
            "method": self.option("otu_dist")
        }
        self.discalc_ref_tool.set_options(options)
        self.discalc_ref_tool.run()

    def discalc_run_query(self):
        """
        运行计算query距离矩阵，比如自定义上传的文件
        :return:
        """
        options = {
            "otutable": self.option('asso_table').prop['path'],
            "method": self.option("asso_dist")
        }
        self.discalc_query_tool.set_options(options)
        self.discalc_query_tool.run()

    def pcoa_run_ref(self):
        options = {
            'dis_matrix': self.discalc_ref_tool.option('dis_matrix'),
            #'group': self.option('group_file').prop['path']
        }
        self.pcoa_ref_tool.set_options(options)
        self.pcoa_ref_tool.run()

    def pcoa_run_query(self):
        options = {
            'dis_matrix': self.discalc_query_tool.option('dis_matrix'),
            #'group': self.option('group_file').prop['path']
        }
        self.pcoa_query_tool.set_options(options)
        self.pcoa_query_tool.run()

    def pca_run_ref(self):
        options = {
            'otutable': self.option('otu_table').prop['path'],
            #'group_table': self.option('group_file').prop['path']
        }
        self.pca_ref_tool.set_options(options)
        self.pca_ref_tool.run()

    def pca_run_query(self):
        options = {
            'otutable': self.option('asso_table').prop['path'],
            #'group_table': self.option('group_file').prop['path']
        }
        self.pca_query_tool.set_options(options)
        self.pca_query_tool.run()

    def procrustes_run(self):
        coord_ref = ''
        coord_query = ''
        if self.method == 'pca':
            coord_ref = os.path.join(self.pca_ref_tool.work_dir, "pca/pca.txt")
            coord_query = os.path.join(self.pca_query_tool.work_dir, "pca/pca.txt")
        elif self.method == 'pcoa':
            coord_ref = os.path.join(self.pcoa_ref_tool.work_dir, "pcoa/pcoa.txt")
            coord_query = os.path.join(self.pcoa_query_tool.work_dir, "pcoa/pcoa.txt")

        options = {
            'coord_ref': coord_ref,
            'coord_query': coord_query
        }
        self.procrustes_tool.set_options(options)
        self.procrustes_tool.run()



    def set_output(self):
        target_files = glob.glob(self.procrustes_tool.output_dir + '/*transformed_reference.txt')
        target_files += glob.glob(self.procrustes_tool.output_dir + '/*transformed_q*.txt')
        target_files += glob.glob(self.procrustes_tool.output_dir + '/procrustes_results.txt')
        for each in target_files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
        self.end()


    def run(self):
        super(ProcrustesModule, self).run()

        '''
        定义依赖关系和执行顺序
        '''

        if self.method == 'pca':
            self.run_tools.append(self.pca_ref_tool)
            self.run_tools.append(self.pca_query_tool)
        elif self.method == 'pcoa':
            self.run_tools.append(self.pcoa_ref_tool)
            self.run_tools.append(self.pcoa_query_tool)
            self.discalc_ref_tool.on('end', self.pcoa_run_ref)
            self.discalc_query_tool.on('end', self.pcoa_run_query)

        if self.method == 'pca':
            self.pca_run_ref()
            self.pca_run_query()
        elif self.method == 'pcoa':
            self.discalc_run_ref()
            self.discalc_run_query()
        self.on_rely(self.run_tools,self.procrustes_run)
        self.procrustes_tool.on("end", self.set_output)


        
    def end(self):
        super(ProcrustesModule, self).end()



class TestFunction(unittest.TestCase):
    '''
    This is test for the module. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liulinmeng/metabolome/package/procrustes/module'
        method = "pcoa"
        data = {
            "id": "procrustes" + str(random.randint(1,10000)),
            "type": "module",
            "name": "metabolome.procrustes",
            "instant": False,
            "options": dict(
                metab_table = test_dir + "/metab.select_table.xls",
                asso_table = test_dir + "/corr.select_table.xls",
                #group_file = test_dir + "/group.txt",
                metab_dist = "bray_curtis",
                asso_dist = "bray_curtis",
                method = method,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

