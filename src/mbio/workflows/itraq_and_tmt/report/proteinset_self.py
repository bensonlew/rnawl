# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import os


class ProteinsetSelfWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ProteinsetSelfWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='proteins', type='string'),
            dict(name='trait_path', type='string'),
            dict(name="name", type="string", default=None),
            dict(name="update_info", type='string'),
            dict(name="proteinset_id", type='string'),
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        self.match_db(self.option("trait_path"), self.option("proteins"))
        self.set_db()

    def set_db(self):
        proteinset_self = self.api.api("itraq_and_tmt.proteinset")
        proteinset_table = os.path.join(self.output_dir, 'proteinset_self.txt')
        proteinset_self.add_proteinset_self(proteinset_output_dir=proteinset_table, main_id=self.option('proteinset_id'), name=self.option('name'))
        self.end()

    def end(self):
        super(ProteinsetSelfWorkflow, self).end()

    def match_db(self, file, proteins):
        with open(proteins) as pr:
            query_dict = {query.strip(): 1 for query in pr if query.strip()}

        with open(self.work_dir + "/proteinset_self.txt", "wb") as f, open(file,'r') as proteinf:
            protein_dict = {protein.strip(): 1 for protein in proteinf if protein.strip()}
            f.write("protein_list\n")
            for protein in protein_dict:
                if query_dict.get(protein):
                    f.write(protein + "\n")

        self.set_output()

    def set_output(self):
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        proteinset_file = self.work_dir + "/" + "proteinset_self.txt"
        os.link(proteinset_file, os.path.join(self.output_dir, "proteinset_self.txt"))
        self.logger.info("设置蛋白集创建结果目录")