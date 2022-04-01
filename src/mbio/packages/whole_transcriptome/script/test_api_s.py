# -*- coding:utf-8 -*-
# __author__ = 'qinjincheng'

import json
import os
import unittest

from biocluster.config import Config


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'smallrna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False

        task_id = 'smallrna'
        project_sn = 'smallrna'

        ################################################################
        api_qc = wf.api.api('whole_transcriptome.qc')
        sample_list = '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/smallrna/rawdata/list.txt'
        group_file = '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/group.txt'
        compare_file = '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/control.txt'
        qc_stat_before = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/HiseqReadsStat/output'
        qc_stat_after = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/HiseqReadsStat1/output'
        qc_result_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/MirnaQc/output/clean_data'

        api_qc.add_sample_info(sample_list=sample_list, library='small')
        group_id, specimen_names, category_names = api_qc.add_sample_group(group_file=group_file, library='small')
        control_id, compare_names = api_qc.add_group_compare(compare_file=compare_file, library='small',
                                                             group_id=group_id)
        qc_id = api_qc.add_qc(fq_type='PE', library='small')
        api_qc.add_qc_detail(qc_id, qc_stat_before, qc_stat_after, 'small', qc_result_dir)
        api_qc.add_qc_graph(qc_id, qc_stat_before, 'small', 'before')
        api_qc.add_qc_graph(qc_id, qc_stat_after, 'small', 'after', qc_result_dir)
        ################################################################

        ################################################################
        database = Config().get_mongo_client(mtype='ref_rna_v2')[Config().get_mongo_dbname('ref_rna_v2')]
        collection = database['sg_genome_db']
        wf.genome_doc = collection.find_one({'genome_id': 'GM0259'})
        wf.db_path = os.path.join(wf.config.SOFTWARE_DIR, 'database/Genome_DB_finish')

        api_genome_info = wf.api.api('whole_transcriptome.genome_info')
        file_path = os.path.join(wf.db_path, wf.genome_doc['gene_stat'])
        species_name = 'Homo_sapiens'
        ref_anno_version = wf.genome_doc['assembly']
        hyperlink = wf.genome_doc['ensemble_web']
        api_genome_info.add_genome_info(file_path, species_name, ref_anno_version, hyperlink)
        ################################################################

        ################################################################
        api_mapping = wf.api.api('whole_transcriptome.mapping')
        stat_file = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/MapperAndStat/output'
        distribution = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/MapperAndStat/output'
        sample_list_file = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/MirnaQc/output/clean_data/list.txt'
        method = 'bowtie'
        params = json.dumps({'task_id': task_id, 'submit_location': 'mapping', 'task_type': 2}, sort_keys=True)
        api_mapping.add_mapping_stat(stat_file=stat_file, library='small', method=method,
                                     sample_list_file=sample_list_file)
        api_mapping.add_chrom_distribution_table(distribution=distribution, params=params, library='small',
                                                 sample_list_file=sample_list_file)
        ################################################################

        ################################################################
        api_srna = wf.api.api('whole_transcriptome.srna')
        mirna_stat = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/Srna/SrnaStat/output/mirna_stat.xls'
        params = json.dumps({'task_id': task_id, 'submit_location': 'mirnastat', 'task_type': 2, 'method': 'mirdeep'},
                            sort_keys=True)
        api_srna.add_mirna_stat(mirna_stat=mirna_stat, task_id=task_id, project_sn=project_sn, params=params)
        ################################################################

        ################################################################
        api_srna = wf.api.api('whole_transcriptome.srna')
        srna_stat = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/Srna/SrnaStat/output/srna_stat.xls'
        srna_stat_for_graph = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/Srna/output/srna_stat/srna_stat_for_graph.xls'
        params = json.dumps({'task_id': task_id, 'submit_location': 'mirnastat', 'task_type': 2, 'method': 'mirdeep'},
                            sort_keys=True)
        api_srna.add_srna_stat(srna_stat=srna_stat, srna_stat_for_graph=srna_stat_for_graph, task_id=task_id,
                               project_sn=project_sn, params=params)
        ################################################################

        ################################################################
        api_diff_exp_stat = wf.api.api('whole_transcriptome.diff_exp_stat')
        map_dict = {
            'control': '/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/whole_transcriptome/demo_human/control.txt',
            's_output_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Smallrna_workflow_4521_6888/DiffExp/output'
        }
        arg_dict = {'program': 'DESeq2', 'fc': 2.0, 'qvalue': 0.05, 'method': 'BH'}
        api_diff_exp_stat.add_diff_exp_stat(map_dict=map_dict, arg_dict=arg_dict, task_id=task_id, project_sn=project_sn, library='small')
        ################################################################


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
