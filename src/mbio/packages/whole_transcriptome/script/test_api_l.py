# -*- coding:utf-8 -*-
# __author__ = 'qinjincheng'

import json
import unittest


class TestFunction(unittest.TestCase):
    '''
    This is test for the api. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.whole_transcriptome.whole_transcriptome_test_api import WholeTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'longrna_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = WholeTranscriptomeTestApiWorkflow(wheet)
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False

        task_id = 'whole_transcriptome'
        project_sn = 'whole_transcriptome'

        ################################################################
        api_qc = wf.api.api('whole_transcriptome.qc')
        sample_list = '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190118060/largeRNA/little/list.txt'
        group_file = '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190118060/largeRNA/example_group.txt'
        compare_file = '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190118060/largeRNA/example_control.txt'
        qc_stat_before = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/HiseqReadsStat/output'
        qc_stat_after = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/HiseqReadsStat1/output'

        api_qc.add_sample_info(sample_list=sample_list, library='long')
        group_id, specimen_names, category_names = api_qc.add_sample_group(group_file=group_file, library='long')
        control_id, compare_names = api_qc.add_group_compare(compare_file=compare_file, library='long',
                                                             group_id=group_id)
        qc_id = api_qc.add_qc(fq_type='PE', library='long')
        api_qc.add_qc_detail(qc_id, qc_stat_before, qc_stat_after, 'long')
        api_qc.add_qc_graph(qc_id, qc_stat_before, 'long', 'before')
        api_qc.add_qc_graph(qc_id, qc_stat_after, 'long', 'after')
        ################################################################

        ################################################################
        api_mapping = wf.api.api('whole_transcriptome.mapping')
        stat_file = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/RnaseqMapping/output/stat'
        satur_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/MapAssessment/output/saturation'
        coverage_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/MapAssessment/output/coverage'
        chrstat_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/MapAssessment/output/chr_stat'
        bamstat_dir = '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/MapAssessment/output/distribution'
        library = 'long'
        method = 'hisat'
        param_dict = json.dumps({'task_id': task_id, 'submit_location': 'mapping', 'task_type': 2}, sort_keys=True)

        api_mapping.add_mapping_stat(stat_file=stat_file, library=library, method=method)
        api_mapping.add_rpkm_table(satur_dir, params=param_dict, detail=True, library=library)
        api_mapping.add_coverage_table(coverage_dir, params=param_dict, detail=True, library=library)
        api_mapping.add_distribution_table(distribution=bamstat_dir, params=param_dict, library=library)
        api_mapping.add_chrom_distribution_table(distribution=chrstat_dir, params=param_dict, library=library)
        ################################################################

        ################################################################
        api_exp_stat = wf.api.api('whole_transcriptome.exp_stat')
        map_dict = {
            't_type': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/filter_by_express/filtered_file/trans_type.xls',
            't_count': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/ExpMake/output/count/T.reads.txt',
            'c_count': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/ExpMake/output/count/C.reads.txt',
            'g_type': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/filter_by_express/filtered_file/gene_type.xls',
            'g_count': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/ExpMake/output/count/G.reads.txt'
        }
        api_exp_stat.add_exp_stat(map_dict=map_dict, task_id=task_id, project_sn=project_sn)
        ################################################################

        ################################################################
        api_annotation = wf.api.api('whole_transcriptome.annotation')
        map_dict = {
            'all_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/all_tran2gene.txt',
            'ref_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/refannot_class/all_tran2gene.txt',
            'new_t2g': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/newannot_class/all_tran2gene.txt',

            'T_go': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/go/go_venn_tran.txt',
            'T_kegg': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/kegg/kegg_venn_tran.txt',
            'T_cog': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/cog/cog_venn_tran.txt',
            'T_nr': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/nr/nr_venn_tran.txt',
            'T_swissprot': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/swissprot/swissprot_venn_tran.txt',
            'T_pfam': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/pfam/pfam_venn_tran.txt',
            'G_go': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/go/go_venn_gene.txt',
            'G_kegg': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/kegg/kegg_venn_gene.txt',
            'G_cog': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/cog/cog_venn_gene.txt',
            'G_nr': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/nr/nr_venn_gene.txt',
            'G_swissprot': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/swissprot/swissprot_venn_gene.txt',
            'G_pfam': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Annotation/output/allannot_class/pfam/pfam_venn_gene.txt',

            'T_count': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/ExpMake/output/count/T.reads.txt',
            'G_count': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/ExpMake/output/count/G.reads.txt'
        }
        api_annotation.add_annotation_stat(map_dict=map_dict, task_id=task_id, project_sn=project_sn)
        ################################################################

        ################################################################
        api_diff_exp_stat = wf.api.api('whole_transcriptome.diff_exp_stat')
        map_dict = {
            'control': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/whole_transcriptome/MJ20190118060/largeRNA/example_control.txt',
            't_type': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/LargeGush/output/filter_by_express/filtered_file/trans_type.xls',
            't_output_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191017/Single_diff_exp_5116_1524/DiffExp/output',
            'c_output_dir': '/mnt/ilustre/users/sanger-dev/workspace/20191017/Single_diff_exp_1583_7158/DiffExp/output',
        }
        api_diff_exp_stat.add_diff_stat(map_dict=map_dict, task_id=task_id, project_sn=project_sn, library='long')
        ################################################################

        ################################################################
        api_assembly = wf.api.api('whole_transcriptome.assembly')
        map_dict = {
            'step': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Assembly/output/step.pk',
            'code': '/mnt/ilustre/users/sanger-dev/workspace/20191016/Longrna_workflow_8323_6543/Assembly/output/code.pk'
        }
        api_assembly.add_assembly(map_dict=map_dict, task_id=task_id, project_sn=project_sn)
        ################################################################


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
