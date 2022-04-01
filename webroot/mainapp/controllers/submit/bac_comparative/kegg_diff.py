# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from cog_diff import CogDiffAction


class KeggDiffAction(CogDiffAction):
    def __init__(self, bind_object=None):
        super(KeggDiffAction, self).__init__(bind_object)
        self.main_table = 'kegg_diff'
        self.functiontype = 'KEGG'


if __name__ == '__main__':
    '''
    test
    '''
    import os
    opt_dict = dict(
        task_id='tsg_36106',
        #group_file="/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/group2.txt",
        gfile='5dc53932e4efe6fcf0548318',
        groups='ALL',
        annotable="/mnt/ilustre/users/sanger-dev/sg-users/xieshichang/coding/cbac/test/a/kegg.annot.xls",
        test='anova',
        testtype='two.sided',
        level='Pathway',
        submit_location='kegg_diff',
        type='group',
        norm='T'
    )

    n = []
    d = []
    [[n.append(a), d.append(b)] for a, b in opt_dict.items()]
    n = ';'.join(n)
    d = ';'.join(d)

    api = 's/bac_comparative/kegg_diff'
    cmd = "python ~/biocluster/bin/webapitest.py post '{}' -n '{}' -d '{}' -c client03 -b http://192.168.12.101:9090 -fr n".format(
        api, n, d)
    # exit(cmd)
    #sg_task_info()
    os.system(cmd)

