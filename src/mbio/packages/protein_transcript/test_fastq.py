import random
from mbio.workflows.single import SingleWorkflow
from biocluster.wsheet import Sheet
import datetime
import os

test_dir='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_rna_v2/demo_Mouse/rawdata'
tmp_list = os.path.join('/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/test_qc_module', datetime.datetime.now().strftime('%H-%M-%S'))
with open(test_dir + '/list.txt', 'r') as list_r, open(tmp_list, 'w') as tmp_w:
    for line in list_r:
        tmp_w.write(test_dir + '/' + line)
data = {
    "id": "wenyao_testfastq" + datetime.datetime.now().strftime('%H-%M-%S'),
    "type": "module",
    "name": "datasplit.fastp",
    "instant": False,
    "options": dict(
        sample_path=tmp_list,
    )
}
wsheet = Sheet(data=data)
wf = SingleWorkflow(wsheet)
wf.run()
