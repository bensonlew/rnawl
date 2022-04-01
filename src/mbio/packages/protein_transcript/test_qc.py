import random
from mbio.workflows.single import SingleWorkflow
from biocluster.wsheet import Sheet
import datetime

test_dir='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_rna_v2/demo_Mouse/rawdata'

data = {
    "id": "wenyao_testqc" + datetime.datetime.now().strftime('%H-%M-%S'),
    "type": "module",
    "name": "ref_rna_v2.hiseq_qc",
    "instant": False,
    "options": dict(
        fastq_dir=test_dir,
        fq_type='PE',
    )
}
wsheet = Sheet(data=data)
wf = SingleWorkflow(wsheet)
wf.run()
