#!/mnt/ilustre/users/sanger/app/Python/bin/python
from mbio.workflows.single import SingleWorkflow
from biocluster.wsheet import Sheet
import time
import sys

b=time.strftime("%Y%m%d",time.localtime(time.time()))
datas = {
  "id": "nr_anno_tool_20180418",
  "name": "wgs.split_fastawgs",
  "type": "tool",
  "options":{
    "fasta": "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/1.Var/0.tool/1.tools.test/2.denovo.fasta/Lands.denovo.scafSeq",
   }
}

wsheet = Sheet(data=datas)
wf = SingleWorkflow(wsheet)
wf.run()