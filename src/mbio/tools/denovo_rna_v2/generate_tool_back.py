# coding=utf-8
import os
from jinja2 import Environment
__author__ = 'gdq'
# tool_template.jinja2 is the tool template. DO NOT CHANGE IT unless necessary.

# 定制以下变量的值，直接运行当前的脚本即可获得定制的tool，然后去修改需要完善的地方，完成tool开发。
tool_name = 'Corr'  # 请定制
raw_tool_name = 'corrxxxx'  # 请定制
tool_description = 'Correlation calculation'  # 请定制

# tool_parent_dir is a directory name in /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/
tool_parent_dir = 'denovo_rna_v2'  # 请定制

# called_script is expected located at: /bioinfo/rna/scripts/
# called_script is expected to be a python script which may wrap outer software written in C/java.
# we expected that all argument names are designed to be same with this tool.
called_script = 'called_script_name'  # 请定制， 最好的情况是该脚本为你对外部程序的python包装

# option list. Very important
type_list = ['int', 'float', 'string', 'infile', 'outfile', ]
option_attributes = ['name', 'type', 'default', 'format', ]
option_list = [  # 请定制，该部分参数名与called_script所需参数名保持一致会减少你的工作量。
    dict(name='express_matrix', type=type_list[3], format="denovo_rna_v2.express_matrix", ),
    dict(name='group_table', type=type_list[3], format="denovo_rna_v2.group_table", ),
    dict(name='method', type=type_list[2], default='pearson', ),
    # dict(name='arg4', type=type_list[0], default=0, ),
]


# ----rendering-------do not change the following codes unless necessary---------
env = Environment()
# env.trim_blocks = True
tool_template = env.from_string(open('tool_template.jinja2').read())
result = tool_template.render(tool_name=tool_name,
                              raw_tool_name=raw_tool_name,
                              tool_parent_dir=tool_parent_dir,
                              called_script=called_script,
                              option_list=option_list,
                              tool_description=tool_description,
                              )
tools_dir = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(tools_dir, raw_tool_name+'.py'), 'w') as f:
    f.write(result)
