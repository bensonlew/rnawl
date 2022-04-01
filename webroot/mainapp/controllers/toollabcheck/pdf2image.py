# -*- coding: utf-8 -*-
# __author__ = 'scp'
import json
import re


def params_check(toollab_params):
    antype = toollab_params["sample_num"]
    if antype != "single":
        # filename = re.search("{(.*)}",str(toollab_params["pdf"])).group(1)
        filename = toollab_params["pdf"]["alias"]
        supported_formats=["tgz","tar.gz","gz","tar","rar","zip"]
        if not any(filename.endswith(s) for s in supported_formats):
            info = {'success': False,
                    'info': "%s  not supported,Please enter the correct compressed format!" % (toollab_params["filename"])}
            return json.dumps(info)
        return None
