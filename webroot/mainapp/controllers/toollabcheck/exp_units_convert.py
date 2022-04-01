# -*- coding: utf-8 -*-
# __author__ = 'scp'
import json


def params_check(toollab_params):
    convert = toollab_params["from_unit"].lower() + "2" + toollab_params["to_unit"].lower()
    if convert not in ['count2cpm', 'fpkm2tpm', 'count2tpm', 'count2fpkm', 'cpm2fpkm', 'cpm2tpm', 'fpkm2cpm']:
        info = {'success': False,
                'info': "%s to %s not supported!" % (toollab_params["from_unit"], toollab_params["to_unit"])}
        return json.dumps(info)
    if convert not in ['count2cpm', 'fpkm2tpm']:
        if not toollab_params["gene_length"] or toollab_params["gene_length"] == "null":
            info = {'success': False,
                    'info': "%s to %s should provie gene length file!" % (toollab_params["from_unit"], toollab_params["to_unit"])}
            return json.dumps(info)
    return None
