# -*- coding: utf-8 -*-
from meta_genomic_v3 import MetaGenomicWorkflow as mg_v3
from meta_genomic_old import MetaGenomicWorkflow as mg_old


def wfwrap(func):
    def _wf(wsheet_object):
        options = wsheet_object.options()
        if "pipeline" in options:
            wf = mg_v3(wsheet_object)
        else:
            wf = mg_old(wsheet_object)
        return wf
    return _wf


@wfwrap
def MetaGenomicWorkflow(wsheet_object):
    pass
