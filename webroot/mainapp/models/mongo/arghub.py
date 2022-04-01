# -*- coding: utf-8 -*-

import types
import json
import random, os, re
from bson import SON
from bson.objectid import ObjectId
from biocluster.config import Config
from mainapp.models.workflow import Workflow
from mainapp.models.mongo.meta import Meta
from biocluster.file import exists,list_dir


class Arghub(Meta):
    """
    获取mongo字段
    """
    def __init__(self, bind_object=None):
        self._bind_object = bind_object
        super(Arghub, self).__init__(self._bind_object)
        self._project_type = 'arghub'

