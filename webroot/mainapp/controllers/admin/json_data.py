# -*- coding: utf-8 -*-
# __author__ = 'guoquan'  
from .admin import Admin
import json


class JsonDataAction(Admin):
    def __init__(self):
        super(JsonDataAction, self).__init__()

    def POST(self):
        key = "json"
        if self.params.type in ["api", "response"]:
            model = self.get_model("apilog")
            if self.params.type == "api":
                key = "data"
            else:
                key = "response"
        else:
            model = self.get_model("workflow")

        data = model.get_by_id(self.params.id, key)
        if data:
            return json.dumps(json.loads(getattr(data, key)), indent=4)
