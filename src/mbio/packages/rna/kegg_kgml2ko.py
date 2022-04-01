# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

import re
import os
import sys
import xml.etree.ElementTree as ET
import lxml.html

class KeggKgml(object):
    def __init__(self, kgml_path):
        self.path = kgml_path

    def get_ko(self):
        map_name = os.path.basename(self.path).split(".")[0]
        ko_name = map_name.replace("map", "ko")

        xml = ET.parse(self.path)
        root = xml.getroot()
        all_entrys =  root.findall('entry')
        ko_set = set()
        for entry in all_entrys:
            kos = entry.get("name")
            for ko in kos.split(" "):
                if ko in ko_set:
                    pass
                else:
                    print "path:{}\t{}".format(map_name, ko)
                    print "path:{}\t{}".format(ko_name, ko)
                    ko_set.add(ko)



if __name__ == "__main__":
    kgml_path = sys.argv[1]
    a = KeggKgml(kgml_path)
    a.get_ko()
