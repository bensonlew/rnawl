# -*- coding: utf-8 -*-
# author fengyitong 2019-01

import sys
sys.path.insert(0, '/mnt/ilustre/users/yitong.feng/scripts/run_commands')

from controller import Controller
import os
import json
import time

class shiyunxing(Controller):
    def __init__(self, test_fa):
        self.fa = test_fa

    def trans_fa(self):
        params = dict(
            cmd = "awk '{if (/>/){printf '\n%s\t',$0}else{printf '%s',$0}}' " + self.fa  + "|sed '1d;s/>//g' > hhh",

        )

