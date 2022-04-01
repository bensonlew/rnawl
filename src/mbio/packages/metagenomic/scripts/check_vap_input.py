# -*- coding: utf-8 -*-
#/usr/bin/env python
import sys
import os

def get_splt(f):
    with open(tax_table, 'r') as tax:
        splt = tax.readline().strip('\r\n').split('\t')[2:]
    return splt 


def env_check(splt, env_file):
    with open(env_file, 'r') as ef:
        envs = ef.readline().strip('\r\n').split('\t')
        env_value_dict = {}
        for l in ef:
            line = l.strip('\r\n').split('\t')
            if line[0] in splt:
                for i in range(1, len(envs)):
                    if envs[i] not in env_value_dict:
                        env_value_dict[envs[i]] = []
                    env_value_dict[envs[i]].append(line[i])
    errs = []
    for k,v in env_value_dict.items():
        if len(set(v)) == 1:
            errs.append(k)
    if os.path.exists('check_vap_input.txt'):
        os.remove('check_vap_input.txt')
    if errs:
        with open('check_vap_input.txt', 'w') as w:
            w.write(str(errs))


if __name__ == '__main__':
    tax_table = sys.argv[1]
    env_file = sys.argv[2]
    
    splt = get_splt(tax_table)
    env_check(splt, env_file)
    
