# -*- coding: utf-8 -*-

import os
import sys
from concurrent.futures import ThreadPoolExecutor

if len(sys.argv) != 3:
    exit('USAGE: %s stat_dir output' %sys.argv[0])

count_dir = sys.argv[1]
if not os.path.exists(count_dir):
    exit('the %s not exists'% count_dir)
person_dict = dict()

def count_file(user):
    if user not in person_dict:
        person_dict[user] = dict()
        person_dict[user]['files_number'] = 0
        person_dict[user]['files_sizes_sum'] = 0
    wf = open(user+'_file_counts.txt', 'w')
    for root, dir, files in os.walk(os.path.join(count_dir, user)):
        for file in files:
            try:
                wf.write(os.path.join(root, file) + '\t' + str(os.path.getsize(os.path.join(root, file))) + '\n')
                person_dict[user]['files_number'] += 1
                person_dict[user]['files_sizes_sum'] += os.path.getsize(os.path.join(root, file))
            except:
                wf.write(os.path.join(root, file) + '\t' + '___' + '\n')
    else:
        wf.write('\n')
        wf.write(user + '\t' + str(person_dict[user]['files_number']) + '\t' + str(person_dict[user]['files_sizes_sum']) + '\n')
        wf.close()

users = [x for x in os.listdir(count_dir) if os.path.isdir(os.path.join(count_dir, x))]
# count_file(users[0])
# print(person_dict)
# count_file(users[1])
# print(person_dict)
# count_file(users[3])
# print(person_dict)
with ThreadPoolExecutor(35) as pool:
    pool.map(count_file, users)
print(person_dict.keys())
with open(sys.argv[2], 'w') as out:
    out.write('user_name\tfiles_count_num\tfiles_sum_sizes\n')
    for user in person_dict:
        print(user)
        out.write(user + '\t' + str(person_dict[user]['files_number']) + '\t' + str(person_dict[user]['files_sizes_sum']) + '\n')
