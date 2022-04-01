# -*- coding: utf-8 -*-
# __author__ = 'shijin'

import os
import re
import sys

task_id_list = []

def import_report_aferter_end():
	with open("import_report_aferter_end","r") as r:
		for line in r:
			line = line.split("\"")
			task_id_list.append(line[1]) 
		
		
target_dir = '/mnt/ilustre/users/sanger/workspace/'
def find(target_dir,task_id_list):
	for oneday in os.listdir(target_dir):
		if oneday == 'tmp':
			continue
		# print('Day:', oneday)
		day_path = os.path.join(target_dir, oneday)
		for one_task in os.listdir(day_path):
			task_info = one_task.split('_')
			if task_info[0] == 'SampleExtract':
				task_id = '_'.join(task_info[1:])
				# print('TaskID: ', task_id)
				if task_id in task_id_list:
					dir_path = os.path.join(day_path,one_task)
					sample_extract_path = os.path.join(dir_path,"SampleExtract")
					info_path = os.path.join(sample_extract_path,"info.txt")
					print task_id + "\t" + dir_path


# find(target_dir,[sys.argv[1]])
dir_list = []
def findall(target_dir):
	for oneday in os.listdir(target_dir):
		if oneday == 'tmp':
			continue
		# print('Day:', oneday)
		day_path = os.path.join(target_dir, oneday)
		for one_task in os.listdir(day_path):
			task_info = one_task.split('_')
			if task_info[0] == 'SampleExtract':
				task_id = '_'.join(task_info[1:])
				dir_path = os.path.join(day_path,one_task)
				sample_extract_path = os.path.join(dir_path,"SampleExtract")
				info_path = os.path.join(sample_extract_path,"info.txt")
				if not os.path.exists(info_path) and len(os.listdir(dir_path)) != 3:
					dir_list.append(dir_path)


findall(target_dir)
for dir in sorted(dir_list):
	print dir + "\n"