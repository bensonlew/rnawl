# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'

import os
import argparse


def get_location(task_id, specimen_id, gene_id):
    collection = db['gene_predict']
    main_collection = collection.find_one({"task_id": task_id})
    main_id = _get_objectid(main_collection['_id'])
    detail = db['gene_predict_detail']
    detail_doc = detail.find_one({"predict_id": main_id, "specimen_id": specimen_id, "gene_id": gene_id})
    return detail_doc['location']

def add_gene_info(file, fasta, sample, split=False):
    if split == "False":
        split = False
    elif split == "True":
        split = True
    global location
    location = {}
    global location_name
    location_name = []
    with open(fasta, 'rb') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                line = line.strip().split(' ')
                gene_name = line[0][1:]
                print(gene_name)
                location_info = line[3].split('_')
                location[gene_name] = location_info[0]
                if location_info[0] not in location_name:
                    location_name.append(location_info[0])
    # if location_name[0].startswith('scaffold') or location_name[0].startswith('Scaffold'):
    #     split = False
    with open(file, 'rb') as f:
        lines = f.readlines()
        line = lines[1].split("\t")
        if line[0] in location:
            add_sample_location(lines, sample, file, split=split)
        else:
            add_sample(lines, sample, file)


def add_sample(lines, sample, out_file):
    new_file = out_file + 'new'
    with open(new_file, 'w')as o:
        ##guanqing.zou 20180917  >>>
        spline0 = lines[0].split('\t')
        s_index = 0
        for i in range(len(spline0)):
            if spline0[i] != 'Sample Name':
                s_index = i
                break
        #o.write("Sample Name\t" + lines[0])
        head = '\t'.join(spline0[s_index:])
        o.write("Sample Name\t" + head)
        ### <<<
        for line in lines[1:]:
            spline = line.split('\t')    ##guanqing.zou 20180917
            tmpline = '\t'.join(spline[s_index:])       ##guanqing.zou 20180917
            #o.write(sample + "\t" + line)
            o.write(sample + "\t" + tmpline)
    os.remove(out_file)
    os.rename(new_file, out_file)


def add_sample_location(lines, sample, out_file, split=False):
    global location
    global location_name
    file_name = {}
    dir = os.path.dirname(out_file)
    old_name = os.path.basename(out_file)
    new_file = out_file + '_new'
    with open(new_file, 'w')as o:
        #head = lines[0].split('\t', 1)
        head = lines[0].split('\t')  #guanqing.zou 20180917
        s_index = 1
        for id in range(1,len(head)):         #guanqing.zou 20180917
            if head[id] != 'Location' and head[id] != 'Sample Name':
                s_index = id
                break
        tmphead = '\t'.join(head[s_index:])
        o.write("\t".join([head[0], "Location", "Sample Name", tmphead]))
        if split:
            for i in location_name:
                l = i
                file_name[i] = open(dir + '/' + sample + '_' + l + old_name.replace('gene_', '_'), 'w')
                file_name[i].write("\t".join([head[0], "Location", "Sample Name", head[1]]))
            for line in lines[1:]:
                line = line.split('\t')
                tmpline = '\t'.join(line[s_index:])      #guanqing.zou 20180917
                w_line = "\t".join([line[0], location[line[0]], sample, tmpline])
                o.write(w_line)
                file_name[location[line[0]]].write(w_line)
            for i in location_name:
                file_name[i].close()
        else:
            for line in lines[1:]:
                line = line.split('\t')
                tmpline = '\t'.join(line[s_index:])     #guanqing.zou 20180917
                o.write("\t".join([line[0], location[line[0]], sample, tmpline]))
    if split:
        os.remove(out_file)
        os.rename(new_file, dir + '/' + sample + old_name.replace('gene_', '_whole_genome_'))
    else:
        os.remove(out_file)   ##20180912 zouguanqing
        os.rename(new_file, out_file)


def _main():
    parser = argparse.ArgumentParser(description='add gene_info in your table ')
    parser.add_argument('-i', '--file', help="file")
    parser.add_argument('-f', '--fasta', help="fasta")
    parser.add_argument('-s', '--sample', help="sample_name")
    parser.add_argument('-d', '--div', default=False, help="div result by location")
    args = parser.parse_args()
    add_gene_info(args.file, args.fasta, args.sample, args.div)


if __name__ == "__main__":
    _main()
