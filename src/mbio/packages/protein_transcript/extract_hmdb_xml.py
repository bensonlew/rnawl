#!/usr/bin/env python
#coding:utf-8
"""
@time    : 2019/3/7 11:04
@file    : extract_hmdb_xml.py
@author  : yitong.feng
@contact: yitong.feng@majorbio.com
"""

import xml.etree.ElementTree as ET
import sys

reload(sys)
sys.setdefaultencoding("utf-8")

def extract_xml(xml, out):
    def get_children_info(block):
        children_info = list()
        sacc = block.getchildren()
        if sacc:
            children_info = [c.text for c in sacc]
        return children_info
    def deal_ontology(ontology):
        term = ''
        for onto in ontology:
            if onto.tag.endswith(('root', 'descendants', 'descendant')):
                tmp = onto
                deal_ontology(tmp)
            if onto.tag.endswith('term'):
                term = onto.text
            if onto.tag.endswith('level'):
                level = onto.text
                if ontology_dict[level] == '_':
                    ontology_dict[level] = term
                else:
                    ontology_dict[level] += ';' + term

    tree = ET.parse(xml)
    root = tree.getroot()
    ow = open(out, 'w')
    ow.write('accession\tsecondary_accessions\tname\tchemical_formula\taverage_molecular_weight\tcas_registry_number\tcellular_locations\tbiospecimen_locations\ttissue_locations\tkegg_id\tontology_level1\tontology_level2\tontology_level3\tontology_level4\tontology_level5\n')
    accession = secondary_accessions = name = chemical_formula = average_molecular_weight = cas_registry_number = cellular_locations = biospecimen_locations = tissue_locations = kegg_id = '_'
    # ontology = {'1':'', '2':'', '3':'', '4':'', '5':''}
    acc2locations = dict()
    for child in root:
        for child1 in child:
            ontology_dict = {'1': '_', '2': '_', '3': '_', '4': '_', '5': '_'}
            if child1.tag.endswith('accession') and 'secondary_' not in child1.tag:
                accession = child1.text if child1.text else '_'
                if accession != '_' and accession not in acc2locations:
                    acc2locations[accession] = {'cellular_locations' : set(), 'biospecimen_locations' : set(),'tissue_locations' : set(),}
                if accession != '_':
                    for key in acc2locations[accession]:
                        if acc2locations[accession][key] == '_':
                            acc2locations[accession][key] = set()
                        elif type(acc2locations[accession][key]) == str:
                            acc2locations[accession][key] = set(acc2locations[accession][key].split(';'))
            if child1.tag.endswith('secondary_accessions'):
                secondary_accessions = get_children_info(child1)
                if secondary_accessions:
                    secondary_accessions = ';'.join(secondary_accessions)
                else:
                    secondary_accessions = '_'
            if child1.tag.endswith('name'):
                name = child1.text if child1.text else '_'
            if child1.tag.endswith('chemical_formula'):
                chemical_formula = child1.text if child1.text else '_'
            if child1.tag.endswith('average_molecular_weight'):
                average_molecular_weight = child1.text if child1.text else '_'
            if child1.tag.endswith('cas_registry_number'):
                cas_registry_number = child1.text if child1.text else '_'
            if child1.tag.endswith('cellular_locations'):
                cellular_locations = get_children_info(child1)
                if accession != '_':
                    # print(acc2locations[accession]['cellular_locations'])
                    acc2locations[accession]['cellular_locations'].update(set(cellular_locations))
            if child1.tag.endswith('biospecimen_locations'):
                biospecimen_locations = get_children_info(child1)
                if accession != '_':
                    acc2locations[accession]['biospecimen_locations'].update(set(biospecimen_locations))
            if child1.tag.endswith('tissue_locations'):
                tissue_locations = get_children_info(child1)
                if accession != '_':
                    acc2locations[accession]['tissue_locations'].update(set(tissue_locations))
            if child1.tag.endswith('biological_properties'):
                for child2 in child1:
                    if child2.tag.endswith('cellular_locations'):
                        cellular_locations = get_children_info(child2)
                        if accession != '_':
                            if acc2locations[accession]['cellular_locations'] == '_':
                                acc2locations[accession]['cellular_locations'] = set()
                            # print(acc2locations[accession]['cellular_locations'])
                            acc2locations[accession]['cellular_locations'].update(set(cellular_locations))
                    if child2.tag.endswith('biospecimen_locations'):
                        biospecimen_locations = get_children_info(child2)
                        if accession != '_':
                            if acc2locations[accession]['biospecimen_locations'] == '_':
                                acc2locations[accession]['biospecimen_locations'] = set()
                            acc2locations[accession]['biospecimen_locations'].update(set(biospecimen_locations))
                    if child2.tag.endswith('tissue_locations'):
                        tissue_locations = get_children_info(child2)
                        if accession != '_':
                            if acc2locations[accession]['tissue_locations'] == '_':
                                acc2locations[accession]['tissue_locations'] = set()
                            acc2locations[accession]['tissue_locations'].update(set(tissue_locations))
            if child1.tag.endswith('kegg_id'):
                kegg_id = child1.text if child1.text else '_'
            if child1.tag.endswith('ontology'):
                deal_ontology(child1)
            if_w = False
            for i in ontology_dict:
                if ontology_dict[i] != '_':
                    if_w = True
                    break
            if if_w:
                # print(acc2locations[accession])
                if acc2locations[accession]['cellular_locations']:
                    cellular_locations = ';'.join(acc2locations[accession]['cellular_locations'])
                else:
                    cellular_locations = '_'
                if acc2locations[accession]['biospecimen_locations']:
                    biospecimen_locations = ';'.join(acc2locations[accession]['biospecimen_locations'])
                else:
                    biospecimen_locations = '_'
                if acc2locations[accession]['tissue_locations']:
                    tissue_locations = ';'.join(acc2locations[accession]['tissue_locations'])
                else:
                    tissue_locations = '_'
                # print(accession,secondary_accessions, name,chemical_formula, average_molecular_weight,cas_registry_number,acc2locations[accession]['cellular_locations'],acc2locations[accession]['biospecimen_locations'],acc2locations[accession]['tissue_locations'],kegg_id,ontology_dict['1'], ontology_dict['2'], ontology_dict['3'],ontology_dict['4'], ontology_dict['5'] )

                ow.write(accession+ '\t' + secondary_accessions+ '\t' + name+ '\t' + chemical_formula+ '\t' + average_molecular_weight+ '\t' + cas_registry_number+ '\t' + cellular_locations + '\t' + biospecimen_locations + '\t' + tissue_locations + '\t' + kegg_id+ '\t' + ontology_dict['1'] + '\t' + ontology_dict['2'] + '\t' + ontology_dict['3'] + '\t' + ontology_dict['4'] + '\t' + ontology_dict['5'] + '\n')
    # print(acc2locations)
            # ow.write(
            #     accession + '\t' + secondary_accessions + '\t' + name + '\t' + chemical_formula + '\t' + average_molecular_weight + '\t' + cas_registry_number + '\t' + cellular_locations + '\t' + biospecimen_locations + '\t' + tissue_locations + '\t' + kegg_id + '\t' +
            #     ontology_dict['1'] + '\t' + ontology_dict['2'] + '\t' + ontology_dict['3'] + '\t' + ontology_dict[
            #         '4'] + '\t' + ontology_dict['5'] + '\n')
    ow.close()
    with open(out, 'r') as o:
        header = o.readline()
        w_info = o.readlines()
    with open(out, 'w') as ow:
        ow.write(header)
        for line in w_info:
            line = line.strip().split('\t')
            if line[0] in acc2locations:
                line[6] = ';'.join(acc2locations[line[0]]['cellular_locations']) if acc2locations[line[0]]['cellular_locations'] else '_'
                line[7] = ';'.join(acc2locations[line[0]]['biospecimen_locations']) if acc2locations[line[0]]['biospecimen_locations'] else '_'
                line[8] = ';'.join(acc2locations[line[0]]['tissue_locations']) if acc2locations[line[0]]['tissue_locations'] else '_'
            ow.write('\t'.join(line) + '\n')



if __name__ == '__main__':

    xml = sys.argv[1]
    out = sys.argv[2]
    extract_xml(xml, out)