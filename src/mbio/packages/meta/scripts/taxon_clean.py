# -*- coding: utf-8 -*-
# __author__ = 'sheng.he'
# 脚本用于规范化注释物种信息格式化，两个参数，第一个是旧的物种注释文件，必须得有 8 层物种， 第二个为新的物种注释文件
# 脚本为了统一未分类物种信息，统一将 'norank', 'uncultured', "incertae_sedis" 改为 norank
# 统一将 "unidentified", "Unclassified", "Unknown" 改为 unclassified
import sys
import re
from collections import defaultdict


check_dict = {}


class TaxonTree(object):
    def __init__(self, name, parent=None):
        self._parent = parent
        self._name = name
        self.new_name = name
        self._children = {}
        self.num = 0
        self.seq_id = []
        pass

    def parse(self, taxs):
        self.num += 1
        if len(taxs) == 2:
            if taxs[0] not in self._children:
                self._children[taxs[0]] = TaxonTree(taxs[0], parent=self).parse_id(taxs[1])
            else:
                self._children[taxs[0]].parse_id(taxs[1])
        else:
            if taxs[0] not in self._children:
                self._children[taxs[0]] = TaxonTree(taxs[0], parent=self).parse(taxs[1:])
            else:
                self._children[taxs[0]].parse(taxs[1:])
        return self

    def parse_id(self, seqid):
        self.num += 1
        self.seq_id.append(seqid)
        return self

    def get_good_parent(self):
        if not self._parent._name:
            raise Exception('存在除了特殊非物种类型外，所有父节点相同情况')
        search = re.search(r'^[dkpcofgs]__(norank|uncultured|Incertae_Sedis|unidentified|Unclassified|Unknown)$', self._parent._name, re.IGNORECASE)
        if search:
            this_tax = re.split(r'__', self._name, maxsplit=1)[1]
            if this_tax in ['norank', 'unclassified'] and search.groups()[0] == this_tax:
                return self._parent.get_good_parent()
            else:
                return self._parent
        else:
            return self._parent

    def all_parent_name(self, taxs):
        if self._parent:
            taxs.insert(0, self.new_name)
            self._parent.all_parent_name(taxs)

    def new_tax(self):
        self.new_tax_str = []
        self.all_parent_name(self.new_tax_str)
        return ';'.join(self.new_tax_str)


def parse_file(fp):
    origin_tree = TaxonTree(0)
    taxon_f = open(fp)
    for line in taxon_f:
        line_sp = line.split('\t')
        taxs = line_sp[1].strip().split(';')
        taxs = correct_tax_unknown(taxs)
        taxs.append(line_sp[0])
        origin_tree.parse(taxs)
    return origin_tree


def correct_tax_unknown(taxs):
    for index, tax in enumerate(taxs):
        temp = re.split('__', tax, maxsplit=1)
        if temp[1].lower() in ['norank', 'uncultured', "incertae_sedis"]:
            taxs[index] = temp[0] + '__norank'
        elif temp[1].lower() in ["unidentified", "unclassified", "unknown"]:
            taxs[index] = temp[0] + '__unclassified'
        else:
            pass
    return taxs


def get_level(tax_tree, level=1):
    level_tax = defaultdict(list)
    if level == 1:
        return dict([(i, [tax_tree._children[i]]) for i in tax_tree._children])
        return tax_tree._children.keys()
    for i1 in tax_tree._children.values():
        if level == 2:
            for i in i1._children:
                level_tax[i].append(i1._children[i])
            continue
        for i2 in i1._children.values():
            if level == 3:
                for i in i2._children:
                    level_tax[i].append(i2._children[i])
                continue
            for i3 in i2._children.values():
                if level == 4:
                    for i in i3._children:
                        level_tax[i].append(i3._children[i])
                    continue
                for i4 in i3._children.values():
                    if level == 5:
                        for i in i4._children:
                            level_tax[i].append(i4._children[i])
                        continue
                    for i5 in i4._children.values():
                        if level == 6:
                            for i in i5._children:
                                level_tax[i].append(i5._children[i])
                            continue
                        for i6 in i5._children.values():
                            if level == 7:
                                for i in i6._children:
                                    level_tax[i].append(i6._children[i])
                                continue
                            for i7 in i6._children.values():
                                for i in i7._children:
                                    level_tax[i].append(i7._children[i])
    return level_tax


def get_good_parent(child, origin):
    try:
        good_parent = child.get_good_parent()
    except Exception as e:
        print "ERROR: " + str(e)
        print origin._name
        return
    origin.new_name = origin.new_name + "_" + good_parent._name
    if origin.new_name in check_dict:
        get_good_parent(good_parent, origin)
    else:
        check_dict[origin.new_name] = 0
    return good_parent


def change_new_name(tree, level):
    for level, children in get_level(tree, level).iteritems():
        if len(children) > 1:
            for child in children:
                get_good_parent(child, child)


def change_all_name(tree):
    change_new_name(tree, 1)
    change_new_name(tree, 2)
    change_new_name(tree, 3)
    change_new_name(tree, 4)
    change_new_name(tree, 5)
    change_new_name(tree, 6)
    change_new_name(tree, 7)
    change_new_name(tree, 8)


def write_new(tree, fp):
    new_tax_file = open(fp, 'w')
    for one, ids in get_level(tree, 8).iteritems():
        for i in ids:
            new_tax = i.new_tax()
            for one_id in i.seq_id:
                new_tax_file.write(one_id + "\t" + new_tax + '\n')

def add_class(taxon):
    """
    为taxon注释信息添加层级关系
    :param taxon:
    :return:
    """
    out_taxon = "taxon_file.xls"
    with open(taxon, 'r') as f, open(out_taxon, "w") as w:
        for line in f:
            line = line.strip().split("\t")
            taxline = line[1]
            taxs = re.split(r';\s*', taxline)
            if len(taxs) != 8:
                pre_listt = []
                for tax in taxs:
                    tax_pre = re.split("__", tax, 1)[0]
                    if tax_pre not in pre_listt:
                        pre_listt.append(tax_pre)
                origin_list = ["d", "k", "p","c", "o", "f", "g", "s"]
                pre_listt = set(pre_listt)
                add_list = list((set(origin_list)).difference(pre_listt))
                for add in add_list:
                    add_index = add_list.index(add) + 1
                    add_name = add + "__norank"
                    taxs.insert(add_index, add_name)
                line[1] = ";".join(taxs)
                new_line = "\t".join(line)
            else:
                new_line = "\t".join(line)
            w.write('{}\n'.format(new_line))
    return out_taxon


if __name__ == "__main__":
    args = sys.argv[1:]
    taxon_file = args[0]
    new_file = args[1]
    new_taxon_file = add_class(taxon_file)
    mytree = parse_file(new_taxon_file)
    change_all_name(mytree)
    write_new(mytree, new_file)
