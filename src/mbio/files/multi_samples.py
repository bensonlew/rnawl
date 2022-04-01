# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import os
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
import json
from biocluster.api.file.remote import RemoteFileManager
import shutil
import re


class MultiSamplesFile(Directory):
    def __init__(self):
        super(MultiSamplesFile, self).__init__()
        self.self_managed = True
        self._sample_config = []
        self._sample_files = {}

    def set_path(self, path):
        try:
            samples = json.loads(path)
        except Exception, e:
            raise FileError("输入格式不正确: %s, %s" % (e, path))
        if not isinstance(samples, list):
            raise FileError("输入样品必须为list")
        self._sample_config = samples
        self._check_reduplicate()
        self._download_samples()
        current_path = os.path.join(self.option.bind_obj.work_dir, "remote_input", self.option.name)
        super(MultiSamplesFile, self).set_path(current_path)
        for sample in self._sample_config:
            sample_path = os.path.join(current_path, sample["sample_name"])
            if not os.path.exists(sample_path):
                os.makedirs(sample_path)
            self.add_files("seq_sample", sample)

    @property
    def logger(self):
        return self.option.bind_obj.logger

    def _download_samples(self):
        current_path = os.path.join(self.option.bind_obj.work_dir, "remote_input", self.option.name)
        if not os.path.exists(current_path):
            os.makedirs(current_path)
        for name, path in self._sample_files.items():
            remote_file = RemoteFileManager(path)
            if remote_file.type != "local":
                self.option.bind_obj.logger.info("发现远程文件%s,开始复制..." % path)
                remote_file.download(os.path.join(current_path, os.path.dirname(name)))
                file_path = remote_file.local_path
                target = os.path.join(current_path, name)
                if file_path == target:
                    pass
                else:
                    if os.path.exists(target):
                        if os.path.islink(target):
                            os.remove(target)
                        elif os.path.isdir(target):
                            shutil.rmtree(target)
                        else:
                            os.remove(target)
                    if os.path.isfile(file_path):
                        os.link(file_path, target)
                    else:
                        os.symlink(file_path, target)
        return current_path

    def _check_reduplicate(self):
        paths = {}
        current_path = os.path.join(self.option.bind_obj.work_dir, "remote_input", self.option.name)
        if not os.path.exists(current_path):
            os.makedirs(current_path)
        for sample in self._sample_config:
            if not isinstance(sample, dict):
                raise FileError("格式错误!")
            if "sample_name" not in sample.keys() or "library" not in sample.keys():
                raise FileError("必须包含样本名和文库信息!")
            if not sample["sample_name"] or not sample["library"]:
                raise FileError("样本名和文库名不能为空!")
            if "sequence_file_ids" in sample.keys():
                if isinstance(sample["sequence_file_ids"], list):
                    h5_count = 0
                    for s in sample["sequence_file_ids"]:
                        file_name = os.path.join(sample["sample_name"], s["alias"])
                        if file_name in self._sample_files.keys():
                            raise FileError("文件名称重复： %s" % file_name)
                        else:
                            self._sample_files[file_name] = s["file_path"]
                        if s["file_path"] in paths.keys():
                            raise FileError("文件路径重复： %s" % s["file_path"])
                        else:
                            paths[s["file_path"]] = 1
                        if re.search('\.bax\.h5$', s["alias"]):
                            h5_count += 1
                    if h5_count > 0 and (h5_count !=3 or len(sample["sequence_file_ids"]) != 3):
                        raise FileError("样本%s 文库 %s 文件输入不正确: ：bax.h5文件一次只能包含3个!" %
                                        (sample["sample_name"], sample["library"]))
                else:
                    raise FileError("样本%s 文库 %s 格式不正确: sequence_file_ids 必须为list数组!" %
                                    (sample["sample_name"], sample["library"]))
            elif "sequence_r1_file_ids" in sample.keys() and "sequence_r2_file_ids" in sample.keys()\
                    and isinstance(sample["sequence_r1_file_ids"], dict)\
                    and isinstance(sample["sequence_r2_file_ids"], dict):
                file_name1 = os.path.join(sample["sample_name"], sample["sequence_r1_file_ids"]["alias"])
                file_name2 = os.path.join(sample["sample_name"], sample["sequence_r2_file_ids"]["alias"])
                if file_name1 in self._sample_files.keys():
                    raise FileError("文件名称重复： %s" % file_name1)
                else:
                    self._sample_files[file_name1] = sample["sequence_r1_file_ids"]["file_path"]
                if file_name2 in self._sample_files.keys():
                    raise FileError("文件名称重复： %s" % file_name2)
                else:
                    self._sample_files[file_name2] = sample["sequence_r2_file_ids"]["file_path"]
                if sample["sequence_r1_file_ids"]["file_path"] in paths.keys():
                    raise FileError("文件路径重复： %s" % sample["sequence_r1_file_ids"]["file_path"])
                else:
                    paths[sample["sequence_r1_file_ids"]["file_path"]] = 1
                if sample["sequence_r2_file_ids"]["file_path"] in paths.keys():
                    raise FileError("文件路径重复： %s" % sample["sequence_r2_file_ids"]["file_path"])
                else:
                    paths[sample["sequence_r2_file_ids"]["file_path"]] = 1
            else:
                raise FileError("样本%s 文库 %s 文件信息不正确！" % (sample["sample_name"], sample["library"]))

    @property
    def samples(self):
        return self.files