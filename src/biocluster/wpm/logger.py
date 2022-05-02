# -*- coding: utf-8 -*-
# __author__ = 'guoquan'
import logging
from ..config import Config
import time
from ..core.function import hostname


class Logger(object):
    """
    日志类
    """

    def __init__(self, log_type="WPM"):
        """
        """
        self.config = Config()
        log_level = {'debug': logging.DEBUG,
                     'info': logging.INFO,
                     'warning': logging.WARNING,
                     'error': logging.ERROR,
                     'critical': logging.CRITICAL}
        self.format = self.config.LOG_FORMAT
        self.formatter = logging.Formatter(self.format, "%Y-%m-%d %H:%M:%S")
        self.level = log_level[self.config.LOG_LEVEL.lower()]
        self.streem_on = self.config.LOG_STREEM
        self._logger_date = time.strftime('%Y%m%d', time.localtime(time.time()))
        self.file_handler = None
        self.log_type = log_type
        if self.log_type == "WPM":
            self.path_dir = self.config.wpm_log_file
        else:
            self.path_dir = self.config.UPDATE_LOG
        if self.streem_on:
            self.stream_handler = logging.StreamHandler()
            self.stream_handler.setLevel(self.level)
            self.stream_handler.setFormatter(self.formatter)
        if log_type == "":
            self.logger = self.get_logger("")
        else:
            self.logger = self.get_logger("%s[%s]" % (self.log_type, hostname))

    def __del__(self):
        self.destroy()

    def destroy(self):
        if self.stream_handler:
            self.logger.removeHandler(self.stream_handler)
        if self.file_handler:
            self.logger.removeHandler(self.file_handler)

    def get_logger(self, name=""):
        """
        返回一个logger对象

        :param name: logger名字
        """
        logger = logging.getLogger(name)
        logger.propagate = 0
        if self.streem_on:
            logger.addHandler(self.stream_handler)
        logger.setLevel(self.level)
        # self._add_handler(logger)
        return logger

    # def _add_handler(self, logger):
    #     """
    #     """
    #     logger.setLevel(self.level)
    #     now_date = time.strftime('%Y%m%d', time.localtime(time.time()))
    #     if now_date != self._logger_date or self.file_handler is None:
    #         self._logger_date = now_date
    #         logger.removeHandler(self.file_handler)
    #         file_path = os.path.join(self.path_dir, now_date + ".log")
    #         self.file_handler = logging.FileHandler(file_path)
    #         self.file_handler.setLevel(self.level)
    #         self.file_handler.setFormatter(self.formatter)
    #         logger.addHandler(self.file_handler)

    def debug(self, *args,  **kargs):
            self.logger.debug(*args, **kargs)

    def info(self, *args, **kargs):
            self.logger.info(*args, **kargs)

    def warning(self, *args, **kargs):
            self.logger.warning(*args, **kargs)

    def error(self, *args, **kargs):
            self.logger.error(*args, **kargs)

    def critical(self, *args, **kargs):
            self.logger.critical(*args, **kargs)
