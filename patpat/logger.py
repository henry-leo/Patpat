"""Logger module 日志记录器模块

The class CoreLogger is used to generate logs and search temporary files,
and is recommended to be used with hub.py and viewer.py.
class CoreLogger用于生成日志和搜索临时文件，建议配合hub.py和viewer.py使用。

    Typical usage example:
    典型用法示例：

    import patpat_env.logger as logger

    l = logger.CoreLogger(uuid)
    l.set_core
    l.set_tmp
"""
import logging
import sys


class CoreLogger:
    """Core Log Module 核心日志模块

    Used to generate logs and search temporary files, recommended to use with hub.py and viewer.py.
    用于生成日志和搜索临时文件，
    """

    def __init__(self, uuid):
        """Initialization Class 初始化类

        Args:
            uuid: The uuid created by class MapperHub, or enter the uuid of each search task by yourself.
                类MapperHub创建的uuid，或自行输入曾经搜索任务的uuid。
        """
        self.id = uuid
        # Create the Logger object for log generation. 创建用于生成日志的Logger对象。
        self.logger = logging.getLogger('core')
        self.fmt = logging.Formatter('%(asctime)s\t%(levelname)s\t%(module)s\t%(funcName)s\t%(message)s')
        self.logger.setLevel(logging.DEBUG)

        # Create the Logger object used to generate the temporary file. 创建用于生成临时的Logger对象
        self.tmp = logging.getLogger('tmp')
        self.tmp.setLevel(logging.INFO)

    def set_core(self):
        """Used to set up two log processors, one sends logs to sys.stderr and the other sends them to the uuid.log file.
           用于设置两个日志处理器，一个将日志发送到sys.stderr，另一个将日志发送到uuid.log文件。
        """
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setFormatter(logging.Formatter('%(message)s'))
        sh.setLevel(logging.INFO)
        fh = logging.FileHandler(filename=f'patpat_env/logs/{self.id}.log', encoding='utf-8')
        fh.setFormatter(self.fmt)

        self.logger.addHandler(sh)
        self.logger.addHandler(fh)

    def set_tmp(self):
        """Processor used to set the log to be sent to the temporary file. 用于设置发送日志至临时文件的处理器。"""
        fh = logging.FileHandler(filename=f'patpat_env/tmp/{self.id}.log', encoding='utf-8')
        fh.setFormatter(logging.Formatter('%(message)s'))
        self.tmp.addHandler(fh)


class TestLogger:
    def __init__(self):
        self.logger = logging.getLogger('core')
        self.logger.setLevel(logging.DEBUG)

        self.sh = logging.StreamHandler(stream=sys.stdout)
        self.sh.setFormatter(logging.Formatter('%(message)s'))

        self.logger.addHandler(self.sh)
