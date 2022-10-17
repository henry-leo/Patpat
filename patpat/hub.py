"""This module provides users with a one-stop shop for using patpat_env classes. 本模块为用户提供一站式使用patpat功能的类。

This module is the core module of patpat_env and allows users to aggregate other classes of patpat_env through the classes
in this module. init function provides the runtime environment required for patpat_env. queryHub class provides querying
of protein metadata and generating peptides to be retrieved, while MapperHub class provides breakpoints for merging
and retrieving multiple proteomic databases . Furthermore, both classes are pluggable, so it is easy to insert
self-built methods or classes into them, as long as the developer adheres to the interface design.
本模块是patpat的核心模块，用户可以通过本模块中的类将patpat的其他类聚合在一起使用。
QueryHub类提供对蛋白质元数据的查询和需要检索肽段的生成功能，MapperHub类则提供多蛋白质组数据库合并检索和检索断点续传的功能。
其次，这两个类都是可拔插的，因此只要开发者遵守接口设计，很容易将自建的方法或类插入其中。

    Typical usage example:
    典型用法示例：

    import patpat_env.hub as hub
    import patpat_env.mapper as mapper

    hub.init()

    identifier_ = 'P05067'

    q = hub.QueryHub()
    q.identifier = identifier_
    q.simple_query()

    conf_ = q.get_query_config()
    mappers_ = [mapper.PrideMapper(), mapper.IProXMapper()]

    m = hub.MapperHub(config=conf_,
                      mappers=mappers_,
                      task=None
                      # task=[your task's uuid]
                      )

    m.mapping()
    result_ = m.export()
"""
import json
import os
import re
import sys
import time
import logging
import uuid

from . import utility
from . import querier
from . import mapper
from . import logger


class QueryHub(object):
    """This class is used to build the configuration information that needs to be retrieved.
       本类用于建立需要检索的配置信息。

    This class uses the UniProt database by default to query protein metadata and generate the peptides to
    be retrieved ,using a local library build method, the process is encapsulated in the simple_query function.
    The developer encourages users to experiment with the pluggable QueryHub class to match specific search
    requirements or to insert new peptide generation methods.
    本类默认使用UniProt数据库查询蛋白质元数据，使用本地建库的方法生成需要检索肽段，这个过程被封装在simple_query函数中。
    开发者鼓励用户尝试使用可拔插的QueryHub类以符合特殊的搜索要求，或插入新的肽段生成方法。
    """

    def __init__(self):
        self.identifier: str = ''  # protein sequence
        self.peptides: list = []  # Filtered peptides for output to external
        self.organism: dict = dict()  # dict, Organism to which the protein belongs.

        self.fasta: dict = dict()  # Structuring protein fasta into a dictionary.
        self.source: str = ''  # Local proteome file address

        self.digestion_params: list = []  # Parameters of protein in-silico enzymatic

        # pluggable
        self._protein_querier = None  # class ProteinQuerier
        self._peptide_querier = None  # class PeptideQuerier

    @property
    def protein_querier(self):
        return self._protein_querier

    @protein_querier.setter
    def protein_querier(self, protein_querier):
        """Set class ProteinQuerier"""
        self._protein_querier = protein_querier

    @property
    def peptide_querier(self):
        return self._peptide_querier

    @peptide_querier.setter
    def peptide_querier(self, peptide_querier):
        """Set class PeptideQuerier"""
        self._peptide_querier = peptide_querier

    def protein_query(self):
        """Call the class ProteinQuerier in itself. 在自身中调用ProteinQuerier类。

        Call the query method and get_properties method of class ProteinQuerier to query protein metadata and export it.
        调用ProteinQuerier类的query方法和get_properties方法，查询蛋白质元数据并导出。

        Returns:
            self.fasta['sequence']: str
            self.organism: dict

        Updates:
            self.protein_querier
            self.identifier
            self.organism
            self.fasta
        """
        self.protein_querier.query()
        self.identifier, self.organism, self.fasta = self.protein_querier.get_properties()

        return self.fasta['sequence'], self.organism

    def peptide_query(self):
        """Call the class PeptideQuerier in itself. 在自身中调用PeptideQuerier类。

        Call the query method and get_properties method of the ProteinQuerier class to generate the sequence of
        peptides to be retrieved and export them.
        调用ProteinQuerier类的query方法和get_properties方法，生成要检索的肽段序列并导出。

        Returns:
            self.digestion_params: dict
            self.source: str
            self.peptides: list

        Updates:
            self.peptide_querier
            self.digestion_params
            self.source
            self.peptides
        """
        self.peptide_querier.query()
        self.digestion_params, self.source, self.peptides = self.peptide_querier.get_properties()

        return self.digestion_params, self.source, self.peptides

    def simple_query(self):
        """Default query methods.

        This function encapsulates the process of using the UniProt database to query protein metadata and
        the local library building method to generate the peptides to be retrieved.
        本函数封装了使用UniProt数据库查询蛋白质元数据，和本地建库的方法生成需要检索肽段的过程。

        Returns:
            None

        Updates:
            self.protein_querier
            self.peptide_querier
        """
        self.protein_querier = querier.UniProtProteinQuerier()
        self.protein_querier.set_params(accession=self.identifier)
        self.protein_query()

        self.peptide_querier = querier.LocalPeptideQuerier()
        self.peptide_querier.set_params(sequence=self.fasta['sequence'],
                                        organism=self.organism)
        self.peptide_query()

    def get_query_config(self):
        """Get query config.

        Returns:
            config: dict, {'identifier': self.identifier,
                           'peptides': self.peptides,
                           'organism': self.organism,
                           'digestion': self.digestion_params,
                           'proteome_source': self.source,
                           'description': self.fasta['description']
                          }
        """
        config = {'identifier': self.identifier,
                  'peptides': self.peptides,
                  'organism': self.organism,
                  'digestion': self.digestion_params,
                  'proteome_source': self.source,
                  'description': self.fasta['description']
                  }

        return config


class MapperHub:
    """This class is used to aggregate multiple proteomic database Mapper, merge and retrieve proteins.
       本类用于聚合多个蛋白质组数据库Mapper，对蛋白质进行合并检索。

    The merge retrieve function of this class only requires the user to insert the Mapper that needs to be retrieved,
    and the developer encourages users to try to build their own Mapper. In addition to the merge search function,
    this class also provides functions such as breakpoint and logging, which users can use by setting the task property.
    本类的合并检索功能只需要用户插入需要检索的Mapper即可，开发者鼓励用户尝试自建Mapper。除了合并搜索功能，
    本类还提供断点续传和日志等功能，用户可通过设置task属性使用。
    """

    def __init__(self, config: dict, mappers: list, task=None):
        # Set working task
        if task:
            self._task = task
        else:
            self._task = str(uuid.uuid4())
        self._tmp = dict()  # temporary data

        # Load the list of classes in the Mapper module, only the classes in this list can be inserted.
        # 载入Mapper模块中的类列表，只有在该列表内的类才可插入。
        self.supported_mappers = mapper.Mapper.__subclasses__()
        self.config = None  # configuration information that needs to be retrieved 用于检索的配置信息
        self.logger = None  # class Logger
        self.mappers = set()  # class Mapper

        self._set_logging()
        self._load_tmp()
        self._set_mappers(mappers)
        self._set_mapping_config(config)
        self._set_task_info()

    def _set_logging(self):
        """Configuring Logger"""
        self.logger = logger.CoreLogger(self._task)
        self.logger.set_core()
        self.logger.set_tmp()

    def _load_tmp(self):
        """Loading temporary data, implementing breakpoint transfer.

        Insert the new Mapper if you want to implement breakpoint transfer, you need to change the function.
        插入新建的Mapper类，若想要实现断点续传，需更改该函数。
        """
        task = self._task
        if os.path.exists(f'patpat_env/tmp/{task}.log'):
            with open(f'patpat_env/tmp/{task}.log', mode='r') as f:
                data = [re.split('\t', i) for i in f.readlines()]
            self._tmp['PRIDE_peptide'] = [[i[1], i[2][:-1]] for i in data if i[0] == 'PRIDE_peptide']
            self._tmp['iProX_peptide'] = [[i[1], i[2][:-1]] for i in data if i[0] == 'iProX_peptide']
        else:
            self._tmp['PRIDE_peptide'] = []
            self._tmp['iProX_peptide'] = []

    def _set_mapping_config(self, config: dict):
        """Update configuration"""
        self.config = config.copy()
        self.config['task'] = self._task
        self.config['mappers'] = []
        for m in self.mappers:
            self.config['mappers'].extend([m.__str__])
        self.config['state'] = 'Preparation'

    def _set_mappers(self, mappers: {list, set}):
        """Set Mappers"""
        for mapper_ in mappers:
            if (type(mapper_) in self.supported_mappers and
                    not type(mapper_) in [type(m) for m in self.mappers]):
                self.mappers.add(mapper_)

    def mapping(self):
        """Core methods of this class

        Insert the new Mapper if you want to implement breakpoint transfer, you need to change the function.
        插入新建的Mapper类，若想要实现断点续传，需更改该函数。
        """
        mappers = set()
        self.config['startTime'] = time.time()
        self.config['state'] = 'Running'
        self._set_task_info()

        try:
            for mapper_ in self.mappers:

                # These three keys (protein UniProt identifier, peptide to be searched, species) are
                # required in the configuration.
                # 这三个键（蛋白UniProt识别符、需要搜索的肽段、物种）是配置中必须的。
                mapper_._identifier = self.config['identifier']
                mapper_._peptides = self.config['peptides']
                mapper_._organism = self.config['organism']

                # Implement breakpoint transfer
                tmp = None
                if mapper_.__class__ is mapper.PrideMapper and self._tmp['PRIDE_peptide']:
                    tmp = self._tmp['PRIDE_peptide']
                elif mapper_.__class__ is mapper.IProXMapper and self._tmp['iProX_peptide']:
                    tmp = self._tmp['iProX_peptide']

                if tmp:
                    tmp_peptides = set([i[0] for i in tmp])
                    mapper_._peptides = list(set(mapper_._peptides) - tmp_peptides)
                    for p in tmp_peptides:
                        logging.getLogger('core').debug(f'Local loading: {p}')

                    tmp = [i for i in tmp if utility.usi_detect(i[1])]  # 去掉空白搜索

                mapper_.mapping(peptide_tmp=tmp)
                mapper_.filtering()
                mappers.add(mapper_)
            self.mappers = mappers
        except:
            self.config['state'] = 'Error'
            logging.getLogger('core').error('RUNNING ERROR')
            self._set_task_info()
            raise BufferError(f'RUNNING ERROR:{sys.exc_info()}')
        else:
            self.config['state'] = 'Success'
            self._set_task_info()

    def _set_task_info(self):
        """Write information about this task to a local file."""
        if os.path.exists('patpat_env/logs/tasks.json'):
            with open('patpat_env/logs/tasks.json', 'r') as fr:
                configs = json.loads(fr.readline())
                config = self.config.copy()
                config['peptides'] = len(config['peptides'])
                configs['tasks'][config['task']] = config
            fr.close()
            with open('patpat_env/logs/tasks.json', 'w') as fw:
                configs_json = json.dumps(configs)
                fw.write(configs_json)
        else:
            config = self.config.copy()
            config['peptides'] = len(config['peptides'])
            configs = dict()
            configs['tasks'] = {config['task']: config}
            configs_json = json.dumps(configs)
            with open('patpat_env/logs/tasks.json', mode='w') as f:
                f.write(configs_json)

    def export(self):
        """Exporting retrieved results.

        Returns:
            output: dir
        """
        output = dict()
        for m in self.mappers:
            output[m.__str__] = m.export()

        if os.path.exists(f'patpat_env/result/{self._task}'):
            pass
        else:
            os.mkdir(f'patpat_env/result/{self._task}')

        with open(f'patpat_env/result/{self._task}/config.json', 'w') as fw:
            config = json.dumps(self.config)
            fw.write(config)

        with open(f'patpat_env/result/{self._task}/result.json', 'w') as fw:
            output_json = json.dumps(output)
            fw.write(output_json)

        with open(f'patpat_env/result/{self._task}/result.tsv', 'w') as fw:
            fw.write('title\tsummary\twebsite\n')
            for projects in output.values():
                for project in projects.values():
                    try:
                        fw.write(f"{project['title']}\t{project['summary']}\t{project['website']}\n")
                    except UnicodeEncodeError:
                        fw.write(
                            f"{project['title']}\t{project['summary'].encode(errors='ignore')}\t{project['website']}\n")

        return output
