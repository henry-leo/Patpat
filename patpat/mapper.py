"""Call the retriever module and implement the retrieve function. 调用retriever模块并实现检索功能

This module builds different classes according to the database oriented, and each class calls the corresponding
retriever to realize the retrieve by itself. The protein/peptide retriever is used to implement the protein/peptide
level search, and the protein-to-project mapping is implemented by the project retriever. In addition, each class
also needs to implement other functions: 1. a filtering function to remove incorrect mappings; 2. an export function
to structure the retrieve results for downstream modules to call.
此模块根据面向数据库的不同构建不同的类，每个类自行调用对应的retriever实现搜索。通过protein/peptide retriever实现蛋白质/肽段水平的检索，
在通过project retriever实现蛋白质对项目的映射。此外，每个类也需要实现其他功能：1.过滤功能，以去除错误的映射；
2.导出功能，结构化检索结果，以供下游模块调用。

    Typical usage example:
    典型用法示例：

    from patpat import hub
    from patpat import mapper
    from patpat import utility

    utility.init()
    utility.initiate_uniprot_proteome_catalog()

    identifier = 'E9PV96'
    q = hub.QueryHub()
    q.identifier = identifier
    q.simple_query()

    conf = q.get_query_config()

    m = mapper.PrideMapper()

    # Add configs
    m._identifier = conf['identifier']
    m._peptides = conf['peptides']
    m._organism = conf['organism']

    m.mapping()
    m.filtering()

    output = m.export()
"""
import time
import re
import logging

import tqdm
import requests

from . import retriever
from . import utility


class Mapper(object):
    """Base class for implementing functions such as retrieve and mapping. 用于实现检索和映射等功能的基类。

    This base class gives preliminary interface definitions to build base classes against databases.
    这个基类给出了初步的接口定义，以便针对数据库建立基类。
    """

    def mapping(self):
        """Functions that implement the retrieve and complete the mapping. 实现搜索和完成映射的函数。

        Parameters can be set according to actual needs.
        可以根据实际需要设置参数。
        """
        raise NotImplementedError

    def filtering(self):
        """Implementing filtering functions for incorrect mappings. 实现对不正确映射的过滤功能。"""
        raise NotImplementedError

    def export(self):
        """Implement export function for downstream module call. 为下游模块的调用实现导出功能。"""
        raise NotImplementedError

    @staticmethod
    def retriever_running(any_retriever, payload: dict = None, overtime: int = 5, sleep: int = 3):
        """Decorator or Closures. Pack the retriever to prevent processes from being interrupted
        during network fluctuations.
        装饰器或闭包。封装检索器，防止因网络波动而中断进程。

        Args:
            any_retriever: classes in retriever.py retriever模块中的类
            payload: parameters of class retriever retriever类的参数
            overtime: The number of retrieve timeouts, the default is to skip the retrieve after five errors.
                    检索超时的次数，默认是在出现五个错误后跳过检索。
            sleep: time.sleep(sleep: int)

        Returns:
            any_retriever: retriever class after retrieve, retrieve may or may not succeed.
                    经检索之后的retriever类，但检索可能成功也可能失败。
            flag: bool. False means retrieve succeed. 布尔数。False意味着检索成功。
        """
        if payload is None:
            pass
        else:
            for k, v in payload.items():
                any_retriever.payload[k] = v

        t = 0
        flag = False
        while t <= overtime:
            t += 1
            try:
                any_retriever.retrieve()
            except requests.exceptions.RequestException as e:
                logging.getLogger('core').warning(f'Error:{e}')
                time.sleep(sleep)

                if t > overtime:
                    flag = True
                    logging.getLogger('core').info(f'{t}/{overtime} times, Stop!')
                else:
                    logging.getLogger('core').info(f'{t}/{overtime} times, Try again!')

            else:
                break

        return any_retriever, flag


class PrideMapper(Mapper):
    """Mapper subclass for PRIDE database. 针对PRIDE数据库的Mapper子类。"""

    def __init__(self, payload=None):
        self.__str__ = 'PRIDE'  # name of database
        self._identifier = None  # str, Protein accession
        self._peptides = None  # list,  Peptide sequences to be retrieved
        self._organism = None  # dict, Organism to which the protein belongs. E.g.: {'name': 'xx', 'accession': 'xx'}
        self._protein_level = None  # list, retrieve results for protein levels
        self._peptides_level = None  # list, retrieve results for peptide levels
        self._project_base = None  # dict, retrieved items
        self.payload = payload
        self.overtime = 5

    def mapping(self, peptide_tmp: list = None):
        """Functions that implement the retrieve and complete the mapping. 实现检索和完成映射的函数。"""
        logging.getLogger('core').info('Select PRIDE to search.')

        # Conduct protein & peptide retrieve
        peptides_level = self.mapping2peptides()
        protein_level = self.mapping2protein()

        # Load temporary files 实现断点续传
        if peptide_tmp:
            peptides_level.extend(peptide_tmp)

        # Mapping to project
        self.mapping2project(protein_level, peptides_level)

    def mapping2protein(self):
        """Implementing protein level retrieve. 实现蛋白质水平的检索。"""
        if self._identifier is None:
            protein_level = None
        else:
            protein_retriever = retriever.PrideProteinRetriever()
            protein_retriever.request_word = self._identifier
            protein_retriever.get_payloads_on_web()

            protein_retriever, _ = self.retriever_running(protein_retriever, payload=self.payload)

            protein_level = self.protein_response_parse(protein_retriever)

        protein_level = [i for i in protein_level if utility.usi_detect(i[1])]
        self._protein_level = protein_level
        return protein_level

    @staticmethod
    def protein_response_parse(protein_retriever):
        """Implementing structured results of protein level retrieve. 实现蛋白质水平检索的结果结构化。"""
        t = [p['peptideevidences'] for p in protein_retriever.response.values()]
        ans_list = []
        for i in t:
            for j in i:
                ans_list.append([j['peptideSequence'], j['_links']['psms']['href']])

        ans_list = [[i[0], utility.url_split(i[1])[1].get('usi')] for i in ans_list
                    if not utility.url_split(i[1])[1].get('usi') is None]
        return ans_list

    def mapping2peptides(self):
        """Implementing peptide level retrieve. 实现肽段水平的检索。"""
        if self._peptides is None:
            peptides_level = None
            logging.getLogger('core').warning('No peptide can be retrieve!')
        elif isinstance(self._peptides, list):
            peptides_level = []
            logging.getLogger('core').info(f'Find {len(self._peptides)} peptides all.')
        else:
            raise TypeError

        for peptide in self._peptides:
            peptide_retriever = retriever.PridePeptideRetriever()
            peptide_retriever.request_word = peptide
            peptide_retriever.get_payloads_on_web()

            peptide_retriever, pass_flag = self.retriever_running(peptide_retriever, payload=self.payload)

            ans_list = self.peptide_response_parse(peptide_retriever)
            if ans_list:
                peptides_level.extend(ans_list)
                for i in ans_list:
                    logging.getLogger('tmp').info(f'PRIDE_peptide\t{i[0]}\t{i[1]}')
            else:
                logging.getLogger('tmp').info(f'PRIDE_peptide\t{peptide}\tNone')  # 把空白搜索加入缓存中以避免重复搜索

            if pass_flag:
                continue

        peptides_level = [i for i in peptides_level if utility.usi_detect(i[1])]
        self._peptides_level = peptides_level
        return peptides_level

    @staticmethod
    def peptide_response_parse(peptide_retriever):
        """Implementing structured results of peptide level retrieve. 实现肽段水平检索的结果结构化。"""
        t = [p['spectraevidences'] for p in peptide_retriever.response.values()]
        ans_list = []
        for i in t:
            for j in i:
                ans_list.append([j.get('peptideSequence'), j.get('usi')])
        return ans_list

    def mapping2project(self, protein_level, peptides_level):
        """Implementing mapping between retrieve result and project. 实现检索结果和项目之间的映射。"""
        logging.getLogger('core').info('Mapper to project...')
        project_base = dict()

        if protein_level:
            protein_level = [[p[0], p[1]] + [i for i in utility.usi_split(p[1]).values()] for p in protein_level]
            self._protein_level = protein_level
        if peptides_level:
            peptides_level = [[p[0], p[1]] + [i for i in utility.usi_split(p[1]).values()] for p in peptides_level]
            self._peptides_level = peptides_level

        projects = set([i[2] for i in peptides_level]).union(set([i[2] for i in protein_level]))

        for project in tqdm.tqdm(projects):
            project_retriever = retriever.PrideProjectRetriever()
            project_retriever.request_word = project

            project_retriever, pass_flag = self.retriever_running(project_retriever, payload=self.payload)

            project_base[project] = project_retriever.response

            if pass_flag:
                continue

        self._project_base = project_base
        logging.getLogger('core').info('PRIDE mapping completed.')

    def filtering(self):
        """Implementing remove for incorrect mappings. 排除不正确的映射。

        Remove Criteria: remove of peptide level mapping that are not from the same species.
        排除标准：排除不是同一物种的肽段水平映射。

        Returns:
                None
        Updates:
            self.project_base: dict, key: project accession, value: project content
            self.peptide_level: n*5 list,
                            5: peptide sequence, usi, project accession, msRun+index, interpretation
        """
        del_projects = []
        for p in self._project_base.values():
            project_org = [i['accession'] for i in p['organisms']]
            if not self._organism['accession'] in project_org:
                del_projects.extend([p['accession']])
        for p in del_projects:
            self._project_base.pop(p)

    def export(self):
        """Implement export function for downstream module call. 为下游模块的调用实现导出功能。

        Returns:
            projects: dict. Structured project information. 结构化的项目信息。
        """
        projects = self._project_base.copy()
        for project in projects:
            if self._peptides_level:
                peptides = [i[1] for i in self._peptides_level if i[2] == project]
            else:
                peptides = None
            if self._protein_level:
                protein = [i[1] for i in self._protein_level if i[2] == project]
            else:
                protein = None

            # These four keys are necessary 这四个key是必要的

            # Used to record the mapping of protein and peptide levels
            projects[project]['peptides'] = peptides
            projects[project]['protein'] = protein

            # Used to record the summary of project
            projects[project]['summary'] = projects[project]['projectDescription']

            # Used to record the website of project
            projects[project]['website'] = f'http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={project}'

        return projects


class IProXMapper(Mapper):
    """Mapper subclass for iProX database. 针对iProX数据库的Mapper子类。"""

    def __init__(self, payload=None):
        self.__str__ = 'iProX'  # name of database
        self._identifier = None  # str, Protein accession
        self._peptides = None  # list,  Peptide sequences to be retrieved
        self._organism = None  # dict, Organism to which the protein belongs. E.g.: {'name': 'xx', 'accession': 'xx'}
        self._protein_level = None  # list, retrieve results for protein levels
        self._peptides_level = None  # list, retrieve results for peptide levels
        self._project_base = None  # dict, retrieved items
        self.payload = payload
        self.overtime = 5

    def mapping(self, peptide_tmp: list = None):
        """Functions that implement the retrieve and complete the mapping. 实现检索和完成映射的函数。"""
        logging.getLogger('core').info('Select iProX to search.')

        # Conduct protein & peptide retrieve
        peptides_level = self.mapping2peptides()
        protein_level = self.mapping2protein()

        # Load temporary files
        if peptide_tmp:
            peptides_level.extend(peptide_tmp)

        # Mapping to project
        self.mapping2project(protein_level, peptides_level)

    def mapping2protein(self):
        """Implementing protein level retrieve. 实现蛋白质水平的检索。"""
        if self._identifier is None:
            protein_level = None
        else:
            protein_retriever = retriever.IProXProteinRetriever()
            protein_retriever.request_word = self._identifier
            protein_retriever.get_payloads_on_web()

            protein_retriever, _ = self.retriever_running(protein_retriever, payload=self.payload)

            protein_level = self.protein_response_parse(protein_retriever)

        protein_level = [i for i in protein_level if utility.usi_detect(i[1])]
        self._protein_level = protein_level
        return protein_level

    @staticmethod
    def protein_response_parse(protein_retriever):
        """Implementing structured results of protein level retrieve. 实现蛋白质水平检索的结果结构化。"""
        ans_list = []

        for ps in protein_retriever.response.values():
            ans_list.extend([[p['peptideSequence'], p['usi']] for p in ps])

        return ans_list

    def mapping2peptides(self):
        """Implementing peptide level retrieve. 实现肽段水平的检索。"""
        if self._peptides is None:
            peptides_level = None
            logging.getLogger('core').warning('No peptide can be retrieve!')
        elif isinstance(self._peptides, list):
            peptides_level = []
            logging.getLogger('core').info(f'Find {len(self._peptides)} peptides all.')
        else:
            raise TypeError

        for n, peptide in enumerate(self._peptides):
            peptide_retriever = retriever.IProXPeptideRetriever()
            peptide_retriever.request_word = peptide
            peptide_retriever.get_payloads_on_web()

            peptide_retriever, pass_flag = self.retriever_running(peptide_retriever, payload=self.payload)

            time.sleep(n / len(self._peptides))

            ans_list = self.peptide_response_parse(peptide_retriever)
            if ans_list:
                peptides_level.extend(ans_list)

                for i in ans_list:
                    logging.getLogger('tmp').info(f'iProX_peptide\t{i[0]}\t{i[1]}')
            else:
                logging.getLogger('tmp').info(f'iProX_peptide\t{peptide}\tNone')
            if pass_flag:
                continue

        peptides_level = [i for i in peptides_level if utility.usi_detect(i[1])]
        self._peptides_level = peptides_level
        return peptides_level

    @staticmethod
    def peptide_response_parse(peptide_retriever):
        """Implementing structured results of peptide level retrieve. 实现肽段水平检索的结果结构化。"""
        ans_list = []
        requests_word = peptide_retriever.request_word
        for ps in peptide_retriever.response.values():
            for p in ps:
                interpretation = utility.usi_split(p['usi'])['interpretation']
                answer_word = re.search('.*(?=/)', re.sub('\\[.*?]', '', interpretation)).group()
                if requests_word == answer_word:
                    ans_list.append([requests_word, p['usi']])

        return ans_list

    def mapping2project(self, protein_level, peptides_level):
        """Implementing mapping between retrieve result and project. 实现检索结果和项目之间的映射。"""
        logging.getLogger('core').info('Mapper to project...')
        project_base = dict()

        if protein_level:
            protein_level = [[p[0], p[1]] + [i for i in utility.usi_split(p[1]).values()] for p in protein_level]
            self._protein_level = protein_level
        if peptides_level:
            peptides_level = [[p[0], p[1]] + [i for i in utility.usi_split(p[1]).values()] for p in peptides_level]
            self._peptides_level = peptides_level

        projects = set([i[2] for i in peptides_level]).union(set([i[2] for i in protein_level]))

        for project in tqdm.tqdm(projects):
            project_retriever = retriever.IProXProjectRetriever()
            project_retriever.request_word = project

            project_retriever, pass_flag = self.retriever_running(project_retriever, payload=self.payload)

            project_base[project] = project_retriever.response

            if pass_flag:
                continue

        self._project_base = project_base
        logging.getLogger('core').info('iProX mapping completed.')

    def filtering(self):
        """Implementing remove for incorrect mappings. 排除不正确的映射。

        Remove Criteria: remove of peptide level mapping that are not from the same species.
        排除标准：排除不是同一物种的肽段水平映射。

        Returns:
            None

        Updates:
            Update self.project_base: dict, key: project accession, value: project content
            Update self.peptide_level: n*5 list,
                                5: peptide sequence, usi, project accession, msRun+index, interpretation
        """
        projects = [
            p['accession']['value'] for p in self._project_base.values()
            if self._organism['accession'] in [
                i['name'] for i in p['species']
                if i['accession'] == 'MS:1001467'
            ]
        ]

        project_base = dict()
        for p in projects:
            project_base[p] = self._project_base[p]
        self._project_base = project_base

    def export(self):
        """Implement export function for downstream module call. 为下游模块的调用实现导出功能。

        Returns:
            projects: dict. Structured project information. 结构化的项目信息。
        """
        projects = self._project_base.copy()
        for project in projects:
            if self._peptides_level:
                peptides = [i[1] for i in self._peptides_level if i[2] == project]
            else:
                peptides = None
            if self._protein_level:
                protein = [i[1] for i in self._protein_level if i[2] == project]
            else:
                protein = None

            # These four keys are necessary 这四个key是必要的

            # Used to record the mapping of protein and peptide levels
            projects[project]['peptides'] = peptides
            projects[project]['protein'] = protein

            # Used to record the summary of project
            projects[project]['summary'] = projects[project]['summary']

            # Used to record the website of project
            projects[project]['website'] = f'http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={project}'

        return projects
