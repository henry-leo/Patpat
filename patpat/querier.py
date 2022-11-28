"""This module is used to build the configuration information that needs to be retrieved.
   该模块用于建立需要检索的配置信息。

The configuration information required to complete the mapping are the protein metadata and the peptide sequences
to be retrieved. The necessary protein metadata are: the protein identifier, the species or tissue to which
the protein belongs, and the protein sequence, where the protein sequence is used to generate the peptide sequences
to be retrieved.
完成映射所需要的配置信息有蛋白质的元数据和需要检索的肽段序列。
必要的蛋白质元数据有：蛋白质识别符、蛋白质所属物种或组织、蛋白质序列，其中蛋白质序列用于产生需要检索的肽段序列。

    Typical usage example:
    典型用法示例：

    from patpat import querier

    p1 = querier.UniProtProteinQuerier()
    p1.identifier = 'E9PV96'
    p1.query()
    _, organism_, fasta_ = p1.get_properties()

    p2 = querier.LocalPeptideQuerier()
    p2.set_params(sequence=fasta_['sequence'],
                  organism=organism_)
    p2.query()
    digestion_params_, source_, filtered_peptides_ = p2.get_properties()

"""

import collections
import re

import requests
import tqdm
from pyteomics import parser, fasta

from . import utility


class ProteinQuerier(object):
    """Base class for implementing the protein metadata search function. 实现蛋白质元数据搜索功能的基类。

    This base class gives preliminary interface definitions to build subclasses against protein databases.
    这个基类给出了初步的接口定义，以便针对蛋白质数据库建立子类。
    """

    def __init__(self):
        # Structuring protein fasta into a dictionary. 将蛋白质fasta结构化为字典。
        self.fasta: dict = {
            'description': '',
            'sequence': ''
        }
        self.organism: dict = dict()  # the species or tissue to which the protein belongs 蛋白质物种或组织信息
        self.identifier: str = ''  # Identifier of protein 蛋白质的识别符
        self._headers: dict = {'Accept': 'application/json'}

    def _request_api(self, *args):
        """Functions that actually interact with the protein database. 实际与蛋白质数据库交互的函数。"""
        raise NotImplementedError

    def set_params(self, *args):
        """Functions for setting query parameters. 设置查询参数的函数。"""
        raise NotImplementedError

    def get_properties(self):
        """Functions to get query results. 获取查询结果的函数。"""
        raise NotImplementedError

    def query(self):
        """Functions for external calls. 供外部调用的函数。"""
        raise NotImplementedError


class PeptideQuerier(object):
    """Generate the base class of peptide sequences to be retrieved. 生成要检索的肽段序列的基类。

    This base class gives preliminary interface definitions to create subclasses for the peptide generation methods.
    这个基类给出了初步的接口定义，以便针对肽段生成方法建立子类。
    """

    def __init__(self):
        self.sequence: str = ''  # protein sequence 蛋白质序列
        self.organism: dict = dict()  # the species or tissue to which the protein belongs 蛋白质物种或组织信息

        # Filtered peptides for output to external 经过筛选的肽段，可输出至外部
        self.filtered_peptides: list = []

        # Unfiltered full peptides generated from protein sequences 根据蛋白质序列产生的未经过滤全部肽段
        self._peptides: list = []

    def set_params(self, *args):
        """Functions for setting query parameters. 设置查询参数的函数。"""
        raise NotImplementedError

    def get_properties(self):
        """Functions to get query results. 获取查询结果的函数。"""
        raise NotImplementedError

    def query(self):
        """Functions for external calls. 供外部调用的函数。"""
        raise NotImplementedError


class UniProtProteinQuerier(ProteinQuerier):
    """ProteinQuerier subclass for querying protein metadata through the UniProt database.
       通过UniProt数据库查询蛋白质元数据的ProteinQuerier子类。
    """

    def __init__(self):
        super().__init__()

    def _request_api(self, accession):
        """Functions that actually interact with the UniProt database. 实际与UniProt数据库交互的函数。

        Args:
            accession: str. UniProt accession of protein 蛋白质的UniProt编号

        Returns:
            protein_fasta: str. fasta of protein 蛋白质fasta

        """
        acc = accession
        url = f"https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query={acc}"
        protein_fasta = requests.get(url,
                                     headers=self._headers).text
        return protein_fasta

    def set_params(self, accession):
        """Functions for setting query parameters. 设置查询参数的函数。

        Args:
            accession: str. UniProt accession of protein 蛋白质的UniProt编号

        Returns:
            None

        Updates:
            self.identifier
        """
        if re.match('[OPQ]\\d[A-Z\\d]{3}\\d|[A-NR-Z]\\d([A-Z][A-Z\\d]{2}\\d){1,2}', accession).group() == accession:
            self.identifier = accession
        else:
            raise TypeError('Can not identify identifier, please input UniProt accession.')

    def query(self):
        """Functions for external calls. 供外部调用的函数。

        Returns:
            None

        Updates:
            self.fasta
            self.organism

        Raises:
            ValueError: self.identifier is empty.
        """
        if not self.identifier:
            raise ValueError('Please run func self.set_params() first. 请先运行self.set_params()函数')

        protein_fasta = self._request_api(self.identifier)
        fasta_sep = re.split('\n', protein_fasta)
        self.fasta['description'] = fasta_sep[0]
        self.fasta['sequence'] = ''.join(fasta_sep[1:])

        self.organism['name'] = re.search('(?<=OS=).*(?=OX=)',
                                          self.fasta['description']).group()[:-1]
        self.organism['accession'] = re.search('(?<=OX=).*(?=GN=)',
                                               self.fasta['description']).group()[:-1]

    def get_properties(self):
        """Functions to get query results. 获取查询结果的函数。

        Returns:
            self.identifier: str. UniProt accession of protein 蛋白质的UniProt编号
            self.organism: dict. UniProt taxonomy of protein 蛋白质所属的UniProt物种分类
            self.fasta: dict.

        Raises:
            ValueError: self.identifier is empty
        """
        if self.identifier:
            return self.identifier, self.organism, self.fasta
        elif self.identifier is None:
            raise ValueError('Please run func self.set_params() first. 请先运行self.set_params()函数')


class LocalPeptideQuerier(PeptideQuerier):
    """The PeptideQuerier subclass used to generate the peptide sequences by local library building.
       通过本地建库的方法产生肽段序列的PeptideQuerier子类。
    """

    def __init__(self):
        super().__init__()

        self.digestion_params = None  # Parameters of protein in-silico enzymatic 蛋白质模拟酶切的参数
        self.source = None  # Local proteome file address 本地蛋白质组文件地址
        self._searched_peptides = None  # Peptide information after local database search 经过本地数据库搜索的肽段信息
        self._base = None  # Local Peptide Database 本地肽段数据库

    def set_params(self,
                   sequence,
                   organism,
                   digestion_params: dict = None,
                   source: str = None):
        """Functions for setting query parameters. 设置查询参数的函数。

        Args:
            sequence: str. protein sequence 蛋白质序列
            organism: dict. the species or tissue to which the protein belongs 蛋白质物种或组织信息
            digestion_params: dict. Parameters of protein in-silico enzymatic 蛋白质模拟酶切的参数
            source: str. Local proteome file address 本地蛋白质组文件地址

        Returns:
            None

        Updates:
            self.sequence
            self.organism
            self.digestion_params
            self.source
        """
        self.sequence = sequence
        self.organism = organism

        if digestion_params is None:
            digestion_params = {
                'rules': 'trypsin',
                'miss': 1,
                'min_length': 7,
                'max_length': 35}

        self.digestion_params = digestion_params
        self.source = source

    def digestion(self):
        """Functions for protein in-silico enzymatic 用于蛋白质模拟酶切的函数

        Returns:
            None

        Updates:
            self._peptides
        """
        self._peptides = list(parser.cleave(self.sequence,
                                            rule=self.digestion_params['rules'],
                                            missed_cleavages=self.digestion_params['miss'],
                                            min_length=self.digestion_params['min_length'],
                                            max_length=self.digestion_params['max_length']))

    def search(self, threshold=1):
        """Functions for local library search. 用于本地数据库搜索的函数

        Args:
            threshold: int. Filtering threshold, which refers to the peptide is shared by how many proteins,
                            greater than the threshold is filtered
                            过滤阈值，指被多少个蛋白质共享的肽段，大于阈值则被过滤

        Returns:
            None

        Updates:
            self._searched_peptides
            self._base
            self.filtered_peptides
        """
        print('Choose local peptide search.')

        source = self._select_local_proteome()

        print(f'Proteome file: {source}')

        f = fasta.read(source)
        input_peps = self._peptides
        try:
            rule = self.digestion_params['rules']
            miss = self.digestion_params['miss']
            min_length = self.digestion_params['min_length']
            max_length = self.digestion_params['max_length']

        except ValueError:
            raise ValueError('digestion() params error!')

        base = collections.defaultdict(list)
        print('Creating database...')
        for protein in tqdm.tqdm(f):
            peptides = parser.cleave(protein.sequence,
                                     rule=rule,
                                     missed_cleavages=miss,
                                     min_length=min_length,
                                     max_length=max_length)
            for peptide in peptides:
                try:
                    base[peptide][0] += 1
                    base[peptide][1].add(protein.description)
                except IndexError:
                    base[peptide] = [1, {protein.description}]

        ans = collections.defaultdict(list)
        for pep in input_peps:
            try:
                ans[pep].append(base[pep][0])
                ans[pep].append(base[pep][1])
            except IndexError:
                ans[pep] = [0, 'loss']
        ans = sorted(zip(ans.values(), ans.keys()))
        ans = [[j for j in utility.flatten(i)] for i in ans]
        errs = [i[2] for i in ans if i[0] == 0]
        if len(errs):
            print(f'There are some peptides that are not in the query database: {errs}')

        self._searched_peptides = ans
        self._base = base
        self._search_filter(threshold=threshold)

    def _select_local_proteome(self):
        """The function is used to automatically search for locally available proteome file
           此函数用于自动化搜索本地可用的蛋白质组文件

        Returns:
            None

        Updates:
            self.source
        """
        source = self.source

        if source is None:
            taxonomy_id = self.organism['accession']
            up_id = utility.uniprot_proteome_id_detect(taxonomy_id)
            file_list = utility.list_local_proteome()
            file_info = utility.proteome_file_info_split(file_list)

            target_files = [i for i in file_info if i[1] == up_id]
            if not target_files:
                print(f"The {self.organism['name']} {up_id} proteome file was not found locally.")

                flag = input('Do you want to download it?(y/n)')
                if flag == 'y':
                    source = utility.download_uniprot_opg(taxonomy_id)
                else:
                    raise FileNotFoundError(f'User cancel operation.')

            else:
                target_file = [i[0] for i in target_files if i[2] == max([j[2] for j in target_files])][0]
                source = f'patpat_env/proteome/uniprot-proteome_{target_file}.fasta'

        else:
            source = f'patpat_env/proteome/{source}'

        self.source = source

        return source

    def _search_filter(self, threshold=1):
        """The function is used to filter peptide sequences that are larger than the threshold
           本函数用于过滤大于阈值的肽段序列
        """
        ans = self._searched_peptides

        ans = [i[2] for i in ans if i[0] <= threshold and i[0] != 0]
        self.filtered_peptides = ans

    def get_properties(self):
        """Functions to get query results. 获取查询结果的函数。

        Returns:
            self.digestion_params: dict. Parameters of protein in-silico enzymatic 蛋白质模拟酶切的参数
            self.source: str. Local proteome file address 本地蛋白质组文件地址
            self.filtered_peptides: list. Filtered peptides for output to external 经过筛选的肽段，可输出至外部

        Raises:
            ValueError: self.sequence or self.organism is empty.
        """
        if not self.sequence or not self.organism:
            raise ValueError('Please run func self.set_params() first. 请先运行self.set_params()函数')

        return self.digestion_params, self.source, self.filtered_peptides

    def query(self):
        """Functions for external calls. 供外部调用的函数。

        Raises:
            ValueError: self.sequence or self.organism is empty.
        """
        if not self.sequence or not self.organism:
            raise ValueError('Please run func self.set_params() first. 请先运行self.set_params()函数')

        self.digestion()
        self.search()


class UniProtPeptideQuerier(PeptideQuerier):
    """Under Development..."""

    def __init__(self,
                 sequence,
                 digestion_params: dict = None,
                 search_args=None):
        super().__init__()
        self.sequence = sequence

        if digestion_params is None:
            digestion_params = {
                'rules': 'trypsin',
                'miss': 1,
                'min_length': 7,
                'max_length': 35}

        self.digestion_params = digestion_params
        self.uniprot_job = dict()
        self.search_args = search_args
        self.source = 'UNIPROT'

        self._searched_peptides = None
        self._digestion_params = None
        self._peptides = None

        self.filtered_peptides = None

        # CALL ERROR
        raise NotImplementedError("This class is under development...")

    def set_params(self):
        pass

    def digestion(self):
        self._peptides = list(parser.cleave(self.sequence,
                                            rule=self.digestion_params['rules'],
                                            missed_cleavages=self.digestion_params['miss'],
                                            min_length=self.digestion_params['min_length'],
                                            max_length=self.digestion_params['max_length']))

    def search(self):
        print('Choose UniProt peptide search.')

    def _uniprot_peptide_search_request(self):
        print('job id={}'.format(self.uniprot_job['job']))

    def _uniprot_peptide_search_check(self):
        # print("This took: {} sec".format(int(end - start)))
        pass

    def _uniprot_peptide_search_answer(self):
        pass

    def _uniprot_search_filter(self, *args, threshold=1):
        pass

    def get_properties(self):
        return self.digestion_params, self.source, self.filtered_peptides

    def query(self):
        self.digestion()
        self.search()


if __name__ == '__main__':
    identifier_ = 'E9PV96'
    p1 = UniProtProteinQuerier()
    p1.set_params(accession=identifier_)
    p1.query()
    _, organism_, fasta_ = p1.get_properties()

    p2 = LocalPeptideQuerier()
    p2.set_params(sequence=fasta_['sequence'],
                  organism=organism_)
    p2.query()
    digestion_params_, source_, filtered_peptides_ = p2.get_properties()
