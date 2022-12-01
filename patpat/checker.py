"""This module is used to check the connectivity of external proteomics databases 此模块用于确认外部蛋白质组学数据库的连通性

    Typical usage example:
    典型用法示例：

    from patpat import checker
    from patpat import retriever

    c = checker.GenericChecker()
    c.peptide_retrievers = [
        retriever.PridePeptideRetriever(),
        retriever.IProXPeptideRetriever()
    ]
    c.check()
"""
import requests
import tqdm

from patpat import retriever


class GenericChecker:
    def check(self):
        raise NotImplementedError


class PrideChecker(GenericChecker):
    def __init__(self, times=3):
        self.retrievers = None
        self.times = times
        self.flag = {
            'project': 0,
            'protein': 0,
            'peptide': 0
        }
        self.message = {
            'project': [],
            'protein': [],
            'peptide': []
        }

    def _implement(self):
        retrievers = retriever.GenericPrideRetriever.__subclasses__()

        self.retrievers = [r() for r in retrievers]

    def _ans(self):
        print("\nPRIDE Connectivity Report")
        print(f"\tProject Retriever: {bool(self.flag['project'])}")
        print(f"\tProtein Retriever: {bool(self.flag['protein'])}")
        print(f"\tPeptide Retriever: {bool(self.flag['peptide'])}")

    def check(self):
        self._implement()

        flag = None
        message = []

        for r in self.retrievers:
            if isinstance(r, retriever.PrideProjectRetriever):
                flag, message = self._project_check(project_retriever=r)
                if flag:
                    self.flag['project'] = 1
                self.message['project'] = message

            elif isinstance(r, retriever.PrideProteinRetriever):
                flag, message = self._protein_check(protein_retriever=r)
                if flag:
                    self.flag['protein'] = 1
                self.message['protein'] = message

            elif isinstance(r, retriever.PridePeptideRetriever):
                flag, message = self._peptide_check(peptide_retriever=r)
                if flag:
                    self.flag['peptide'] = 1
                self.message['peptide'] = message
        self._ans()

    def _project_check(self, project_retriever):
        t = project_retriever
        flag = False
        message = []

        print(f"\nCheck the connectivity of the PRIDE Project API")
        for _ in tqdm.tqdm(range(self.times)):
            try:
                c = requests.get(t.api, params=t.example)
            except requests.exceptions.RequestException as e:
                message.extend([e])
                continue
            else:
                if c.ok:
                    flag = True
                    break
                else:
                    message.extend([c.text])
        if flag:
            return flag, message
        else:
            return flag, message

    def _protein_check(self, protein_retriever):
        t = protein_retriever
        flag = False
        message = []

        print(f"\nCheck the connectivity of the PRIDE Protein API")
        for _ in tqdm.tqdm(range(self.times)):
            try:
                payloads = {'proteinAccession': t.example,
                            'pageSize': 5}
                c = requests.get(t.api, params=payloads)
            except requests.exceptions.RequestException as e:
                message.extend([e])
                continue
            else:
                if c.ok:
                    flag = True
                    break
                else:
                    message.extend([c.text])
        if flag:
            return flag, message
        else:
            return flag, message

    def _peptide_check(self, peptide_retriever):
        t = peptide_retriever
        flag = False
        message = []

        print(f"\nCheck the connectivity of the PRIDE Peptide API")
        for _ in tqdm.tqdm(range(self.times)):
            try:
                payloads = {'peptideSequence': t.example,
                            'pageSize': 5}
                c = requests.get(t.api, params=payloads)
            except requests.exceptions.RequestException as e:
                message.extend([e])
                continue
            else:
                if c.ok:
                    flag = True
                    break
                else:
                    message.extend([c.text])
        if flag:
            return flag, message
        else:
            return flag, message


class IProXChecker(GenericChecker):
    def __init__(self, times=3):
        self.retrievers = None
        self.times = times
        self.flag = {
            'project': 0,
            'protein': 0,
            'peptide': 0
        }
        self.message = {
            'project': [],
            'protein': [],
            'peptide': []
        }

    def _implement(self):
        retrievers = retriever.GenericIProXRetriever.__subclasses__()

        self.retrievers = [r() for r in retrievers]

    def _ans(self):
        print("\niProX Connectivity Report")
        print(f"\tProject Retriever: {bool(self.flag['project'])}")
        print(f"\tProtein Retriever: {bool(self.flag['protein'])}")
        print(f"\tPeptide Retriever: {bool(self.flag['peptide'])}")

    def check(self):
        self._implement()

        flag = None
        message = []

        for r in self.retrievers:
            if isinstance(r, retriever.IProXProjectRetriever):
                flag, message = self._project_check(project_retriever=r)
                if flag:
                    self.flag['project'] = 1
                self.message['project'] = message

            elif isinstance(r, retriever.IProXProteinRetriever):
                flag, message = self._protein_check(protein_retriever=r)
                if flag:
                    self.flag['protein'] = 1
                self.message['protein'] = message

            elif isinstance(r, retriever.IProXPeptideRetriever):
                flag, message = self._peptide_check(peptide_retriever=r)
                if flag:
                    self.flag['peptide'] = 1
                self.message['peptide'] = message
        self._ans()

    def _project_check(self, project_retriever):
        t = project_retriever
        flag = False
        message = []

        print(f"\nCheck the connectivity of the iProX Project API")
        for _ in tqdm.tqdm(range(self.times)):
            try:
                c = requests.get(t.api, params=t.example)
            except requests.exceptions.RequestException as e:
                message.extend([e])
                continue
            else:
                if c.ok:
                    flag = True
                    break
                else:
                    message.extend([c.text])
        if flag:
            return flag, message
        else:
            return flag, message

    def _protein_check(self, protein_retriever):
        t = protein_retriever
        flag = False
        message = []

        print(f"\nCheck the connectivity of the iProX Protein API")
        for _ in tqdm.tqdm(range(self.times)):
            try:
                payloads = {'proteinAccession': t.example,
                            'pageNumber': '1',
                            'pageSize': '5',
                            'resultType': 'compact'}
                c = requests.get(t.api, params=payloads)
            except requests.exceptions.RequestException as e:
                message.extend([e])
                continue
            else:
                if c.ok:
                    flag = True
                    break
                else:
                    message.extend([c.text])
        if flag:
            return flag, message
        else:
            return flag, message

    def _peptide_check(self, peptide_retriever):
        t = peptide_retriever
        flag = False
        message = []

        print(f"\nCheck the connectivity of the PRIDE Peptide API")
        for _ in tqdm.tqdm(range(self.times)):
            try:
                payloads = {'peptideSequence': t.example,
                            'pageNumber': '1',
                            'pageSize': '200',
                            'resultType': 'compact'}
                c = requests.get(t.api, params=payloads)
            except requests.exceptions.RequestException as e:
                message.extend([e])
                continue
            else:
                if c.ok:
                    flag = True
                    break
                else:
                    message.extend([c.text])
        if flag:
            return flag, message
        else:
            return flag, message