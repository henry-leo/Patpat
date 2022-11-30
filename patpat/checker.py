"""This module is used to check the connectivity of external proteomics databases 此模块用于确认外部蛋白质组学数据库的连通性

    Typical usage example:
    典型用法示例：

    from patpat import checker
    from patpat import retriever

    c = checker.Checker()
    c.peptide_retrievers = [
        retriever.PridePeptideRetriever(),
        retriever.IProXPeptideRetriever()
    ]
    c.check()
"""
import re

import requests


class Checker:
    def __init__(self,
                 threshold=3
                 ):
        self._peptide_retrievers = None
        self.threshold = threshold

    @property
    def peptide_retrievers(self):
        return self.peptide_retrievers

    @peptide_retrievers.setter
    def peptide_retrievers(self, peptide_retrievers: list):
        self._peptide_retrievers = peptide_retrievers

    def check(self):
        threshold = self.threshold

        for retriever in self._peptide_retrievers:
            retriever.peptide_retriever = 'TCVADESAENCDK'
            n = 1
            while True:
                database = re.search('.*(?=Peptide)', retriever.__class__.__name__).group().upper()
                print(f"Check {database} connectivity: {n}/{threshold}")
                try:
                    request = retriever.get_payloads_on_web()
                except requests.exceptions.RequestException as e:
                    print(f'Error:{e}')

                else:
                    if request:
                        print(f'Connectable ({retriever.url})')
                        break
                    else:
                        print(f'Status Code: {request.status_code}')
                finally:
                    if n >= threshold:
                        print(f"Not connectable, Please check API. ({retriever.url})")
                        break
                    n += 1
