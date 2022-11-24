"""This module contains classes for interacting with external proteomics databases.此模块包含与外部蛋白质组学数据库交互的类。

The classes in this module use three levels of inheritance. First,(1) all subclasses need to inherit from
the Retriever base class, second,(2) depending on the database to be targeted, a base class  needs to be
created for that database, and finally, the subclass that collects the specified data from that database needs
to inherit from that database's base class (3). For example, the base class of PRIDE database is named
GenericPrideRetriever (2), which inherits Retriever (1) base class, and the subclass that collects protein
information from PRIDE is named PrideProteinRetriever (3).
本模块的类使用了三层继承。首先，所有的子类都需要继承Retriever基类（1），其次，根据要针对的数据库，需要为该数据库创建一个基类（2），
最后，从该数据库收集指定数据的子类（3）需要继承该数据库的基类。例如，PRIDE数据库的基类名为GenericPrideRetriever（2），
它继承了Retriever基类（1），从PRIDE收集蛋白质信息的子类名为PrideProteinRetriever（3）

    Typical usage example:
    典型用法示例：

    import patpat_env.retriever as retriever

    r = retriever.PrideProteinRetriever()
    r.request_word = "Q9CWY9"
    r.get_payloads_on_web()
    r.payloads["pageSize"] = 200 # Optional: Parameters can be modified 可选：参数可修改
    r.retrieve()
"""

import re
import time
import logging
import random

import requests

from . import utility


class Retriever(object):
    """Base class for accessing proteomics database URL. 用于访问蛋白质组数据库网址的基类。

    Common proteomics databases are PRIDE, PeptideAtlas, MassIVE, jPOST, iProx, Panorama Public.
    This base class gives preliminary interface definitions to build base classes against databases.
    常见的蛋白质组学数据库有 PRIDE、PeptideAtlas、MassIVE、jPOST、iProx、Panorama Public。
    这个基类给出了初步的接口定义，以便针对数据库建立基类。

    """

    def retrieve(self):
        """The call makes the class Retriever  run and structures the incoming data.
           该调用使Retriever类运行，并对传入的数据进行结构化。

        This is the most important method of the class, so please design it according to your actual needs!
        It is worth noting that this method is required whether the class is called directly,
        or subsequently through another class.
        这是该类中最重要的方法，所以请根据你的实际需要来设计它。值得注意的是，不管是直接调用该类，还是随后通过其他类调用，都需要这个方法。
        """
        raise NotImplementedError

    def url_requester(self):
        """For direct access to the database URL. 用于直接访问数据库URL。

        Access parameters can be set according to actual needs.
        可以根据实际需要设置访问参数。
        """
        raise NotImplementedError

    def get_payloads_on_web(self):
        """Used to get the list of optional parameters from the server of the database.用于从数据库的服务器中获取可选参数列表。

        Optionally, according to the actual needs of the database using.
        可选的，根据数据库实际需要采用。
        """
        raise NotImplementedError


class GenericPrideRetriever(Retriever):
    """Generic base class for PRIDE database. 针对PRIDE数据库的通用基类。

    """

    def __init__(self):
        self._request_word = ''  # Keywords to be requested 需要请求的关键词
        self.headers = None  # HTTP protocol headers HTTP协议的标头
        self.payloads = None  # parameters. 参数列表
        self.url = ''  # URL 网址

    def retrieve(self):
        """Call the methods run by the class PrideRetriever. 调用PrideRetriever类运行的方法。"""
        raise NotImplementedError

    def url_requester(self):
        """Direct access to PRIDE database API. 直接访问PRIDE数据库API的方法。"""
        raise NotImplementedError

    def get_payloads_on_web(self):
        """Access PRIDE API and get the list of parameters of the URL through redirection.
           访问PRIDE API并通过重定向获得URL的参数列表。
        """
        raise NotImplementedError

    def check_end_page(self, url, payloads):
        """Determine if the page is the end page. 确定该页是否为结束页。

        Determine if the page is an end page by the presence of the returned request.json()['_embedded'].
        通过回传的request.json()['_embedded']是否存在，判断该页是否为尾页。

        Args:
            url: URL without parameters 没有参数的网址
            payloads: Parameters to be entered 需要输入的参数

        Returns:
            Bool

        Raises:
            KeyError: Catching this error means that request.json()['_embedded'] does not exist,
             i.e. the page is the end page. 捕捉到这个错误意味着request.json()['_embedded']不存在，也就是说，该页面是结束页。
        """
        page = payloads['page']

        # Almost all of the if-else judgments in this function are used to determine which page number should be
        # written to the log.
        # 这个函数中几乎所有的if-else判断都是用来确定哪一个页码应该被写入日志的。
        if int(page) == 0:
            logging.getLogger('core').debug("----------------------------------------")
            logging.getLogger('core').debug("\tTest if the {} page is the end page.".format(str(int(page))))
        else:
            logging.getLogger('core').debug("----------------------------------------")
            logging.getLogger('core').debug("\tTest if the {} page is the end page.".format(str(int(page) - 1)))

        n = 0
        end = 3
        for n in range(0, end, 1):
            d = requests.get(url, params=payloads)
            time.sleep(random.random())  # Take a break and avoid frequent access. 歇歇，避免频繁访问。
            if d.ok:
                try:
                    _ = d.json()['_embedded']
                except KeyError:
                    pass
                else:
                    logging.getLogger('core').debug("\tthe {} page is not the end page.".format(str(int(page) - 1)))
                    logging.getLogger('core').debug("----------------------------------------")
                    return False
        if n == end - 1 and int(page) != 0:
            logging.getLogger('core').debug("\tthe {} page is the end page.".format(str(int(page) - 1)))
            logging.getLogger('core').debug("----------------------------------------")
            return True
        elif n == end - 1 and int(page) == 0:
            logging.getLogger('core').debug("\tthe {} page is the end page.".format(str(int(page))))
            logging.getLogger('core').debug("----------------------------------------")
            return True


class PrideProjectRetriever(GenericPrideRetriever):
    """Used to collect project information in PRIDE database. 用于收集PRIDE数据库中的项目信息。"""

    def __init__(self):
        super().__init__()
        self.response = dict()  # Collect the returned data through this property. 通过这个属性收集返回的数据。

    @property
    def request_word(self):
        return self._request_word

    @request_word.setter
    def request_word(self, request_word):
        """Used to set request keyword. 用于设置请求关键词。

        Only allows searching for items assigned a number by ProteomeXchange.
        仅允许搜索被ProteomeXchange分配编号的项目。

        Args:
            request_word: keyword 请求关键词
        """
        if re.search("pxd", request_word, re.I) is None:
            raise TypeError('Input error! 输入错误！')
        else:
            self._request_word = request_word

    def retrieve(self):
        """Used to request project data from PRIDE database. 用于向PRIDE数据库请求项目数据。

        Returns:
            None

        Updates:
             self.response

        """
        url_response = self.url_requester()
        logging.getLogger('core').debug('query:{}'.format(self._request_word))

        if url_response.ok:
            url_response = url_response.json()
        else:
            logging.error('No correct return transmission received! 未收到正确回传！')
            url_response = dict()

        time.sleep(.3)
        self.response = url_response
        return self.response

    def url_requester(self):
        """API for direct access to PRIDE database project data. 用于直接访问PRIDE数据库项目数据的API。

        Returns:
            Accessing the API returns a requests object. 访问API返回一个requests对象。

        """
        if self.headers is None:
            self.headers = {'Accept': 'application/json'}

        self.url = "/".join(('https://www.ebi.ac.uk/pride/ws/archive/v2/projects', self._request_word))

        url_response = requests.get(self.url, headers=self.headers, params=self.payloads,
                                    timeout=120  # Avoid blocking
                                    )

        return url_response

    def get_payloads_on_web(self):
        """Occupy a seat 占位"""
        raise NotImplementedError("Don\'t need that in this class. 该类不需要这个！")


class PridePeptideRetriever(GenericPrideRetriever):
    """Used to collect peptide information in PRIDE database. 用于收集PRIDE数据库中的肽段信息。"""

    def __init__(self):
        super().__init__()
        self.response = dict()  # Collect the returned data through this property. 通过这个属性收集返回的数据。

    @property
    def request_word(self):
        return self._request_word

    @request_word.setter
    def request_word(self, request_word):
        """Used to set request keyword. 用于设置请求关键词。

        Args:
            request_word: keyword 请求关键词

        Returns:
            None

        Updates:
            self._request_word

        """
        self._request_word = request_word

    def retrieve(self):
        """Used to request peptide data from PRIDE database. 用于向PRIDE数据库请求肽段数据。

        Returns:
            Updated self.response

        """
        url = self.url
        payloads = self.payloads
        headers = self.headers
        response = dict()
        flag = False

        logging.getLogger('core').info('query:{}'.format(self._request_word))
        while flag is False:
            page = payloads["page"]
            url_response = requests.get(url, headers=headers, params=payloads,
                                        timeout=60  # Avoid blocking
                                        ).json()
            try:
                response["page:{}".format(page)] = {'spectraevidences': url_response['_embedded']['spectraevidences']}
            except KeyError:
                flag = True
                # flag = self.check_end_page(url, payloads)
            else:
                response["page:{}".format(page)]['link'] = url_response["_links"]["self"]["href"]
                logging.getLogger('core').info("page {} is completed.".format(page))
            finally:
                payloads["page"] = str(int(page) + 1)
        self.response = response
        return self.response

    def url_requester(self):
        """API for direct access to PRIDE database peptide data. 用于直接访问PRIDE数据库肽段数据的API。

        Returns:
            Accessing the API returns a requests object. 访问API返回一个requests对象。

        """
        if self.headers is None:
            self.headers = {'Accept': 'application/json'}
        if self.payloads is None:
            self.payloads = {"peptideSequence": self._request_word,
                             "pageSize": 200}

        self.url = "https://www.ebi.ac.uk/pride/ws/archive/v2/spectra"

        url_response = requests.get(self.url, headers=self.headers, params=self.payloads)
        return url_response

    def get_payloads_on_web(self):
        """Access the URL and get the list of parameters of the URL through redirection.
           访问URL并通过重定向获得URL的参数列表。
        """
        url_response = self.url_requester()

        if url_response.ok:
            url_response = url_response.json()
            self.url, self.payloads = utility.url_split(url_response["_links"]['self']['href'])

        return url_response


class PrideProteinRetriever(GenericPrideRetriever):
    """Used to collect protein information in PRIDE database. 用于收集PRIDE数据库中的蛋白质信息。"""

    def __init__(self):
        super().__init__()
        self.response = dict()

    @property
    def request_word(self):
        return self._request_word

    @request_word.setter
    def request_word(self, request_word):
        """Used to set request keyword. 用于设置请求关键词。

        Only UniProt numbers are allowed.
        仅允许UniProt编号

        Args:
            request_word: keyword 请求关键词

        Returns:
            None

        Updates:
            self._request_word

        """
        if re.match('[OPQ]\\d[A-Z\\d]{3}\\d|[A-NR-Z]\\d([A-Z][A-Z\\d]{2}\\d){1,2}',
                    request_word).group() == request_word:
            self._request_word = request_word
        else:
            raise TypeError('Input error! 输入错误！')

    def retrieve(self):
        """Used to request peptide data from PRIDE database. 用于向PRIDE数据库请求蛋白质数据。

        Returns:
            None

        Updates:
            self.response
        """
        url = self.url
        payloads = self.payloads
        headers = self.headers
        response = dict()
        flag = False

        logging.getLogger('core').info('query:{}'.format(self._request_word))
        while flag is False:
            try:
                page = payloads['page']
            except KeyError as e:
                logging.getLogger('core').error(e)
                logging.getLogger('core').error('Reply 1 time')
                self.get_payloads_on_web()
                page = payloads['page']

            url_response = requests.get(url, headers=headers, params=payloads,
                                        timeout=120  # Avoid blocking
                                        ).json()
            try:
                response['page:{}'.format(page)] = {'peptideevidences': url_response['_embedded']['peptideevidences']}
            except KeyError:
                flag = True
                # flag = self.check_end_page(url, payloads)
            else:
                response['page:{}'.format(page)]['link'] = url_response["_links"]["self"]["href"]
                logging.getLogger('core').info('page {} is completed.'.format(page))
            finally:
                payloads['page'] = str(int(page) + 1)
        self.response = response

    def url_requester(self):
        """API for direct access to PRIDE database protein data. 用于直接访问PRIDE数据库蛋白质数据的API。

        Returns:
            Accessing the API returns a requests object. 访问API返回一个requests对象。

        """
        if self.headers is None:
            self.headers = {'Accept': 'application/json'}
        if self.payloads is None:
            self.payloads = {'proteinAccession': self._request_word,
                             'pageSize': 200}

        self.url = "https://www.ebi.ac.uk/pride/ws/archive/v2/peptideevidences"

        url_response = requests.get(self.url, headers=self.headers, params=self.payloads)
        return url_response

    def get_payloads_on_web(self):
        """Access the URL and get the list of parameters of the URL through redirection
           访问URL并通过重定向获得URL的参数列表
        """
        url_response = self.url_requester()

        if url_response.ok:
            url_response = url_response.json()
            self.url, self.payloads = utility.url_split(url_response["_links"]['self']['href'])


class GenericIProXRetriever(Retriever):
    """Generic base class for iProX database. 针对iProX数据库的通用基类。"""

    def __init__(self):
        self._request_word = ''  # Keywords to be requested 需要请求的关键词
        self.headers = None  # HTTP protocol headers HTTP协议的标头
        self.payloads = None  # parameters. 参数列表
        self.url = ''  # URL 网址

    def retrieve(self):
        """Call the methods run by the class IProXRetriever. 调用IProXRetriever类运行的方法。"""
        raise NotImplementedError

    def url_requester(self):
        """Direct access to iProX database API. 直接访问iProX数据库API的方法。"""
        raise NotImplementedError

    def get_payloads_on_web(self):
        """Access iProX API and get the list of parameters of the URL through redirection.
           访问iProX API并通过重定向获得URL的参数列表。
        """
        raise NotImplementedError


class IProXProjectRetriever(GenericIProXRetriever):
    """Used to collect project information in iProX database. 用于收集iProX数据库中的项目信息。"""

    def __init__(self):
        super().__init__()
        self.response = dict()

    @property
    def request_word(self):
        return self._request_word

    @request_word.setter
    def request_word(self, request_word):
        """Used to set request keyword. 用于设置请求关键词。

        Only allows searching for items assigned a number by ProteomeXchange.
        仅允许搜索被ProteomeXchange分配编号的项目。

        Args:
            request_word: keyword 请求关键词
        """
        if re.search("pxd", request_word, re.I) is None:
            raise TypeError("输入错误")
        else:
            self._request_word = request_word

    def retrieve(self):
        """Used to request project data from iProX database. 用于向iProX数据库请求项目数据。

        Returns:
            None

        Updates:
            self.response
        """
        url_response = self.url_requester()
        logging.getLogger('core').debug('query:{}'.format(self._request_word))

        if url_response.ok:
            url_response = url_response.json()
        else:
            logging.error('No correct return transmission received! 未收到正确回传！')
            url_response = dict()
        time.sleep(.3)
        self.response = url_response

    def url_requester(self):
        """API for direct access to iProX database project data. 用于直接访问iProX数据库项目数据的API。

        Returns:
            Accessing the API returns a requests object. 访问API返回一个requests对象。

        """
        if self.headers is None:
            self.headers = {'Accept': 'application/json'}

        self.url = "/".join(('https://www.iprox.cn/proxi/datasets', self._request_word))

        url_response = requests.get(self.url, headers=self.headers, params=self.payloads,
                                    timeout=120  # Avoid blocking
                                    )
        return url_response

    def get_payloads_on_web(self):
        """Occupy a seat 占位"""
        raise NotImplementedError("Don\'t need that in this class. 该类不需要这个！")


class IProXPeptideRetriever(GenericIProXRetriever):
    """Used to collect peptide information in iProX database. 用于收集iProX数据库中的肽段信息。"""

    def __init__(self):
        super().__init__()
        self.response = dict()

    @property
    def request_word(self):
        return self._request_word

    @request_word.setter
    def request_word(self, request_word):
        """Used to set request keyword. 用于设置请求关键词。

        Args:
            request_word: keyword 请求关键词

        Returns:
            None

        Updates:
            self._request_word
        """
        if len(request_word) < 7:
            raise TypeError("输入错误")
        else:
            self._request_word = request_word

    def retrieve(self):
        """Used to request peptide data from iProX database. 用于向iProX数据库请求肽段数据。

        Returns:
            None

        Updates:
            self.response
        """
        url = self.url
        payloads = self.payloads
        headers = self.headers
        response = dict()
        flag = False

        logging.getLogger('core').info(f'query:{self._request_word}')
        while flag is False:
            page = payloads['pageNumber']
            url_response = requests.get(url, headers=headers, params=payloads,
                                        timeout=60  # Avoid blocking
                                        )

            if url_response.ok:
                if url_response.json():
                    response[f'page:{page}'] = url_response.json()
                    logging.getLogger('core').info("page {} is completed.".format(page))
                else:
                    flag = True
            else:
                logging.getLogger('core').error(url_response.json())
                raise requests.exceptions.RequestException(url_response.json())
            payloads['pageNumber'] = str(int(page) + 1)

        self.response = response

    def url_requester(self):
        """API for direct access to iProX database peptide data. 用于直接访问iProX数据库肽段数据的API。

        Returns:
            Accessing the API returns a requests object. 访问API返回一个requests对象。

        """
        if self.headers is None:
            self.headers = {"Accept": "application/json"}
        if self.payloads is None:
            self.payloads = {'peptideSequence': self._request_word,
                             'pageNumber': '1',
                             'pageSize': '200',
                             'resultType': 'compact'}

        self.url = "https://www.iprox.cn/proxi/spectra"

        url_response = requests.get(self.url, headers=self.headers, params=self.payloads)
        return url_response

    def get_payloads_on_web(self):
        """Access the URL and get the list of parameters of the URL through redirection.
           访问URL并通过重定向获得URL的参数列表。
        """
        url_response = self.url_requester()

        if url_response.ok:
            url_response = url_response
            self.url, self.payloads = utility.url_split(url_response.url)

        return url_response


class IProXProteinRetriever(GenericIProXRetriever):
    """Used to collect protein information in iProX database. 用于收集iProX数据库中的蛋白质信息。"""

    def __init__(self):
        super().__init__()
        self.response = dict()

    @property
    def request_word(self):
        return self._request_word

    @request_word.setter
    def request_word(self, request_word):
        """Used to set request keyword. 用于设置请求关键词。

        Only UniProt numbers are allowed.
        仅允许UniProt编号

        Args:
            request_word: keyword 请求关键词

        Returns:
            None

        Updates:
            self._request_word
        """
        if re.match('[OPQ]\\d[A-Z\\d]{3}\\d|[A-NR-Z]\\d([A-Z][A-Z\\d]{2}\\d){1,2}',
                    request_word).group() == request_word:
            self._request_word = request_word
        else:
            raise TypeError('Input error! 输入错误！')

    def retrieve(self):
        """Used to request peptide data from iProX database. 用于向iProX数据库请求蛋白质数据。

        Returns:
            Updated self.response

        """
        url = self.url
        payloads = self.payloads
        headers = self.headers
        response = dict()
        flag = False

        logging.getLogger('core').info(f'query:{self._request_word}')

        while flag is False:
            page = payloads['pageNumber']
            url_response = requests.get(url, headers=headers, params=payloads,
                                        timeout=120  # Avoid blocking
                                        )

            if url_response.ok:
                if url_response.json():
                    response[f'page:{page}'] = url_response.json()
                    logging.getLogger('core').info("page {} is completed.".format(page))
                else:
                    flag = True
            else:
                logging.getLogger('core').error(url_response.json())
                raise requests.exceptions.RequestException(url_response.json())
            payloads['pageNumber'] = str(int(page) + 1)

        self.response = response

    def url_requester(self):
        """API for direct access to iProX database protein data. 用于直接访问iProX数据库蛋白质数据的API。

        Returns:
            Accessing the API returns a requests object. 访问API返回一个requests对象。

        """
        if self.headers is None:
            self.headers = {'Accept': 'application/json'}
        if self.payloads is None:
            self.payloads = {'proteinAccession': self._request_word,
                             'pageNumber': '1',
                             'pageSize': '200',
                             'resultType': 'compact'}

        self.url = 'https://www.iprox.cn/proxi/psms'

        url_response = requests.get(self.url, headers=self.headers, params=self.payloads)
        return url_response

    def get_payloads_on_web(self):
        """Access the URL and get the list of parameters of the URL through redirection
           访问URL并通过重定向获得URL的参数列表
        """
        url_response = self.url_requester()

        if url_response.ok:
            url_response = url_response
            self.url, self.payloads = utility.url_split(url_response.url)
