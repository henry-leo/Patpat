"""This module contains patpat_env's self-built methods. 本模块包含patpat自建的方法。"""

import collections
import os
import re
import time
import json
import gzip
from urllib import parse
from ftplib import FTP

import requests
import tqdm
import pandas as pd


def init():
    """Used to create package runtime environment. 用于创建包运行环境。

    Directory structure of the runtime environment:
        patpat_env/
            |-- logs/
            |-- tmp/
            |-- result/
            |-- proteome/
    """

    def init_subdir():
        """Used to create subdirectory. 用于创建子目录"""
        dir_list = ['logs', 'tmp', 'result', 'proteome']
        for dir_ in dir_list:
            try:
                os.makedirs(f'patpat_env/{dir_}')
            except FileExistsError:
                print(f'\033[0;31mWarning! This directory already exists. 警告！该目录已经存在。[patpat_env/{dir_}]\033[0m')

    if os.path.exists('patpat_env/'):
        print(f'\033[0;31mWarning! This directory already exists. 警告！该目录已经存在。[patpat_env/]\033[0m')
        init_subdir()
    else:
        os.mkdir('patpat_env/')
        init_subdir()
    print('Runtime environment completed.')


def flatten(item, ignore_types=(str, bytes, set)):
    for x in item:
        if isinstance(x, collections.abc.Iterable) and not isinstance(x, ignore_types):
            yield from flatten(x)
        else:
            yield x


def initiate_uniprot_proteome_catalog():
    """Initiate UniProt's proteome catalog

    Returns:
        patpat_env/
            |-- proteome/
                |-- UP_README_yyyy-mm-dd
    """

    ftp = FTP('ftp.uniprot.org')
    ftp.login()
    ftp.cwd('./pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/')
    download_time = time.strftime('%Y-%m-%d', time.localtime())
    with open(f'patpat_env/proteome/UP_README_{download_time}', mode='wb') as f:
        ftp.retrbinary(f'RETR README', f.write)
    ftp.close()
    print('Uniprot_proteome_list has been initiated.')


def create_list_uniprot_proteome():
    """Create UniProt Proteome Catalog"""
    os.chdir('patpat_env/proteome/')
    file = [i for i in os.listdir() if re.search('UP_README', i)][0]
    with open(file, mode='r') as f:
        flag = ''
        while flag != 'Proteome_ID	Tax_ID	OSCODE	SUPERREGNUM	#(1)	#(2)	#(3)	Species Name\n':
            flag = f.readline()
        r = [re.split('\t', flag[:-1])]
        while True:
            flag = f.readline()
            try:
                flag[1]
            except IndexError:
                break
            else:
                r.append(re.split('\t', flag[:-1]))
        df = pd.DataFrame(r[1:], columns=r[0])
        os.chdir('../..')
    return df


def list_local_proteome():
    """List local proteome fasta"""
    proteome = os.listdir('patpat_env/proteome')
    return [re.search('(?<=uniprot-proteome_).*(?=.fasta)', p).group()
            for p in proteome if re.search('(?<=uniprot-proteome_).*(?=.fasta)', p)]


def download_uniprot_opg(taxonomy_id):
    """Download OPG Proteome fasta from UniProt"""
    df = create_list_uniprot_proteome()
    target = df.loc[df['Tax_ID'] == str(taxonomy_id)]
    tax_division = target['SUPERREGNUM'].values[0].capitalize()
    up_id = target['Proteome_ID'].values[0]

    ftp = FTP('ftp.uniprot.org')
    ftp.login()
    ftp.cwd(f'./pub/databases/uniprot/current_release/'
            f'knowledgebase/reference_proteomes/{tax_division}/{up_id}')

    file_name = up_id + '_' + target['Tax_ID'].values[0] + '.fasta'
    download_time = time.strftime('%Y-%m-%d', time.localtime())

    with open(f'patpat_env/proteome/{file_name}.gz',
              mode='wb') as f_gz:
        ftp.retrbinary(f'RETR {file_name}.gz', f_gz.write)
    ftp.close()

    gz = gzip.GzipFile(f'patpat_env/proteome/{file_name}.gz', mode='rb')
    with open(f'patpat_env/proteome/uniprot-proteome_{up_id}_OGP_{download_time}.fasta', mode='wb') as fo:
        fo.write(gz.read())
    gz.close()

    print(f'{file_name} download is complete.')
    os.remove(f'patpat_env/proteome/{file_name}.gz')
    return f'patpat_env/proteome/uniprot-proteome_{up_id}_OGP_{download_time}.fasta'


def download_uniprot_proteome(identifier: str):
    """Download proteome fasta from UniProt database

    Args:
        identifier: str.
                The 'proteome' filter value has invalid format, it should match the regular expression UP[0-9]{9}
    """
    url = f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3A{identifier}%29%29'

    protein_fasta = requests.get(url).text

    # FASTA name format: uniprot-proteome_UPxxxxxxxxx_yyyy-mm-dd.fasta
    download_time = time.strftime('%Y-%m-%d', time.localtime())
    fasta_name = '_'.join(['patpat_env/proteome/uniprot-proteome',
                           identifier,
                           download_time]) + '.fasta'
    with open(fasta_name, mode='w') as fo:
        fo.write(protein_fasta)
    print('Download complete.')


def pagination_download_uniprot_proteome(identifier: str, reviewed: bool = True, size: int = 500):
    """Under Development..."""
    if reviewed:
        reviewed = 'true'
    else:
        reviewed = 'false'
    session = requests.Session()
    base_url = f'https://rest.uniprot.org/uniprotkb/search?&' \
               f'query=proteome:{identifier} AND (reviewed:{reviewed})&format=fasta&size={size}'
    response = session.get(base_url)
    total = int(response.headers['x-total-results'])
    fasta = [response.text]

    for _ in tqdm.tqdm(range(size, total, size)):
        try:
            response = session.get(response.links['next']['url'])
            fasta.extend([response.text])
        except KeyError:
            print(1)
            break

    fasta = '\n'.join(fasta)
    return fasta


def url_split(url_with_params):
    """Extract parameters from url 从网址中提取参数

    Args:
        url_with_params: URL with parameters 带参数的网址

    Returns:
        base_url: URL without parameters 没有参数的网址
        base_params: Separated parameters 被分离的参数
    """
    base_url = url_with_params[:re.search("\\?", url_with_params).start()]
    params = re.split("&", url_with_params[re.search("\\?", url_with_params).end():])
    params = [re.split("=", i) for i in params]
    base_params = dict()
    for parameter in params:
        base_params[parameter[0]] = parse.unquote(parameter[1])
    return base_url, base_params


def usi_split(usi):
    """Splitting USI"""
    component = dict()

    collection = re.search('(?<=mzspec:)\\S*?(?=:)', usi)

    component['collection'] = collection.group()

    split_point1 = collection.end() + 1
    split_point2 = re.search('(index|scan|trace|nativeId):\\d*', usi).end()

    component['msRun+index'] = usi[split_point1:split_point2]
    component['interpretation'] = usi[split_point2 + 1:]

    return component


def usi_detect(usi):
    """Detect if the string is USI"""
    try:
        c = usi_split(usi)
    except AttributeError:
        return False
    else:
        if re.search("pxd", c['collection'], re.I) and c['msRun+index'] and c['interpretation']:
            return True
        else:
            return False


def pride_usi_split(usi):
    """Splitting PRIDE USI"""
    component = dict()

    collection = re.search('\\S*?(?=:)', usi)

    component['collection'] = collection.group()

    split_point1 = collection.end() + 1
    split_point2 = re.search('[OPQ]\\d[A-Z\\d]{3}\\d|[A-NR-Z]\\d([A-Z][A-Z\\d]{2}\\d){1,2}:\\d*', usi).end()

    component['msRun+index'] = usi[split_point1:split_point2 - 1]
    component['interpretation'] = usi[split_point2:]

    return component


def pride_usi_detect(usi):
    """Detect if the string is PRIDE USI"""
    try:
        c = pride_usi_split(usi)
    except AttributeError:
        return False

    else:
        if re.search("pxd", c['collection'], re.I) and c['msRun+index'] and c['interpretation']:
            return True
        else:
            return False


def proteome_file_info_split(file_list: list):
    file_info = [[i,
                  re.match('UP\\d*', i).group(),
                  re.search('(?<=_)\\d{4}-\\d{2}-\\d{2}', i).group()
                  ] for i in file_list]
    file_info = [[i[0], i[1], time.strptime(i[2], '%Y-%m-%d')] for i in file_info]
    return file_info


def uniprot_proteome_id_detect(taxonomy_id):
    """Find proteome files based on UniProt taxonomy ID"""
    df = create_list_uniprot_proteome()
    target = df.loc[df['Tax_ID'] == str(taxonomy_id)]
    try:
        up_id = target['Proteome_ID'].values[0]
    except IndexError:
        raise IndexError('Taxonomy is not found.')
    return up_id


def proteome_file_detect(taxonomy_id):
    """Detection of the local presence of the corresponding proteomic fasta"""
    up_id = uniprot_proteome_id_detect(taxonomy_id)
    file_list = list_local_proteome()
    file_info = proteome_file_info_split(file_list)

    target_files = [i for i in file_info if i[1] == up_id]
    if not target_files:
        return False
    else:
        return True


def singleton(cls):
    """Singleton Pattern"""
    _instance = {}

    def inner():
        if cls not in _instance:
            _instance[cls] = cls()
        return _instance[cls]

    return inner


def get_result_from_file(task=None):
    """Get Patpat search result from patpat_env/result/<task>/result.json

    Args:
        task: uuid

    Returns:
        dict, project metadata
    """
    if task:
        with open(f'patpat_env/result/{task}/result.json') as f:
            data_json = ''.join(f.readlines())
            data_dict = json.loads(data_json)
        f.close()

        return data_dict

    else:
        print('Need to input task uuid.')
