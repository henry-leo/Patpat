# Patpat
Patpat means Proteomics Aiders Telescope, it is a public proteomics dataset search framework that
only requires protein identifiers to be passed in to search for relevant datasets.

## Quickly Use
Load Patpat package and create runtime environment:

```Python
import patpat.hub as hub
import patpat.mapper as mapper
import patpat.utility as utility

utility.init()
utility.initiate_uniprot_proteome_catalog()
```
Directory structure of the runtime environment is as follows:
```
patpat_env/
    |-- logs/
    |-- tmp/
    |-- result/
    |-- proteome/
        |-- UP_README_yyyy-mm-dd
```
As an example, take the mouse protein P23950 (UniProt) and search for the peptide to be searched by QueryHub
```python
identifier_ = 'P23950'
q = hub.QueryHub()
q.identifier = identifier_
q.simple_query()
```
Having checked that the corresponding FASTA file for *Mus musculus* does not exist locally, consider obtaining from UniProt:
```
Choose local peptide search.
The Mus musculus UP000000589 proteome file was not found locally.
Do you want to download it?(y/n)
```
Get the search configs:
```python
conf_ = q.get_query_config()
```
Set up Mappers based on connectivity and add MapperHub's configuration, search and get results.:
```python
mappers_ = hub.CheckerHub().checker()

m = hub.MapperHub(config=conf_,
                  mappers=mappers_,
                  )
m.mapping()

result_ = m.export()
```

Result files store in ```patpat_envs/result/<task_uuid>```, you can find ```<task_uuid>``` by ```m.config```

In its current version, Patpat supports both PRIDE and iProX databases. In addition,
Patpat is an extensible framework and users are encouraged to extend it with databases of interest to Patpat or
to build their processes.

For more information, see the [Wiki](https://github.com/henry-leo/Patpat/wiki).
