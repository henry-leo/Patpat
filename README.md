# Patpat

Patpat is public proteomics dataset search framework that
only requires protein identifiers to be passed in to search for relevant datasets

## Base Usage

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
patpat_env
    |-- logs
    |-- tmp
    |-- result
    |-- proteome
```

Get the search configs via the QueryHub:
```python
identifier_ = 'P23950'
q = hub.QueryHub()
q.identifier = identifier_
q.simple_query()

conf_ = q.get_query_config()
```
Set up Mappers for MapperHub, search and get results:
```python
mappers_ = [mapper.PrideMapper(), mapper.IProXMapper()]

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


For more information, see the Wiki (https://github.com/henry-leo/Patpat/wiki). 

(But the Wiki is currently under development... (T⌓T)




