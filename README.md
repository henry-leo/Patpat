# Patpat

Patpat is public proteomics dataset search framework that
only requires protein identifiers to be passed in to search for relevant datasets

## Base Usage

Load Patpat package and create running envs.

```Python
import patpat.hub as hub
import patpat.mapper as mapper
import patpat.querier as querier

hub.init()
```
Get the search configs via the Querier.
```python
identifier_ = 'P23950'
q = hub.QueryHub()
q.simple_query()
conf_ = q.get_query_config()
```
Set up mappers for Hub module, search and get results
```python
mappers_ = [mapper.PrideMapper(), mapper.IProXMapper()]

m = hub.MapperHub(config=conf_,
                  mappers=mappers_,
                  task=None
                  # task=[your task's uuid]
                  )
m.mapping()

result_ = m.export()
```
In its current version, Patpat supports both PRIDE and iProX databases. In addition, 
Patpat is an extensible framework and users are encouraged to extend it with databases of interest to Patpat or
to build their processes. 

For more information, see the Wiki. 

But the Wiki is currently under development... (T⌓T)




