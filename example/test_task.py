import patpat.hub as hub
import patpat.mapper as mapper
import patpat.utility as utility

utility.init()
utility.initiate_uniprot_proteome_catalog()

identifier_ = 'P23950'
q = hub.QueryHub()
q.identifier = identifier_
q.simple_query()
conf_ = q.get_query_config()

mappers_ = [mapper.PrideMapper(), mapper.IProXMapper()]

m = hub.MapperHub(config=conf_,
                  mappers=mappers_,
                  )
m.mapping()

result_ = m.export()
