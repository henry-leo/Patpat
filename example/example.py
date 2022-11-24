"""Basic Usage Example"""

# Load Patpat package
import patpat.hub as hub
import patpat.mapper as mapper
import patpat.utility as utility

# Create runtime environment and initiate UniProt's proteome catalog.
utility.init()
utility.initiate_uniprot_proteome_catalog()

# Get the search configs via the QueryHub
identifier_ = 'P23950'
q = hub.QueryHub()
q.identifier = identifier_
q.simple_query()

# Set up Mappers and the configs for MapperHub, search and get results.
conf_ = q.get_query_config()
mappers_ = [mapper.PrideMapper()]

m = hub.MapperHub(config=conf_,
                  mappers=mappers_,
                  )
m.mapping()

result_ = m.export()
