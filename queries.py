#!/usr/bin/env python

from string import Template
from SPARQLWrapper import SPARQLWrapper, JSON
import numpy, terminaltables, csv, string

PREFIXES = """\
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX foaf: <http://xmlns.com/foaf/0.1/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX cco: <http://rdf.ebi.ac.uk/terms/chembl#>
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX cheminf: <http://semanticscience.org/resource/>
PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/>
PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX efo: <http://www.ebi.ac.uk/efo/>
PREFIX biosd-terms: <http://rdf.ebi.ac.uk/terms/biosd/>
PREFIX pav: <http://purl.org/pav/2.0/>
PREFIX prov: <http://www.w3.org/ns/prov#>
PREFIX oac: <http://www.openannotation.org/ns/>
"""

# SPARQL query template
QUERY = Template('$prefixes\n$query')


def retrieve_compound_labels(chemblid):
    """Retrieve all the compound names and idetifiers for the specified compound. Only
    unique names are returned.

    Args:
        chemblId (str): Specifies chemblId of the compound.
    
    Returns:
        list: A numpy array of compound names.
    """
    
    query = Template("""\
    SELECT DISTINCT ?compound
    WHERE {
        ?molecule <http://rdf.ebi.ac.uk/terms/chembl#chemblId> "$chemblid" .
        {
            ?molecule sio:SIO_000008 ?sid .
            {
                ?sid a cheminf:CHEMINF_000059 . # get inchikey
            }
            UNION
            {
                ?sid a cheminf:CHEMINF_000113 . # get inchi
            }
            UNION
            {
                ?sid a cheminf:CHEMINF_000018 . # get smiles
            }
            ?sid sio:SIO_000300 ?compound . # take the values of these compound identifiers
        }
        UNION
        {
            ?molecule skos:altLabel ?compound . # get all compound alternative labels
        }
        UNION
        {
            ?molecule skos:prefLabel ?compound . # get all compound preffered labels
        }
        UNION
        {
            ?molecule rdfs:label ?compound . # get all compound labels
        }
    }
    LIMIT $limit
    OFFSET $offset
    """)
    # I noticed a strange effect here. If I use the limit or offset on the ChEMBL public
    # SPARQL output on the web it didn't work. However it works if I use it from python.
    
    # Retrieve the results from ChEMBL
    limit = 50
    offset = 0
    step = 50
    names = []
    sparql = SPARQLWrapper("https://www.ebi.ac.uk/rdf/services/chembl/sparql")
    while True:
        sparql.setQuery(QUERY.substitute(
            prefixes = PREFIXES, query=query.substitute(chemblid = chemblid,
            limit = limit, offset = offset)))
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        if len(results['results']['bindings']) == 0:
            break
        elif len(results['results']['bindings']) > 0:
            for i in results['results']['bindings']:
                names.append(i['compound']['value'])
        offset += step
    
    # Create a list of compounds from the results
    compounds = numpy.unique(numpy.char.lower(numpy.array(names)))
    
    return compounds


def retrieve_samples_and_labels_for_compound(compounds):
    """Retrieve samples and labels for compound and it's alternative names from ChEMBL.
    
    Initially we used the VALUES parameter in SPARQL to define all possible name values
    for a compound and noted that it was was too slow. Instead we use one-by-one queries
    where each compound name is retrieved independently. This also allows us to monitor
    which compound names might cause problems for the queries.
    
    Args:
        compounds (list): A list of possible compound names / identifiers.
    
    Returns:
        dict: A dictionary containing the samples associated with a specific compound
        name.
    """
    
    sparql = SPARQLWrapper("https://www.ebi.ac.uk/rdf/services/biosamples/sparql")
    query = Template("""SELECT DISTINCT ?sample ?lab
    WHERE {
        ?attribute atlasterms:propertyValue ?lab .
        FILTER(CONTAINS(LCASE(STR(?lab)), LCASE(STR("$name")))) .
        ?sample biosd-terms:has-sample-attribute ?attribute .
        ?sample a biosd-terms:Sample .
    }
    LIMIT $limit
    OFFSET $offset
    """)
    
    results_dict = {}
    limit = 50
    step = 50
    for compound in numpy.nditer(compounds):
        offset = 0
        compound = str(compound)
        results_dict[compound] = []
        while True:
            print("Compound name: " + compound)
            print("\n"+QUERY.substitute(
                prefixes=PREFIXES, query=query.substitute(name=compound,
                limit = limit, offset = offset)))
            sparql.setQuery(QUERY.substitute(
                prefixes = PREFIXES, query=query.substitute(name=compound,
                limit = limit, offset = offset)))
            sparql.setReturnFormat(JSON)
            results = sparql.query().convert()
            if len(results['results']['bindings']) == 0:
                break
            elif len(results['results']['bindings']) > 0:
                for i in results['results']['bindings']:
                    # The first element in the array is the label and the second is the sample
                    results_dict[compound].append([i['lab']['value'], i['sample']['value']])
                offset += step
    
    return results_dict


def retrieve_samples_by_chebiid(chebiid):
    """Retrieve samples by chebiid
    
    Args:
        chebiid (str): Specifies ChEBI id of the compound.
    
    Returns:
        list: A list of samples related to the specified ChEBI id.
    """
    chebiid = chebiid.replace(':', '_')
    sparql = SPARQLWrapper("https://www.ebi.ac.uk/rdf/services/biosamples/sparql")
    query = Template("""SELECT DISTINCT ?sample ?attribute
    WHERE {
        ?attribute ?b <http://purl.obolibrary.org/obo/$chebiid> .
        ?sample biosd-terms:has-sample-attribute ?attribute .
        ?sample a biosd-terms:Sample .
    }
    LIMIT $limit
    OFFSET $offset
    """)
    
    limit = 50
    step = 50
    offset = 0
    samples = []
    
    while True:
        print("\n"+QUERY.substitute(
            prefixes=PREFIXES, query=query.substitute(chebiid = chebiid,
            limit = limit, offset = offset)))
        sparql.setQuery(QUERY.substitute(
            prefixes = PREFIXES, query=query.substitute(chebiid = chebiid,
            limit = limit, offset = offset)))
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        if len(results['results']['bindings']) == 0:
            break
        elif len(results['results']['bindings']) > 0:
            for i in results['results']['bindings']:
                # The first element in the array is the label and the second is the sample
                samples.append([i['sample']['value'], i['attribute']['value']])
            offset += step
    
    return samples


def retrieve_samples_by_chebiid_and_molar(chebiid):
    """Retrieve samples by chebiid and molarity unit
    
    Args:
        chebiid (str): Specifies ChEBI id of the compound.
    
    Returns:
        list: A list of samples related to the specified ChEBI id
        and molar unit.
    """
    
    # glucosamine - http://purl.obolibrary.org/obo/CHEBI_5417
    chebiid = chebiid.replace(':', '_')
    sparql = SPARQLWrapper("https://www.ebi.ac.uk/rdf/services/biosamples/sparql")
    query = Template("""SELECT DISTINCT ?sample ?compound
    WHERE {
        ?subunit rdfs:subClassOf obo:UO_0000061 .
        FILTER ( ?subunit != efo:EFO_0002902 ) # exclude this unit
        ?molarity a ?subunit .
        ?compound a <http://purl.obolibrary.org/obo/$chebiid> .
        ?sample biosd-terms:has-sample-attribute ?compound .
        ?sample a biosd-terms:Sample .
        ?sample biosd-terms:has-sample-attribute ?compound, ?molarity .
    }
    LIMIT $limit
    OFFSET $offset
    """)
    
    limit = 50
    step = 50
    offset = 0
    samples = []
    
    while True:
        print("\n"+QUERY.substitute(
            prefixes=PREFIXES, query=query.substitute(chebiid = chebiid,
            limit = limit, offset = offset)))
        sparql.setQuery(QUERY.substitute(
            prefixes = PREFIXES, query=query.substitute(chebiid = chebiid,
            limit = limit, offset = offset)))
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        if len(results['results']['bindings']) == 0:
            break
        elif len(results['results']['bindings']) > 0:
            for i in results['results']['bindings']:
                # The first element in the array is the label and the second is the sample
                samples.append([i['sample']['value'], i['compound']['value']])
            offset += step
    
    return samples


def retrieve_samples_by_chebiid_and_celline(chebiid):
    """Retrieve samples by chebiid and cellline ontology
    
    Args:
        chebiid (str): Specifies ChEBI id of the compound.
    
    Returns:
        list: A list of samples related to the specified ChEBI id
        and cell-line.
    """
    
    # glucosamine - http://purl.obolibrary.org/obo/CHEBI_5417
    chebiid = chebiid.replace(':', '_')
    sparql = SPARQLWrapper("https://www.ebi.ac.uk/rdf/services/biosamples/sparql")
    query = Template("""SELECT DISTINCT ?sample ?cellline
    WHERE {
        ?compound a <http://purl.obolibrary.org/obo/$chebiid> .
        ?subline rdfs:subClassOf* <http://www.ebi.ac.uk/efo/EFO_0000322> . #celline
        ?cellline a ?subline .
        ?sample a biosd-terms:Sample .
        ?sample biosd-terms:has-sample-attribute ?compound, ?cellline .
    }
    LIMIT $limit
    OFFSET $offset
    """)
    
    limit = 50
    step = 50
    offset = 0
    samples = []
    
    while True:
        print("\n"+QUERY.substitute(
            prefixes=PREFIXES, query=query.substitute(chebiid = chebiid,
            limit = limit, offset = offset)))
        sparql.setQuery(QUERY.substitute(
            prefixes = PREFIXES, query=query.substitute(chebiid = chebiid,
            limit = limit, offset = offset)))
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        if len(results['results']['bindings']) == 0:
            break
        elif len(results['results']['bindings']) > 0:
            for i in results['results']['bindings']:
                # The first element in the array is the label and the second is the sample
                samples.append([i['sample']['value'], i['cellline']['value']])
            offset += step
    
    return samples


def retrieve_samples_by_chebiid_and_strain(chebiid):
    """Retrieve samples by chebiid and strain ontology
    
    Args:
        chebiid (str): Specifies ChEBI id of the compound.
    
    Returns:
        list: A list of samples related to the specified ChEBI id
        and strain ontology.
    """
    
    # glucosamine - http://purl.obolibrary.org/obo/CHEBI_5417
    chebiid = chebiid.replace(':', '_')
    sparql = SPARQLWrapper("https://www.ebi.ac.uk/rdf/services/biosamples/sparql")
    query = Template("""SELECT DISTINCT ?sample ?cellline
    WHERE {
        ?compound a <http://purl.obolibrary.org/obo/$chebiid> .
        ?strain rdfs:subClassOf* <http://www.ebi.ac.uk/efo/EFO_0001329> . #strain
        ?sample a biosd-terms:Sample .
        ?sample biosd-terms:has-sample-attribute ?compound, ?strain .
    }
    LIMIT $limit
    OFFSET $offset
    """)
    
    limit = 50
    step = 50
    offset = 0
    samples = []
    
    while True:
        print("\n"+QUERY.substitute(
            prefixes=PREFIXES, query=query.substitute(chebiid = chebiid,
            limit = limit, offset = offset)))
        sparql.setQuery(QUERY.substitute(
            prefixes = PREFIXES, query=query.substitute(chebiid = chebiid,
            limit = limit, offset = offset)))
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        if len(results['results']['bindings']) == 0:
            break
        elif len(results['results']['bindings']) > 0:
            for i in results['results']['bindings']:
                # The first element in the array is the label and the second is the sample
                samples.append([i['sample']['value'], i['cellline']['value']])
            offset += step
    
    return samples


def samples_associated_with_molarity(samplesandlabels):
    """Retrieve all samples associated with a molarity unit
    
    Args:
        samplesandlabels (dict): A dictionary of all the compound
        labels and associated samples of the compound
    
    Returns:
        list: A list of samples related to all compound labels
    """
    samples = []
    for compound in samplesandlabels:
        allsamples = samplesandlabels[compound]
        for sample in allsamples:
            samples.append(sample[1])


def no_of_samples(samplesandlabels):
    """Calculate the number of samples for the given compound for all its labels
    """
    tabledata = [['Compound name', 'No of samples']]
    for compound in samplesandlabels.keys():
        tabledata.append([compound, len(samplesandlabels[compound])])
    return tabledata


def no_of_samples_beautiful(samplesandlabels):
    """Prints out the number of samples of a compound and all its labels in a table
    format
    """
    table = terminaltables.AsciiTable(no_of_samples(samplesandlabels))
    print(table.table)


def check_duplicity_of_samples(samplesandlabels):
    """Creates a intersection table between compound labels and the number of samples
    that overlap between the labels
    """
    compoundrows = []
    headers = ['Compound Name']
    for compound1 in samplesandlabels:
        headers.append(compound1)
        compoundrow = []
        compoundrow.append(compound1)
        for compound2 in samplesandlabels:
            if len(samplesandlabels[compound1]) == 0 or len(samplesandlabels[compound2]) == 0:
                compoundrow.append(0)
            else:
                samples1 = [x[1] for x in samplesandlabels[compound1]]
                samples2 = [x[1] for x in samplesandlabels[compound2]]
                overlap = len(set(samples1).intersection(samples2))
                compoundrow.append(overlap)
        compoundrows.append(compoundrow)
    table = []
    table.append(headers)
    table.extend(compoundrows)
    return table


def check_chebisamples_duplicity(chebisamples):
    """Checks if there is no duplicity in the samples retrieved using ChEBI id
    """
    samples = [x[0] for x in chebisamples]
    return(samples)


def write2table(header = None, data = None, filename = None):
    """Write the data in a table format. The data is a list of lists where each sublist
    represents a row.
    
    Args:
        header (list): A list of header values
        data (list): A list of rows where each row is a list of elements
        filename (str): Specifies the filename and location
    
    Returns:
        list: A list of samples related to the specified ChEBI id.
    """
    output = []
    if header != None:
        output.append(header)
    output.extend(data)
    with open(filename, 'w') as f:
        csv.writer(f).writerows(output)


if __name__ == "__main__":
    print ("Program started")
    # rosiglitazone
    rosiglitazone = retrieve_compound_labels("CHEMBL121")
    rosiglitazone_samplesandlabels = retrieve_samples_and_labels_for_compound(rosiglitazone)
    rosiglitazone_no_of_samples = no_of_samples(rosiglitazone_samplesandlabels)
    no_of_samples_beautiful(rosiglitazone_samplesandlabels)
    rosiglitazone_duplicity = check_duplicity_of_samples(rosiglitazone_samplesandlabels)
    rosiglitazone_chebisamples = retrieve_samples_by_chebiid('CHEBI:50122')
    retrieve_samples_by_chebiid_and_molar('CHEBI:50122')
    rosiglitazone_celllines = retrieve_samples_by_chebiid_and_celline('CHEBI:50122')
    rosiglitazone_strains = retrieve_samples_by_chebiid_and_strain('CHEBI:50122')
    # aspirin
    aspirin = retrieve_compound_labels("CHEMBL25")
    aspirin_samplesandlabels = retrieve_samples_and_labels_for_compound(aspirin)
    aspirin_no_of_samples = no_of_samples(aspirin_samplesandlabels)
    no_of_samples_beautiful(aspirin_samplesandlabels)
    aspirin_duplicity = check_duplicity_of_samples(aspirin_samplesandlabels)
    aspirin_chebisamples = retrieve_samples_by_chebiid('CHEBI:15365')
    retrieve_samples_by_chebiid_and_molar('CHEBI:15365')
    aspirin_celllines = retrieve_samples_by_chebiid_and_celline('CHEBI:15365')
    aspirin_strains = retrieve_samples_by_chebiid_and_strain('CHEBI:15365')
    # valproic acid
    valproate = retrieve_compound_labels("CHEMBL109")
    valproate_samplesandlabels = retrieve_samples_and_labels_for_compound(valproate)
    valproate_no_of_samples = no_of_samples(valproate_samplesandlabels)
    no_of_samples_beautiful(valproate_samplesandlabels)
    valproate_duplicity = check_duplicity_of_samples(valproate_samplesandlabels)
    valproate_chebisamples = retrieve_samples_by_chebiid('CHEBI:39867')
    retrieve_samples_by_chebiid_and_molar('CHEBI:39867')
    valproate_celllines = retrieve_samples_by_chebiid_and_celline('CHEBI:39867')
    valproate_strains = retrieve_samples_by_chebiid_and_strain('CHEBI:39867')
    # glucosamine
    retrieve_samples_by_chebiid_and_molar('CHEBI:5417')
    glucosamine_celllines = retrieve_samples_by_chebiid_and_celline('CHEBI:5417')
    glucosamine_strains = retrieve_samples_by_chebiid_and_strain('CHEBI:5417')
