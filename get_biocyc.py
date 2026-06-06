import requests
import sys
import argparse 
import time
import pandas as pd
import numpy as np
import networkx as nx
import xml.etree.ElementTree as ET
import itertools

session = requests.Session()
session.post('https://websvc.biocyc.org/credentials/login/', data={'email':'bh523@ic.ac.uk', 'password':'oRvKrA.%yU'})

BASE_URL_BIOCYC = "https://websvc.biocyc.org/apixml?"
BASE_URL_UNIPROT = "https://rest.uniprot.org/uniprotkb/"
ORG = "ECOLI"
#target_dir = "C:/Users/bella/OneDrive - Imperial College London/FYP/Code/Pathways"

def get_pathways(session, baseURL, organism): #get a list of all pathways for a particular
    r = session.get(
        f"{baseURL}",
         params={
             "fn": "get-class-all-instances",
             "id": f"{organism}:Pathways",
             "detail": "none"
            },
    )   
    r.raise_for_status()

    #print(r.text)
    root = ET.fromstring(r.text)

    pathways = []
    
    for pathway in root.findall("Pathway"):
        pathways.append(pathway.get("frameid"))

    return pathways

    

def get_genes(session, baseURL, organism, pathwayID): #get an np array of all the genes in a given pathway
    
    r = session.get(
        f"{baseURL}fn=genes-of-pathway",
        params={
            "id": f"{organism}:{pathwayID}", #in the format ORGANISM:PATHWAY
            "detail": "none" 
        },
    )
    r.raise_for_status()

    #print(r.text)
    root = ET.fromstring(r.text) #convert XML to JSON

    genes = []
    for gene in root.findall("Gene"): #iterate through XML to find genes only
        genes.append(gene.get("frameid"))

    return genes



def get_uniprot(session, baseURL, gene_id, organism):
    r = session.get(
        f"{baseURL}fn=all-products-of-gene", 
        params={
            "id": f"{organism}:{gene_id}",
            "detail": "full"
        }
    )
    r.raise_for_status()

    root = ET.fromstring(r.text)

    uniprots = []
    for dblink in root.findall('.//dblink'):
        db = dblink.find('dblink-db')
        oid = dblink.find('dblink-oid')
        
        # Check if this is a UNIPROT dblink
        if db is not None and db.text == 'UNIPROT':
            if oid is not None:
                uniprots.append(oid.text)
    
    return uniprots



def get_sequence(baseURL, uniprotID):
    
    params = {
    "fields": [
        "sequence"
    ]
    }
    headers = {
    "accept": "text/plain;format=tsv"
    }
    url = f"{baseURL}/{uniprotID}"

    response = requests.get(url, headers=headers, params=params)
    
    if not response.ok:
        response.raise_for_status()

    sequence = response.text


    return sequence[10:-1]



def get_compounds(session, baseURL_biocyc, organism, pathway):

    EXCLUDE_COMPOUNDS = { ## CREATED MANUALLY
    "PROTON",
    "WATER",
    "ATP",
    "ADP",
    "NAD",
    "NADH",
    "NADP",
    "NADPH",
    "Pi",
    "PPI",
    "CO2",
    "AMP",
    "Protein-Histidines",
    "Protein-pi-phospho-L-histidines",
    "Hpr-Histidine",
    "Hpr-pi-phospho-L-histidines",
    "CO-A",
    "AMMONIUM",
    "NADH-P-OR-NOP",
    "NAD-P-OR-NOP",
    "CARBON-DIOXIDE",
    "FMN",
    "FMNH2",
    "OXYGEN-MOLECULE",
    "Ubiquinones",
    "Ubiquinols",
    "GTP",
    "FORMATE",
    "AMMONIA",
    "GDP",
    "GMP",
    "ACETALD",
    "NITRATE",
    "NITRITE",
    "FE+2",
    "P3I",
    "Alcohols",
    "HYDROGEN-MOLECULE",
    "ETOH",
    "HCO3",
    "CMP",
    "CDP",
    "CTP",
    "UMP",
    "UDP",
    "UTP",
    "UREA",
    "SO3",
    "HS",
    "Menaquinones",
    "Menaquinols",
    "Aldehydes",
    "Alkanesulfonates",
    "TMP",
    "TDP",
    "TTP",
    "IMP",
    "Beta-D-Glucuronides",
    "HMP",
    "NI+2",
    "DUMP",
    "CH4",
    "SEPO3",
    "METOH",
    "ACETOL",
    "CARBON-MONOXIDE",
    "FAD",
    "H2SO3",
    "HCN",
}
    
    r = session.get(
        f"{baseURL_biocyc}fn=compounds-of-pathway", 
        params={
            "id": f"{organism}:{pathway}",
            "detail": "none"
        }
    )
    r.raise_for_status()

    root = ET.fromstring(r.text)

    compounds = []

    for compound in root.findall(".//Compound"):
        
            compound_id = compound.get("frameid")

            if compound_id not in EXCLUDE_COMPOUNDS:

                compounds.append(compound_id)

    return compounds



def get_compound_reactions(session, baseURL_biocyc, organism, compound):

    compound_reactions = []

    r = session.get(
    f"{baseURL_biocyc}fn=reactions-of-compound", 
    params={
        "id": f"{organism}:{compound}",
        "detail": "none"
    }
    )

    r.raise_for_status()
    root = ET.fromstring(r.text)

    for reaction in root.findall(".//Reaction"):

        reaction_id = reaction.get("frameid")

        compound_reactions.append(reaction_id)

    return compound_reactions



def get_gene_reactions(session, baseURL_biocyc, organism, gene):

    gene_reactions = []

    r = session.get(
    f"{baseURL_biocyc}fn=reactions-of-gene", 
    params={
        "id": f"{organism}:{gene}",
        "detail": "none"
    }
    )

    r.raise_for_status()
    root = ET.fromstring(r.text)

    for reaction in root.findall(".//Reaction"):

        reaction_id = reaction.get("frameid")

        gene_reactions.append(reaction_id)


    return gene_reactions




def get_network(session, baseURL_biocyc, organism):

    pathways = get_pathways(session, baseURL_biocyc, organism)

    G = nx.Graph()

    all_genes = []

    for pathway in pathways:

        genes = get_genes(session, baseURL_biocyc, organism, pathway)

        compounds = get_compounds(session, baseURL_biocyc, organism, pathway)
        
        all_genes.extend(genes)
        G.add_nodes_from(genes)

        gene_rxns = {}
        compound_rxns = {}
        
        for compound in compounds:

            reactionsC = get_compound_reactions(session, baseURL_biocyc, organism, compound)
            
            compound_rxns[compound] = reactionsC

        for gene in genes:

            reactionsG = get_gene_reactions(session, baseURL_biocyc, organism, gene)

            gene_rxns[gene] = reactionsG

        print(gene_rxns, compound_rxns)
        print(compound_rxns)
        
        gene_compounds = {}

        for gene, reactionG in gene_rxns.items():

            for compound, reactionC in compound_rxns.items():

                if set(reactionG) & set(reactionC):

                    gene_compounds.setdefault(gene, []).append(compound)
        

        for gene_a, gene_b in itertools.combinations(gene_compounds.keys(), 2):

            shared = set(gene_compounds[gene_a]) & set(gene_compounds[gene_b])

            if shared:

                G.add_edge(gene_a, gene_b, compounds=list(shared))


    all_genes = list(dict.fromkeys(all_genes))    
    adj = nx.adjacency_matrix(G, nodelist=all_genes)
    adj_df = pd.DataFrame(adj.todense(), index=all_genes, columns=all_genes)

    adj_df.to_csv(f"C:/Users/bella/OneDrive - Imperial College London/FYP/Code/biocyc_adj/{organism}_all_adj.csv")
    
    print(f"saved {organism}_all_adj.csv to C:/Users/bella/OneDrive - Imperial College London/FYP/Code/biocyc_adj")

    return adj_df



###      MAIN LOOP       ###

def main(session, baseURL_uniprot, baseURL_biocyc, organism):

    pathways = get_pathways(session, baseURL_biocyc, organism) #get a list of all pathways for the organism
    pathways = pathways[76:]

    for pathway in pathways:

        print(f"processing {pathway}")

        genes_biocyc = get_genes(session, baseURL_biocyc, organism, pathwayID=pathway) #get an array of all genes in a pa

        pathwayData = []

        for gene in genes_biocyc:

            uniprotIDs = get_uniprot(session, baseURL_biocyc, gene, organism)
            
            if uniprotIDs == []:
                print(f"no uniprot ID for {gene}")
                continue

            else:
                uniprotID = uniprotIDs[0] #extract string
                
                sequence = get_sequence(baseURL_uniprot, uniprotID)
                concat = [uniprotID, sequence]

            pathwayData.append([uniprotID, sequence])

        df = pd.DataFrame(pathwayData,
                          columns = ["uniprot_ID", "sequence"]
        )

        df.to_hdf(f"{pathway}.h5", key="df", index=False)       


