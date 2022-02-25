#Bibliotecas a usar en el examen
import re
import csv
import itertools
import pandas as pd
from Bio import Entrez

#Funcion de busqeda con el keyword
def download_pubmed(keyword): 
    """En esta funcion se realiza una busqeuda en PubMed sin la necesidad de descargar un archivo, para ello se usa Entrez que permite hacer la busqueda directa con el servidor."""
    #"Ecuador genomics [Title/Abstract]"
    Entrez.email = "jessica.quinonez@est.ikiam.edu.ec"
    handle = Entrez.esearch (db = "pubmed",
                             term = keyword,
                             usehistory = "y")
    record = Entrez.read (handle)
    id_list = record ["IdList"]
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    handle = Entrez.efetch(db="pubmed", 
                           rettype="medline", 
                           retmode="text", 
                           retstart=0, 
                           retmax = 1500, 
                           webenv = webenv, 
                           query_key = query_key)
    #Se lee los archivos "descargados" del servidor
    textEcGen1 = handle.read()
    #Se cambia los saltos de linea y espacios por un solo salto al final
    textEcGen2 = re.sub(r'\n\s{6}', ' ', textEcGen1)
    #El valor que retorna para la funcion siguiente
    return (textEcGen2)


def mining_pubs(tipo):
    """Docstring mining_pubs"""
    if tipo == "AD":
        
    
    return 

    