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

def mining_pubs(tipo,textEcGen2):
    """Docstring mining_pubs
    La función mining_pubs tendra dos argumentos que seran los correspondientes al que proviene del archivo y el tipo de busuqeda, que puede ser: DP, AU o AD.
           Si el tipo es "DP" recupera el año de publicación del artículo. El retorno es un dataframe con el PMID y el DP_year.
           Si el tipo es "AU" recupera el número de autores por PMID. El retorno es un dataframe con el PMID y el num_auth.
           Si el tipo es "AD" recupera el conteo de autores por país. El retorno es un dataframe con el country y el num_auth."""
    #Primero se extraen los datos de PMID y para ello se usa splitlines para serar en una nueva lista y se usa starwith para unicamente obtener las que empiezan con PMID, este proceso se hace para todos los tipos de extracion de datos, ya sean de tipos PD, AU o AD 
    #Para esto se crea una lista onde se almacenaran los datos, usando append
    #Si el tipo es "DP" recupera el año de publicación del artículo. El retorno es un dataframe con el PMID y el DP_year
    #Listas donde se guardaran los datos despues de buscarlos en el keyword
    PMID1 = []
    PMID2 = []
    DP_year1 = []
    DP_year2 = []
    #Se hace una busqueda linea por liena de lo seleccionado
    for line in textEcGen2.splitlines():
        #Si empieza con PMID
        if line.startswith("PMID-"):
            PMID1.append(line[:])
            PMID4 = line[:]
            PMID3 = re.findall(r'\d{8}$', PMID4)
            PMID2.append(PMID3)
        #Si es la fecha de publicacion
        if line.startswith("DP  -"):
            DP_year1.append(line[:])
            DP_year4 = line[:]
            DP_year3 = re.findall(r'\d{4}', DP_year4)
            DP_year2.append(DP_year3)
    #Los datos se los ubica en forma iretrada
    PMID2 = list(itertools.chain.from_iterable(PMID2))
    DP_year2 = list(itertools.chain.from_iterable(DP_year2))
    #Se ubican los datos al dataframe
    tablaAnios = pd.DataFrame()
    tablaAnios['PMID'] = PMID2
    tablaAnios['DP_year'] = DP_year2
    
    #Si el tipo es "AU" recupera el número de autores por PMID. El retorno es un dataframe con el PMID y el num_auth
    #Se busca el patron en especifico y luego se crea una lista para poder contar el numero de autores
    AU = textEcGen2.split("PMID- ")
    AU.pop(0) #Se borran espacios vacios
    num_AU = []
    #Bucle que permite contar los autores
    for line in range(len(AU)):
        cont = re.findall("AU -", AU[line])
        cant = (len(cont))
        num_AU.append(cant)
    #Se ubican los datos al dataframe
    tablaAutores = pd.DataFrame()
    tablaAutores["PMID"] = PMID2 
    tablaAutores["num_auth"] = num_AU
    
    #Si el tipo es "AD" recupera el conteo de autores por país. El retorno es un dataframe con el country y el num_auth
    #Listas a usar en esta sección
    AD = []
    paisesAlFinal1 = []
    paisesAlFinal2 = []
    paisesAlFinal3 = []
    paisesAlFinalSinEA1 = []
    paisesAlFinalSinEA2 = []
    paisesAlFinalSinEA3 = []
    paisesAlFinalConEA1 = []
    paisesAlFinalConEA2 = []
    paisesAlFinalConEA3 = []
    paisesAlFinalConNum1 = []
    paisesAlFinalConNum2 = []
    paisesAlFinalConNum3 = []
    paisesTodosPatrones = []
    #Se busca la informacion esperada
    for line in textEcGen2.splitlines():
        if line.startswith("AD  -"):
            AD.append(line[:])
    for line in textEcGen2.splitlines():
        if line.startswith("AD  -"):
            AD = line[:]
            #Se ingresan los patrones mas comunes por la dimension del archivo, es dificil leer todos
            #Pais al final de una palabra
            paisAlFinal1 = re.findall(r'\,\s(\w{2,16})\.', AD)
            paisesAlFinal1.append(paisAlFinal1)
            #Pais al final de dos palabras
            paisAlFinal2 = re.findall(r'\,\s(\w{2,16}[^0-9\,]\s\w{2,16}[^0-9])\.', AD)
            paisesAlFinal2.append(paisAlFinal2)
            #Pais al final de tres palabras
            paisAlFinal3 = re.findall(r'\,\s(\w{3,16}[^0-9\,]\s\w{2,3}[^0-9\,]\s\w{3,16}[^0-9\,])\.', AD)
            paisesAlFinal3.append(paisAlFinal3)
            #Pais al final de una palabra y correo sin EA
            paisAlFinalSinEA1 = re.findall(r'\,\s(\w{2,16})\.\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}', AD)
            paisesAlFinalSinEA1.append(paisAlFinalSinEA1)
            #Pais al final de dos palabras y correo sin EA
            paisAlFinalSinEA2 = re.findall(r'\,\s(\w{2,16}[^0-9\,]\s\w{2,16}[^0-9])\.\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}', AD)
            paisesAlFinalSinEA2.append(paisAlFinalSinEA2)
            #Pais al final de tres palabras y correo sin EA
            paisAlFinalSinEA3 = re.findall(r'\,\s(\w{3,16}[^0-9\,]\s\w{2,3}[^0-9\,]\s\w{3,16}[^0-9\,])\.\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}', AD)
            paisesAlFinalSinEA3.append(paisAlFinalSinEA3)
            #Pais al final de una palabra y correo con EA
            paisAlFinalConEA1 = re.findall(r'\,\s(\w{2,16})\. Electronic address:\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}\.', AD)
            paisesAlFinalConEA1.append(paisAlFinalConEA1)
            #Pais al final de dos palabras y correo con EA
            paisAlFinalConEA2 = re.findall(r'\,\s(\w{2,16}[^0-9\,]\s\w{2,16}[^0-9])\. Electronic address:\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}\.', AD)
            paisesAlFinalConEA2.append(paisAlFinalConEA2)
            #Pais al final de tres palabras y correo con EA
            paisAlFinalConEA3 = re.findall(r'\,\s(\w{3,16}[^0-9\,]\s\w{2,3}[^0-9\,]\s\w{3,16}[^0-9\,])\. Electronic address:\s[a-z0-9_\.-]+@[\da-z\.-]+\.[a-z\.]{2,6}\.', AD)
            paisesAlFinalConEA3.append(paisAlFinalConEA3)
            #Pais al final de una palabra y numero antes
            paisAlFinalConNum1 = re.findall(r'\,\s\w{3,9}[0-9\-]\,\s(\w{2,16})\.', AD)
            paisesAlFinalConNum1.append(paisAlFinalConNum1)
            #Pais al final de dos palabras y numero antes
            paisAlFinalConNum2 = re.findall(r'\,\s\w{3,9}[0-9\-]\,\s(\w{2,16}[^0-9\,]\s\w{2,16}[^0-9])\.', AD)
            paisesAlFinalConNum2.append(paisAlFinalConNum2)
            #Pais al final de tres palabras y numero antes
            paisAlFinalConNum3 = re.findall(r'\,\s\w{3,9}[0-9\-]\,\s(\w{3,16}[^0-9\,]\s\w{2,3}[^0-9\,]\s\w{3,16}[^0-9\,])\.', AD)
            paisesAlFinalConNum3.append(paisAlFinalConNum3)

            #Al existir errores por la cantidad de datos no se usa | para unir en un solo bucle todos los patrones
            #En alternativa se une todos en una sola lista
            paisesTodosPatrones=paisesAlFinal1+paisesAlFinal2+paisesAlFinal3+paisesAlFinalSinEA1+paisesAlFinalSinEA2+paisesAlFinalSinEA3+paisesAlFinalConEA1+paisesAlFinalConEA2+paisesAlFinalConEA3+paisesAlFinalConNum1+paisesAlFinalConNum2+paisesAlFinalConNum3
        #print(paises)
    #A continuacion se creara unan lista ordenada de los paises, sin repetirlos
    paisesTodosPatrones= list(itertools.chain.from_iterable(paisesTodosPatrones))
    len(paisesTodosPatrones)
    #10 primeros paises
    unique_paisesTodosPatrones = list(set(paisesTodosPatrones))
    unique_paisesTodosPatrones.sort()
    len(unique_paisesTodosPatrones)
    #Se importa la liberia capaz de abrir y leer las columas y filas de la tabala de coordenadas del mundo para asiganarlas a un diccionario
    import csv
    #Declarar diccionario vacio
    coordenadas = {}
    #Abrir el archivo de coordenadas
    with open('Data/coordenadas.txt') as f:
        csvr = csv.DictReader(f)
        for row in csvr:
            coordenadas[row['Name']] = [row['Latitude'], row['Longitude']]
            #print(coordenadas)
    #Se crea una lista por cada elemento de unico y se compara con el documento de coordenadas para asi graficarlo y contarlo
    pais = []
    longitud = []
    latitud = []
    contador = []
    for z in unique_paisesTodosPatrones:
        if z in coordenadas.keys():
            pais.append(z)
            latitud.append(float(coordenadas[z][0]))
            longitud.append(float(coordenadas[z][1]))
            contador.append(paisesTodosPatrones.count(z))
    tablaPaises = pd.DataFrame()
    tablaPaises["country"] = pais 
    tablaPaises["num_auth"] = contador
    
    #Dato que retorna del tipo
    if tipo == 'AD':
        return tablaAnios
    if tipo == 'AU':
        return tablaAutores
    if tipo == 'PD':
        return tablaPaises

    