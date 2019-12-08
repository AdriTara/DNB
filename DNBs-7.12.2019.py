#####################################################
### @author: Adrián Tarazona Sánchez              ###
###                                               ###
###                                               ###
### sys.argv[1] = normalized data file            ###
### sys.argv[2] = probe.spot-Nt.ID file           ###
### sys.argv[3] = Nt.ID-At.ID file                ###
### sys.argv[4] = interaction-network file        ###
#####################################################
"""
This script takes as arguments:
    1-An expression file separated by tabs
    2-A file determining the correspondence between the spotID and NTgenID
    3-A file determining the correspondence between the NtgenID and AtgenID
    4-A network interaction file with the interactor1 and interactor2 separated by 2 tabs
First, it gives lists of biphasic genes ordered from the gene having the highest
expression jump to the gene having the lower expression jump.
"""

import sys
import numpy as np
import os
import pandas as pd
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import multiprocessing as mp
import matplotlib as mpl
import datetime 
mpl.use('Agg')

#%% Functions to be used
nlines = 0
for line in open(sys.argv[1]): #'datos_trabajo/GSE99838_normalized_data_with_probename.txt' #sys.argv[1]
    nlines += 1    
def phaseanalysis(lista):
    FCs = np.array(lista).astype(np.float)
    media = np.mean(FCs)
    phases = 1
    for j in range(0, len(FCs)):
        if j == 0:
            FCi = FCs[j]
            if FCi > media:
                marcador1 = 1
            elif FCi < media:
                marcador1 = 0
        else:
            FC = FCs[j]
            if FC > media:
                marcador2 = 1
            elif FC < media:
                marcador2 = 0
            if marcador1 != marcador2:
                phases += 1
                marcador1 = marcador2
                delta = FC-FCi
                FCi = FC
                salto = j
    return phases, delta, salto

def biplot(list1, list2):
    plt.plot(list1,list2)                                             
    plt.xlabel('Strain')                                                 
    plt.ylabel('Log2(FC)')                                               
    plt.hlines(np.mean(list2),0,len(list1)-1, linestyles='dotted')           
    plt.xticks(rotation=70)
    plt.tight_layout()
    plt.show()
    
def unico(lista):
    unica = []
    for element in lista:
        if element not in unica:
            unica.append(element)
    return unica
def ajuste(lista):
    minimo = float(min(lista))
    normalizada = []
    for valor in lista:
        valor = float(valor)/minimo
        normalizada.append(valor)
    return normalizada

def toma_datos(d):
    sprobes = []
    smedias = []
    sexpre = []
    ncores = mp.cpu_count()
    z = 0
    for line in open(sys.argv[1]): #'datos_trabajo/GSE99838_normalized_data_with_probename.txt' #sys.argv[1]
        line = line.strip('\n')
        fields = line.split('\t')
        z += 1
        line = line.strip('\n')
        media_fila = []
        expre_fila = []
        if z > 1:
            if z > (d*nlines/ncores) and z <= ((d+1)*nlines/ncores):
                sprobes.append(line.split('\t')[0])
                for j in range(2,len(fields)-3,3):
                    exp = np.array(fields[j:j+3]).astype(np.float)
                    if j >= 20:
                        exp = np.array(fields[j:j+4]).astype(np.float)
                    media_fila.append(np.mean(exp))
                    expre_fila.append(exp)
                smedias.append(media_fila)
                sexpre.append(expre_fila)
    return [smedias, sprobes, sexpre]
def list_to_txt(lista, file):
    out = open(file, 'w')
    for item in lista:
        out.write(item+'\n')
#%% Data processing
probe_gen = {}
i = 0
for line in open(sys.argv[2]): #'datos_trabajo/021113_D_GeneList_20130122.txt' #sys.argv[2]
    line = line.strip('\n')
    if i >= 1:
        try:
            NT = line.split('\t')[1]
            probe = line.split('\t')[0]
            if NT != '':
                probe_gen[probe] = NT
        except:
            continue
    i += 1
print('probe-gen se crea')
#crear diccionario de ortologos
NT_AT = {}
for line in open(sys.argv[3]): #'datos_trabajo/ortologos-Nt_At.txt' #sys.argv[3]
    line = line.strip('\n')
    NT = line.split('\t')[0]
    AT = line.split('\t')[1]
    NT_AT[NT] = AT
print('el diccionario de ortólogos se crea', datetime.datetime.now())

i = 0
num_cepa = {}
cepas = []
for line in open(sys.argv[1]): #'datos_trabajo/GSE99838_normalized_data_with_probename.txt' #sys.argv[1]
    line = line.strip('\n')
    i += 1
    for j in range(2,len(line.split('\t'))-8,3):
        if i == 1:
            cepa = line.split('\t')[j].split('_')[0]
            num_cepa[j] = cepa
            cepas.append(cepa)
    for j in range(len(line.split('\t'))-7,len(line.split('\t'))-3,3):
        if i == 1:
            cepa = line.split('\t')[j].split('_')[0]
            num_cepa[j-1] = cepa
            cepas.append(cepa)
    if i > 1:
        break
#%% TOMA DE DATOS

pool = mp.Pool(mp.cpu_count())
resultado = pool.map(toma_datos, [j for j in range(mp.cpu_count())])
pool.close()

lmedias = []
lexpre = []
probes = []
for core in resultado:
    lmedias += core[0]
    lexpre += core[2]
    probes += core[1]
print('datos de expresion tomados', datetime.datetime.now())
ntID = []
atID = []
for probe in probes:
    try:
        ntID.append(probe_gen[probe])
        try:
            atID.append(NT_AT[probe_gen[probe]])
        except:
            atID.append('Sin ortologo')
    except:
        ntID.append('Sin ntID')
        atID.append('Sin ntID')
print('ntID y atID creados', datetime.datetime.now())
#%%
medias = pd.DataFrame(lmedias, columns = cepas)
medias.insert(0,'probeID',probes, allow_duplicates=True)
medias.insert(0,'ntID', ntID, allow_duplicates=True)
medias.insert(0,'atID', atID, allow_duplicates=True)
expresion = pd.DataFrame(lexpre, columns = cepas)
expresion.insert(0,'probeID',probes, allow_duplicates=True)
expresion.insert(0,'ntID', ntID, allow_duplicates=True)
expresion.insert(0,'atID', atID, allow_duplicates=True)
print('cálculo de medias---OK\nvectores de expresion---OK')

FCs = medias[medias.columns[:-1]].copy() #pd.DataFrame(columns = expresion.columns[:-1])
for column in FCs.columns[3:]:
    FCs.loc[:,column] = np.log2((medias[column]/medias['Negative']).astype(np.float64))

#%% Bifasic analysis
print('empieza el análisis de bifásicos', datetime.datetime.now())
bifasicos = pd.DataFrame(columns = FCs.columns)
for index in FCs.index:
    exp = FCs.iloc[index,3:]
    if phaseanalysis(exp)[0] == 2:
#        biplot(FCs.columns,exp) #UNCOMMENT THIS LINE IF BIPHASIC PLOT IS DESIRED
        bifasicos.loc[index,:] = FCs.loc[index,:]
        bifasicos.loc[index, 'DeltaFC'] = phaseanalysis(exp)[1]
        bifasicos.loc[index, 'jump'] = int(phaseanalysis(exp)[2])

#ordenar los bifásicos según su DeltaFC        
bifasicos = bifasicos.reindex(bifasicos['DeltaFC'].abs().sort_values(ascending = False).index)
bifasicos = bifasicos.set_index('atID')

#separar los bifásicos según si aumenta o disminuye su expresión
bifasicos_up = bifasicos[bifasicos['DeltaFC'] > 0]
bifasicos_down = bifasicos[bifasicos['DeltaFC'] < 0]
print('análisis de bifasicos---OK', datetime.datetime.now())

#%% Module obtention
print('empieza la obtención de modulos', datetime.datetime.now())
interacciones = []
i = 0
for line in open(sys.argv[4]):#'datos_trabajo/protein-network.txt' #sys.argv[4]
    line = line.strip('\n')
    if i > 0:
        gen1 = line.split('\t')[0]
        gen2 = line.split('\t')[2]
        interaccion = [gen1,gen2]
        interacciones.append(interaccion)
    i += 1  

i = 0
submodulos = {} #key = gen bifasico, value = {key=gen de interacción, value = linea de interaccion}
submoduls = []
for gen in bifasicos.index.tolist():
    i += 1
    if i % 1000 == 0:
        print('avanza la creación de submodulos', gen)
    for interaccion in interacciones:
        if (gen == interaccion[0]) or (gen == interaccion[1]):
            #tenemos un nuevo módulo de la red
            if gen in submodulos:
                if interaccion[0] not in submodulos[gen]:
                    submodulos[gen].append(interaccion[0])
                if interaccion[1] not in submodulos[gen]:
                    submodulos[gen].append(interaccion[1])
            elif gen not in submodulos:
                #checkeamos que el gen no pertenezca ya a un módulo
                submodulos[gen] = [interaccion[0],interaccion[1]]
                submoduls.append(gen)

modulos = {}
i = 0 #contador de submodulos
for submodulo1 in sorted(submodulos):
    if i % 100 == 0:
        print('avanza la creacion de modulos', submodulo1)
    try:
        i += 1    
        modulo = sorted(submodulos[submodulo1])
        submoduls.remove(submodulo1)
        for submodulo2 in sorted(submoduls):
            for gen in sorted(submodulos[submodulo2]):
                if gen in modulo:
                    modulo += sorted(submodulos[submodulo2])
                    submoduls.remove(submodulo2)
                    break
    except:
        continue
    modulos['modulo'+str(i)] = unico(modulo)

print('los modulos se crean', datetime.datetime.now())
#%% CI calculation
print('Empiezan los cálculos de CI', datetime.datetime.now())
expresionFC = expresion.copy()
expresionFC = expresionFC.drop_duplicates(subset='atID', keep = 'last')
for gen in expresionFC.index:
    expresionFC.loc[gen, 'Negative'] = np.mean(expresion.loc[gen,'Negative'])
for column in expresionFC.columns[3:-1]:
    for gen in expresionFC.index:
        exp = []
        for valor in expresionFC.loc[gen,column]:
            exp.append(np.log2(float(valor)/expresionFC.loc[gen,'Negative']))
        expresionFC.loc[gen,column] = np.array(exp).astype(np.float)
expresionFC = expresionFC.set_index('atID')
print('calculo de expresionFC---OK ', datetime.datetime.now())
dfs_modulos = {}
for modulo in modulos:
    dfs_modulos[modulo] = pd.DataFrame(columns = expresionFC.columns[2:-1]).copy()
    for gen in modulos[modulo]:
        try:
            dfs_modulos[modulo].loc[gen,:] = expresionFC.loc[gen, :'WT']
        except:
            continue
print('los dataframes de los modulos se crean', datetime.datetime.now())

def CI_calc(columna, df):
    sSD = 0
    sPCCin = 0
    sPCCout = 0
    calcsd = 0
    calcin = 0
    calcout = 0
    expresionFC2 = expresionFC.drop(df.index, axis=0)
    for exp1 in df[columna]:
        exp1 = np.array(exp1)
        sSD += np.std(exp1)
        calcsd += 1
        for exp2 in df[calcsd:][columna]:
            exp2 = np.array(exp2)
            sPCCin += abs(pearsonr(exp1,exp2)[0])
            calcin += 1
        for exp2 in expresionFC2[calcsd:][columna]:
            sPCCout += abs(pearsonr(exp1,exp2)[0])
            calcout += 1
    mSDin = sSD/calcsd
    mPCCin = sPCCin/calcin
    mPCCout = sPCCout/calcout
    ci = mSDin*mPCCin/mPCCout
    #return [mSDin, mPCCin, mPCCout,ci]
    return ci

pool = mp.Pool(mp.cpu_count())
for modulo in dfs_modulos:
    print('Cálculo CIs '+modulo)
    df = dfs_modulos[modulo]
    if len(df) <= 1:
        df.loc['CI',:] = np.nan
    else:
        df.loc['CI',:] = pool.starmap(CI_calc, [(columna, dfs_modulos[modulo]) for columna in dfs_modulos[modulo]])

pool.close()
print('Cálculos de CI---OK', datetime.datetime.now())

#%% OUTPUTS
print('Empieza la generación de outputs', datetime.datetime.now())
os.system('mkdir biphasic-outputs')
expresion.to_csv('biphasic-outputs/expresion-vector.csv')
list_to_txt(bifasicos.index, 'biphasic-outputs/ath-idlist-bifasicos.txt')
bifasicos.iloc[:,:-1].to_csv('biphasic-outputs/bifasicos.csv')
bifasicos_up.iloc[:,:-1].to_csv('biphasic-outputs/bifasicos_up.csv')
list_to_txt(bifasicos_up.index, 'biphasic-outputs/ath-idlist-bifasicos_up.txt')
bifasicos_down.iloc[:,:-1].to_csv('biphasic-outputs/bifasicos_down.csv')
list_to_txt(bifasicos_down.index, 'biphasic-outputs/ath-idlist-bifasicos_down.txt')

for j in range(0,len(cepas)-2):
    pareja = '-'.join(cepas[j:j+2])
    bifasicos[bifasicos['jump']==1].to_csv('biphasic-outputs/bifasicos-'+str(pareja)+'.csv')
    list_to_txt(bifasicos[bifasicos['jump']== j+1].index, 'biphasic-outputs/idlist-bifasicos'+str(pareja)+'.txt')

os.system('mkdir dnb-outputs')
i = 0
j = 0
for modulo in dfs_modulos:
    CI = dfs_modulos[modulo].loc['CI',:]
    if CI[1] == np.nan:
        break
    CI = ajuste(CI)
    dfs_modulos[modulo].loc['CI',:] = CI
    dfs_modulos[modulo].to_csv(r'dnb-outputs/'+modulo+'.csv')
    list_to_txt(dfs_modulos[modulo].index[:-1], 'dnb-outputs/'+modulo+'-ath-idlist.txt')
    #media = np.mean(CI)
    #SD = np.std(CI)
    genenum = str(len(dfs_modulos[modulo]))
    #for valor in CI:
        #if (valor > (media+2*SD)):
    plt.plot(cepas[:-1], CI)
    plt.xticks(rotation=70)
    plt.title(modulo)
    plt.xlabel('Strain')
    plt.ylabel('CI')
    plt.tight_layout()
    plt.text(5, max(CI), genenum+' genes')
    #plt.show()
    plt.savefig('dnb-outputs/CI-'+modulo+'.png')
    plt.clf()
            #break
print('Los outputs se han generado', datetime.datetime.now())