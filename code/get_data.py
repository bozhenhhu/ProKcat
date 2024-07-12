import requests
import numpy as np
import pandas as pd
from pubchempy import Compound, get_compounds
from matplotlib import pyplot as plt
import seaborn as sns
from urllib import request
from brendapyrser import BRENDA
import html
import pickle
from math import exp
from feature_functions import *
import random
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs


brenda = BRENDA('../brenda.txt')


def get_smiles(substrate):
    try :
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/CanonicalSMILES/TXT'%substrate
        req = requests.get(url)
        if req.status_code != 200:
            smiles = 'NaN'
        else:
            smiles = req.content.splitlines()[0].decode()
    except :
        smiles = 'NaN'
    return smiles

def get_seq(ID):
    url = "https://www.uniprot.org/uniprot/%s.fasta" % ID
    try :
        data = requests.get(url)
        if data.status_code != 200:
            seq = 'NaN'
        else:
            seq =  "".join(data.text.split("\n")[1:])
    except :
        seq = 'NaN'
    return seq

def check_mutations(seq, mut_list):
    no_error = True
    for mut in mut_list:
        ind = int(mut[1:-1])-1
        old = mut[0].upper()
        if (ind > len(seq)-1) or (seq[ind] != old):
            no_error = False
            break
    return no_error

def apply_mutations(seq, mut_list):
    mut_seq = seq
    for mut in mut_list:
        ind = int(mut[1:-1])-1
        new = mut[-1].upper()
        temp_list = list(mut_seq)
        temp_list[ind] = new
        mut_seq = ''.join(temp_list)
    return mut_seq


brenda_ec_list = []
for rxn in brenda.reactions:
    brenda_ec_list.append( rxn.ec_number )
brenda_ec_list = list(set(brenda_ec_list))
print(len(brenda_ec_list))


QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'

with open('../data/enzyme.dat', 'r') as outfile :
    lines = outfile.readlines()

ec_list = []
for line in lines :
    if line.startswith('ID') :
        ec = line.strip().split('  ')[1]
        ec_list.append(ec.strip())
# print(ec_list)
print(len(ec_list))


def get_entry_kcat( ec ):
    r = brenda.reactions.get_by_id(ec)
    all_data = r.Kcatvalues
    result = []
    for sub in all_data.keys():
        sub_data = all_data[sub]
        for entry in sub_data:
            if ('Â°C' not in entry['meta'] ) or ( '#' not in entry['meta']) \
                        or (';' in entry['meta']) or ('specified' in entry['meta'] ):
                continue
            else:
                value = entry['value']
                temperature = html.unescape( entry['meta'].split('Â°C')[0] ) [-2:]
                try :
                    temperature = float(temperature)  
                except:
                    continue
                    
                if 'mutant' not in entry['meta']:
                    enz_type = 'WT'
                    mutation = 'None'
                else:
                    mut4 = re.findall('[A-Z][0-9][0-9][0-9][0-9][A-Z]',entry['meta'])
                    mut3 = re.findall('[A-Z][0-9][0-9][0-9][A-Z]',entry['meta'])
                    mut2 = re.findall('[A-Z][0-9][0-9][A-Z]',entry['meta'])
                    mut1 = re.findall('[A-Z][0-9][A-Z]',entry['meta'])
                    mut_list = mut4 + mut3 + mut2 + mut1
                    if len(mut_list) < 1:
                        continue
                    else:
                        enz_type = 'MUT'
                        mutation = '/'.join(mut_list) 
                    
                p_ref = entry['meta'].split('#')[1].strip()
                if ',' in p_ref:
                    p_ref_list = p_ref.split(',')
                else:
                    p_ref_list = [ p_ref ]
                p_ids = []
                for ref in p_ref_list:
                    p_ids.append( r.proteins[ref]['proteinID']  )
                    
                for p_id in p_ids:
                    if p_id == '':
                        continue
                    else:
                        result.append( {'EC':ec,'temperature':float(temperature),'sub': sub, 
                                'UniProtID':p_id,'EnzymeType':enz_type,'Mutation':mutation,'kcat': float(value) } )
    return result


result = []
idx = 0
for ec in brenda_ec_list:
    if idx % 500 == 0:
        print(str(idx) + ' done')
    result += get_entry_kcat( ec )
    idx+=1

rawdata_brenda = pd.DataFrame(result)
rawdata_brenda = (rawdata_brenda[rawdata_brenda['kcat']>0]).dropna().reset_index().drop(['index'],axis=1)
proteinIDs = []
for i in range(len(rawdata_brenda['UniProtID'])):
    ID = list( rawdata_brenda['UniProtID'] )[i]
    proteinIDs.append( ID.split(' ')[0] )
rawdata_brenda['UniProtID'] =  proteinIDs  

kcat_brenda = []
total = len(rawdata_brenda['sub'])

for i in range(len(rawdata_brenda['sub'])):
    ec, T, sub, pid, enz_type, muts, kcat = rawdata_brenda.iloc[i]
    data={'EC':ec,'Temp':T,'sub':sub,'ProtID':pid,'EnzymeType':enz_type,'Mutation':muts,'kcat':kcat}
    data['smiles']=get_smiles( sub )
    if data['smiles']  == 'NaN' or data['smiles'] == '':
        continue
    temp_seq = get_seq( pid )
    if temp_seq == 'NaN' or temp_seq == '':
        continue
    if enz_type == 'WT':
        data['seq'] = temp_seq
    else:
        mut_list = muts.split('/')
        if check_mutations(temp_seq, mut_list):
            temp_seq = apply_mutations(temp_seq, mut_list)
            data['seq'] = temp_seq
        else:
            continue
              
    kcat_brenda.append(data)
    if i%2000 == 0:
        print(str(i/total)+'% done')

raw_kcat_brenda = pd.DataFrame( kcat_brenda )
raw_kcat_brenda = raw_kcat_brenda.dropna().reset_index().drop(['index'],axis=1)

raw_kcat_brenda.to_csv('../data/raw_kcat_brenda.csv',index = None)

