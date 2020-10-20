#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:30:11 2020

@author: vivekmodi
"""

from PyIgClassify_App import db, perCDR, PDBrows
#from flask_sqlalchemy import SQLAlchemy
import sys
import pandas as pd
from datetime import datetime

filename=sys.argv[1]
today=str(datetime.now())[0:10].strip()

df=pd.read_csv(filename,sep=',',header='infer')
db.drop_all()
db.create_all()
for i in df.index:
    df.at[i,'PDB_Chain_CDR']=f"{df.at[i,'PDB']}_{df.at[i,'Chain']}_{df.at[i,'CDR']}"
    if df.at[i,'CDR seqid']=='-':
        df.at[i,'CDR seqid']=999.0
    pdb1=perCDR(pdb_chain_cdr=df.at[i,'PDB_Chain_CDR'],\
                 cdr=df.at[i,'CDR'],\
                 cdr_length=df.at[i,'CDR Length'],\
                 pdb=df.at[i,'PDB'],\
                 chain=df.at[i,'Chain'],\
                 aho_resnum=df.at[i,'AHO Resnum'],\
                 author_resnum=df.at[i,'Author Resnum'],\
                 sequence=df.at[i,'Sequence'],\
                 germline_sequence=df.at[i,'Germline Sequence'],\
                 gene=df.at[i,'Gene'],\
                 pdb_species=df.at[i,'PDB Species'],\
                 cluster=df.at[i,'Cluster'],\
                 distance=df.at[i,'Distance'],\
                 cdr_germline=df.at[i,'CDR Germline'],\
                 cdr_seqid=df.at[i,'CDR seqid'],\
                 rama4=df.at[i,'Rama4'],\
                 beta_turns=df.at[i,'Beta Turns'],\
                 minimum_edia=df.at[i,'Minimum EDIA'],\
                 keywords=df.at[i,'Keywords'])

    db.session.add(pdb1)
    db.session.commit()

df=pd.read_csv('PDB_wise_data.csv',sep='\t',header='infer')
for i in df.index:
    pdb1=PDBrows(pdb=df.at[i,'PDB'],\
                 light_chain=df.at[i,'LightChain'],\
                 heavy_chain=df.at[i,'HeavyChain'],\
                 h1cluster=df.at[i,'H1Cluster'],\
                 l1cluster=df.at[i,'L1Cluster'],\
                 h2cluster=df.at[i,'H2Cluster'],\
                 l2cluster=df.at[i,'L2Cluster'],\
                 h3cluster=df.at[i,'H3Cluster'],\
                 l3cluster=df.at[i,'L3Cluster'])
    db.session.add(pdb1)
    db.session.commit()
