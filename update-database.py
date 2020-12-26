#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:30:11 2020

@author: vivekmodi
"""

from PyIgClassify_App import db, perVRegion, perCDR
#from flask_sqlalchemy import SQLAlchemy
import sys
import pandas as pd
from datetime import datetime

#filename=sys.argv[1]
today=str(datetime.now())[0:10].strip()

df=pd.read_csv('PerVRegion_df.csv',sep=',',header='infer')
#df['Cluster'].fillna('None',inplace=True)
db.drop_all()
db.create_all()
for i in df.index:
    df.at[i,'PDB']=df.at[i,'PDB'].upper()
    df.at[i,'PDB_heavy_light_chain']=f"{df.at[i,'PDB']}_{df.at[i,'VH Chain']}_{df.at[i,'VL Chain']}"
    
    v_region=perVRegion(pdb_heavy_light_chain=df.at[i,'PDB_heavy_light_chain'],\
                        pdb=df.at[i,'PDB'],\
                        vh_chain=df.at[i,'VH Chain'],\
                        vh_framework=df.at[i,'VH Framework'],\
                        vh_framework_seqid=df.at[i,'VH Framework Seqid'],\
                        h1_cluster=df.at[i,'H1 Cluster'],\
                        h1_cluster_distance=df.at[i,'H1 Cluster Distance'],\
                        h1_seqid=df.at[i,'H1 Seqid'],\
                        h2_cluster=df.at[i,'H2 Cluster'],\
                        h2_cluster_distance=df.at[i,'H2 Cluster Distance'],\
                        h2_seqid=df.at[i,'H2 Seqid'],\
                        h3_cluster=df.at[i,'H3 Cluster'],\
                        h3_cluster_distance=df.at[i,'H3 Cluster Distance'],\
                        vl_chain=df.at[i,'VL Chain'],\
                        vl_framework=df.at[i,'VL Framework'],\
                        vl_framework_seqid=df.at[i,'VL Framework Seqid'],\
                        l1_cluster=df.at[i,'L1 Cluster'],\
                        l1_cluster_distance=df.at[i,'L1 Cluster Distance'],\
                        l1_seqid=df.at[i,'L1 Seqid'],\
                        l2_cluster=df.at[i,'L2 Cluster'],\
                        l2_cluster_distance=df.at[i,'L2 Cluster Distance'],\
                        l2_seqid=df.at[i,'L2 Seqid'],\
                        l3_cluster=df.at[i,'L3 Cluster'],\
                        l3_cluster_distance=df.at[i,'L3 Cluster Distance'])

    db.session.add(v_region)
    db.session.commit()

df=pd.read_csv('PerCDR_df.csv',sep=',',header='infer')
for i in df.index:
    df.at[i,'PDB']=df.at[i,'PDB'].upper()
    df.at[i,'PDB_chain_cdr']=f"{df.at[i,'PDB']}_{df.at[i,'Chain']}_{df.at[i,'CDR']}"
    cdr_row=perCDR(pdb_chain_cdr=df.at[i,'PDB_chain_cdr'],\
                   cdr=df.at[i,'CDR'],\
                   cdr_length=df.at[i,'CDR Length'],\
                   pdb=df.at[i,'PDB'],\
                   chain=df.at[i,'Chain'],\
                   resolution=df.at[i,'Resolution'],\
                   aho_resnum=df.at[i,'AHO Resnum'],\
                   author_resnum=df.at[i,'Author Resnum'],\
                   sequence=df.at[i,'Sequence'],\
                   germline_sequence=df.at[i,'Germline Sequence'],\
                   gene=df.at[i,'Gene'],\
                   pdb_species=df.at[i,'PDB Species'],\
                   frame_germline=df.at[i,'Frame Germline'],\
                   cluster=df.at[i,'Cluster'],\
                   distance=df.at[i,'Distance'],\
                   cdr_germline=df.at[i,'CDR Germline'],\
                   cdr_seqid=df.at[i,'CDR seqid'],\
                   rama4=df.at[i,'Rama4'],\
                   beta_turns=df.at[i,'Beta Turns'],\
                   minimum_edia=df.at[i,'Minimum EDIA'],\
                   keywords=df.at[i,'Keywords'])
    
    db.session.add(cdr_row)
    db.session.commit()
