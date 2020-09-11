#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:30:11 2020

@author: vivekmodi
"""

from PyIgClassify_App import db, Cluster
#from flask_sqlalchemy import SQLAlchemy
import sys
import pandas as pd
from datetime import datetime

filename=sys.argv[1]
today=str(datetime.now())[0:10].strip()

df=pd.read_csv(filename,sep=',',header='infer',nrows=3550)
db.drop_all()
db.create_all()
for i in df.index:
    df.at[i,'PDB_Chain_CDR']=f"{df.at[i,'PDB']}_{df.at[i,'original_chain']}_{df.at[i,'CDR']}"
    df.at[i,'CDR_length']=f"{df.at[i,'CDR']}-{df.at[i,'length']}"
    pdb1=Cluster(pdb_chain_cdr=df.at[i,'PDB_Chain_CDR'],\
                 datatag=df.at[i,'datatag'],\
                 pdb=df.at[i,'PDB'],\
                 original_chain=df.at[i,'original_chain'],\
                 CDR=df.at[i,'CDR'],\
                 length=int(df.at[i,'length']),\
                 CDR_length=df.at[i,'CDR_length'],\
                 cluster=df.at[i,'cluster'],\
                 length_type=df.at[i,'length_type'],\
                 fullcluster=df.at[i,'fullcluster'],\
                 center=df.at[i,'center'],\
                 seq=df.at[i,'seq'],\
                 dis=df.at[i,'dis'],\
                 normDis=df.at[i,'normDis'],\
                 DistDegree=round(float(df.at[i,'DistDegree']),2),\
                 bb_rmsd_cdr_align=df.at[i,'bb_rmsd_cdr_align'],\
                 bb_rmsd_stem_align=df.at[i,'bb_rmsd_stem_align'],\
                 ss=df.at[i,'ss'],\
                 rama=df.at[i,'rama'],\
                 dihedrals=df.at[i,'dihedrals'],\
                 gene=df.at[i,'gene'],\
                 species=df.at[i,'species'],\
                 method=df.at[i,'method'],\
                 resolution=df.at[i,'resolution'],\
                 rfactor=df.at[i,'rfactor'],\
                 SeqStart=int(df.at[i,'SeqStart']),\
                 SeqEnd=int(df.at[i,'SeqEnd']),\
                 IsRep=int(df.at[i,'IsRep']),\
                 GSpecies=df.at[i,'GSpecies'],\
                 IG=df.at[i,'IG'],\
                 Germ=df.at[i,'Germ'],\
                 Gpercent=df.at[i,'Gpercent'],\
                 WebCluster=df.at[i,'WebCluster'],\
                 WebDistance=df.at[i,'WebDistance'])
    
    db.session.add(pdb1)
    db.session.commit()
