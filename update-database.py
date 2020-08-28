
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
    pdb1=Cluster(pdb_chain_cdr=df.at[i,'PDB_Chain_CDR'],\
                 datatag=df.at[i,'datatag'],\
                 pdb=df.at[i,'PDB'],\
                 original_chain=df.at[i,'original_chain'],\
                 CDR=df.at[i,'CDR'],\
                 length=int(df.at[i,'length']),\
                 cluster=df.at[i,'cluster'])


    db.session.add(pdb1)
    db.session.commit()
