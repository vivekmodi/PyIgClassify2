#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:30:11 2020

@author: vivekmodi
"""

import os, subprocess, sys, glob, random, numpy as np
from datetime import datetime
from flask import Flask, render_template, url_for, redirect, session, request
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired
from werkzeug.utils import secure_filename
from sqlalchemy import func
from collections import defaultdict
from flask import Markup
from sqlalchemy import desc,asc
import pandas as pd
from natsort import natsorted

pwd=os.getcwd()
sys.path.append(pwd+'/scripts')

UPLOAD_FOLDER = (pwd+'/server/uploads')
ALLOWED_EXTENSIONS = {'txt', 'pdb','cif'}

app=Flask(__name__)
app.config['SECRET_KEY'] = 'mysecretkey'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['ALLOWED_EXTENSIONS'] = ALLOWED_EXTENSIONS
app.jinja_env.add_extension('jinja2.ext.loopcontrols')

############SQL Databse and Models############################
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///cdr_data.sqlite'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

db = SQLAlchemy(app)
Migrate(app,db)

class Cluster(db.Model):
    pdb_chain_cdr=db.Column(db.Text,primary_key=True)
    datatag=db.Column(db.Text)
    pdb=db.Column(db.Text)
    original_chain=db.Column(db.Text)
    CDR=db.Column(db.Text)
    length=db.Column(db.Integer)
    CDR_length=db.Column(db.Text)
    cluster=db.Column(db.Text)
    length_type=db.Column(db.Text)
    fullcluster=db.Column(db.Text)
    center=db.Column(db.Integer)
    seq=db.Column(db.Text)
    dis=db.Column(db.REAL)
    normDis=db.Column(db.REAL)
    DistDegree=db.Column(db.REAL)
    bb_rmsd_cdr_align=db.Column(db.REAL)
    bb_rmsd_stem_align=db.Column(db.REAL)
    ss=db.Column(db.Text)
    rama=db.Column(db.Text)
    dihedrals=db.Column(db.Text)
    gene=db.Column(db.Text)
    species=db.Column(db.Text)
    method=db.Column(db.Text)
    resolution=db.Column(db.REAL)
    rfactor=db.Column(db.REAL)
    SeqStart=db.Column(db.Integer)
    SeqEnd=db.Column(db.Integer)
    IsRep=db.Column(db.Integer)
    GSpecies=db.Column(db.Text)
    IG=db.Column(db.Text)
    Germ=db.Column(db.Text)
    Gpercent=db.Column(db.REAL)
    WebCluster=db.Column(db.Text)
    WebDistance=db.Column(db.REAL)
    
def __init__(self,pdb_chain_cdr,datatag,pdb,original_chain,CDR,length,CDR_length,cluster,length_type,fullcluster,center,seq,dis,normDis,DistDegree,bb_rmsd_cdr_align,\
             bb_rmsd_stem_align,ss,rama,dihedrals,gene,species,method,resolution,rfactor,SeqStart,SeqEnd,IsRep,GSpecies,IG,Germ,Gpercent,WebCluster,WebDistance):
    
    self.pdb_chain_cdr=pdb_chain_cdr;self.datatag=datatag;self.pdb=pdb;self.original_chain=original_chain;self.CDR=CDR;self.length=length;\
    self.CDR_length=CDR_length;self.cluster=cluster;self.length_type=length_type;self.fullcluster=fullcluster;self.center=center,self.seq=seq;self.dis=dis;\
    self.normDis=normDis;self.DistDegree=DistDegree;self.bb_rmsd_cdr_align=bb_rmsd_cdr_align;self.bb_rmsd_stem_align=bb_rmsd_stem_align;\
    self.ss=ss;self.rama=rama;self.dihedrals=dihedrals;self.gene=gene;self.species=species;self.method=method;self.resolution=resolution;\
    self.rfactor=rfactor;self.SeqStart=SeqStart;self.SeqEnd=SeqEnd;self.IsRep=IsRep;self.GSpecies=GSpecies;self.IG=IG;self.Germ=Germ;\
    self.Gpercent=Gpercent;self.WebCluster=WebCluster;self.WebDistance=WebDistance
    
def __repr__(self):
    return f'{self.pdb_chain_cdr} {self.datatag} {self.pdb} {self.original_chain} {self.CDR} {self.length} {self.CDR_length} {self.cluster} {self.length_type}\
             {self.fullcluster} {self.center} {self.seq} {self.dis} {self.normDis} {self.DistDegree} {self.bb_rmsd_cdr_align} {self.bb_rmsd_stem_align}\
             {self.ss} {self.rama} {self.dihedrals} {self.gene} {self.species} {self.method} {self.resolution} {self.rfactor} {self.SeqStart}\
             {self.SeqEnd} {self.IsRep} {self.GSpecies} {self.IG} {self.Germ} {self.Gpercent} {self.WebCluster} {self.WebDistance}'

def create_lists():
    pdbListDb=list();chainListDb=list();clusterListDb=list()
    for pdbs in Cluster.query.with_entities(Cluster.pdb_chain_cdr):
        pdbid,chainid,cdr=pdbs[0].split('_')
        pdbListDb.append(pdbid);chainListDb.append(pdbid+chainid);
    for fullcluster in Cluster.query.with_entities(Cluster.fullcluster):
        clusterListDb.append(fullcluster[0])
    clusterListDb=list(set(sorted(clusterListDb)))
    return pdbListDb,chainListDb,clusterListDb
    
@app.route('/index')
@app.route('/home')
@app.route('/')
def home():
    return render_template('home.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/browse')
def browse():
    retrieve_all=Cluster.query.all()
    return render_template('browse.html',retrieve_all=retrieve_all)

@app.route('/statistics')
def statistics():
    species_unique=dict();gene_unique=dict();entries=dict();chains=dict();
    (pdbListDb,chainListDb,clusterListDb)=create_lists()
    for fullcluster in clusterListDb:
        #cluster_list=Cluster.query.filter(Cluster.fullcluster.contains(queryname)).all()
        entries[fullcluster]=Cluster.query.filter(Cluster.fullcluster.contains(fullcluster)).count()
        species_list=Cluster.query.filter(Cluster.fullcluster.contains(fullcluster)).with_entities('species')
        species_unique[fullcluster]=pd.Series(names[0] for names in species_list).unique()
        gene_list=Cluster.query.filter(Cluster.fullcluster.contains(fullcluster)).with_entities('gene')
        gene_unique[fullcluster]=pd.Series(names[0] for names in gene_list).unique()
        
        
    return render_template('statistics.html',clusterListDb=natsorted(clusterListDb,key=str),entries=entries,species_unique=species_unique,gene_unique=gene_unique)

@app.route('/formSearch', methods=['GET','POST'])
def formSearch():
    (pdbListDb,chainListDb,clusterListDb)=create_lists()
    
    if request.method=='POST':
        inputString=request.form['inputString'].upper()

        if inputString in pdbListDb:    #match without chain
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
        if inputString in chainListDb:  #match with chain
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
        else:
            return render_template('nomatch.html')
        
    return render_template('search.html')

    
@app.route('/webserver', methods=['GET','POST'])
def webserver():
    render_template('webserver.html')

@app.route('/download', methods=['GET','POST'])
def download():
    render_template('download.html')
    
@app.route('/help', methods=['GET','POST'])
def help():
    render_template('help.html')
    
@app.route('/contact', methods=['GET','POST'])
def contact():
    render_template('contact.html')

@app.route('/dunbrackLab')
def dunbrackLab():
    return redirect("http://dunbrack.fccc.edu")

@app.route('/<settings>/<queryname>')
def uniqueQuery(settings,queryname):
    if settings=='PDB':
        queryname=queryname.upper()
        if len(queryname)==5:
            queryname=queryname[0:-1]           #Remove chain so that only PDB id is always searched
        pdb_list=Cluster.query.filter(Cluster.pdb.contains(queryname)).all()
        pdb_count=Cluster.query.filter(Cluster.pdb.contains(queryname)).count()
        
        for items in pdb_list:
            resolution=items.resolution;
            species=items.species
    
#        pymolSession=f'downloads/pymolSessions/{pdbGroup}_{pdbGene}_{queryname}.pse.zip'
#        pymolScript=f'downloads/pymolSessionScripts/{pdbGroup}_{pdbGene}_{queryname}.zip'
#        coordinateFiles=f'downloads/coordinateFiles/{pdbGroup}_{pdbGene}_{queryname}'

        return render_template('pdbs.html',queryname=queryname,pdb_list=pdb_list,resolution=resolution,species=species,pdb_count=pdb_count)

    if settings=='fullcluster':
        cluster_list=Cluster.query.filter(Cluster.fullcluster.contains(queryname)).all()
        cluster_count=Cluster.query.filter(Cluster.fullcluster.contains(queryname)).count()
        species_list=Cluster.query.filter(Cluster.fullcluster.contains(queryname)).with_entities('species')
        species_unique=pd.Series(names[0] for names in species_list).unique()
        gene_list=Cluster.query.filter(Cluster.fullcluster.contains(queryname)).with_entities('gene')
        gene_unique=pd.Series(names[0] for names in gene_list).unique()
        return render_template('fullcluster.html',queryname=queryname,cluster_list=cluster_list,cluster_count=cluster_count,species_unique=species_unique,gene_unique=gene_unique)
    
    if settings=='CDR':
        cdr_list=Cluster.query.filter(Cluster.CDR.contains(queryname)).all()
        cdr_count=Cluster.query.filter(Cluster.CDR.contains(queryname)).count()
        species_list=Cluster.query.filter(Cluster.CDR.contains(queryname)).with_entities('species')
        species_unique=pd.Series(names[0] for names in species_list).unique()
        gene_list=Cluster.query.filter(Cluster.CDR.contains(queryname)).with_entities('gene')
        gene_unique=pd.Series(names[0] for names in gene_list).unique()
        return render_template('cdr.html',queryname=queryname,cdr_list=cdr_list,cdr_count=cdr_count,species_unique=species_unique,gene_unique=gene_unique)
    
    if settings=='CDR_length':
        cdr_length_list=Cluster.query.filter(Cluster.CDR_length.contains(queryname)).all()
        cdr_length_count=Cluster.query.filter(Cluster.CDR_length.contains(queryname)).count()
        species_list=Cluster.query.filter(Cluster.CDR_length.contains(queryname)).with_entities('species')
        species_unique=pd.Series(names[0] for names in species_list).unique()
        gene_list=Cluster.query.filter(Cluster.CDR_length.contains(queryname)).with_entities('gene')
        gene_unique=pd.Series(names[0] for names in gene_list).unique()
        return render_template('cdr_length.html',queryname=queryname,cdr_length_list=cdr_length_list,cdr_length_count=cdr_length_count,species_unique=species_unique,gene_unique=gene_unique)
    
#@app.route('/multipleQuery/<L1_cluster>/<H1_cluster>/<L2_cluster>/<H2_cluster>/<L3-cluster>/<H3_cluster>')
#def multipleQuery(L1_cluster,H1_cluster,L2_cluster,H2_cluster,L3_cluster,H3_cluster):
#    retrieve_cluster=Cluster.query.filter(Cluster.fullcluster.contains(L1_cluster)).all()
#                total_count[tabs]=Cluster.query.filter(Cluster.ligand_type.contains(ligTypeSelect),Cluster.ligand_type.notlike(dontmatch)).count()



if __name__ == '__main__':
    app.run(debug=True)
