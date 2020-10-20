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

class perCDR(db.Model):
    pdb_chain_cdr=db.Column(db.Text,primary_key=True)    #Might not be unique, check with Roland
    cdr=db.Column(db.Text)
    cdr_length=db.Column(db.Text)
    pdb=db.Column(db.Text)
    chain=db.Column(db.Text)
    aho_resnum=db.Column(db.Text)
    author_resnum=db.Column(db.Text)
    sequence=db.Column(db.Text)
    germline_sequence=db.Column(db.Text)
    gene=db.Column(db.Text)
    pdb_species=db.Column(db.Text)
    cluster=db.Column(db.Text)
    distance=db.Column(db.REAL)
    cdr_germline=db.Column(db.Text)
    cdr_seqid=db.Column(db.REAL)
    rama4=db.Column(db.Text)
    beta_turns=db.Column(db.Text)
    minimum_edia=db.Column(db.REAL)
    keywords=db.Column(db.Text)


    def __init__(self,pdb_chain_cdr,cdr,cdr_length,pdb,chain,aho_resnum,author_resnum,sequence,germline_sequence,gene,pdb_species,cluster,distance,cdr_germline,cdr_seqid,\
                rama4,beta_turns,minimum_edia,keywords):
                self.pdb_chain_cdr=pdb_chain_cdr;self.cdr=cdr;self.cdr_length=cdr_length;self.pdb=pdb;self.chain=chain;self.aho_resnum=aho_resnum;self.author_resnum=author_resnum;\
                self.sequence=sequence;self.germline_sequence=germline_sequence;self.gene=gene;self.pdb_species=pdb_species;self.cluster=cluster;self.distance=distance;\
                self.cdr_germline=cdr_germline;self.cdr_seqid=cdr_seqid;self.rama4=rama4;self.beta_turns=beta_turns;self.minimum_edia=minimum_edia;self.keywords=keywords


    def __repr__(self):
        return f'{self.pdb_chain_cdr} {self.cdr} {self.cdr_length} {self.pdb} {self.chain} {self.aho_resnum} {self.author_resnum} {self.sequence} {self.germline_sequence}\
                 {self.gene} {self.pdb_species} {self.cluster} {self.distance} {self.cdr_germline} {self.cdr_seqid} {self.rama4} {self.beta_turns}\
                 {self.minimum_edia} {self.keywords}'

class PDBrows(db.Model):
    pdb=db.Column(db.Text,primary_key=True)
    heavy_chain=db.Column(db.Text)
    light_chain=db.Column(db.Text)
    h1cluster=db.Column(db.Text)
    l1cluster=db.Column(db.Text)
    h2cluster=db.Column(db.Text)
    l2cluster=db.Column(db.Text)
    h3cluster=db.Column(db.Text)
    l3cluster=db.Column(db.Text)

    def __init__(self,pdb,heavy_chain,light_chain,h1cluster,l1cluster,h2cluster,l2cluster,h3cluster,l3cluster):
        self.pdb=pdb;self.heavy_chain=heavy_chain;self.light_chain=light_chain;self.h1cluster=h1cluster;self.l1cluster=l1cluster;
        self.h2cluster=h2cluster;self.l2cluster=l2cluster;self.h3cluster=h3cluster;self.l3cluster=l3cluster

    def __repr__(self):
        return f'{self.pdb} {self.heavy_chain} {self.light_chain} {self.h1cluster} {self.l1cluster} {self.h2cluster} {self.l2cluster}\
            {self.h3cluster} {self.l3cluster}'

def create_lists():
    pdbListDb=list();chainListDb=list();clusterListDb=list()
    for pdbs in perCDR.query.with_entities(perCDR.pdb_chain_cdr):
        pdbid,chainid,cdr=pdbs[0].split('_')
        pdbListDb.append(pdbid);chainListDb.append(pdbid+chainid);
    for cluster in perCDR.query.with_entities(perCDR.cluster):
        clusterListDb.append(cluster[0])
    clusterListDb=list(set(sorted(clusterListDb)))
    return pdbListDb,chainListDb,clusterListDb

def percent_loop_length(chain_count,queryname):
    (cdr,length,cluster)=queryname.split('-')
    cdr_length=cdr+'_'+length
    cdr_length_count=perCDR.query.filter(perCDR.cdr_length.contains(cdr_length)).count()
    per_loop=round((chain_count/cdr_length_count)*100,2)
    return per_loop

def most_common_rama(rama_list):
    rama_dict=dict()
    for item in rama_list:
        rama_dict[item]=0
    for item in rama_list:
        rama_dict[item]+=1
    most_common_rama_count=max(rama_dict.values())
    most_common_rama_string=max(rama_dict, key=lambda key: rama_dict[key]).rama4
    return most_common_rama_count,most_common_rama_string    
    
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
    retrieve_all=PDBrows.query.all()
    return render_template('browse.html',retrieve_all=retrieve_all)

@app.route('/statistics')
def statistics():
    pdb_species_unique=dict();gene_unique=dict();entries=dict()
    (pdbListDb,chainListDb,clusterListDb)=create_lists()
    for cluster in clusterListDb:
        #cluster_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).all()
        entries[cluster]=perCDR.query.filter(perCDR.cluster.contains(cluster)).count()
        pdb_species_list=perCDR.query.filter(perCDR.cluster.contains(cluster)).with_entities('pdb_species')
        pdb_species_unique[cluster]=pd.Series(names[0] for names in pdb_species_list).unique()
        gene_list=perCDR.query.filter(perCDR.cluster.contains(cluster)).with_entities('gene')
        gene_unique[cluster]=pd.Series(names[0] for names in gene_list).unique()


    return render_template('statistics.html',clusterListDb=natsorted(clusterListDb,key=str),entries=entries,pdb_species_unique=pdb_species_unique,gene_unique=gene_unique)

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

@app.route('/formSearchMultiple',methods=['GET','POST'])
def formSearchMultiple():
    if request.method=='POST':
        h1Select=request.form['h1Select']
        l1Select=request.form['l1Select']
        h2Select=request.form['h2Select']
        l2Select=request.form['l2Select']
        h3Select=request.form['h3Select']
        l3Select=request.form['l3Select']
        return redirect(url_for('multipleQuery',h1Select=h1Select,l1Select=l1Select,h2Select=h2Select,l2Select=l2Select,h3Select=h3Select,l3Select=l3Select))
    return render_template ('search.html')

@app.route('/formSearchMultipleCDR',methods=['GET','POST'])
def formSearchMultipleCDR():
    if request.method=='POST':
        h1CDRSelect=request.form['h1CDRSelect']
        l1CDRSelect=request.form['l1CDRSelect']
        h2CDRSelect=request.form['h2CDRSelect']
        l2CDRSelect=request.form['l2CDRSelect']
        h3CDRSelect=request.form['h3CDRSelect']
        l3CDRSelect=request.form['l3CDRSelect']
        return redirect(url_for('multipleQueryCDR',h1CDRSelect=h1CDRSelect,l1CDRSelect=l1CDRSelect,h2CDRSelect=h2CDRSelect,l2CDRSelect=l2CDRSelect,h3CDRSelect=h3CDRSelect,\
        l3CDRSelect=l3CDRSelect))
    return render_template ('search.html')


@app.route('/webserver', methods=['GET','POST'])
def webserver():
    return render_template('webserver.html')

@app.route('/download')
def download():
    return render_template('download.html')

@app.route('/help')
def help():
    return render_template('help.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/dunbrackLab')
def dunbrackLab():
    return redirect("http://dunbrack.fccc.edu")


@app.route('/<settings>/<queryname>')
def uniqueQuery(settings,queryname):
    if settings=='PDB':
        queryname=queryname.upper()
        if len(queryname)==5:
            queryname=queryname[0:-1]           #Remove chain so that only PDB id is always searched
        pdb_list=perCDR.query.filter(perCDR.pdb.contains(queryname)).all()
        pdb_count=perCDR.query.filter(perCDR.pdb.contains(queryname)).count()
        pdb_resolution=perCDR.query.filter(perCDR.pdb.contains(queryname)).with_entities('resolution').first()[0]
        pdb_rfactor=perCDR.query.filter(perCDR.pdb.contains(queryname)).with_entities('rfactor').first()[0]
        vh_gene=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('H')).with_entities('gene').first()[0]
        vl_gene=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('L')).with_entities('gene').first()[0]
        vh_chainid=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('H')).with_entities('chain').first()[0]   #Might not be always true; what if there are two heavy chains in the PDB?
        vl_chainid=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('L')).with_entities('chain').first()[0]
        vh_species=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('H')).with_entities('pdb_species').first()[0]
        vl_species=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('L')).with_entities('pdb_species').first()[0]


        return render_template('pdbs.html',queryname=queryname,pdb_list=pdb_list,pdb_count=pdb_count,pdb_resolution=pdb_resolution,pdb_rfactor=pdb_rfactor,\
                               vh_gene=vh_gene,vl_gene=vl_gene,vh_chainid=vh_chainid,vl_chainid=vl_chainid,vh_species=vh_species,vl_species=vl_species)

    if settings=='cluster':
        cluster_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).all()
        cluster_cdr_length=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('cdr_length').first()[0]
        cluster_cdr=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('cdr').first()[0]
        pdb_count=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('pdb')
        pdb_unique_count=len(pd.Series(names[0] for names in pdb_count).unique())
        chain_count=perCDR.query.filter(perCDR.cluster.contains(queryname)).count()    #This is number of chains with a given cluster
        per_loop=percent_loop_length(chain_count,queryname)
        seq_in_cluster=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('sequence')
        seq_unique_count=len(pd.Series(names[0] for names in seq_in_cluster).unique())
        pdb_species_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('pdb_species')
        pdb_species_unique=pd.Series(names[0] for names in pdb_species_list).unique()
        gene_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('gene')
        gene_unique=pd.Series(names[0] for names in gene_list).unique()
        rama_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('rama4')
        (most_common_rama_count,most_common_rama_string )=most_common_rama(rama_list)
        return render_template('cluster.html',queryname=queryname,cluster_list=cluster_list,cluster_cdr_length=cluster_cdr_length,cluster_cdr=cluster_cdr,pdb_unique_count=pdb_unique_count,chain_count=chain_count,per_loop=per_loop,seq_unique_count=seq_unique_count,pdb_species_unique=pdb_species_unique,gene_unique=gene_unique,most_common_rama_count=most_common_rama_count,most_common_rama_string=most_common_rama_string)

    if settings=='cdr':
        cdr_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).all()
        pdb_count=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('pdb')
        pdb_unique_count=len(pd.Series(names[0] for names in pdb_count).unique())
        chain_count=perCDR.query.filter(perCDR.cdr.contains(queryname)).count()    #This is number of chains with a given cluster
        pdb_species_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('pdb_species')
        pdb_species_unique=pd.Series(names[0] for names in pdb_species_list).unique()
        gene_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('gene')
        gene_unique=pd.Series(names[0] for names in gene_list).unique()
        return render_template('cdr.html',queryname=queryname,cdr_list=cdr_list,pdb_unique_count=pdb_unique_count,chain_count=chain_count,pdb_species_unique=pdb_species_unique,gene_unique=gene_unique)

    if settings=='cdr_length':
        cdr_length_list=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).all()
        cdr_cdr_length=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).with_entities('cdr').first()[0]
        pdb_count=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).with_entities('pdb')
        pdb_unique_count=len(pd.Series(names[0] for names in pdb_count).unique())
        chain_count=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).count()    #This is number of chains with a given cluster
        pdb_species_list=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).with_entities('pdb_species')
        pdb_species_unique=pd.Series(names[0] for names in pdb_species_list).unique()
        gene_list=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).with_entities('gene')
        gene_unique=pd.Series(names[0] for names in gene_list).unique()
        return render_template('cdr_length.html',queryname=queryname,cdr_cdr_length=cdr_cdr_length,cdr_length_list=cdr_length_list,pdb_unique_count=pdb_unique_count,chain_count=chain_count,pdb_species_unique=pdb_species_unique,gene_unique=gene_unique)

@app.route('/multipleQuery/<h1Select>/<l1Select>/<h2Select>/<l2Select>/<h3Select>/<l3Select>')
def multipleQuery(h1Select,l1Select,h2Select,l2Select,h3Select,l3Select):
        if h1Select=='All':
            h1Select=''    #Does not match None, so skips the rows in which h1cluster value is missing - change in future
        if l1Select=='All':
            l1Select=''
        if h2Select=='All':
            h2Select=''
        if l2Select=='All':
            l2Select=''
        if h3Select=='All':
            h3Select=''
        if l3Select=='All':
            l3Select=''

        cluster_list=PDBrows.query.filter(PDBrows.h1cluster.contains(h1Select),PDBrows.l1cluster.contains(l1Select),\
                                          PDBrows.h2cluster.contains(h2Select),PDBrows.l2cluster.contains(l2Select),\
                                          PDBrows.h3cluster.contains(h3Select),PDBrows.l3cluster.contains(l3Select),).all()

        return render_template('clusterquery.html',cluster_list=cluster_list)

@app.route('/multipleQueryCDR/<h1CDRSelect>/<l1CDRSelect>/<h2CDRSelect>/<l2CDRSelect>/<h3CDRSelect>/<l3CDRSelect>')
def multipleQueryCDR(h1CDRSelect,l1CDRSelect,h2CDRSelect,l2CDRSelect,h3CDRSelect,l3CDRSelect):
        if h1CDRSelect=='All':
            h1CDRSelect=''    #Does not match None, so skips the rows in which h1cluster value is missing - change in future
        if l1CDRSelect=='All':
            l1CDRSelect=''
        if h2CDRSelect=='All':
            h2CDRSelect=''
        if l2CDRSelect=='All':
            l2CDRSelect=''
        if h3CDRSelect=='All':
            h3CDRSelect=''
        if l3CDRSelect=='All':
            l3CDRSelect=''

        cluster_list=PDBrows.query.filter(PDBrows.h1cluster.contains(h1Select),PDBrows.l1cluster.contains(l1Select),\
                                          PDBrows.h2cluster.contains(h2Select),PDBrows.l2cluster.contains(l2Select),\
                                          PDBrows.h3cluster.contains(h3Select),PDBrows.l3cluster.contains(l3Select),).all()

        return render_template('clusterquery.html',cluster_list=cluster_list)


if __name__ == '__main__':
    app.run(debug=True)
