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

class perVRegion(db.Model):
    pdb_heavy_light_chain=db.Column(db.Text,primary_key=True)
    pdb=db.Column(db.Text)
    vh_chain=db.Column(db.Text)
    vh_framework=db.Column(db.Text)
    vh_framework_seqid=db.Column(db.REAL)
    h1_cluster=db.Column(db.Text)
    h1_cluster_distance=db.Column(db.REAL)
    h1_seqid=db.Column(db.REAL)
    h2_cluster=db.Column(db.Text)
    h2_cluster_distance=db.Column(db.REAL)
    h2_seqid=db.Column(db.REAL)
    h3_cluster=db.Column(db.Text)
    h3_cluster_distance=db.Column(db.REAL)
#   h3_seqid=db.Column(db.REAL)     #No seqid for H3
    vl_chain=db.Column(db.Text)
    vl_framework=db.Column(db.Text)
    vl_framework_seqid=db.Column(db.REAL)
    l1_cluster=db.Column(db.Text)
    l1_cluster_distance=db.Column(db.REAL)
    l1_seqid=db.Column(db.REAL)
    l2_cluster=db.Column(db.Text)
    l2_cluster_distance=db.Column(db.REAL)
    l2_seqid=db.Column(db.REAL)
    l3_cluster=db.Column(db.Text)
    l3_cluster_distance=db.Column(db.REAL)
#   l3_seqid=db.Column(db.REAL)     #No seqid for L3
    
    def __init__(self,pdb_heavy_light_chain,pdb,vh_chain,vh_framework,vh_framework_seqid,h1_cluster,h1_cluster_distance,h1_seqid,\
                 h2_cluster,h2_cluster_distance,h2_seqid,h3_cluster,h3_cluster_distance,vl_chain,vl_framework,vl_framework_seqid,\
                 l1_cluster,l1_cluster_distance,l1_seqid,l2_cluster,l2_cluster_distance,l2_seqid,l3_cluster,l3_cluster_distance):
        
        self.pdb_heavy_light_chain=pdb_heavy_light_chain;self.pdb=pdb;self.vh_chain=vh_chain;self.vh_framework=vh_framework;self.vh_framework_seqid=vh_framework_seqid;
        self.h1_cluster=h1_cluster;self.h1_cluster_distance=h1_cluster_distance;self.h1_seqid=h1_seqid;
        self.h2_cluster=h2_cluster;self.h2_cluster_distance=h2_cluster_distance;self.h2_seqid=h2_seqid;self.h3_cluster=h3_cluster;self.h3_cluster_distance=h3_cluster_distance;
        self.vl_chain=vl_chain;self.vl_framework=vl_framework;self.vl_framework_seqid=vl_framework_seqid;
        self.l1_cluster=l1_cluster;self.l1_cluster_distance=l1_cluster_distance;self.l1_seqid=l1_seqid;self.l2_cluster=l2_cluster;self.l2_cluster_distance=l2_cluster_distance;
        self.l2_seqid=l2_seqid;self.l3_cluster=l3_cluster;self.l3_cluster_distance=l3_cluster_distance


    def __repr__(self):
        return f'{self.pdb_heavy_light_chain} {self.pdb} {self.vh_chain} {self.vh_framework} {self.vh_framework_seqid} {self.h1_cluster} {self.h1_cluster_distance}\
            {self.h1_seqid} {self.h2_cluster} {self.h2_cluster_distance} {self.h2_seqid} {self.h3_cluster} {self.h3_cluster_distance} {self.vl_chain} {self.vl_framework}\
            {self.vl_framework_seqid} {self.l1_cluster} {self.l1_cluster_distance} {self.l1_seqid} {self.l2_cluster} {self.l2_cluster_distance} {self.l2_seqid}\
            {self.l3_cluster} {self.l3_cluster_distance}'

class perCDR(db.Model):
    pdb_chain_cdr=db.Column(db.Text,primary_key=True)
    cdr=db.Column(db.Text)
    cdr_length=db.Column(db.Text)
    pdb=db.Column(db.Text)
    chain=db.Column(db.Text)
    resolution=db.Column(db.REAL)
    aho_resnum=db.Column(db.Text)
    author_resnum=db.Column(db.Text)
    sequence=db.Column(db.Text)
    germline_sequence=db.Column(db.Text)
    gene=db.Column(db.Text)
    pdb_species=db.Column(db.Text)
    frame_germline=db.Column(db.Text)
    cluster=db.Column(db.Text)
    distance=db.Column(db.REAL)
    cdr_germline=db.Column(db.Text)
    cdr_seqid=db.Column(db.Text)   #There are - in records so can not use REAL e.g. 1a2y
    rama4=db.Column(db.Text)
    beta_turns=db.Column(db.Text)
    minimum_edia=db.Column(db.Text)
    keywords=db.Column(db.Text)
    
    def __init__(self,pdb_chain_cdr,cdr,cdr_length,pdb,chain,resolution,aho_resnum,author_resnum,sequence,germline_sequence,gene,pdb_species,\
                 frame_germline,cluster,distance,cdr_germline,cdr_seqid,rama4,beta_turns,minimum_edia,keywords):
        
        self.pdb_chain_cdr=pdb_chain_cdr;self.cdr=cdr;self.cdr_length=cdr_length;self.pdb=pdb;self.chain=chain;self.resolution=resolution;
        self.aho_resnum=aho_resnum;self.author_resnum=author_resnum;self.sequence=sequence;self.germline_sequence=germline_sequence;
        self.gene=gene;self.pdb_species=pdb_species;self.frame_germline=frame_germline;self.cluster=cluster;self.distance=distance;
        self.cdr_germline=cdr_germline;self.cdr_seqid=cdr_seqid;self.rama4=rama4;self.beta_turns=beta_turns;self.minimum_edia=minimum_edia;
        self.keywords=keywords
    
    def __repr__(self):
        return f'{self.pdb_chain_cdr} {self.cdr} {self.cdr_length} {self.pdb} {self.chain} {self.resolution} {self.aho_resnum} {self.author_resnum}\
            {self.sequence} {self.germline_sequence} {self.gene} {self.pdb_species} {self.frame_germline} {self.cluster} {self.distance} {self.cdr_germline}\
            {self.cdr_seqid} {self.rama4} {self.beta_turns} {self.minimum_edia} {self.keywords}'
    

def create_lists():
    pdbListDb=list();chainListDb=list();clusterListDb=list();clusterChainCount=dict();cluster_cdr_list=list()
    for pdbs in perCDR.query.with_entities(perCDR.pdb_chain_cdr):
        pdbid,chainid,cdr=pdbs[0].split('_')
        pdbListDb.append(pdbid);chainListDb.append(pdbid+chainid);
    for cluster in perCDR.query.with_entities(perCDR.cluster):
        clusterListDb.append(cluster[0])
        cluster_cdr_list.append(cluster[0])
        cdr_length='-'.join(cluster[0].split('-')[0:2])  #Get the first two element after the split like H1-10 from H1-10-1
        cluster_cdr_list.append(cdr_length)     #This list contains both clusters and cdr-lengths
        
    clusterListDb=list(natsorted(set(sorted(clusterListDb))))
    cluster_cdr_list=list(natsorted(set(sorted(cluster_cdr_list))))
    pdbListDb=list(natsorted(set(sorted(pdbListDb))))
    
    for cluster in cluster_cdr_list:
            clusterChainCount[cluster]=perCDR.query.filter(perCDR.cluster.contains(cluster)).count()  #Incorrect count for cdr-length due to substring mismatch
        
    clusterChainCount['H1-All']=perCDR.query.filter(perCDR.cluster.contains('H1')).count()
    clusterChainCount['H2-All']=perCDR.query.filter(perCDR.cluster.contains('H2')).count()
    clusterChainCount['H3-All']=perCDR.query.filter(perCDR.cluster.contains('H3')).count()
    clusterChainCount['L1-All']=perCDR.query.filter(perCDR.cluster.contains('L1')).count()
    clusterChainCount['L2-All']=perCDR.query.filter(perCDR.cluster.contains('L2')).count()
    clusterChainCount['L3-All']=perCDR.query.filter(perCDR.cluster.contains('L3')).count()
    return pdbListDb,chainListDb,clusterListDb, clusterChainCount, cluster_cdr_list

def percent_loop_length(chain_count,queryname):
    (cdr,length,cluster)=queryname.split('-')      #How to split in clusters with 'cis' in the name? 
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

def write_text_file(sublist,tsvFile,cdr_format):
    print(cdr_format)
    fhandle_textFile=open(f'{pwd}/static/{tsvFile}','w')
    if cdr_format:     #Header for perCDR files
        fhandle_textFile.write('CDR\tCDR Length\tPDB\tChain\tResolution\tAHO Resnum\tAuthor Resnum\tSequence\tGermline Sequence\tGene\tPDB Species\tFrame Germline\tCluster\tDistance\tCDR Germline\tCDR Seqid\tRama4\tBeta Turns\tMinimum EDIA\n')
        for item in sublist:
            fhandle_textFile.write(f'{item.cdr}\t{item.cdr_length}\t{item.pdb}\t{item.chain}\t{item.resolution}\t{item.aho_resnum}\t\
                                   {item.author_resnum}\t{item.sequence}\t{item.germline_sequence}\t{item.gene}\t{item.pdb_species}\t\
                                   {item.frame_germline}\t{item.cluster}\t{item.distance}\t{item.cdr_germline}\t{item.cdr_seqid}\t\
                                   {item.rama4}\t{item.beta_turns}\t{item.minimum_edia}\n')
    
    else:    #Header for perVRegion files
        fhandle_textFile.write('PDB\tVH Chain\tVH Framework\tVH Framework Seqid\tH1 Cluster\tH1 Cluster Distance\tH1 Seqid\tH2 Cluster\tH2 Cluster Distance\tH2 Seqid\tH3 Cluster\tH3 Cluster Distance\tVL Chain\tVL Framework\tVL Framework Seqid\tL1 Cluster\tL1 Cluster Distance\tL1 Seqid\tL2 Cluster\tL2 Cluster Distance\tL2 Seqid\tL3 Cluster\tL3 Cluster Distance\n')
        for item in sublist:
            fhandle_textFile.write(f'{item.pdb}\t{item.vh_chain}\t{item.vh_framework}\t{item.vh_framework_seqid}\t{item.h1_cluster}\t{item.h1_cluster_distance}\t{item.h1_seqid}\t\
                 {item.h2_cluster}\t{item.h2_cluster_distance}\t{item.h2_seqid}\t{item.h3_cluster}\t{item.h3_cluster_distance}\t{item.vl_chain}\t{item.vl_framework}\t{item.vl_framework_seqid}\t\
                 {item.l1_cluster}\t{item.l1_cluster_distance}\t{item.l1_seqid}\t{item.l2_cluster}\t{item.l2_cluster_distance}\t{item.l2_seqid}\t{item.l3_cluster}\t{item.l3_cluster_distance}\n')
    fhandle_textFile.close()
    
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
    retrieve_all=perVRegion.query.all()
    return render_template('browse.html',retrieve_all=retrieve_all)

@app.route('/statistics')
def statistics():
    pdb_species_unique=dict();gene_unique=dict();entries=dict()
    (pdbListDb,chainListDb,clusterListDb,clusterChainCount,cluster_cdr_list)=create_lists()
    for cluster in clusterListDb:
        #cluster_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).all()
        entries[cluster]=perCDR.query.filter(perCDR.cluster.contains(cluster)).count()
        pdb_species_list=perCDR.query.filter(perCDR.cluster.contains(cluster)).with_entities('pdb_species')
        pdb_species_unique[cluster]=pd.Series(names[0] for names in pdb_species_list).unique()
        gene_list=perCDR.query.filter(perCDR.cluster.contains(cluster)).with_entities('gene')
        gene_unique[cluster]=pd.Series(names[0] for names in gene_list).unique()


    return render_template('statistics.html',clusterListDb=clusterListDb,entries=entries,pdb_species_unique=pdb_species_unique,gene_unique=gene_unique)

@app.route('/formSearch', methods=['GET','POST'])
def formSearch():
    (pdbListDb,chainListDb,clusterListDb,clusterChainCount,cluster_cdr_list)=create_lists()

    if request.method=='POST':
        inputString=request.form['pdb_select'].upper()

        if inputString in pdbListDb:    #match without chain
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
        #if inputString in chainListDb:  #match with chain
        #    return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
        else:
            return render_template('nomatch.html')

    return render_template('search.html', pdbListDb=pdbListDb, clusterListDb=clusterListDb, clusterChainCount=clusterChainCount, cluster_cdr_list=cluster_cdr_list)

@app.route('/formSearchMultiple',methods=['GET','POST'])
def formSearchMultiple():
    (pdbListDb,chainListDb,clusterListDb,clusterChainCount,cluster_cdr_list)=create_lists()
    if request.method=='POST':
        cdr_select=request.form['cdr_select']
        
        if cdr_select=='All':
            return redirect(url_for('browse'))
        elif 'All' in cdr_select:
            cdr_select=cdr_select[0:-4]
            return redirect(url_for('uniqueQuery',settings='cdr',queryname=cdr_select))
        elif cdr_select in clusterListDb:
            return redirect(url_for('uniqueQuery',settings='cluster',queryname=cdr_select))
        else:
            return redirect(url_for('uniqueQuery',settings='cdr_length',queryname=cdr_select))
        
    return render_template ('search.html', pdbListDb=pdbListDb, clusterListDb=clusterListDb, clusterChainCount=clusterChainCount, cluster_cdr_list=cluster_cdr_list)


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
        
        vh_gene=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('H')).with_entities('gene').first()[0]
        vl_gene=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('L')).with_entities('gene').first()[0]
        vh_chainid=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('H')).with_entities('chain').first()[0]   #Might not be always true; what if there are two heavy chains in the PDB?
        vl_chainid=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('L')).with_entities('chain').first()[0]
        vh_species=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('H')).with_entities('pdb_species').first()[0]
        vl_species=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('L')).with_entities('pdb_species').first()[0]
        tsvFile=f'downloads/text-files/{queryname}.tab'
        write_text_file(pdb_list,tsvFile,True)

        return render_template('pdbs.html',queryname=queryname,pdb_list=pdb_list,pdb_count=pdb_count,pdb_resolution=pdb_resolution,\
                               vh_gene=vh_gene,vl_gene=vl_gene,vh_chainid=vh_chainid,vl_chainid=vl_chainid,vh_species=vh_species,vl_species=vl_species,tsvFile=tsvFile)

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
        pdb_species_unique=pd.Series(names[0] for names in pdb_species_list).sort_values().unique()
        gene_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('gene')
        gene_unique=pd.Series(names[0] for names in gene_list).unique()
        rama_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('rama4')
        (most_common_rama_count,most_common_rama_string )=most_common_rama(rama_list)
        tsvFile=f'downloads/text-files/{queryname}.tab'
        write_text_file(cluster_list,tsvFile,True)
        return render_template('cluster.html',queryname=queryname,cluster_list=cluster_list,cluster_cdr_length=cluster_cdr_length,\
                               cluster_cdr=cluster_cdr,pdb_unique_count=pdb_unique_count,chain_count=chain_count,per_loop=per_loop,\
                                   seq_unique_count=seq_unique_count,pdb_species_unique=pdb_species_unique,gene_unique=gene_unique,\
                                       most_common_rama_count=most_common_rama_count,most_common_rama_string=most_common_rama_string,tsvFile=tsvFile)

    if settings=='cdr':
        cdr_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).all()
        pdb_count=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('pdb')
        pdb_unique_count=len(pd.Series(names[0] for names in pdb_count).unique())
        chain_count=perCDR.query.filter(perCDR.cdr.contains(queryname)).count()    #This is number of chains with a given cluster
        pdb_species_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('pdb_species')
        pdb_species_unique=pd.Series(names[0] for names in pdb_species_list).unique()
        gene_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('gene')
        gene_unique=pd.Series(names[0] for names in gene_list).unique()
        tsvFile=f'downloads/text-files/{queryname}.tab'
        write_text_file(cdr_list,tsvFile,True)
        return render_template('cdr.html',queryname=queryname,cdr_list=cdr_list,pdb_unique_count=pdb_unique_count,chain_count=chain_count,\
                               pdb_species_unique=pdb_species_unique,gene_unique=gene_unique,tsvFile=tsvFile)

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
        tsvFile=f'downloads/text-files/{queryname}.tab'
        write_text_file(cdr_length_list,tsvFile,True)
        return render_template('cdr_length.html',queryname=queryname,cdr_cdr_length=cdr_cdr_length,cdr_length_list=cdr_length_list,\
                               pdb_unique_count=pdb_unique_count,chain_count=chain_count,pdb_species_unique=pdb_species_unique,\
                                   gene_unique=gene_unique,tsvFile=tsvFile)

    if settings=='germline':
        queryname=queryname.split('*')[0]
        germline_list=perVRegion.query.filter(perVRegion.vh_framework.contains(queryname)).all()
        
        if germline_list:    #First check for VH germline and then for VL
            tsvFile=f'downloads/text-files/{queryname}.tab'
            write_text_file(germline_list,tsvFile,False)
            return render_template('germline.html',queryname=queryname,germline_list=germline_list,tsvFile=tsvFile)
        else:
            germline_list=perVRegion.query.filter(perVRegion.vl_framework.contains(queryname)).all()
            tsvFile=f'downloads/text-files/{queryname}.tab'
            write_text_file(germline_list,tsvFile,False)
            return render_template('germline.html',queryname=queryname,germline_list=germline_list,tsvFile=tsvFile)
        
    if settings=='frame_germline':
        queryname=queryname.split('*')[0]
        germline_list=perCDR.query.filter(perCDR.frame_germline.contains(queryname)).all()
        tsvFile=f'downloads/text-files/{queryname}.tab'
        write_text_file(germline_list,tsvFile,True)
           
        return render_template('germline.html',queryname=queryname,germline_list=germline_list,tsvFile=tsvFile)
    
    if settings=='cdr_germline':
        queryname=queryname.split('*')[0]
        germline_list=perCDR.query.filter(perCDR.cdr_germline.contains(queryname)).all()
        tsvFile=f'downloads/text-files/{queryname}.tab'
        write_text_file(germline_list,tsvFile,True)
           
        return render_template('germline.html',queryname=queryname,germline_list=germline_list,tsvFile=tsvFile)
            

if __name__ == '__main__':
    app.run(debug=True)
