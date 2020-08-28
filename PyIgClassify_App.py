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
    cluster=db.Column(db.Text)
    
def __init__(self,pdb_chain_cdr,datatag,pdb,original_chain,CDR,length,cluster):
    self.pdb_chain_cdr=pdb_chain_cdr;self.datatag=datatag;self.pdb=pdb;self.original_chain=original_chain
    self.CDR=CDR;self.length=length;self.cluster=cluster
    
def __repr__(self):
    return f'{self.pdb_chain_cdr} {self.datatag} {self.pdb} {self.original_chain} {self.CDR} {self.length} {self.cluster}'

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

@app.route('/formSearch', methods=['GET','POST'])
def formSearch():
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

if __name__ == '__main__':
    app.run(debug=True)
