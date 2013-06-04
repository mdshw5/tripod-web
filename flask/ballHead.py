import os
import random
import re
import string
import time
import ConfigParser
from flask import Flask, flash, request, redirect, url_for, render_template, session, escape, Response
from werkzeug import secure_filename
from celery.result import AsyncResult
from celery import Celery
from tasks import *

celery = Celery()
celery.config_from_object('celeryconfig')

config = ConfigParser.RawConfigParser()
config.read(os.path.join('/opt/tripod/flask','default.cfg'))

perl = config.get('paths', 'perl')
tripodPath = os.path.join(config.get('paths','install'),'triPOD.pl')

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = config.get('paths', 'upload')
app.config['MAX_CONTENT_LENGTH'] = config.get('param', 'max_upload')
app.secret_key = config.get('param', 'secret_key')

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in set(['txt', 'csv', 'tsv'])

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['GET', 'POST'])
def upload():
    """ Display upload page template. Handle file upload and 
    triPOD parameter validation. 
    """
    if request.method == 'POST':
        file = request.files['file']

        if file and allowed_file(file.filename):
            salt = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(12))
            os.mkdir(os.path.join(app.config['UPLOAD_FOLDER'], salt))
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], salt, filename))
            filename = os.path.join(app.config['UPLOAD_FOLDER'], salt, filename)
            flash(u"{0} was uploaded successfully!".format(filename))
            session['output'] = os.path.join(app.config['UPLOAD_FOLDER'], salt)
            session['build'] = request.form['build']
            command = [perl, tripodPath,
                       '--gender=' + request.form['gender'],
                       '--graph=png',
                       '--alpha=' + request.form['alpha'],
                       '--build=' + os.path.join(config.get('paths','install'),request.form['build']),
                       '--' + request.form['pod'], 
                       '--' + request.form['podhd'], 
                       '--' + request.form['podmi1'], 
                       '--' + request.form['podcr'], 
                       '--out=' + session['output'],
                       filename
            ]
            p = run.delay(command)
            session['celeryid'] = p.id
            
            return redirect(url_for('status', id=session['celeryid']))

        elif not allowed_file(file.filename):
            flash(u"File type must be .txt .csv or .tsv", 'error')
            return redirect(url_for('upload'))

    return render_template('upload.html')

@app.route('/progress/<id>')
def progress(id):
    return render_template('progress.html', id=id)

@app.route('/status/<id>')
def status(id):
    """ Wait until the job is finished and report success."""
    result = AsyncResult(id, app=celery)
    while not result.ready():
        time.sleep(2)
    return redirect(url_for('results'))
    
@app.route('/results')
def results():
    """ Format the final results page and return template."""
    outdir = os.listdir(session['output'])
    png = re.compile("png$")
    txt = re.compile("txt$")
    bed = re.compile("bed$")
    images = []
    for file in outdir:
        if re.search(png, file):
            images.append(file)
        elif re.search(txt, file):
            txtResults = file
        elif re.search(bed, file):
            bedResults = item
        else:
            continue
    table = open(session['output'] + '/' + txtResults, 'r+')
    txtTable = table.read().splitlines()
    del txtTable[0:6] # remove before table
    table.close()
    for _ in range(4):
        del txtTable[-1] #remove after table
    splitTable = [line.split() for line in txtTable]
    headerTable = splitTable.pop(0)
    chroms = list([int(i[1]) for i in splitTable])
    uniqChroms = list(set(chroms))
    freqs = [(i, chroms.count(i)) for i in uniqChroms]
    # Find the largest anomaly detected, so we can pass coords to UCSC
    starts = list([int(i[2]) for i in splitTable])
    ends = list([int(i[3]) for i in splitTable])
    bp = list([int(i[8]) for i in splitTable])
    maxBp = max(bp)
    maxIndex = bp.index(maxBp)
    maxStart = starts[maxIndex]
    maxEnd = ends[maxIndex]
    maxChr = chroms[maxIndex]
    ucscPos = 'chr' + str(maxChr) + ':' + str(maxStart) + '-' + str(maxEnd)
    ucscURL = outdir + '/' + bedResults

    return render_template('results.html', filename, txtResults, bedResults, build, ucscPos, ucscURL)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
