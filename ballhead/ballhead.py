import os
import random
import re
import string
import time
from PIL import Image
from collections import OrderedDict
import ConfigParser
from flask import Flask, flash, request, redirect, url_for, render_template, session, escape, Response, send_from_directory
from werkzeug import secure_filename
from celery.result import AsyncResult
from celery import Celery
from tasks import *

configPath = '/opt/ballhead/ballhead/default.cfg'
config = ConfigParser.RawConfigParser()
config.read(os.path.join(configPath))

perl = config.get('paths', 'perl')
tripodPath = os.path.join(config.get('paths', 'install'),
                                       'tripod', 'triPOD.pl')
celery = Celery()
celery.config_from_object('celeryconfig')

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.join(config.get('paths', 'install'), 'upload')
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
        salt = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(12))
        os.mkdir(os.path.join(app.config['UPLOAD_FOLDER'], salt))
        session['out'] = os.path.join(app.config['UPLOAD_FOLDER'], salt)
        file = request.files['file']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], salt, filename))
            session['filename'] = os.path.join(app.config['UPLOAD_FOLDER'], salt, filename)
        elif not allowed_file(file.filename):
            flash(u"File type must be .txt", 'error')
            return redirect(url_for('upload'))
        flash(u"{0} was uploaded successfully!".format(os.path.basename(session['filename'])))
        session['build'] = request.form['build']
        command = [perl, tripodPath,
                   '--gender=' + request.form['gender'],
                   '--graph=png',
                   '--alpha=' + request.form['alpha'],
                   '--build=' + os.path.join(config.get('paths','install'),
                                             'ballhead', 'static', 'genome_build',
                                             request.form['build']),
                   '--' + request.form['pod'], 
                   '--' + request.form['podhd'], 
                   '--' + request.form['podmi1'], 
                   '--' + request.form['podcr'], 
                   '--out=' + session['out'],
                   session['filename']
        ]
        p = run.delay(command)
        session['celeryid'] = p.id

        return redirect(url_for('progress', id=session['celeryid']))

    return render_template('upload.html')

@app.route('/progress/<id>')
def progress(id):
    return render_template('progress.html', id=id)

@app.route('/status/<id>')
def status(id):
    """ Wait until the job is finished and report success."""
    result = AsyncResult(id, app=celery)
    while not result.ready():
        if result.failed():
            flash(u"error running triPOD. please check {0}".format(session['filename']),'error')
            return redirect(url_for('upload'))
        time.sleep(5)
    return url_for('results')
    
@app.route('/results')
def results():
    """ Format the final results page and return template."""
    outdir = session['out']
    outdirList = os.listdir(outdir)
    png = re.compile("png$")
    txt = re.compile("triPOD_Results.txt$")
    bed = re.compile("bed$")
    images = []
    for file in outdirList:
        if re.search(png, file):
            images.append(file)
        elif re.search(txt, file):
            session['txt'] = file
            with open(os.path.join(outdir, file), 'r') as txtfile:
                table = extract_table(txtfile)
        elif re.search(bed, file):
            session['bed'] = file
        else:
            continue
    

    return render_template('results.html', 
                           filename=os.path.basename(session['filename']),
                           ucsc=os.path.basename(session['out']),
                           images=reversed(images),
                           table=table,
                           tablerange=range(0,len(table['Sample']) + 1))

@app.route('/data/<file>')
def data(file):
    """ Return requested file to results page """
    if file == 'txt':
        return send_from_directory(session['out'],session['txt'])
    elif file == 'bed':
        return send_from_directory(session['out'],session['bed'])
    else:
        return send_from_directory(session['out'],file)

@app.route('/external/<id>/<filename>')
def external(id, filename):
    """ Return results independently of session data
    id = upload path directory """
    outdir = os.path.join(app.config['UPLOAD_FOLDER'],id)
    outdirList = os.listdir(outdir)
    png = re.compile("png$")
    txt = re.compile("triPOD_Results.txt$")
    bed = re.compile("bed$")
    images = []
    for file in outdirList:
        if re.search(png, file):
            images.append(file)
        elif re.search(txt, file):
            txt = file
        elif re.search(bed, file):
            bed = file
        else:
            continue

    bulkResize(outdir, 22)

    if filename == 'txt':
        return send_from_directory(os.path.join(app.config['UPLOAD_FOLDER'],id),txt)
    elif filename == 'bed':
        return send_from_directory(os.path.join(app.config['UPLOAD_FOLDER'],id),bed)
    
def extract_table(txt):
    """ Extract the useful parts of the text table from triPOD output """
    for line in txt:
        line = line.split()
        if line == []:
            continue
        elif line[0] == "Sample":
            table = OrderedDict((k,[]) for k in line)
            break
    for line in txt:
        line = line.split()
        if line == []:
            break
        else:
            for k, v in zip(table.keys(), line):
                table[k].append(v)
    return table
    
def resize(folder, fileName, factor):
    filePath = os.path.join(folder, fileName)
    im = Image.open(filePath)
    w, h  = im.size
    newIm = im.resize((int(w*factor), int(h*factor)))
    newIm.save(filePath)

def bulkResize(imageFolder, factor):
    imgExts = ["png", "bmp", "jpg"]
    for path, dirs, files in os.walk(imageFolder):
        for fileName in files:
            ext = fileName[-3:].lower()
            if ext not in imgExts:
                continue

            resize(path, fileName, factor)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
