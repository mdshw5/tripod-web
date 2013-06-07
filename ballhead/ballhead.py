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

installpath = config.get('paths', 'install')
perl = config.get('paths', 'perl')
tripod = os.path.join(installpath,
                      'tripod', 'triPOD.pl')
celery = Celery()
celery.config_from_object('celeryconfig')

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.join(installpath, 'upload')
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
        outdir = os.path.join(app.config['UPLOAD_FOLDER'], salt)
        os.mkdir(outdir)
        file = request.files['file']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], salt, filename)
            file.save(filepath)
            basename = os.path.basename(filepath)
            flash(u"{0} was uploaded successfully!".format(basename))
        elif not allowed_file(file.filename):
            flash(u"File type must be .txt", 'error')
            return redirect(url_for('upload'))

        command = OrderedDict([('bin',perl),
                   ('script',tripod),
                   ('cores','--cores=2'),
                   ('gender','--gender=' + request.form['gender']),
                   ('graphics','--graph=png'),
                   ('alpha','--alpha=' + request.form['alpha']),
                   ('build','--build=' + os.path.join(installpath, 
                                                     'ballhead', 
                                                     'static', 
                                                     request.form['build'])),
                   ('pod','--' + request.form['pod']), 
                   ('podhd','--' + request.form['podhd']),
                   ('podmi1','--' + request.form['podmi1']), 
                   ('podcr','--' + request.form['podcr']),
                   ('out','--out=' + outdir),
                   ('filepath',filepath)
               ])

        p = run.delay(command)
        celeryid = p.id

        return redirect(url_for('progress', id=celeryid, name=basename))

    return render_template('upload.html')

@app.route('/progress/<id>')
def progress(id):
              return render_template('progress.html', id=id, name=request.args.get('name'))

@app.route('/status/<id>')
def status(id):
    """ Wait until the job is finished and report success."""
    result = AsyncResult(id, app=celery)
    while not result.ready():
        if result.failed():
            flash(u"error running triPOD. please check input file", 'error')
            return redirect(url_for('upload'))
        time.sleep(5)
    return url_for('results', id=id)
    
@app.route('/results/<id>')
def results(id):
    result = AsyncResult(id, app=celery)
    """ Format the final results page and return template."""
    command, exitstatus, stdout, stderr = result.get()
    if exitstatus == '3':
        flash(u"Please check your input file: {0}".format(stdout), 'error')
    outdir = command['out'].split('=')[-1]

    if not any([re.search('.resize.png', file) for file in os.listdir(outdir)]):
        bulkResize(outdir, width=640, height=480)

    png = re.compile("png$")
    txt = re.compile("triPOD_Results.txt$")
    bed = re.compile("bed$")

    images = []
    for path, dirs, files in os.walk(outdir):
        for file in files:
            if re.search(png, file):
                images.append(file)
            elif re.search(txt, file):
                with open(os.path.join(outdir, file), 'r') as f:
                    table = extract_table(f)
                txtfile = file
            elif re.search(bed, file):
                bedfile = file
            else:
                continue

    r = re.compile('resize.png$')
    thumbnails = filter(r.search, images)

    buildfile = command['build'].split('=')[-1]
    build = os.path.basename(buildfile).split('_')[0]
    
    return render_template('results.html',
                           id=os.path.basename(outdir),
                           name=os.path.basename(command['filepath']),
                           build=build,
                           txtfile=txtfile,
                           bedfile=bedfile,
                           images=reversed(thumbnails),
                           table=table,
                           tablerange=range(0,len(table['Sample']) + 1))

@app.route('/data/<id>/<file>')
def data(id, file):
    """ Return requested file to results page """
    outdir = os.path.join(app.config['UPLOAD_FOLDER'],id)
    return send_from_directory(outdir,file)
    
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
    
def resize(folder, fileName, width, height):
    name, ext = os.path.splitext(fileName)
    filePath = os.path.join(folder, fileName)
    im1 = Image.open(filePath)
    im = im1.copy()
    size = width, height
    im.thumbnail(size, Image.ANTIALIAS)
    im.save(os.path.join(folder, name) + '.resize.png')

def bulkResize(imageFolder, width, height):
    imgExts = ["png", "bmp", "jpg"]
    for path, dirs, files in os.walk(imageFolder):
        for fileName in files:
            ext = fileName[-3:].lower()
            if ext not in imgExts:
                continue

            resize(path, fileName, width, height)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
