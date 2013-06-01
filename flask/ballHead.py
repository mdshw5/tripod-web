import os
import subprocess
import random
import string
import time
import ConfigParser
from flask import Flask, flash, request, redirect, url_for, render_template, session, escape, Response
from werkzeug import secure_filename

config = ConfigParser.RawConfigParser()
config.read('default.cfg')

perlPath = ' '.join([config.get('paths', 'perl'), os.path.join(os.getcwd(),'triPOD.pl')])

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
            flash(u"{0} was uploaded successfully!".format(filename))
            session['gender'] = request.form['gender']
            session['alpha'] = request.form['alpha']
            session['build'] = request.form['build']
            session['pod'] = request.form['pod']
            session['podhd'] = request.form['podhd']
            session['podmi1'] = request.form['podmi1']
            session['podcr'] = request.form['podcr']
            session['filename'] = os.path.join(app.config['UPLOAD_FOLDER'], salt, filename)
            session['output'] = os.path.join(app.config['UPLOAD_FOLDER'], salt)
            return redirect(url_for('progress', status='success'))

        elif not allowed_file(file.filename):
            return redirect(url_for('progress', status='type_not_allowed'))

        
    return render_template('upload.html')

@app.route('/error')
def upload_error():
    return render_template('error.html')

@app.route('/success')
def success():
    return render_template('success.html')

@app.route('/progress/<status>')
def progress(status):
    if status == 'type_not_allowed':
        flash(u"File type must be .txt .csv or .tsv", 'error')
        return redirect(url_for('upload'))

    elif status == 'success':
        command = [perlPath,
                   '--gender=' + session['gender'],
                   '--graph=png',
                   '--alpha=' + session['alpha'],
                   '--build=' + os.path.join(os.getcwd(),session['build']),
                   '--' + session['pod'], 
                   '--' + session['podhd'], 
                   '--' + session['podmi1'], 
                   '--' + session['podcr'], 
                   '--out=' + session['output'],
                   session['filename']
        ]
        tripod = subprocess.Popen(' '.join(command),
                                  shell=True, 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE
                                  )
        flash(tripod.stdout.readlines())
        return redirect(url_for('progress', status='results'))
    
    elif status == 'results':
        return render_template('progress.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
