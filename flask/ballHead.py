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
pid_mask = [int(n) for n in config.get('param', 'pid_mask').split(',')]
print pid_mask

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
            command = [perlPath,
                       '--gender=' + request.form['gender'],
                       '--graph=png',
                       '--alpha=' + request.form['alpha'],
                       '--build=' + os.path.join(os.getcwd(),request.form['build']),
                       '--' + request.form['pod'], 
                       '--' + request.form['podhd'], 
                       '--' + request.form['podmi1'], 
                       '--' + request.form['podcr'], 
                       '--out=' + request.form['output'],
                       request.form['filename']
            ]

            tripod = subprocess.Popen(' '.join(command),
                                      shell=True, 
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE
                                      )

            session['filename'] = os.path.join(app.config['UPLOAD_FOLDER'], salt, filename)
            session['output'] = os.path.join(app.config['UPLOAD_FOLDER'], salt)
            session['pid'] = tripod.pid
            return redirect(url_for('progress', pid=mask_pid(tripod.pid)))

        elif not allowed_file(file.filename):
            flash(u"File type must be .txt .csv or .tsv", 'error')
            return redirect(url_for('upload'))

    return render_template('upload.html')

@app.route('/progress/<pid>')
def progress(pid):
    return render_template('progress.html', pid=pid)

@app.route('/status/<pid>')
def status(pid):
    """ Wait until the PID does not exist, and report success."""
    print '/proc/' + str(unmask_pid(int(pid)))
    while os.path.exists('/proc/' + str(unmask_pid(int(pid)))):
        time.sleep(2)
        print os.path.exists('/proc/' + str(unmask_pid(int(pid))))
    return 'success'
    
@app.route('/results')
def results():
    """ Format the final results page and return template."""
    return render_template('results.html')

def mask_pid(pid):
    """ Mask the PID before passing it around """
    pid = pid * pid_mask[0]
    pid = pid + pid_mask[1]
    return pid

def unmask_pid(pid):
    """ Return actual PID from masked PID """
    pid = pid - pid_mask[1]
    pid = pid / pid_mask[0]
    return pid

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
