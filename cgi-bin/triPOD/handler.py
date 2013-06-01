#!/usr/bin/env python26
""" This script is meant to receive form input form index.php, process file
uploads, create analysis and results directories, and launch triPOD.pl 
and generate_results.py
"""
import cgi
import os
import sys
import subprocess
import cgitb
import time
import hashlib
cgitb.enable()

# Import config file
configFile = 'config.py'
execfile(configFile)

# Generate our hashes so that we have unique files and directories
salt = time.ctime()
md5Hash = hashlib.md5(salt).hexdigest()

# Get form input from index.php
form = cgi.FieldStorage()
fileItem = form['file']
gender = form.getfirst('gender', 'NA')
alpha = form.getfirst('alpha', '0.1')
build = form.getfirst('build', 'hg18_centromeres.txt')
buildFile = os.getcwd() + '/' + build
pod = '--' + form.getfirst('pod', 'pod')
podHD = '--' + form.getfirst('podhd', 'nohd')
podMI1 = '--' + form.getfirst('podmi1', 'nomi1')
podCR = '--' + form.getfirst('podcr', 'nopodcr')
build = build.replace('_centromeres.txt', '')
sampleData = form.getfirst('sampledata', '')

# Make sure we use at least one detection method
if (pod is 'nopod' and podHD is 'nohd' 
    and podMI1 is 'nomi1' and podCR is 'nopodcr'):
    print 'Content-Type: text/html\n\n'
    print 'Please specify at least one detection method'
    sys.exit()

def fbuffer(f, chunk_size):
    '''Generator to buffer file chunks'''  
    while True:
        chunk = f.read(chunk_size)      
        if not chunk: break
        yield chunk

uploadDir = 'upload/' + md5Hash
os.mkdir(uploadDir)
fileName = os.path.basename(fileItem.filename) # base filename
uploadedFile = uploadDir + '/' + fileName # path to uploaded file

# Run generator and move file, report success
if fileItem.filename:
    file = open(uploadedFile, 'w', 10000)
    for chunk in fbuffer(fileItem.file, 10000):
        file.write(chunk)
    file.close()
    message = '''
The file %s was uploaded successfully. \n
This page will automatically refresh every 10 seconds until results are ready. \n 
Do not close this window.
Processing...''' % fileName
elif sampleData:
    uploadedFile = os.getcwd() + '/' + sampleData
    message = '''
Using sample data.
This page will automatically refresh every 10 seconds until results are ready. \n 
Do not close this window.
Processing...'''
else:
    message = 'No file was uploaded'

# Create output directory
resultsDir = apacheRoot + md5Hash
os.mkdir(resultsDir)
resultsURL = authURL + 'triPOD/results/' +  md5Hash

# Execute triPOD.pl
logFile = resultsDir + '/progress'
command = '%s triPOD.pl --cores %s --gender %s --graph=png --alpha %s --build %s %s %s %s %s --verbose --out %s %s > %s' % (perlPath, cores, gender, alpha, buildFile, pod, podHD, podMI1, podCR, resultsDir, uploadedFile, logFile)

tripod = subprocess.Popen(command,
                 shell=True,
                 stdout=subprocess.PIPE,
                 stderr=subprocess.PIPE)


results = open(resultsDir + '/results.php', 'w')
# Write the intermediate results page
results.write('''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html;charset=utf-8">
<script>
self_refresh=window.setTimeout(function()
{window.location.href=window.location.href},10000);
</script>
<noscript>
<meta http-equiv="refresh" content="10">
</noscript>
<title>triPOD results for sample %s</title>
<link rel="stylesheet" type="text/css" href="../../triPOD.css">
<?php include '../../header.php'; ?>
</head>
<body>
<pre>
%s
</pre>
</body></html> ''' % (fileName, message))
results.close()

# Call script that waits until triPOD is finished, formats the results page
command = ([pythonPath, 'generate_results.py', 
            str(tripod.pid), resultsDir,
             uploadedFile, str(md5Hash), fileName,
             resultsURL, message, build, logFile, sampleData])
generateResults = subprocess.Popen(command,
                 stdout=subprocess.PIPE,
                 stderr=subprocess.PIPE)

# Redirect user to the intermediate results page
redirectURL = (apacheURL + 'triPOD/results/' + md5Hash + '/results.php')
print('''Location:%s''' % redirectURL)
print '' # to end the CGI response headers.
