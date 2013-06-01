#!/usr/bin/env python26
"""This script is meant to wait for handler.py to finish,
and then format the results from triPOD"""

import os
import sys
import time
import re
import string

pid = sys.argv[1]
resultsDir = sys.argv[2]
uploadedFile = sys.argv[3]
md5Hash = sys.argv[4]
fileName = sys.argv[5]
resultsURL = sys.argv[6]
message = sys.argv[7]
build = sys.argv[8]
logFile = sys.argv[9]
config = 'config.py'
execfile(config)
log = open(resultsDir + '/log', 'w')
origStderr = sys.stderr
sys.stderr = log

log.write('initialized\n')
log.write('%s\n' % resultsDir)
log.write('%s\n' % uploadedFile)

time.sleep(5)

timer = 0
progress = open(logFile, 'r')
while os.path.exists("/proc/" + pid):
    time.sleep(1)
    timer += 1
    if timer is 10:
        timer = 0
        progress.seek(0)
        progressText = progress.read().splitlines()
        progressText = [i.replace(apacheRoot, '') for i in progressText]
        progressText = "\n".join(map(str, progressText))
        results = open(resultsDir + '/results.php', 'w')
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
<title>Processing sample %s...</title>
<link rel="stylesheet" type="text/css" href="../../triPOD.css">
<?php include '../../header.php'; ?>
</head>
<body>
<pre>
%s<p>
%s
</pre>
</body></html> ''' % (fileName, message, progressText))
        results.close()
progress.close()
resultsDirListing = os.listdir(resultsDir)

png = re.compile("png$")
txt = re.compile("txt$")
bed = re.compile("bed$")
images = []
for item in resultsDirListing:
    if re.search(png, item):
        images.append(item)
    elif re.search(txt, item):
        txtResults = item
    elif re.search(bed, item):
        bedResults = item
    else:
        continue

table = open(resultsDir + '/' + txtResults, 'r+')
txtTable = table.read().splitlines()
del txtTable[0:6] # remove before table
txtTable[-1] = txtTable[-1].replace(cgiRoot,'')
txtTable[-1] = txtTable[-1].replace(htmlRoot,'')
table.seek(0)
for line in txtTable:
    table.write("%s\n" % line) # write sanitized table
table.truncate()
table.close()
for _ in range(4):
    del txtTable[-1] #remove after table
splitTable = [line.split() for line in txtTable]
headerTable = splitTable.pop(0)
chroms = list([int(i[1]) for i in splitTable])
uniqChroms = list(set(chroms))
freqs = [(i, chroms.count(i)) for i in uniqChroms]
results = open(resultsDir + '/results.php', 'w')
results.write('''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html;charset=utf-8">
<title>triPOD results for sample %s</title>
<link rel="stylesheet" type="text/css" href="../../triPOD.css">
<style type="text/css">
table
{border-collapse:collapse;}
table, td, th
{border:1px solid black;}
</style>
</head>
<body>
<?php include '../../header.php'; ?>
''' % fileName)

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
ucscURL = resultsURL + '/' + bedResults
results.write('''<p><center><table id="top">
<tr><td rowspan="3">Regions detected in sample %s</td>
<td><a href='%s'>Download results table as text file.</a></td></tr>
<tr><td><a href='%s'>Download results table as bed file.</a></td></tr>
<tr><td><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s&org=human&hgt.customText=%s" target="_blank">View results in UCSC Genome Browser</a></td></tr></table>''' % (fileName, txtResults, bedResults, build, ucscPos, ucscURL))

results.write('<table><tr><th>Image</th>')
for header in headerTable[1:9]:
    results.write('''<th>%s</th>''' % header)

for c in uniqChroms:
    indices = [i for i, x in enumerate(chroms) if x is c]
    results.write('''
<tr><td rowspan="%s"><a href="#%s">View Chromosome %s</a>''' % (str(len(indices)), c, c))
    count = 0
    for n in indices:
        if count > 0:
            results.write('<tr>')
        for cell in splitTable[n][1:9]:
            results.write('<td>%s</td>' % cell)
        count += 1
    results.write('</tr>')

# Write images in to the table
for c in uniqChroms:
    indices = [i for i, x in enumerate(chroms) if x is c]
    r = '\D' + str(c) + '.png$'
    childImage = str([x for x in images if re.search(r, x)]).strip('[]')
    results.write('''<tr><td colspan=9><img width="640" id="%s" src=%s><img></td></tr>
<tr><td colspan=9><a href="#top">Back to top</a></td></tr>''' % (c, childImage))
results.write('</table>')

results.write('''</center></body></html>''' )
results.close()
sys.stderr = origStderr
log.close()
if sampleData is None:
    os.remove(uploadedFile)
