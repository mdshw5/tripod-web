import subprocess
import time
from celery import Celery
from collections import OrderedDict

celery = Celery()
celery.config_from_object('celeryconfig')

@celery.task
def run(command):
    p = subprocess.Popen(command.values(),
                         shell=False,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    startTime = time.time()
    stdout, stderr = p.communicate()
    endTime = time.time()
    return (command, p.returncode, stdout, stderr, round(endTime - startTime))
