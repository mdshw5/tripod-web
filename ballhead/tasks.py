import subprocess
from celery import Celery
from collections import OrderedDict

celery = Celery()
celery.config_from_object('celeryconfig')

@celery.task
def run(command.values()):
    p = subprocess.Popen(command,
                         shell=False,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return (command, p.returncode, stdout, stderr)
