import subprocess
from celery import Celery

celery = Celery()
celery.config_from_object('celeryconfig')

@celery.task
def run(command):
    p = subprocess.Popen(command,
                         shell=False,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return (stdout, stderr)

