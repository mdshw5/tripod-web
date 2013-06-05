# Overview

Ballhead is a Flask application	for running triPOD analysis. It	utilizes the Celery distributed	task queue for process management.

# Requirements

- Linux	
- Perl 5.8.8
    - Algorithm::Cluster; Time::HiRes; Tree::Interval
- Python 2.7
    - Flask 0.9
    - Celery 3.0
- R 2.14
    - shape, TTR, Rscript
- Apache2
- RabbitMQ

# Installation

## Celery

Install	Celery:

    sudo pip install Celery

After installing Celery, copy the init.d script	to your init.d directory:	`cp celery/celeryd /etc/init.d/celeryd`.
Also `mv celery/defaults.celeryd /etc/defaults/celeryd`. Make sure to edit default values to match installation	directories.

## RabbitMQ

Install	RabbitMQ:

    sudo apt-get install rabbitmq-server

## Configure Apache

