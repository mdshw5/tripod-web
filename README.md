#Overview

This document contains information relevant to installing and running
the triPOD web interface.

#System requirements

- Linux or MacOS
- Perl > 5.8.8
-- Modules: Algorithm::Cluster; Time::HiRes; Tree::Interval
- Python > 2.6
-- Modules: cgi, os, subprocess, cgitb, time, hashlib, re
- R > 2.14
-- Packages: shape, TTR
- Rscript
- Apache
- PHP5

#Installation

Install cgi-bin and html directories to webserver root, then modify
file paths in the header of each script.