"""
Utility for generating PDF reports

© 2023 University of Cambridge
SPDX-FileCopyrightText: 2023 University of Cambridge
SPDX-License-Identifier: GPL-3.0-or-later
"""

import json
import http.server
import socketserver
import requests
import sys
import os
import webbrowser
import threading
import time
from os.path import join
from shutil import copyfile
import tempfile
from pathlib import Path

TMPDIR = tempfile.gettempdir()
pedigree = join(TMPDIR, 'pedigree.txt')


class Thread_With_Trace(threading.Thread):

    def __init__(self, *args, **keywords):
        threading.Thread.__init__(self, *args, **keywords)
        self.killed = False

    def start(self):
        self.__run_backup = self.run
        self.run = self.__run
        threading.Thread.start(self)

    def __run(self):
        sys.settrace(self.globaltrace)
        self.__run_backup()
        self.run = self.__run_backup

    def globaltrace(self, frame, event, arg):
        if event == 'call':
            return self.localtrace
        else:
            return None

    def localtrace(self, frame, event, arg):
        if self.killed:
            if event == 'line':
                raise SystemExit()
        return self.localtrace

    def kill(self):
        self.killed = True


class QuietServer(http.server.SimpleHTTPRequestHandler):

    def log_message(self, _fmt, *args):
        pass


class HttpServer:
    cwd = os.getcwd()
    server_thread = None
    base = join(TMPDIR, 'base.html')

    def run_server(self):
        PORT = 8081
        os.chdir(TMPDIR)
        try:
            with socketserver.TCPServer(("", PORT), QuietServer) as httpd:
                print("serving at port", PORT)
                httpd.serve_forever()
        except OSError as e:
            print(e)
            sys.exit(1)

    def start_www(self, url):
        self.base_html_setup(url)
        HttpServer.server_thread = Thread_With_Trace(target=HttpServer().run_server)
        self.server_thread.start()

    def stop_www(self):
        self.server_thread.kill()
        os.remove(HttpServer.base)

    def base_html_setup(self, url):
        url = url[:-1] if url.endswith("/") else url
        f = open(join(os.path.dirname(os.path.realpath(__file__)), 'base.html'), "r")
        new_content = ""
        for line in f:
            new_content += line.replace("{URL}", url)
        f.close()

        bf = open(HttpServer.base, "w")
        bf.write(new_content)
        bf.close()

        # add CanRisk icon
        imgDir = join(TMPDIR, 'static', 'img')
        Path(imgDir).mkdir(parents=True, exist_ok=True)
        img_data = requests.get(url+"/static/img/CanRisk250x83.png").content
        img = open(join(imgDir, 'CanRisk250x83.png'), 'wb')
        img.write(img_data)
        img.close()
        
        # add favicon
        fav = requests.get(url+"/static/favicon.ico").content
        ico = open(join(TMPDIR, 'favicon.ico'), 'wb')
        ico.write(fav)
        ico.close()
        

def rm_file(fname):
    ''' Remove file if exists. '''
    if os.path.isfile(fname):
        os.remove(fname)
        print(f"REMOVED {fname}")


def wait_for_pdf_download(fname="canrisk_report.pdf", max_time=10, time_interval=0.25, rename=None):
    ''' Wait for file to exist on download. '''
    fname = join(Path.home(), "Downloads", fname)
    for _i in range(int(max_time/time_interval)):
        if os.path.isfile(fname) and os.stat(fname).st_size != 0:
            # print(str(_i*time_interval)+"s wait for file save")
            if rename:
                print("PDF FILE: "+rename)
                os.rename(fname, rename)
                rm_file(pedigree)
            return
        time.sleep(time_interval)
    print("************************* "+(rename if rename is not None else fname)+" NOT SAVED! *************************")


def create_pdf(url, token, ows_result, bws_result, bwa, cwd):
    rm_file(pedigree)
    copyfile(bwa, pedigree)

    # send to server to get web-page and write output.html
    data = {"ows_result": json.dumps(ows_result.json(), separators=(',', ':')),
            "bws_result": json.dumps(bws_result.json(), separators=(',', ':'))}
    r = requests.post(url+'combine/', data=data, headers={'Authorization': "Token "+token})
    if r.status_code == 200:
        f = open(join(TMPDIR, 'output.html'), 'wb')
        f.write(r.content)
        f.close()
    else:
        sys.stderr.write("Web-services error status: "+str(r.status_code))
        sys.stderr.write(r.text)
        exit(1)

    # open browser to generate results page and PDF
    rm_file(join(Path.home(), 'Downloads', "canrisk_report.pdf"))
    webbrowser.open('http://localhost:8081/base.html')
    _dir, fname = os.path.split(bwa)
    wait_for_pdf_download(rename=join(Path.home(), 'Desktop', f"{fname}.pdf"))
