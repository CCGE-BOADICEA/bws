import json
import http.server
import socketserver
import requests
import sys
import os
import webbrowser
import threading
import time
from os.path import expanduser, join
from shutil import copyfile

DIRECTORY = "/tmp"


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


def rm_report(fname="canrisk_report.pdf"):
    ''' Wait for file to exist on download. '''
    home = expanduser("~")
    fname = join(home, 'Downloads', fname)
    if os.path.isfile(fname):
        os.remove(fname)


def wait_for_statfile(fname="canrisk_report.pdf", max_time=10, time_interval=0.25, rename=None):
    ''' Wait for file to exist on download. '''
    home = expanduser("~")
    fname = join(home, 'Downloads', fname)
    for _i in range(int(max_time/time_interval)):
        if os.path.isfile(fname) and os.stat(fname).st_size != 0:
            print(str(_i*time_interval)+"s wait for file save")
            if rename:
                os.rename(fname, rename)
            return
        time.sleep(time_interval)
    print(fname+" not saved")


def run_server():
    PORT = 8081
    os.chdir(DIRECTORY)
    with socketserver.TCPServer(("", PORT), http.server.SimpleHTTPRequestHandler) as httpd:
        print("serving at port", PORT)
        httpd.serve_forever()


def create_pdf(url, token, ows_result, bws_result, bwa, pdf_name):
    cwd = os.getcwd()
    base = join(DIRECTORY, 'base.html')
    pedigree = join(DIRECTORY, 'pedigree.txt')
    copyfile(join(os.path.dirname(os.path.realpath(__file__)), 'base.html'), base)
    copyfile(bwa, pedigree)
    data = {"ows_result": json.dumps(ows_result.json(), separators=(',', ':')),
            "bws_result": json.dumps(bws_result.json(), separators=(',', ':'))}
    r = requests.post(url+'combine/', data=data, headers={'Authorization': "Token "+token})
    if r.status_code == 200:
        f = open('/tmp/output.html', 'wb')
        f.write(r.content)
        f.close()
    else:
        sys.stderr.write("Web-services error status: "+str(r.status_code))
        sys.stderr.write(r.text)
        exit(1)

    t = Thread_With_Trace(target=run_server)
    t.start()
    rm_report()
    webbrowser.open('http://0.0.0.0:8081/base.html', new=2)
    wait_for_statfile(rename=join(cwd, pdf_name))
    os.remove(base)
    t.kill()
