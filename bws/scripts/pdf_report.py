# Utility for generating PDF reports
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

TMPDIR = "/tmp"
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


class HttpServer:
    cwd = os.getcwd()
    server_thread = None
    base = join(TMPDIR, 'base.html')

    def run_server(self):
        PORT = 8081
        os.chdir(TMPDIR)
        try:
            with socketserver.TCPServer(("", PORT), http.server.SimpleHTTPRequestHandler) as httpd:
                print("serving at port", PORT)
                httpd.serve_forever()
        except OSError as e:
            print(e)
            sys.exit(1)

    def start_www(self, url):
        self.base_html_setup(url)
        copyfile("/home/tim/workspace/boadicea/boadicea/local_apps/fh/static/js/pedigreejs.v2.1.0-rc2.min.js",
                 join(TMPDIR, 'pedigreejs.v2.1.0-rc2.min.js'))
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


def rm_file(fname):
    ''' Remove file if exists. '''
    if os.path.isfile(fname):
        os.remove(fname)
        print(f"REMOVED {fname}")


def wait_for_pdf_download(fname="canrisk_report.pdf", max_time=10, time_interval=0.2, rename=None):
    ''' Wait for file to exist on download. '''
    home = expanduser("~")
    fname = join(home, 'Downloads', fname)
    for _i in range(int(max_time/time_interval)):
        if os.path.isfile(fname) and os.stat(fname).st_size != 0:
            # print(str(_i*time_interval)+"s wait for file save")
            if rename:
                os.rename(fname, rename)
                rm_file(pedigree)
            return
        time.sleep(time_interval)
    print(fname+" not saved")


def create_pdf(url, token, ows_result, bws_result, bwa, cwd):
    copyfile(bwa, pedigree)

    # send to server to get web-page and write output.html
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

    # open browser to generate results page and PDF
    rm_file(join(expanduser("~"), 'Downloads', "canrisk_report.pdf"))
    webbrowser.open('http://0.0.0.0:8081/base.html')
    _dir, fname = os.path.split(bwa)
    wait_for_pdf_download(rename=join(cwd, f"{fname}.pdf"))
