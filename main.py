import urllib.request
import os
from datetime import date
from pathlib import Path

def save_file(text, filename):
    file = open(filename, 'w')
    file.write(text)
    file.close()


def download_non_redundant_set():
    resolution = '3.0A'
    base_url = 'http://rna.bgsu.edu/rna3dhub/nrlist/download'
    release = 'current'
    url = '/'.join([base_url, release, resolution])
    response = urllib.request.urlopen(url)
    try:
        info = response.info()['Content-Disposition']
        filename = info.split('filename=',1)[1]
    except:
        today = date.today()
        timestamp = today.strftime("%d-%m-%Y")
        filename = 'output_d' + timestamp + '.csv'
    data = response.read()
    text = data.decode('utf-8')

    path_to_location = Path('./RNA_SETS')
    save_file(text, path_to_location / filename)


if __name__ == "__main__":
    if not os.path.exists(Path('./RNA_SETS')):
        os.makedirs(Path('RNA_SETS'))
    download_non_redundant_set()
