import urllib.request
import os
import csv
from datetime import date
from pathlib import Path
import difflib

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


def pick_file(n):
    base_path = Path('./RNA_SETS')
    list_of_files = list(base_path.glob('*.csv'))
    list_of_files_sorted = sorted(list_of_files, key=os.path.getmtime)
    picked_file = list_of_files_sorted[-n]
    return picked_file

def parse_output_file(filename):
    list_of_structures = []
    path_to_file = Path(filename)
    with open(path_to_file, newline='') as csvfile:
        csv_file_content = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in csv_file_content:
            list_of_structures.append(row[1][1:5])
    return list_of_structures

def find_difference(list_of_structures_new, list_of_structures_old):
    f = open('files_to_update.txt', 'w')
    for line in difflib.unified_diff(list_of_structures_old, list_of_structures_new, fromfile='file1', tofile='file2', lineterm='', n=0):
        for prefix in ('---', '+++', '@@'):
            if line.startswith(prefix):
                break
        else:
            f.write(line +'\n')
    f.close()


if __name__ == "__main__":
    if not os.path.exists(Path('./RNA_SETS')):
        os.makedirs(Path('RNA_SETS'))
    download_non_redundant_set()
    base_path = Path('./RNA_SETS')
    if len(list(base_path.glob('*.csv'))) > 1:
        list_of_structures_new = parse_output_file(pick_file(1))
        list_of_structures_old = parse_output_file(pick_file(2))
        find_difference(list_of_structures_new, list_of_structures_old)
    else:
        save_file('\n'.join(parse_output_file(pick_file(1))), 'init_set.txt')

