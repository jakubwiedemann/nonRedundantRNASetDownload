import urllib.request


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
        filename = 'output.csv'
    data = response.read()
    text = data.decode('utf-8')

    save_file(text, filename)


if __name__ == "__main__":
    download_non_redundant_set()