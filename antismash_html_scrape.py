#!/miniconda/bin/python
#created by R. Murphy at the SSEG-UCPH
from bs4 import BeautifulSoup
import argparse
from pathlib import Path
import os
import re
import pandas as pd
from shutil import copy

parser = argparse.ArgumentParser(description = ("Script to scrape data from antismash html output into a csv format."))"

parser.add_argument("-p", "--path", help = "give path/to/directory containing antismash outputs", required = True)

args = parser.parse_args()

Path(args.path + '/html').mkdir(parents = True, exist_ok = True)
Path(args.path + '/csv').mkdir(parents = True, exist_ok = True)
Path('tmp').mkdir(parents = True, exist_ok = True)

def file_search(filename):
    seq_name = filename.parent.name
    copy(filename, args.path+'/html/'+seq_name+'.html')

for filename in Path(args.path).rglob('index.html'):
    file_search(filename)

# function to scrape the html files for the BGC cluster closest know cluster, cluster type and similarity to closest known.
def scraper(filename):
    df = pd.DataFrame(columns = ['classification', 'type', 'similar', 'similarity']) # generating the blank dataframe
    soup = BeautifulSoup(open(filename), 'html.parser') # reading in the html file

    with open("tmp/"+filename.stem+".txt", "w") as tmpfile: # opening a tmp file to save the converted html to text output
        # kill all script and style elements
        for script in soup(["script", "style"]):
            script.extract()    # rip it out
        # get text
        text = soup.get_text()
        # break into lines and remove leading and trailing space on each
        lines = (line.strip() for line in text.splitlines())
        # break multi-headlines into a line each
        chunks = (phrase.strip() for line in lines for phrase in line.split("  "))
        # drop blank lines
        tmpfile.write('\n'.join(chunk for chunk in chunks if chunk))
    marker = []
    with open("tmp/"+filename.stem+".txt", "r") as infile: # opening the tmp file for scraping the data
        readfile = infile.readlines() #reads the infile lne by line and returns a list containing the lines
        for i, line in enumerate(readfile[1:], 1): # looping over all the lines in the file from position 1 (so skipping 0) to avoid circular feedback
            if 'Overview' in line:
                start = i
            if 'Identified secondary metabolite regions using strictness' in line:
                    seq_name = readfile[i + 1].split('_')[0]
                    end = i
                    marker = list(map(lambda s: s.strip('\n'), readfile[start + 1:end])) # stripping the '\n' off every element in the list. map executes a function for each element in a sequence
        for location in marker:
            counter = 0
            for i, line in enumerate(readfile[1:], 1):
                if 'Region&nbsp' + location in line and counter < 1 and '%' in readfile[i + 6]: # looking at all line containing the region marker with a % 6 lines below them
                    counter += 1 # since the text file is in duplicate for some reason to ensure our csv is not we have added a counter here to ensure only one per marker is read
                    df = df.append({'classification': 'predicted', 'type': readfile[i + 1].rstrip().strip(), 'similar': readfile[i + 4].rstrip().strip(), 'similarity': readfile[i + 6].rstrip().strip()}, ignore_index = True)
                elif 'Region&nbsp' + location in line and counter < 1 and '%' not in readfile[i + 6]:
                    df = df.append({'classification': 'novel', 'type': readfile[i + 1].rstrip().strip(), 'similar': 'N/A', 'similarity': 'N/A'}, ignore_index = True)
                    counter += 1
        df.to_csv(args.path+"/csv/"+filename.stem+".csv", sep = ',', index = False) # saving the dataframe to a csv file
        infile.close()
        os.remove("tmp/"+filename.stem+".txt")

def source(script, update=1):
    pipe = Popen(". %s; env" % script, stdout=PIPE, shell=True)
    data = pipe.communicate()[0]

    env = dict((line.split("=", 1) for line in data.splitlines()))
    if update:
        environ.update(env)

    return env

for filename in Path(args.path + '/html/').glob("*.html"):
    scraper(filename)

Path('tmp').rmdir()
