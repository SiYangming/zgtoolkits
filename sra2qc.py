#! /bin/python3

# you need preinstall sratoolkit fastqc and multiqc and make sure these were included in the PATH
# This script can only run in the unix-like system

from urllib import request
import re
import os
import subprocess
import glob

# use geo accession to find the ftp address of sra files
# geo query website :https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=
def get_sra_address(geo=None):
    common = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
    query = common + geo
    rq = request.Request(query)
    response = request.urlopen(rq)
    page_content = response.read().decode('utf-8')
    pattern = re.compile('.*<a href="(.*SRP.*?)">\(ftp\)</a>')
    sra_address = re.findall(string=page_content, pattern=pattern)[0]
    return sra_address

# download sra files

def sra2fq(sra_address, geo):
    os.makedirs(geo)
    os.chdir(geo)
    cmd = 'wget -A .sra -nd -r 2 -4 ' + sra_address
    subprocess.run(cmd, shell=True)
    for sra_file in glob.glob('*.sra'):
        cmd = 'fastq-dump –split-3 ' + sra_file
        subprocess.run(cmd, shell=True)

def multi_fastqc():
    for fq_file in glob.glob('*.fastq'):
        cmd = 'fastqc -t 4 ' + fq_file
        subprocess.run(cmd, shell=True)
    subprocess.run('multiqc *fastqc.zip --pdf', shell=True)

def main():
    import sys
    geo = sys.argv[1]
    sra_address = get_sra_address(geo)
    sra2fq(sra_address, geo)
    multi_fastqc()

if __name__ == '__main__':
    main()
