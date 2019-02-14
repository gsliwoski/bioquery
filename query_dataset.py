import sys
import re
import os
import gzip
import pandas as pd

try:
    term = sys.argv[1].strip('"')
except:
    sys.exit("python query_dataset.py query_term dataset(uniprot,medgen, ...addingmore)")

try:
    dataset = sys.argv[2]
except:
    dataset = 'uniprot'

if dataset=='uniprot':
    datafile = "/dors/capra_lab/data/uniprot/2019-02-13/uniprot_sprot_human.dat"
elif dataset=='medgen':
    datafile = "/dors/capra_lab/data/medgen/2019-02-13/MGCONSO.csv.gz"
elif dataset=='clinvar':
    datafile = "/dors/capra_lab/data/clinvar/2019-02-13/GRCh38/clinvar_20190211.vcf.gz"

print "querying {} for genes with {}".format(datafile,term)    

outfilename = "{}_{}.results".format(re.sub(r"\s+","_",term),dataset)
if not os.path.isfile(outfilename) and dataset!="medgen":
    with open(outfilename,"w") as outfile:
        if dataset=="uniprot":
            outfile.write("ID\tAC\tGN\tmatch\tln\n")
#        elif dataset=="medgen":
#            outfile.write("CUI\tAUI\tSAUI\tSCUI\tSDUI\tCODE\tSTR\n")

def log_records(cid,cac,cgn,r):
    cid = re.sub(r"\s+","_",cid)
    cac = re.sub(r"\s+","_",cac)
    cgn = re.sub(r"\s+","_",cgn)
    print "Found term for {}".format(cid)
    with open(outfilename,"a") as outfile:
        for record in r:
            line, ln = record
            outfile.write("\t".join([cid,cac,cgn,line,str(ln)])+"\n")           
            
ln = 0
if dataset=='uniprot':
    records = list()
    current_id = None
    current_ac = None
    current_gn = None
    with open(datafile) as infile:
        for line in infile:
            ln += 1
            if line[:2] == "ID":
                if len(records)>0:
                    log_records(current_id, current_ac, current_gn, records)
                records = list()
                current_id = line.split()[1]
                current_ac = ""
                current_gn = ""
            elif line[:2] == "AC":
                current_ac += "".join(line.split()[1:])
            elif line[:2] == "GN":
                if "Name" not in line and "Synon" not in line: continue
                gn = line[2:].strip().split(";")
                gn = [x.strip().split("=")[1] for x in gn if "=" in x]
            elif len(re.findall(r"\b{}\b".format(term.lower()),line.lower()))>0:
                records.append([re.sub(r"\s+","_",line),ln])

elif dataset=='medgen':
    df = pd.read_csv(datafile,sep=",",quotechar='"')
    matches = df[df.STR.str.contains(r"\b{}\b".format(term),case=False)]
    print "{} matches found".format(matches.shape[0])
    for x in matches.columns:
        matches.loc[:,x].replace(r"\s+","_",regex=True,inplace=True)
    matches[["CUI","AUI","SAUI","SCUI","CODE","STR"]].to_csv(outfilename,sep="\t",header=True,index=False)

elif dataset=='clinvar':
    begin = False
    clinrec = list()
    with gzip.open(datafile) as infile:
        for line in infile:
            if line[:6] == "#CHROM": begin = True
            if not begin: continue
            if "{}".format(term.lower()) in line.lower() and "CLINSIG=Benign" not in line and "CLINSIG=Likely_benign" not in line:
                clinrec.append(line)
    print "found {} matches".format(len(clinrec))
    if len(clinrec)>0:
        with open(outfilename,'w') as outfile:
            outfile.write("".join(clinrec))                
