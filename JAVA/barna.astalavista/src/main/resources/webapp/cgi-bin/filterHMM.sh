#!/bin/bash
## expects input stream with hmm file, 
## and $1 a list of pfam IDs (one per line)
## ignores version suffix (in query and db) 
## and compares case-insensitive
## writes to stdout
# HMMER3/b [3.0 | March 2010]
# NAME  1-cysPrx_C
# ACC   PF10417.4
# DESC  C-terminal domain of 1-Cys peroxiredoxin

awk -v f=$1 'BEGIN{while(getline< f > 0){split($1,a,".");x[toupper(a[1])]=$1}y["NAME"]="";y["ACC"]="";y["DESC"]=""}
($0~/^HMMER/){if(b==1){if(buf!=""){print buf}buf="";b=0}}
($1=="ACC"){split($2,a,".");id=toupper(a[1]);for(i in x){if(i==id){b=1;break}}}
{buf=buf""$0"\n"}
END{if(b==1){print buf}}'

