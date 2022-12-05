#!/bin/bash
awk -F "\t" '{n=split($9,a," ");for(i=1;i<=n;++i){
if(a[i]=="gene_id"){++i;print substr(a[i],2,length(a[i])-3)}
if(a[i]=="transcript_id"){++i;print substr(a[i],2,length(a[i])-3)  > "/dev/stderr"}}}'
