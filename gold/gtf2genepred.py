#!/usr/bin/python

import sys
import os
################################################################################

if len(sys.argv) >= 2:
    gtf_filename = sys.argv[1]

else:
    print("usage: python GTF2refFlat.py GTFfilename")
    print("or ./GTF2refFlat.py GTFfilename")
    sys.exit(1)

################################################################################
print "start: gtf file"
gtf = open(gtf_filename, 'r')
gtf_dict = {}
for line in gtf:
    if line[0] == "#":
        continue
    line_list = line.strip().split("\t")
    chr_name = line_list[0]
    feature = line_list[2]
    line_list[3]= int(line_list[3])
    line_list[4]= int(line_list[4])
    score =line_list[5]
    strand = line_list[6]
    frame = line_list[7]
    attribute_list = line_list[8].strip(";").split("; ")
    attribute_dict = {}
    for item in attribute_list:
        item_list = item.split(" ")
        attribute_dict[item_list[0]]=item_list[1].strip("\"")

    transcript_id = attribute_dict["transcript_id"]
    if attribute_dict.has_key("oId"):
        transcript_id = transcript_id + "|" + attribute_dict["oId"]

          
    if not gtf_dict.has_key(transcript_id):
        gtf_dict[transcript_id] = {}
    if not gtf_dict[transcript_id].has_key(feature):
        gtf_dict[transcript_id][feature] = []
    temp_list = line_list[0:8]
    temp_list.append(attribute_dict)
    gtf_dict[transcript_id][feature].append( temp_list )
          
gtf.close()
print "end: gtf file"

feature = "exon"
output_filename = gtf_filename + "_refFlat.txt"
output=open(output_filename,'w')
for transcript_id in gtf_dict:
    start_list=[]
    end_list=[]
    strand=set()
    gene_id = set()
    chr_name = set()
    classcode= set()
    for line_list in gtf_dict[transcript_id][feature]:
        start_list.append(line_list[3]-1)
        end_list.append(line_list[4])
        strand.add(line_list[6])
        gene = line_list[8]["gene_id"] 
        if line_list[8].has_key("gene_name"):
            gene = gene + "|" + line_list[8]["gene_name"]
        gene_id.add(gene)
        chr_name.add(line_list[0])
        classcode.add(line_list[8]["class_code"])
    start_list.sort()
    end_list.sort()
          
    if not gtf_dict[transcript_id].has_key("CDS"):
        CDS_start = str(start_list[0])
        CDS_end = str(end_list[-1])
    else:
        CDS_start_list=[]
        CDS_end_list=[]
        for line_list in gtf_dict[transcript_id]["CDS"]:
            CDS_start_list.append(line_list[3]-1)
            CDS_end_list.append(line_list[4])
        CDS_start_list.sort()
        CDS_end_list.sort()
        CDS_start = str(CDS_start_list[0])
        CDS_end = str(CDS_end_list[-1])

        
    if not len(strand)*len(gene_id)*len(chr_name)*len(classcode)==1:
        print "warnings: differenet gene id's or strands or chr_name for the transcript: " + transcript_id
    else:
        start_list.sort()
        end_list.sort()
        temp_start=""
        for start in start_list:
            temp_start=temp_start+str(start)+","
        temp_end=""
        for end in end_list:
            temp_end=temp_end+str(end)+","

        Nexon = str(len(start_list))
        gene_id = list(gene_id)[0]
        strand = list(strand)[0]
        chr_name = list(chr_name)[0]
        classcode = list(classcode)[0]
        output.write("\t".join([gene_id,transcript_id+"#"+classcode,chr_name,strand,str(start_list[0]),str(end_list[-1]),CDS_start,CDS_end, Nexon,temp_start,temp_end]) + "\n")
output.close()    
