#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3:
    refFlat_filename = sys.argv[1]
    source = sys.argv[2]
else:
    print("usage: python genephed2gtf.py refFlat.txt source]")
    print("or ./genephed2gtf.py refFlat.txt source")
    sys.exit(1)

def convertframe(CDS_len_list):
    CDS_frame_list = []
    ref_CDS_frame=0
    CDS_frame_list.append(ref_CDS_frame)
    del CDS_len_list[-1]
    for CDS_len in CDS_len_list:
        ref_CDS_frame = (ref_CDS_frame + ( 2 - (CDS_len%3) ) )%3
        CDS_frame_list.append(ref_CDS_frame)
    return CDS_frame_list

    
################################################################################
    
ref=open(refFlat_filename,'r')
output = open(refFlat_filename+".gtf",'w')
print "loading refFlat file: ", refFlat_filename

for refline in ref:
    refline_list=refline.split()
    Exon_start=int(refline_list[4])
    Exon_end=int(refline_list[5])
    CDS_start=int(refline_list[6])
    CDS_end=int(refline_list[7])
    Exon_start_list=refline_list[9].strip(",").split(',')
    Exon_end_list=refline_list[10].strip(",").split(',')
    strand=refline_list[3]
    i=0
    CDS_indicator=0
    CDS_len_list = []
    left_str_list = []
    frame_list = []
    right_str = "\t" + "gene_id \""+refline_list[0]+"\"; transcript_id \""+refline_list[1]+"\";" + " \n" 
    k=0
    for start in Exon_start_list:
        estart= int(start)
        eend = int(Exon_end_list[i])
        if CDS_start+5>=CDS_end and not CDS_start==CDS_end:
            if strand == "+":
                feature1 =  "start_codon"
                feature2 =  "stop_codon"
            elif strand =="-":
                feature2 =  "start_codon"
                feature1 =  "stop_codon"
            left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature1 + "\t" + str(CDS_start + 1 ) + "\t" + str(CDS_start + 3) + "\t" + "0.000000" + "\t" + strand + "\t" )
            frame_list.append(".")
            left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature2 + "\t" + str(max( CDS_start+1, CDS_end-2)  ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t" )
            frame_list.append(".")   
            left_str_list.append( refline_list[2] + "\t" + source + "\t" + "exon" + "\t" + str(estart + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
            frame_list.append(".")     

        elif CDS_start >= estart and CDS_start <= eend and not (CDS_end-2 > estart and CDS_end-2 <= eend):
            CDS_indicator=1
            if strand == "+":
                feature = "start_codon"
                if CDS_start + 3 <= eend:

                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_start + 1 ) + "\t" + str(CDS_start + 3) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    
                else:
                    d = 3 - (eend-CDS_start)
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_start + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(int(Exon_start_list[i+1]) + 1 ) + "\t" + str(int(Exon_start_list[i+1]) + d ) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                left_str_list.append( refline_list[2] + "\t" + source + "\t" + "CDS" + "\t" + str(CDS_start + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
                frame_list.append("")
                CDS_len_list.append(eend-CDS_start - 1 )

            elif strand == "-":
                feature = "stop_codon"
                if CDS_start + 3 <= eend:
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_start + 1 ) + "\t" + str(CDS_start + 3) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    if CDS_start + 4 <= eend:
                        left_str_list.append( refline_list[2] + "\t" + source + "\t" + "CDS" + "\t" + str(CDS_start + 4 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t")
                        frame_list.append("")
                        CDS_len_list.append(eend-CDS_start - 4)
                else:
                    d = 3 - (eend-CDS_start)
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_start + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")                                                                                                       
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(int(Exon_start_list[i+1]) + 1 ) + "\t" + str(int(Exon_start_list[i+1]) + d ) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    k=d


            left_str_list.append(  refline_list[2] + "\t" + source + "\t" + "exon" + "\t" + str(estart + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
            frame_list.append(".")

#######################
                                                                                                       
        elif CDS_end-2 > estart and CDS_end-2 <= eend and not (CDS_start >= estart and CDS_start <= eend):
            CDS_indicator=0
            if strand == "-":
                left_str_list.append(  refline_list[2] + "\t" + source + "\t" + "CDS" + "\t" + str(estart + k + 1 ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t" )
                frame_list.append("")
                CDS_len_list.append(CDS_end-estart - 1)

                feature = "start_codon"
                if CDS_end <= eend:
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_end-2 ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                else:
                    d = 3-(3+eend - CDS_end)
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_end-2 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(int(Exon_start_list[i+1]) + 1 ) + "\t" + str(int(Exon_start_list[i+1]) + d ) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                                                                                                           


            elif strand == "+":


                feature = "stop_codon"
                if CDS_end  <= eend:
                    if estart + 1 <= CDS_end-3:
                        left_str_list.append( refline_list[2] + "\t" + source + "\t" + "CDS" + "\t" + str(estart + 1 ) + "\t" + str(CDS_end-3) + "\t" + "0.000000" + "\t" + strand + "\t" )
                        frame_list.append("")
                    CDS_len_list.append(CDS_end- estart - 4)
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_end-2 ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")

                else:
                    d = 3-(3+eend - CDS_end)
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_end-2 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(int(Exon_start_list[i+1]) + 1 ) + "\t" + str(int(Exon_start_list[i+1]) + d ) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    k=d


            left_str_list.append( refline_list[2] + "\t" + source + "\t" + "exon" + "\t" + str(estart + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
            frame_list.append(".")


#######################
                                                                                                       
        elif CDS_start >= estart and CDS_start <= eend and CDS_end-2 > estart and CDS_end-2 <= eend:
            CDS_indicator=0
            if not CDS_start+5>=CDS_end:
            
                if strand == "+":
                    feature = "start_codon"
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_start + 1 ) + "\t" + str(CDS_start + 3) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + "CDS" + "\t" + str(CDS_start + 1 ) + "\t" + str(CDS_end-3) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append("")
                    CDS_len_list.append(CDS_end-CDS_start - 4)
                    feature = "stop_codon"
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_end-2 ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    
                elif strand == "-":
                    feature = "stop_codon"
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_start + 1 ) + "\t" + str(CDS_start + 3) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + "CDS" + "\t" + str(CDS_start + 4 ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t")
                    frame_list.append("")                
                    CDS_len_list.append(CDS_end-CDS_start - 4)
                    feature = "start_codon"
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature+ "\t" + str(CDS_end-2 ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append(".")
            left_str_list.append( refline_list[2] + "\t" + source + "\t" + "exon" + "\t" + str(estart + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
            frame_list.append(".")

#######################
                                                                                                                       
        else:
            
            if CDS_end>=estart and CDS_end<=eend:
                d= 3 - (CDS_end-estart)
                if strand == "+":
                    feature = "stop_codon"
                    left_str_list[-2] = refline_list[2] + "\t" + source + "\t" + "CDS" + "\t" + str(int(Exon_start_list[i-1]) + k + 1 ) + "\t" + str(int(Exon_end_list[i-1]) - d) + "\t" + "0.000000" + "\t" + strand + "\t" 
                    CDS_indicator=0
                elif strand == "-":
                    feature = "start_codon"
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + "CDS"  + "\t" + str(estart + 1 ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append("")   
                    CDS_len_list.append(CDS_end-estart - 1)
                left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature + "\t" + str(int(Exon_end_list[i-1]) - d + 1 ) + "\t" + Exon_end_list[i-1] + "\t" + "0.000000" + "\t" + strand + "\t" )
                frame_list.append(".")
                left_str_list.append( refline_list[2] + "\t" + source + "\t" + feature  + "\t" + str(estart + 1 ) + "\t" + str(CDS_end) + "\t" + "0.000000" + "\t" + strand + "\t" )
                frame_list.append(".")
              
            else:
                if CDS_indicator==1:
                    left_str_list.append( refline_list[2] + "\t" + source + "\t" + "CDS" + "\t" + str(estart + k + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
                    frame_list.append("")
                    CDS_len_list.append(eend-estart - 1 )
                    k=0

            left_str_list.append( refline_list[2] + "\t" + source + "\t" + "exon" + "\t" + str(estart + 1 ) + "\t" + str(eend) + "\t" + "0.000000" + "\t" + strand + "\t" )
            frame_list.append(".")
        i+=1


    if not CDS_start+5>=CDS_end:
        if strand == "-":
            CDS_len_list.reverse()
        CDS_frame_list = convertframe(CDS_len_list)
        if strand == "-":
            CDS_frame_list.reverse()


    i=0
    j=0
    for left_str in left_str_list:
        frame=frame_list[i]
        if frame == ".":
            output.write(left_str+frame+right_str)
        else:
            output.write(left_str+frame+str(CDS_frame_list[j])+right_str)
            j+=1
        i+=1
 
ref.close()
output.close()
################################################################################
