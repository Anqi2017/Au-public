#!/usr/bin/python

import sys
import math


### generate_transcript_list
############################
def generate_transcript_list(gpd_file, transcript_list):
    
    for line in gpd_file:
        
        if (line[0] == '#'):
            continue
        
        fields = line.split()
        num_exons = int(fields[8])

        start_pos_list = fields[9].split(',')
        end_pos_list = fields[10].split(',')
        
        exon_pos = [0] * num_exons
        for i in range(num_exons):
           exon_pos[i] = [start_pos_list[i], end_pos_list[i]]
        
        transcript_list.append([fields[0], fields[1], fields[2], fields[3], exon_pos])
        

### generate_RPKM_dict
#######################
def generate_RPKM_dict(RPKM_file, RPKM_dict):
    
    for line in RPKM_file:
        fields = line.split()
        RPKM_dict[fields[0]] = fields[1]
        
    

### generate_gpd_format
#######################
def generate_gtf_format(gtf_file, transcript_list, RPKM_dict, tool):

    
    for line in transcript_list:
        exon_pos = line[4]
        # transcript line
        
        # chr name 
        gtf_file.write(line[2] + '\t' + "IDP" + '\t' + "transcript" + '\t')
        # start-end pos, score
        gtf_file.write(exon_pos[0][0] + '\t' + exon_pos[-1][1] + '\t' + '*' + '\t')
        # Direction
        gtf_file.write(line[3] + '\t' + '.' + '\t')
        
        if (RPKM_dict.has_key( line[1]) ):
            RPKM = RPKM_dict[line[1]]
        else:
            RPKM = '*'
        attribute_1 = 'gene_id "' + line[0] + '"; transcript_id "' + line[1] + '"; '
        attribute_2 = ('FPKM "' + RPKM + '"; frac "' + '*' + '"; conf_lo "' + '*' + '"; ' +
                       'conf_hi "' + '*' + '"; cov "' + '*' + '";\n')
        
        gtf_file.write(attribute_1)
        gtf_file.write(attribute_2)
        
        num_exons = len(exon_pos)
        for i in range(num_exons):
            # chr name 
            gtf_file.write(line[2] + '\t' + tool + '\t' + "exon" + '\t')
            # start-end pos, score
            gtf_file.write(exon_pos[i][0] + '\t' + exon_pos[i][1] + '\t' + '*' + '\t')
            # Direction
            gtf_file.write(line[3] + '\t' + '.' + '\t')
            gtf_file.write(attribute_1)
            gtf_file.write('exon_number "' +  str(i) + '"; ')
            gtf_file.write(attribute_2)



### Main
########
def main():
    gpd_file = open(sys.argv[1], 'r')
    RPKM_file = open(sys.argv[2], 'r')
    gtf_file = open(sys.argv[3], 'w')
    
    if (len(sys.argv) > 4):
        tool = sys.argv[4]
    else:
        tool = "IDP"
    transcript_list = []
    RPKM_dict = dict()
    generate_transcript_list(gpd_file, transcript_list)
    generate_RPKM_dict(RPKM_file, RPKM_dict)
    generate_gtf_format(gtf_file, transcript_list, RPKM_dict, tool)
    
    gpd_file.close()
    gtf_file.close()


if __name__ == '__main__':
    main()
