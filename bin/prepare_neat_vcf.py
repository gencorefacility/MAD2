#!/usr/bin/env python3

import sys
import random
import numpy as np
import pysam
import os
from pyfaidx import Fasta

# give it a vcf file
in_vcf = sys.argv[1]

# give it a golden vcf file
golden_vcf = sys.argv[2]

# give it a ref fa
fa = sys.argv[3]

# give it an AF
# expects a float or the string 'random'
in_af = sys.argv[4]

# give it an output prefix
out_prefix = sys.argv[5]

# give it a seed, need this to create
# the same 'random' AFs between replicates
seed = int(sys.argv[6])

vcf = pysam.VariantFile(in_vcf, "r")
assembly = os.path.basename(fa)
contigs = Fasta(fa)

# parse golden vcf file and make index (list) 
# of snps using "<chrom>-<pos>" as the key
golden_snps = []
with open(golden_vcf, 'r') as vcf:
        for line in vcf:
                if line.startswith('#'):
                        continue
                columns = line.split('\t')
                chrom = columns[0]
                pos = columns[1]
                golden_snps.append('{}-{}'.format(chrom, pos))

def get_af(af, golden_snp_index=0):
    if af == 'random' and False: #old way
        vals = np.arange(0.01, 1.01, 0.01).tolist()
        p_under_50 = 3
        p_over_50 = 1
        l1 = [p_under_50] * 50
        l2 = [p_over_50] * 50   
        weights = l1 + l2
        random_af = random.choices(vals, weights=weights, k=1)[0]
        return random_af
    elif af == 'random':
        x = 0
        # in case np.random.poisson() returns 0, try again
        while x == 0:
            # use seed to get the same value for 
            # 'random' AFs between replicates
            np.random.seed(seed + golden_snp_index)
            x = np.random.poisson(10)
        return x / 100
    elif af == 'pcr_mut':
        x = 0
        # in case np.random.poisson() returns 0, try again
        while x == 0:
            np.random.seed()
            x = np.random.poisson(2)
        return x / 100
    else:
        af = float(af)
        if af > 1:
            print("AF cannot be greater than 1! AF = ", af)
            sys.exit()
        return float(af)

contigs_added = False   
with open(in_vcf, 'r') as vcf, open(out_prefix+"_golden.vcf", "w") as golden_vcf, open(out_prefix+'.vcf', 'w') as out_vcf:
    golden_snp_index = 0
    for line in vcf:
        if line.startswith('#'):
            if line.startswith("#CHROM\tPOS"):
                golden_vcf.write(line.strip() + "\tFORMAT\tsample\n")
                out_vcf.write(line.strip() + "\tFORMAT\tsample\n")
            else:
                golden_vcf.write(line)
                out_vcf.write(line)

            ## Add Contigs Here (needed by ReSeq)
            if not contigs_added:
                for key in contigs.keys():
                    golden_vcf.write("##contig=<ID={},length={},assembly={}\n".format(key,len(contigs[key]),assembly))
                    out_vcf.write("##contig=<ID={},length={},assembly={}\n".format(key,len(contigs[key]),assembly))
                contigs_added = True
        else:
            columns = line.split('\t')
            pos = columns[1]
            ref = columns[3]
            alt = columns[4]
            chrom = columns[0]
                        
            # if this snp is in the golden snps, set the desired AF
            # else, it is a 'PCR artifact', assign the AF as such
            if '{}-{}'.format(chrom, pos) in golden_snps:
                af = get_af(in_af, golden_snp_index)
                golden_snp_index+=1
                print("I: , AF: ", in_af, golden_snp_index, af)
            else:
                af = get_af('pcr_mut')
            
            new_gt_field = 'GT:AD:DP:GQ:PL\t1'
            # - 1 in the range function because we've
            # already added '1' in the line above
            for i in range(int(af*float(100)) - 1):
                new_gt_field+= '/1'
            # compare_afs.py will use ad[1]/dp to calculate the AF
            # set it appropriately here. Using dp = 100 for simplicity
            dp = 100
            new_ad_1 = round(dp*af)
            new_ad_0 = dp - new_ad_1
            new_gt_field+= ':{},{}:{}:1:1'.format(new_ad_0, new_ad_1, dp)
            out_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\tAF={}\t{}\n"\
                        .format(chrom, pos, '.', ref, alt, '.', '.', af, new_gt_field)
            if '{}-{}'.format(chrom, pos) in golden_snps:
                golden_vcf.write(out_line)
                out_vcf.write(out_line)
            else:
                out_vcf.write(out_line)
