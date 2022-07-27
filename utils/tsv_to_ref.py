#!/usr/bin/env python
import sys
from glob import glob
import textwrap
import operator

strain = sys.argv[1]
vcf_out = "{}.vcf".format(strain)
fa_out = "{}.fa".format(strain)
gtf_out = "{}.gtf".format(strain)

tsvs = glob("{}_*.tsv".format(strain))

with open(vcf_out, 'w') as vcf:
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("##reference=file://{}.fa\n".format(strain))
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for tsv in sorted(tsvs):
        ref_seq = ""
        print("Working with {}".format(tsv))
        with open(tsv, 'r') as infile:
            segment = tsv.replace("{}_".format(strain),'').replace(".tsv",'')
            infile.readline() #skip the header
            for line in infile:
                columns = line.split('\t')
                pos = columns[0]
                if segment == "MP":
                    if strain in ["H3N2", "H1N1"]:
                        pos = int(pos) - 25
                    elif strain == "Vic":
                        pos = int(pos) - 7
                elif segment == "NS":
                    if strain in ["H3N2", "H1N1"]:
                        pos = int(pos) - 26
                    elif strain == "Vic":
                        pos = int(pos) - 27
                
                ref = columns[2]
                alleles = {
                    "A": int(columns[3]),
                    "T": int(columns[4]),
                    "G": int(columns[5]),
                    "C": int(columns[6]),
                }
                
                # Do some validation before adding to VCF and FASTA
                if ref == "N" or (ref == "-" and pos != "N/A"):
                    ref = max(alleles.items(), key=operator.itemgetter(1))[0]
                    print("ref is N. changing to: {}. Strain: {}. Segment: {}. Line: {}" \
                            .format(ref, strain, segment, line))
                elif ref not in alleles.keys():
                    print("Line skipped. ref {} not in A,T,C,G. Strain: {}. Segment: {}. Pos: {}." \
                            .format(ref, strain, segment, pos))
                    continue

                # Build the FASTA
                ref_seq += ref

                # If the ref count = # sequences, there are no
                # variants at this position, skip it
                if alleles[ref] == columns[8].strip():
                    continue
                else:
                    alts_list = []
                    for k in [k for k,v in alleles.items() if k is not ref and v > 0]:
                        alts_list.append(k)
                    # Sometimes the # of ref alleles != # of sequences (columns[8])
                    # But there are no alt alleles, in this case, skip it
                    if len(alts_list) == 0:
                        continue
                    alts = ",".join(alts_list)
                        
                    vcf.write("{}_{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n" \
                            .format(strain, segment, pos, '.', ref, alts, '.', '.', '.'))
            # Because we are appending, if you run this script twice it will append to
            # the existing fasta. Make it initialize an empty .fa file at the start
            with open(fa_out, "a") as fa:
                fa.write('>{}_{}\n'.format(strain, segment))
                broken = textwrap.wrap(ref_seq, 70)
                fa.write('\n'.join(broken))
                fa.write('\n')
            with open(gtf_out, "a") as gtf:
                gtf.write("{}_{}\t.\tCDS\t1\t{}\t.\t.\t.\tgene_id \"{}_{}\"; transcript_id \"{}_{}.1\";\n".format(strain, segment, len(ref_seq), strain, segment, strain, segment))
