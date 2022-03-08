#!/usr/bin/env python3
import sys
import pysam

vcf_i = sys.argv[1]
vcf_2 = sys.argv[2]

def get_data(vcf):
    # assumes only 1 variant at this position, takes first one
    #var = list(vcf.fetch(chrom, pos - 1, pos))[0]
    vcf = pysam.VariantFile(vcf, "r")
    data = {}
    for var in vcf:
        dp = None
        ad = None
        for sample in var.samples:
            ad = var.samples[sample]['AD']
            dp = int(var.samples[sample]['DP'])
        if ad is not None and dp is not None:
            if type(ad) is tuple or type(ad) is list:
                ad = ad[1]
            else:
                ad = ad
            af_workflow = int(ad) / dp

        # lofreq doesn't have sample (GT) field
        # parse data from INFO field
        # lofreq has DP, AF, and DP4 fields
        # DP4 is ref + count, ref - count, alt + count, alt - count
        tool_list = ['lofreq']
        if any(x in vcf.filename.decode() for x in tool_list):
            dp = var.info["DP"]
            info_af = var.info["AF"]
            if type(info_af) is tuple:
                af_workflow = float(var.info["AF"][0])
            elif type(info_af) in [float, str]:
                af_workflow = float(var.info["AF"])
            dp4 = var.info["DP4"]
            ad = "{},{}".format(sum(dp4[0:2]), sum(dp4[2:]))

        if ad is None or af_workflow is None or dp is None:
            print("Something went wrong. AF: {}; AD: {}; DP: {}".format(af_workflow, ad, dp))
            break
        else:
            data["{}-{}".format(var.contig, var.pos)] = {"af": round(af_workflow,6), "ad": ad, "dp": dp}
    
    return data

vcf_i_dict = get_data(vcf_i)
vcf_2_dict = get_data(vcf_2)

with open(vcf_i, 'r') as vcf, open('out.vcf', 'w') as out:
    gt_keys = "AD:DP:AF"
    for line in vcf:
        if line.startswith('#'):
            if line.startswith("#CHROM\tPOS"):
                # Add all sample (GT) fields for lofreq
                # Add just a new sample for the rest
                if 'lofreq' in vcf_i:
                    out.write(line.strip() + "\tFORMAT\tsample\tm2\n")
                else:
                    out.write(line.strip() + "\tm2\n")
            else:
                out.write(line)
        else:
            columns = line.split('\t')
            pos = columns[1]
            ref = columns[3]
            alt = columns[4]
            chrom = columns[0]
            vcf_i_data = vcf_i_dict["{}-{}".format(chrom, pos)]
            vcf_2_data = vcf_2_dict["{}-{}".format(chrom, pos)]
            af = round((vcf_i_data['af'] + vcf_2_data['af']) / 2, 6)
            dp = round((vcf_i_data['dp'] + vcf_2_data['dp']) / 2)
            vcf_i_gt = "{}:{}:{}".format(vcf_i_data['ad'], vcf_i_data['dp'], vcf_i_data['af'])
            vcf_2_gt = "{}:{}:{}".format(vcf_2_data['ad'], vcf_2_data['dp'], vcf_2_data['af'])
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\tAF={};DP={}\t{}\t{}\t{}\n".format(chrom, pos, '.', ref, alt, '.', '.', af, dp, gt_keys, vcf_i_gt, vcf_2_gt))
