#!/usr/bin/env python3
# requires: pysam/intel/python3.6/0.14.1
# Expects you to have run the command
# `bcftools isec $golden_vcf $workflow_vcf -p isec`
# Expects 3 files 000[123].vcf in a dir called isec

import pysam
import sys
import json
import csv
import argparse

def main():
    ## Retrieve command line arguments
    in_arg = get_input_args()
    golden_vcf = in_arg.golden_vcf
    workflow_vcf = in_arg.workflow_vcf

    # Output from bcftools isec
    false_negatives_vfc = 'isec/0000.vcf'
    false_positives_vcf = 'isec/0001.vcf'
    true_positives_vcf = 'isec/0002.vcf'

    # These are lists which hold the positions of the fn, fp, and tp vars
    false_negatives = get_positions_from_vcf(false_negatives_vfc)
    false_positives = get_positions_from_vcf(false_positives_vcf)
    true_positives = get_positions_from_vcf(true_positives_vcf)

    results = []

    vcf_golden = pysam.VariantFile(golden_vcf)
    vcf_workflow = pysam.VariantFile(workflow_vcf)
    
    sample_id = workflow_vcf.split(".vcf")[0]

    for x in false_negatives:
        data = get_data_golden(vcf_golden, x)
        results.append({
            'sample_id': sample_id, 
            'chrom': x[0],
            'pos': x[1], 
            'af_golden': data['af'], 
            'af_workflow': 0, 
            'ref': data['ref'], 
            'alt': data['alt'], 
            'dp': data['dp']
        })

    for x in false_positives:
        data = get_data_workflow(vcf_workflow, x)
        results.append({
            'sample_id': sample_id, 
            'chrom': x[0],
            'pos': x[1], 
            'af_golden': 0, 
            'af_workflow': data['af'], 
            'ref': data['ref'], 
            'alt': data['alt'], 
            'dp': data['dp']
        })

    for x in true_positives:
        gold_data = get_data_golden(vcf_golden, x)
        workflow_data = get_data_workflow(vcf_workflow, x)
        results.append({
            'sample_id': sample_id, 
            'chrom': x[0],  
            'pos': x[1], 
            'af_golden': gold_data['af'], 
            'af_workflow': workflow_data['af'], 
            'ref': workflow_data['ref'], 
            'alt': workflow_data['alt'], 
            'dp': workflow_data['dp']
        })

    with open("{}.csv".format(in_arg.out), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames = results[0].keys())
        w.writeheader()
        w.writerows(results)
    with open("{}_snp_count.txt".format(in_arg.out), 'w') as f:
            f.write("{}\t{} ".format(in_arg.out, len(true_positives) + len(false_negatives)))

def get_input_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--golden_vcf', type = str, required = True)
    parser.add_argument('--workflow_vcf', type = str, required = True)
    parser.add_argument('--out', type = str, required = True)
    in_arg = parser.parse_args()
    return in_arg

def get_positions_from_vcf(vcf):
    l = []
    vcf = pysam.VariantFile(vcf, "r")
    for var in vcf: 
        l.append([var.contig, var.pos])
    return l

def get_data_workflow(vcf, x):
    # Parse data from variant caller VCFs
    # assumes only 1 variant at this position, takes first one
    chrom = x[0]
    pos = x[1]
    variants = list(vcf.fetch(chrom, pos - 1, pos))
    for v in variants:
        if v.pos == pos:
            var = v
            break
    # Set dp = None as default
    dp = None
    #for sample in var.samples:
    #    ad = var.samples[sample]['AD']
    #    dp = int(var.samples[sample]['DP'])

    # Get AF directly from INFO AF field for some tools if possible
    # Else, compute using AD/DP
    # Use this AF for the _mx_ VCFs because this is the 
    # average AF, what we want in this case (alternatively, 
    # compute it using AD and DP from both samples). 
    use_af_tool_list = ['cliquesnv', 'lofreq', '_mx_']
    if any(x in vcf.filename.decode() for x in use_af_tool_list):
        info_af = var.info["AF"]
        if type(info_af) is tuple:
            af_workflow = float(var.info["AF"][0])
        elif type(info_af) in [float, str]:
            af_workflow = float(var.info["AF"])
        if dp is None and "DP" in var.info:
                dp = var.info["DP"]
    #elif 'varscan' in vcf.filename.decode():
    #    for sample in var.samples:
    #            freq = var.samples[sample]['FREQ']
    #            af_workflow = float(freq.replace('%',''))/100
    else:
        for sample in var.samples:
            ad = var.samples[sample]['AD']
            dp = int(var.samples[sample]['DP'])
                
            if (type(ad) is tuple or type(ad) is list):
                if len(ad) > 1:
                    ad = ad[1]
                elif len(ad) == 1:
                    ad = ad[0]
                else:
                    print("Error: Unexpected AD field", ad)
                    return False
            else:
                ad = ad

        af_workflow = int(ad) / dp
    
    return {"af": round(af_workflow, 6), "dp": dp, "ref": var.ref, "alt": var.alts[0]}

def get_data_golden(vcf, x):
    chrom = x[0]
    pos = x[1]
    var = list(vcf.fetch(chrom, pos - 1, pos))[0]
    for sample in var.samples:
        dp = int(var.samples[sample]['DP'])
    af_golden = round(float(var.info["AF"][0]),2)
    return {"af": af_golden, "dp": dp, "ref": var.ref, "alt": var.alts[0]}

if __name__ == "__main__":
    main()

