#! /usr/bin/env python3

import argparse
import vcfpy


def ParseVCF(input_vcf, region):
    chrom = region.split(":")[0]
    start = int(region.split(":")[1].split("-")[0])
    stop = int(region.split(":")[1].split("-")[1])

    ps = []
    reader = vcfpy.Reader.from_path(input_vcf)
    for record in reader.fetch(chrom, begin=start, end=stop):
        phaseset = record.calls[0].data.get('PS')  # Assuming single sample VCF
        if phaseset:
            ps.append(phaseset)
    if len(set(ps)) == 1:
        return list(set(ps))[0]
    else:  # if not unique ps is detected in roi, pick most prevelant
        count_dic = {}
        total_count = 0
        for ps_group in ps:
            if ps_group not in count_dic:
                count_dic[ps_group] = 1
            else:
                count_dic[ps_group] += 1
            total_count += 1
        for ps_group in count_dic:
            if count_dic[ps_group]/total_count > args.freq:
                return ps_group
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_vcf', help='full path to input VCF')
    parser.add_argument('roi', help='region of interest [chr:start-stop]')
    parser.add_argument(
        '--freq',
        type=float,
        default=0.7,
        help='threshold to determine most likely PS if multiple are detected in roi [default = 0.7]'
    )
    args = parser.parse_args()
    phaseset = ParseVCF(args.input_vcf, args.roi)
    print(phaseset)
