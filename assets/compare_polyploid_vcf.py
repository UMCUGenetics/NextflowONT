#! /usr/bin/env python3

import argparse
import vcf


def ParseIllumina(illumina_vcf, min_ad_illumina):
    vcf_reader = vcf.Reader(open(illumina_vcf, 'rb'))
    region_chr = "chr5"
    region_start = 71300000
    region_end = 71500000
    sampleid = vcf_reader.samples[0]  ## Assume single sample VCF here!

    illumina_dic = {}
    for record in vcf_reader.fetch(region_chr, region_start, region_end):
        if int(record.genotype(sampleid)['DP']) >= min_ad_illumina and record.is_snp:
            chrom = record.CHROM
            pos = record.POS
            gt = record.genotype(sampleid)['GT']
            ad = record.genotype(sampleid)['AD']
            af = ad[1] / (ad[0] + ad[1]) * 100
            if f"{chrom}_{pos}" not in illumina_dic:
                illumina_dic[f"{chrom}_{pos}"] = {
                    "gt": gt,
                    "af":af,
                    "ishet":record.genotype(sampleid).is_het,
                    "issnp": record.is_snp,
                    "het": record.heterozygosity,
                    "filter": record.FILTER
                }

    return illumina_dic


def CompareONT(ont_vcf, illumina_dic, high_conf, low_conf):
    vcf_reader = vcf.Reader(open(ont_vcf, 'rb'))
    sampleid = vcf_reader.samples[0]  ## Assume single sample VCF here!

    tp_high = tp_low = fn = unknown = 0

    print(f"Conclusion\tIllumina_Chr\tIllumina_Position\tIllumina_GT\tIllumina_AF\tONT_Chr\tONT_Position\tONT_GT\tONT_AF")
    for region in illumina_dic:
        ill_chrom, ill_pos = region.split("_")
        ill_gt = illumina_dic[region]["gt"]
        ill_af = illumina_dic[region]["af"]
        ill_filter = illumina_dic[region]["filter"]
        chrom = pos = gt = ad = af = ''
        for record in vcf_reader.fetch(ill_chrom, int(ill_pos)-1, int(ill_pos)):
            chrom = record.CHROM
            pos = record.POS
            gt = record.genotype(sampleid)['GT']
            ad = record.genotype(sampleid)['AD']
            af = ad[1] / (ad[0] + ad[1]) * 100
            if not illumina_dic[region]["ishet"] and not record.genotype(sampleid).is_het:  # Both Illumina and ONT are homozygous
                print(f"TP_homozygous\t\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{chrom}\t{pos}\t{gt}\t{af}") 
                tp_high += 1
            else:
                if af > (ill_af - high_conf) and af < (ill_af + high_conf):
                    print(f"TP_heterozygous_{high_conf}\t\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{chrom}\t{pos}\t{gt}\t{af}")
                    tp_high += 1
                elif af > (ill_af - low_conf) and af < (ill_af + low_conf):
                    print(f"TP_heterozygous_{low_conf}\t\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{chrom}\t{pos}\t{gt}\t{af}")
                    tp_low += 1
                else:
                    if ill_filter:
                        ill_filters  = "_".join(ill_filter)
                        print(f"Possible_Filter_{ill_filters}\t\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{chrom}\t{pos}\t{gt}\t{af}")
                        unknown += 1
                    else:
                        print(f"Unknown\t\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{chrom}\t{pos}\t{gt}\t{af}")
                        unknown += 1

        if not pos:  # Variant not detected
            if ill_filter:
                ill_filters  = "_".join(ill_filter)
                print(f"FN_ONT_{ill_filters}\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\tn/a\tn/a\tn/a\tn/a")
                unknown += 1
            else:
                print(f"FN_ONT\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\tn/a\tn/a\tn/a\tn/a")
                fn += 1
    print(f"TP_high\tTP_low\tTP_total\tFN\tUnknown\tTotal")
    print(f"{tp_high}\t{tp_low}\t{tp_high + tp_low}\t{fn}\t{unknown}\t{tp_high + tp_low + fn + unknown}")
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('illumina_vcf', help='Full path to Illumina VCF')
    parser.add_argument('ont_vcf', help='Full path to ONT polyploid VCF')
    parser.add_argument('--min_ad_illumina', default= 30, type=int, help='Minimum Depth Illumina to be considered usefull variant')
    parser.add_argument('--high_conf', default= 10, type=int, help='Minimum Depth Illumina to be considered usefull variant')
    parser.add_argument('--low_conf', default= 20, type=int, help='Minimum Depth Illumina to be considered usefull variant')
    args = parser.parse_args()

    illumina_dic = ParseIllumina(args.illumina_vcf, args.min_ad_illumina)
    CompareONT(args.ont_vcf, illumina_dic, args.high_conf, args.low_conf)


