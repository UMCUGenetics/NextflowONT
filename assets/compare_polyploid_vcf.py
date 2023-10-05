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
                    "filter": record.FILTER,
                    "qd": record.INFO['QD'],
                    "mq": record.INFO['MQ'],
                    "fs": record.INFO['FS'],
                }

    return illumina_dic


def CompareONT(ont_vcf, illumina_dic, high_conf, low_conf):
    print(f"Conclusion\tILL_Chr\tILL_Position\tILL_GT\tILL_AF\tILL_QD\tILL_MQ\tILL_FS\tONT_Chr\tONT_Position\tONT_GT\tONT_AF\tONT_QD\tONT_MQ\tONT_FS")
    vcf_reader = vcf.Reader(open(ont_vcf, 'rb'))
    if len(vcf_reader.samples) > 1:
        sys.exit("for single sample VCFs only. Exiting")
    sampleid = vcf_reader.samples[0]  ## Assume single sample VCF here!
    tp_high = tp_low = fn = fn_homopolymer = unknown = 0
    for region in illumina_dic:
        ill_chrom, ill_pos = region.split("_")
        ill_gt = illumina_dic[region]["gt"]
        ill_af = illumina_dic[region]["af"]
        ill_filter = illumina_dic[region]["filter"]
        ill_qd = illumina_dic[region]["qd"]
        ill_mq = illumina_dic[region]["mq"]
        ill_fs = illumina_dic[region]["fs"]

        chrom = pos = gt = ad = af = ''
        for record in vcf_reader.fetch(ill_chrom, int(ill_pos)-1, int(ill_pos)):
            chrom = record.CHROM
            pos = record.POS
            gt = record.genotype(sampleid)['GT']
            ad = record.genotype(sampleid)['AD']
            af = ad[1] / (ad[0] + ad[1]) * 100
            qd = record.INFO['QD']
            mq = record.INFO['MQ']
            fs = record.INFO['FS']

            if not illumina_dic[region]["ishet"] and not record.genotype(sampleid).is_het:  # Both Illumina and ONT are homozygous
                print(f"TP_homozygous\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{ill_qd}\t{ill_mq}\t{ill_fs}\t{chrom}\t{pos}\t{gt}\t{af:.2f}\t{qd}\t{mq}\t{fs}") 
                tp_high += 1
            else:
                if af > (ill_af - high_conf) and af < (ill_af + high_conf):
                    print(f"TP_heterozygous_high_conf\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{ill_qd}\t{ill_mq}\t{ill_fs}\t{chrom}\t{pos}\t{gt}\t{af:.2f}\t{qd}\t{mq}\t{fs}")
                    tp_high += 1
                elif af > (ill_af - low_conf) and af < (ill_af + low_conf):
                    print(f"TP_heterozygous_low_conf\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{ill_qd}\t{ill_mq}\t{ill_fs}\t{chrom}\t{pos}\t{gt}\t{af:.2f}\t{qd}\t{mq}\t{fs}")
                    tp_low += 1
                else:
                    if ill_filter:
                        ill_filters  = "_".join(ill_filter)
                        print(f"Possible_Filter_{ill_filters}\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{ill_qd}\t{ill_mq}\t{ill_fs}\t{chrom}\t{pos}\t{gt}\t{af:.2f}\t{qd}\t{mq}\t{fs}")
                        unknown += 1
                    else:
                        print(f"Unknown\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{ill_qd}\t{ill_mq}\t{ill_fs}\t{chrom}\t{pos}\t{gt}\t{af:.2f}\t{qd}\t{mq}\t{fs}")
                        unknown += 1

        if not pos:  # Variant not detected
            if ill_filter:
                ill_filters  = "_".join(ill_filter)
                print(f"FN_ONT_{ill_filters}\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{ill_qd}\t{ill_mq}\t{ill_fs}\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\t/na")
                if "Homopolymer" in ill_filters:
                    fn_homopolymer += 1
                else:
                    unknown += 1
            else:
                print(f"FN_ONT\t{ill_chrom}\t{ill_pos}\t{ill_gt}\t{ill_af:.2f}\t{ill_qd}\t{ill_mq}\t{ill_fs}\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a\t/na")
                fn += 1
    print(f"TP_high\tTP_low\tTP_total\tFN\tFN_due_to_homopolymer\tUnknown\tTotal")
    print(f"{tp_high}\t{tp_low}\t{tp_high + tp_low}\t{fn}\t{fn_homopolymer}\t{unknown}\t{tp_high + tp_low + fn + fn_homopolymer + unknown}")
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('illumina_vcf', help='Full path to Illumina VCF')
    parser.add_argument('ont_vcf', help='Full path to ONT polyploid VCF')
    parser.add_argument('--min_ad_illumina', default= 30, type=int, help='Minimum Depth Illumina to be considered usefull variant')
    parser.add_argument('--high_conf', default= 10, type=int, help='percentage range for high confident af')
    parser.add_argument('--low_conf', default= 20, type=int, help='percentage range for lower confident af')
    parser.add_argument('--ploidy', type=int, help='ploidy of sample. Is used to calculate high and low confident af')
    args = parser.parse_args()

    high_conf = args.high_conf
    low_conf = args.low_conf
    if args.ploidy:
        high_conf = 100 / args.ploidy / 2
        low_conf = 100 / args.ploidy 

    illumina_dic = ParseIllumina(args.illumina_vcf, args.min_ad_illumina)
    CompareONT(args.ont_vcf, illumina_dic, high_conf, low_conf)


