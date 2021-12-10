#!/usr/bin/env python

import io
import os
import sys
from collections import OrderedDict
import numpy
from datetime import datetime
import re
import math

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def parse_MuseSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample_normal = info[9].split(":")[1]
			ad_sample_tumor = info[10].split(":")[1]
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad_normal']=ad_sample_normal
			snvs[chrid]['ad_tumor']=ad_sample_tumor			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_MutectSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample_normal = info[9].split(":")[1]
			ad_sample_tumor = info[10].split(":")[1]
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad_normal']=ad_sample_normal
			snvs[chrid]['ad_tumor']=ad_sample_tumor			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_StrelkaSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample_normal = "."
			ad_sample_tumor = "."
			qual = info[5]
			filt = info[6]
			somatic_evs = info[7]
			somatic_evs = somatic_evs.split(";")[12]
			somatic_evs = somatic_evs.split("=")[1]
			snvs[chrid] = {}
			snvs[chrid]['ad_normal']=ad_sample_normal
			snvs[chrid]['ad_tumor']=ad_sample_tumor			
			snvs[chrid]['qual']=somatic_evs
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_FreeBayesSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample_normal = info[9].split(":")[1]
			ad_sample_tumor = info[10].split(":")[1]
			qual = info[5]
			filt = info[6]
			snvs[chrid] = {}
			snvs[chrid]['ad_normal']=ad_sample_normal
			snvs[chrid]['ad_tumor']=ad_sample_tumor			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_SomaticSniperSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			if len(info[9].split(":")) >= 11 :
				chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
				DP4_sample_normal = info[9].split(":")[3]
				DP4_sample_normal = DP4_sample_normal.split(",")
				ref_count_normal = int(DP4_sample_normal[0]) + int(DP4_sample_normal[1])
				alt_count_normal = int(DP4_sample_normal[2]) + int(DP4_sample_normal[3])
				ad_sample_normal = str(ref_count_normal) + "," + str(alt_count_normal)
				DP4_sample_tumor = info[10].split(":")[3]
				DP4_sample_tumor = DP4_sample_tumor.split(",")
				ref_count_tumor = int(DP4_sample_tumor[0]) + int(DP4_sample_tumor[1])
				alt_count_tumor = int(DP4_sample_tumor[2]) + int(DP4_sample_tumor[3])
				ad_sample_tumor = str(ref_count_tumor) + "," + str(alt_count_tumor)
				qual = info[10].split(":")[12]
				filt = info[6]
				snvs[chrid] = {}
				snvs[chrid]['ad_normal']=ad_sample_normal
				snvs[chrid]['ad_tumor']=ad_sample_tumor			
				snvs[chrid]['qual']=qual
				snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_VarScan2SNV(vcf):
	snvs = {}
	datacolumn = {}
	ss_dict = {"0":"Reference", "1":"Germline", "2":"Somatic", "3":"LOH", "5":"Unknown"}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			rd_normal = info[9].split(":")[3]
			ad_normal = info[9].split(":")[4]
			ad_sample_normal = rd_normal + "," + ad_normal
			rd_tumor = info[10].split(":")[3]
			ad_tumor = info[10].split(":")[4]
			ad_sample_tumor = rd_tumor + "," + ad_tumor
			somatic_score = info[7].split(";")[4]
			somatic_score = somatic_score.split("=")[1]
			statut = info[7].split(";")[3]
			statut = statut.split("=")[1]
			statut = ss_dict[statut]
			qual = somatic_score
			filt = statut
			snvs[chrid] = {}
			snvs[chrid]['ad_normal']=ad_sample_normal
			snvs[chrid]['ad_tumor']=ad_sample_tumor	
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_VarDictSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample_normal = info[9].split(":")[1]
			ad_sample_tumor = info[10].split(":")[1]
			status = info[7]
			status = status.split(";")[10]
			status = status.split("=")[1]
			qual = info[5]
			filt = status
			snvs[chrid] = {}
			snvs[chrid]['ad_normal']=ad_sample_normal
			snvs[chrid]['ad_tumor']=ad_sample_tumor			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_Mutect2SNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			ad_sample_normal = info[9].split(":")[1]
			ad_sample_tumor = info[10].split(":")[1]
			qual = info[5]
			filt = info[6]
			val = info[7]
			tlod = val.split(";")[-1]
			tlod = tlod.split("=")[1]
			qual = tlod
			snvs[chrid] = {}
			snvs[chrid]['ad_normal']=ad_sample_normal
			snvs[chrid]['ad_tumor']=ad_sample_tumor			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def parse_LoFreqSNV(vcf):
	snvs = {}
	datacolumn = {}
	for line in open(vcf, 'r'):
		line=line.strip()
		if not line.startswith("#"):
			info=line.split("\t")
			chrid = info[0] + '\t' + info[1] + '\t' + info[3] + '\t' + info[4]
			qual = info[5]
			filt = info[6]
			val = info[7]
			DP4 = val.split(";")[2]
			DP4 = DP4.split("=")[1]
			DP4 = DP4.split(",")
			ref_count = int(DP4[0]) + int(DP4[1])
			alt_count = int(DP4[2]) + int(DP4[3])
			ad_sample = str(ref_count) + "," + str(alt_count)
			snvs[chrid] = {}
			snvs[chrid]['ad_tumor']=ad_sample
			snvs[chrid]['ad_normal']="."			
			snvs[chrid]['qual']=qual
			snvs[chrid]['filter']=filt

	return {'snvs':snvs}

def get_af(ad):
	try :
		ref_count = float(ad.split(",")[0])
		alt_count = float(ad.split(",")[1])
		af = alt_count/(alt_count+ref_count)
	except ValueError :
		af = numpy.nan
	except ZeroDivisionError :
		af = numpy.nan
	return af

def mergeSNV(freebayes_snv, lofreq_snv, muse_snv, mutect_snv, mutect2_snv, sniper_snv, strelka_snv, vardict_snv, varscan2_snv, output):
	SAMPLE = os.path.basename(sys.argv[1])
	filter1=re.compile('(.*).snp.*')
	SAMPLE=filter1.search(sys.argv[1]).group(1)
	sf = open(output,"w")
	sf.write("%s\n" %("##fileformat=VCFv4.2"))
	sf.write("%s%s\n" %("##date=",str(datetime.now())))
	sf.write("%s\n" %("##source=MergeCaller"))
	sf.write("%s\n" %("##FILTER=<ID=FreeBayes,Description=\"Called by FreeBayes\""))
	sf.write("%s\n" %("##FILTER=<ID=LoFreq,Description=\"Called by LoFreq\""))
	sf.write("%s\n" %("##FILTER=<ID=Muse,Description=\"Called by Muse\""))
	sf.write("%s\n" %("##FILTER=<ID=Mutect,Description=\"Called by Mutect\""))
	sf.write("%s\n" %("##FILTER=<ID=Mutect2,Description=\"Called by Mutect2\""))
	sf.write("%s\n" %("##FILTER=<ID=SomaticSniper,Description=\"Called by SomaticSniper\""))
	sf.write("%s\n" %("##FILTER=<ID=Strelka,Description=\"Called by Strelka\""))
	sf.write("%s\n" %("##FILTER=<ID=Vardict,Description=\"Called by Vardict\""))
	sf.write("%s\n" %("##FILTER=<ID=Varscan2,Description=\"Called by Varscan2\""))
	sf.write("%s\n" %("##INFO=<ID=VAF_NORMAL,Number=1,Type=Float,Description=\"Median vaf between callers in normal\""))
	sf.write("%s\n" %("##INFO=<ID=VAF_TUMOR,Number=1,Type=Float,Description=\"Median vaf between callers in tumor\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADP1,Number=R,Type=Integer,Description=\"Allelic depths reported by FreeBayes for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADLF,Number=R,Type=Integer,Description=\"Allelic depths reported by LoFreq for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADMU,Number=R,Type=Integer,Description=\"Allelic depths reported by Muse for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADM,Number=R,Type=Integer,Description=\"Allelic depths reported by Mutect for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADM2,Number=R,Type=Integer,Description=\"Allelic depths reported by Mutect2 for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADSN,Number=R,Type=Integer,Description=\"Allelic depths reported by SomaticSniper for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADSK,Number=R,Type=Integer,Description=\"Allelic depths reported by Strelka for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADVC,Number=R,Type=Integer,Description=\"Allelic depths reported by Vardict for the ref and alt alleles in the order listed\""))
	sf.write("%s\n" %("##FORMAT=<ID=ADVS2,Number=R,Type=Integer,Description=\"Allelic depths reported by Varscan2 for the ref and alt alleles in the order listed\""))
	sf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('#CHROM', 'POS','ID', 'REF', 'ALT','QUAL', 'FILTER', 'INFO','FORMAT', "NORMAL", "TUMOR"))
	all_snvs = sorted(set( list(freebayes_snv['snvs'].keys()) +  list(lofreq_snv['snvs'].keys()) + list(muse_snv['snvs'].keys()) + list(mutect_snv['snvs'].keys()) + list(mutect2_snv['snvs'].keys()) + list(sniper_snv['snvs'].keys()) + list(strelka_snv['snvs'].keys()) + list(vardict_snv['snvs'].keys()) + list(varscan2_snv['snvs'].keys()) ))
	for snv in all_snvs :
		vcfinfo = {}
		if snv in freebayes_snv['snvs'] :
			vcfinfo['freebayes']=snv
		if snv in lofreq_snv['snvs'] :
			vcfinfo['lofreq']=snv
		if snv in muse_snv['snvs']:
			vcfinfo['muse']=snv
		if snv in mutect_snv['snvs'] :
			vcfinfo['mutect']=snv
		if snv in mutect2_snv['snvs'] :
			vcfinfo['mutect2']=snv
		if snv in sniper_snv['snvs'] :
			vcfinfo['sniper']=snv
		if snv in strelka_snv['snvs'] :
			vcfinfo['strelka']=snv
		if snv in vardict_snv['snvs'] :
			vcfinfo['vardict']=snv
		if snv in varscan2_snv['snvs'] :
			vcfinfo['varscan2']=snv
		called_by = list(vcfinfo.keys())
		if all(value == vcfinfo[called_by[0]] for value in vcfinfo.values()):
			format=''
			gf_normal=''
			gf_tumor=''
			callers=''
			nb_callers_pass=0
			af_normal = []
			af_tumor = []
			for c in called_by :
				if c=='freebayes':
					format=format+'ADP1:'
					gf_normal=gf_normal+freebayes_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(freebayes_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+freebayes_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(freebayes_snv['snvs'][snv]['ad_tumor']))
					if freebayes_snv['snvs'][snv]['filter'] != "REJECT" :
						nb_callers_pass += 1
						callers=callers+'FreeBayes|'
				elif c=='lofreq':
					#filter1=re.compile('min_dp_10')
					#filter2=re.compile('sb_fdr')
					#filter3=re.compile('min_snvqual_76')
					#filter4=re.compile('min_indelqual_54')
					#f1=filter1.search(lofreq_snv['snvs'][snv]['filter'])
					#f2=filter2.search(lofreq_snv['snvs'][snv]['filter'])
					#f3=filter3.search(lofreq_snv['snvs'][snv]['filter'])
					#f4=filter4.search(lofreq_snv['snvs'][snv]['filter'])
					format=format+'ADLF:'
					gf_normal=gf_normal+lofreq_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(lofreq_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+lofreq_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(lofreq_snv['snvs'][snv]['ad_tumor']))
					nb_callers_pass +=1
					callers=callers+'Lofreq|'
				elif c=="muse":
					format=format+'ADMU:'
					gf_normal=gf_normal+muse_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(muse_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+muse_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(muse_snv['snvs'][snv]['ad_tumor']))
					filter1 = re.compile('Tier1')
					filter2 = re.compile('Tier2')
					filter3 = re.compile('Tier3')
					filter4 = re.compile('Tier4')
					filter5 = re.compile('Tier5')
					f1=filter1.search(muse_snv['snvs'][snv]['filter'])
					f2=filter2.search(muse_snv['snvs'][snv]['filter'])
					f3=filter3.search(muse_snv['snvs'][snv]['filter'])
					f4=filter4.search(muse_snv['snvs'][snv]['filter'])
					f5=filter5.search(muse_snv['snvs'][snv]['filter'])
					if not (f1 or f2 or f3 or f4 or f5) :
						nb_callers_pass += 1
						callers=callers+'Muse|'
				elif c=="mutect":
					format=format+'ADMU:'
					gf_normal=gf_normal+mutect_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(mutect_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+mutect_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(mutect_snv['snvs'][snv]['ad_tumor']))
					if mutect_snv['snvs'][snv]['filter'] == "PASS" :
						nb_callers_pass += 1
						callers=callers+'Mutect|'
				elif c=='mutect2':
					filter1=re.compile('alt_allele_in_normal')
					filter2=re.compile('clustered_events')
					filter3=re.compile('germline_risk')
					filter4=re.compile('homologous_mapping_event')
					filter5=re.compile('multi_event_alt_allele_in_normal')
					filter6=re.compile('panel_of_normals')
					filter7=re.compile('str_contraction')
					filter8=re.compile('t_lod_fstar')
					filter9=re.compile('triallelic_site')
					filter10=re.compile('germline')
					f1=filter1.search(mutect2_snv['snvs'][snv]['filter'])
					f2=filter2.search(mutect2_snv['snvs'][snv]['filter'])
					f3=filter3.search(mutect2_snv['snvs'][snv]['filter'])
					f4=filter4.search(mutect2_snv['snvs'][snv]['filter'])
					f5=filter5.search(mutect2_snv['snvs'][snv]['filter'])
					f6=filter6.search(mutect2_snv['snvs'][snv]['filter'])
					f7=filter7.search(mutect2_snv['snvs'][snv]['filter'])
					f8=filter8.search(mutect2_snv['snvs'][snv]['filter'])
					f9=filter9.search(mutect2_snv['snvs'][snv]['filter'])
					f10=filter10.search(mutect2_snv['snvs'][snv]['filter'])
					format=format+'ADM2:'
					gf_normal=gf_normal+mutect2_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(mutect2_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+mutect2_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(mutect2_snv['snvs'][snv]['ad_tumor']))
					if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9 or f10) :
						nb_callers_pass += 1
						callers=callers+'Mutect2|'
				elif c=='sniper' :
					format=format+'ADSN:'
					gf_normal=gf_normal+sniper_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(sniper_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+sniper_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(sniper_snv['snvs'][snv]['ad_tumor']))
					qual = float(sniper_snv['snvs'][snv]['qual'])
					if qual > 15 :
						nb_callers_pass += 1
						callers=callers+'SomaticSniper|'
				elif c=='strelka' :
					format=format+'ADSK:'
					gf_normal=gf_normal+strelka_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(strelka_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+strelka_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(strelka_snv['snvs'][snv]['ad_tumor']))
					if strelka_snv['snvs'][snv]['filter'] == "PASS" :
						nb_callers_pass += 1
						callers=callers+'Strelka|'
				elif c=='vardict':
					#filter1=re.compile('AMPBIAS')
					#filter2=re.compile('Bias')
					#filter3=re.compile('Cluster0bp')
					#filter4=re.compile('InGap')
					#filter5=re.compile('InIns')
					#filter6=re.compile('LongMSI')
					#filter7=re.compile('MSI12')
					#filter8=re.compile('NM5.25')
					#filter9=re.compile('Q10')
					#filter10=re.compile('SN1.5')
					#filter11=re.compile('d3')
					#filter12=re.compile('f0.02')
					#filter13=re.compile('p8')
					#filter14=re.compile('pSTD')
					#filter15=re.compile('q22.5')
					#filter16=re.compile('v2')
					#f1=filter1.search(vardict_snv['snvs'][snv]['filter'])
					#f2=filter2.search(vardict_snv['snvs'][snv]['filter'])
					#f3=filter3.search(vardict_snv['snvs'][snv]['filter'])
					#f4=filter4.search(vardict_snv['snvs'][snv]['filter'])
					#f5=filter5.search(vardict_snv['snvs'][snv]['filter'])
					#f6=filter6.search(vardict_snv['snvs'][snv]['filter'])
					#f7=filter7.search(vardict_snv['snvs'][snv]['filter'])
					#f8=filter8.search(vardict_snv['snvs'][snv]['filter'])
					#f9=filter9.search(vardict_snv['snvs'][snv]['filter'])
					#f10=filter10.search(vardict_snv['snvs'][snv]['filter'])
					#f11=filter11.search(vardict_snv['snvs'][snv]['filter'])
					#f12=filter12.search(vardict_snv['snvs'][snv]['filter'])
					#f13=filter13.search(vardict_snv['snvs'][snv]['filter'])
					#f14=filter14.search(vardict_snv['snvs'][snv]['filter'])
					#f15=filter15.search(vardict_snv['snvs'][snv]['filter'])
					#f16=filter16.search(vardict_snv['snvs'][snv]['filter'])
					format=format+'ADVC:'
					gf_normal=gf_normal+vardict_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(vardict_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+vardict_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(vardict_snv['snvs'][snv]['ad_tumor']))
					#if not (f1 or f2 or f3 or f4 or f5 or f6 or f7 or f8 or f9 or f10 or f11 or f12 or f13 or f14 or f15 or f16) :
					if vardict_snv['snvs'][snv]['filter'] == 'StrongSomatic' or vardict_snv['snvs'][snv]['filter'] == 'LikelySomatic' :
						callers=callers+'Vardict|'
						nb_callers_pass += 1
				elif c=='varscan2':
					format=format+'ADVS2:'
					gf_normal=gf_normal+varscan2_snv['snvs'][snv]['ad_normal']+':'
					af_normal.append(get_af(varscan2_snv['snvs'][snv]['ad_normal']))
					gf_tumor=gf_tumor+varscan2_snv['snvs'][snv]['ad_tumor']+':'
					af_tumor.append(get_af(varscan2_snv['snvs'][snv]['ad_tumor']))
					if float(varscan2_snv['snvs'][snv]['qual']) > 30 and varscan2_snv['snvs'][snv]['filter'] == "Somatic" :
						nb_callers_pass += 1
						callers=callers+'Varscan2|'

			if nb_callers_pass > 0 : 
				vaf_tumor = round(numpy.nanmedian(af_tumor),4)
				vaf_normal = round(numpy.nanmedian(af_normal),4)
				callers = callers[:-1]
				if vaf_tumor < 0.05 :
					filt = "LowVariantFreq|"+callers
				else :
					filt = callers
				if nb_callers_pass > 4 :
					filt =  "CONCORDANT|"+filt
				else :
					filt = "DISCORDANT|"+filt
				info = "VAF_NORMAL="+str(vaf_normal)+";"+"VAF_TUMOR="+str(vaf_tumor)
				format = format[:-1]
				gf_normal = gf_normal[:-1]
				gf_tumor = gf_tumor[:-1]
				qual=nb_callers_pass
				vcfinfolist=vcfinfo[called_by[0]].split("\t")
				baseinfo=vcfinfolist[0]+'\t'+vcfinfolist[1]+'\t.\t'+vcfinfolist[2]+'\t'+vcfinfolist[3]
				sf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(baseinfo,qual, filt, info, format, gf_normal, gf_tumor))
		else :
			print("Conflict in ref and alt alleles between callers at pos "+snv)

freebayes_snv = parse_FreeBayesSNV(sys.argv[1])
lofreq_snv = parse_LoFreqSNV(sys.argv[2])
muse_snv = parse_MuseSNV(sys.argv[3])
mutect_snv = parse_MutectSNV(sys.argv[4])
mutect2_snv = parse_Mutect2SNV(sys.argv[5])
sniper_snv = parse_SomaticSniperSNV(sys.argv[6])
strelka_snv = parse_StrelkaSNV(sys.argv[7])
vardict_snv = parse_VarDictSNV(sys.argv[8])
varscan2_snv = parse_VarScan2SNV(sys.argv[9])
output = sys.argv[10]

mergeSNV(freebayes_snv, lofreq_snv, muse_snv, mutect_snv, mutect2_snv, sniper_snv, strelka_snv, vardict_snv, varscan2_snv, output)
