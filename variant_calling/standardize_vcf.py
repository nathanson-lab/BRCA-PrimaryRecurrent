#############################
#Brad Wubbenhorst
#bwubb@pennmedicine.upenn.edu
#Mar. 2021

#This code is used to standardize tumor-normal vcf files
#Strelka2, Mutect2, VarDictJava, VarScan2.

import os
import vcfpy
import argparse
from copy import copy

def _tumor_normal_allele_freqs(ref,alt,calls):
    """
    Retreive AD, AF fields for Strelka vcf files
    Makes best guess attempt
    """
    def _snv_freq(call):
        try:
            used_depth=call.data['DP']-call.data['FDP']#-call.data['SDP']
            _aad=call.data[f'{alt[0].value}U'][0]
            _aaf=float(_aad)/float(used_depth)
            _rad=used_depth-_aad
        except (ZeroDivisionError,TypeError) as e:
            print(f'{e}: {ref} {alt[0].value} {call}')
            _rad=0
            _aad=0
            _aaf=0.00
        return [_rad,_aad],float(f"{_aaf:.2f}")
    #DUP,DEL,INS,SNV,MNV
    def _indel_freq(call):
        try:
            used_depth=call.data['DP']
            _aad=call.data['TIR'][0]
            _aaf=float(_aad)/float(used_depth)
            _rad=used_depth-_aad
        except (ZeroDivisionError,TypeError) as e:
            print(f'{e}: {ref} {alt} {call}')
            _rad=0
            _aad=0
            _aaf=0.00
        return [_rad,_aad],float(f"{_aaf:.2f}")
    
    for call in calls:
        if alt[0].type=='SNV':
            _ad,_aaf=_snv_freq(call)
        elif alt[0].type in ['DEL','DUP','INS']:
            _ad,_aaf=_indel_freq(call)
        else:
            print(f'{alt[0].type} type UNKNOWN')
        call.data['AD']=_ad
        call.data['AAF']=_aaf
    return calls

def _tumor_normal_genotypes(ref,alt,info):
    """Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.
    Normal -- NT field (ref, het, hom, conflict)
    Tumor -- SGT field
      - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
        can also represent more complex alleles in which case we set at heterozygotes
        pending longer term inclusion of genotypes in Strelka2 directly
        (https://github.com/Illumina/strelka/issues/16)
      - For indels, uses the ref, het, hom convention
      
    Modified from https://github.com/bcbio/bcbio-nextgen strelka.py
    """
    known_names=set(["het","hom","ref","conflict"])
    def name_to_gt(val):
        if val.lower()=="het":
            return "0/1"
        elif val.lower()=="hom":
            return "1/1"
        elif val.lower() in set(["ref","confict"]):
            return "0/0"
        else:
            return "0/1"
    def alleles_to_gt(val):
        gt_indices = {gt.upper(): i for i, gt in enumerate([ref]+[a.value for a in alt])} #combo list of ref and alt, enumerated for GT integers
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts)==0:
                tumor_gt="0/0"
            elif 0 in tumor_gts:
                tumor_gt=f"0/{min([x for x in tumor_gts if x > 0])}" 
            else:
                tumor_gt=f"{min(tumor_gts)}/{max(tumor_gts)}"
        else:
            tumor_gt=name_to_gt(val)
        return tumor_gt
    nt_val=info['NT']
    normal_gt=name_to_gt(nt_val)
    sgt_val=info.get('SGT',False)
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val=sgt_val.split("->")[-1]
        tumor_gt=alleles_to_gt(sgt_val)
    return tumor_gt,normal_gt

def _sample_names(names,tumor,normal):
    D={'NORMAL':normal,'TUMOR':tumor}
    new_names=[]
    for name in names:
        if name in D:
            new_names.append(D[name])
        else:
            new_names.append(name)
    else:
        return new_names

def check_FORMAT(FORMAT):
    add_format=[f for f in ['GT','AD','AAF'] if f not in FORMAT]
    return add_format+FORMAT

def modify_outheader(outheader):
    info_lines={'CALLER':vcfpy.OrderedDict([('ID','CALLER'),('Number','1'),('Type','String'),('Description','Variant Call method')]),
    'SS':vcfpy.OrderedDict([('ID','SS'),('Number','1'),('Type','String'),('Description','Somatic Status from respective Call Method')]),
    'OFS':vcfpy.OrderedDict([('ID','OFS'),('Number','.'),('Type','String'),('Description','Original FILTER State')])}
    #
    format_lines={'GT':vcfpy.OrderedDict([('ID','GT'),('Number','1'),('Type','String'),('Description','Genotype')]),
    'AD':vcfpy.OrderedDict([('ID','AD'),('Number','R'),('Type','Integer'),('Description','Alt Allele Depth')]),
    'AAF':vcfpy.OrderedDict([('ID','AAF'),('Number','1'),('Type','Float'),('Description','Alt Allele Frequency')])}
    #
    for i in info_lines.keys():
        if i not in outheader.info_ids():
            outheader.add_info_line(info_lines[i])
        else:
            outheader.get_info_field_info(i).mapping.update(info_lines[i])
    for f in format_lines.keys():
        if f not in outheader.format_ids():
            outheader.add_format_line(format_lines[f])
        else:
            outheader.get_format_field_info(f).mapping.update(format_lines[f])
    return outheader

def make_outfile(infile):
    dir,file=os.path.split(infile)
    outfile="{0}.std.{1}.{2}".format(*file.rsplit('.',2))
    return f"{dir}/{outfile}"

def std_mutect2(infile,tumor,normal):
    myvcf=vcfpy.Reader.from_path(infile)
    myvcf.header.samples.names=_sample_names(myvcf.header.samples.names,tumor,normal)
    myvcf.header.samples.name_to_idx={k:v for v,k in enumerate(myvcf.header.samples.names)}
    #
    outheader=modify_outheader(myvcf.header.copy())
    outheader.samples.names=[tumor,normal]
    outfile=make_outfile(infile)
    writer=vcfpy.Writer.from_path(outfile,outheader)
    #
    for record in myvcf:
        outrecord=copy(record)
        E={'CALLER':'Mutect2'}
        E['OFS']=record.FILTER
        if outrecord.FILTER!=['PASS']:
            E['SS']='REJECT'
        else:
            E['SS']='SOMATIC'
        outrecord.INFO.update(E)
        outrecord.FORMAT=check_FORMAT(outrecord.FORMAT)
        ref=outrecord.REF
        alt=outrecord.ALT
        info=outrecord.INFO
        for call in outrecord.calls:
            try:
                used_depth=call.data['DP']
                _aad=call.data['AD'][1]
                _aaf=float(_aad)/float(used_depth)
            except (ZeroDivisionError,TypeError) as e:
                print(f'{e}: {ref} {alt[0].value} {call}')
                _aaf=0.00
            call.data['AAF']=float(f"{_aaf:.2f}")
        writer.write_record(outrecord)
    writer.close()
    myvcf.close()

def std_strelka2(infile,tumor,normal):
    myvcf=vcfpy.Reader.from_path(infile)
    myvcf.header.samples.names=_sample_names(myvcf.header.samples.names,tumor,normal)
    myvcf.header.samples.name_to_idx={k:v for v,k in enumerate(myvcf.header.samples.names)}
    #
    outheader=modify_outheader(myvcf.header.copy())
    outheader.samples.names=[tumor,normal]
    outfile=make_outfile(infile)
    writer=vcfpy.Writer.from_path(outfile,outheader)
    #
    for record in myvcf:
        outrecord=copy(record)
        outrecord.FORMAT=check_FORMAT(outrecord.FORMAT)#Check if these items exist
        E={'CALLER':'Strelka2'}
        E['OFS']=record.FILTER
        if record.INFO.get('SOMATIC',False):
            E['SS']='SOMATIC'
        else:
            E['SS']='REJECT'
        outrecord.INFO.update(E)
        ref=outrecord.REF
        alt=outrecord.ALT
        info=outrecord.INFO
        tumor_gt,normal_gt=_tumor_normal_genotypes(ref, alt, info)
        outrecord.calls[0].data['GT']=normal_gt
        outrecord.calls[1].data['GT']=tumor_gt
        outrecord.calls=_tumor_normal_allele_freqs(ref,alt,outrecord.calls)
        writer.write_record(outrecord)
    writer.close()
    myvcf.close()

def std_vardict(infile,tumor,normal):
    myvcf=vcfpy.Reader.from_path(infile)
    myvcf.header.samples.names=_sample_names(myvcf.header.samples.names,tumor,normal)
    myvcf.header.samples.name_to_idx={k:v for v,k in enumerate(myvcf.header.samples.names)}
    #
    outheader=modify_outheader(myvcf.header.copy())
    outheader.samples.names=[tumor,normal]
    outfile=make_outfile(infile)
    writer=vcfpy.Writer.from_path(outfile,outheader)
    #
    for record in myvcf:
        outrecord=copy(record)
        E={'CALLER':'VarDict'}
        E['OFS']=record.FILTER
        if not outrecord.INFO.get('STATUS','UNKNOWN') in ['LikelySomatic','Somatic']:
            E['SS']=outrecord.INFO.get('STATUS','UNKNOWN').upper()
        else:
            E['SS']='SOMATIC'
        outrecord.INFO.update(E)
        outrecord.FORMAT=check_FORMAT(outrecord.FORMAT)
        #check DP
        #check AD
        ref=outrecord.REF
        alt=outrecord.ALT
        info=outrecord.INFO
        for call in outrecord.calls:
            try:
                used_depth=call.data['DP']
                _aad=call.data.get('VD',0)#VD not AD for vardict
                _aaf=float(_aad)/float(used_depth)
                _rad=used_depth-_aad
            except (ZeroDivisionError,TypeError) as e:
                print(f'{e}: {ref} {alt[0].value} {call}')
                _aad=0
                _aaf=0.00
                _rad=0
            call.data['AAF']=float(f"{_aaf:.2f}")
            call.data['AD']=[_rad,_aad]
        writer.write_record(outrecord)
    writer.close()
    myvcf.close()

def std_varscan2(infile,tumor,normal):
    myvcf=vcfpy.Reader.from_path(infile)
    myvcf.header.samples.names=_sample_names(myvcf.header.samples.names,tumor,normal)
    myvcf.header.samples.name_to_idx={k:v for v,k in enumerate(myvcf.header.samples.names)}
    #
    outheader=modify_outheader(myvcf.header.copy())
    outheader.samples.names=[tumor,normal]
    outfile=make_outfile(infile)
    writer=vcfpy.Writer.from_path(outfile,outheader)
    #
    for record in myvcf:
        outrecord=copy(record)
        if not outrecord.INFO.get('SOMATIC',False):
            outrecord.INFO['SS']={'0':'REF','1':'GERMLINE','2':'SOMATIC','3':'LOH','5':'UNKNOWN'}[outrecord.INFO.get('SS','5')]
        else:
            outrecord.INFO['SS']='SOMATIC'
        E={'CALLER':'VarScan2'}
        E['OFS']=record.FILTER
        outrecord.INFO.update(E)
        outrecord.FORMAT=check_FORMAT(outrecord.FORMAT)
        ref=outrecord.REF
        alt=outrecord.ALT
        info=outrecord.INFO
        for call in outrecord.calls:
            try:
                used_depth=call.data['DP']
                _aad=call.data['AD']
                _aaf=float(_aad)/float(used_depth)
                _rad=used_depth-_aad
            except (ZeroDivisionError,TypeError) as e:
                print(f'{e}: {ref} {alt[0].value} {call}')
                _aad=0
                _aaf=0.00
                _rad=0
            call.data['AAF']=float(f"{_aaf:.2f}")
            call.data['AD']=[_rad,_aad]
        writer.write_record(outrecord)
    writer.close()
    myvcf.close()

def parse_arguments():
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input_fp',required=True,help='VCF input filepath')
    p.add_argument('-T','--tumor',required=True,help='Tumor ID')
    p.add_argument('-N','--normal',required=True,help='Normal ID')
    p.add_argument('--mode',choices=['mutect2','strelka2','vardict','varscan2'],required=True,help='Variant Caller used.')
    return vars(p.parse_args())

def main(argv=None):
    RUN={'mutect2':std_mutect2,'strelka2':std_strelka2,'vardict':std_vardict,'varscan2':std_varscan2}
    RUN[argv['mode']](argv['input_fp'],argv['tumor'],argv['normal'])

if __name__=='__main__':
    main(parse_arguments())

