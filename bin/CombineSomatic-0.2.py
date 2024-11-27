# Author: Arnab Ghosh, National Institute of Biomedical Genomics
import pandas as pd
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_args():
    import argparse
    parser = argparse.ArgumentParser(description = 'combine vcf')
    parser.add_argument("--config", help = 'configureation file', required = True)
    parser.add_argument("--sample", help = 'sample name', required = True)
    parser.add_argument("--variantCalls", help = 'VCF files. E.g. Mutect2:muect2.vcf+Strelka2:strelka2.vcf', required = True)
    parser.add_argument("--out", help = 'output TSV file', required = True)
    args = parser.parse_args()
    _dict = {}
    _dict['sample name'] = args.sample
    _dict['config file'] = args.config
    _dict['variant calls'] = {l.split(':')[0]:l.split(':')[1] for l in args.variantCalls.split('+')}
    _dict['output file'] = args.out
    return _dict

def getSomatic(__input_vcf, __var_caller):
    global _config
    global records
    if __var_caller in _config.get('variant callers').keys():            
        _metadata = []
        _header = []
        _data_dump = []
        with open(__input_vcf, 'r') as f:
            for _line in f:
                if _line.startswith('#') == True:
                    if _line.startswith('##') == True:
                        _metadata.append(_line.replace('\n', ''))
                    else:
                        _header = [k for k in _line.replace('\n', '').split('\t')]
                else:
                    _data_dump.append([k for k in _line.replace('\n', '').split('\t')])
        if _config.get('variant callers').get(__var_caller).get('sample name').get('type') == 'vcf embeded':
            __normal = [l for l in _metadata if l.startswith('##' + _config.get('variant callers').get(__var_caller).get('sample name').get('normal prefix'))][0].split('=')[1]
            __tumor = [l for l in _metadata if l.startswith('##' + _config.get('variant callers').get(__var_caller).get('sample name').get('tumor prefix'))][0].split('=')[1]
        elif _config.get('variant callers').get(__var_caller).get('sample name').get('type') == 'fixed':
            __normal = _config.get('variant callers').get(__var_caller).get('sample name').get('normal name')
            __tumor = _config.get('variant callers').get(__var_caller).get('sample name').get('tumor name')
        elif _config.get('variant callers').get(__var_caller).get('sample name').get('type') == 'ordered':
            __normal = [l for l in _header if l not in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']][int(_config.get('variant callers').get(__var_caller).get('sample name').get('normal order'))]
            __tumor = [l for l in _header if l not in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']][int(_config.get('variant callers').get(__var_caller).get('sample name').get('tumor order'))]
        _fields_FORMAT_raw = [l.split(',')[0].split('=')[2] for l in _metadata if l.startswith('##FORMAT')]
        _desired_FORMAT_fields = _config.get('variant callers').get(__var_caller).get('desired FORMAT')
        _fields_FORMAT = [l for l in _fields_FORMAT_raw if l in _desired_FORMAT_fields]
        _final_fields = _config.get('default fields')
        for _sample_type in ['normal', 'tumor']:
            for _element in _fields_FORMAT:
                _final_fields = _final_fields + [_sample_type + ':' + l for l in _config.get('variant callers').get(__var_caller).get('field formats').get(_element)]
        _res_dict = {}
        p = 0
        for _element in _data_dump:
            _dict_1 = {}
            for _field in _config.get('default fields'):
                _dict_1[_field] = _element[_header.index(_field)]
            for _sample_type in ['normal', 'tumor']:
                if _sample_type == 'normal':
                    _samp = __normal
                else:
                    _samp = __tumor
                for _field in _fields_FORMAT:
                    i = 0
                    for _sub_field in _config.get('variant callers').get(__var_caller).get('field formats').get(_field):
                        _dict_1[_sample_type + ':' + _sub_field] = _element[_header.index(_samp)].split(':')[_element[_header.index('FORMAT')].split(':').index(_field)].split(',')[i]
                        i = i + 1
            # oxo-g detection. only available for Mutect2 calls.
            if _config.get('oxog details').get('status') == 'report':
                _oxog_dict = _config.get('oxog details').get('motif')
                if __var_caller == 'Mutect2':
                    if (_dict_1.get('REF') == 'C') & (_dict_1.get('ALT') == 'A'):
                        if (int(_dict_1.get('tumor:alt_F2R1'))+int(_dict_1.get('tumor:alt_F1R2'))) > 0:
                            if (int(_dict_1.get('tumor:alt_F2R1'))/(int(_dict_1.get('tumor:alt_F2R1'))+int(_dict_1.get('tumor:alt_F1R2')))) >= _config.get('oxog details').get('foxog cutoff'):
                                if str(records[_dict_1.get('#CHROM')].seq)[int(_dict_1.get('POS'))-1:int(_dict_1.get('POS'))+2] in _oxog_dict.get('C>A'):
                                    _dict_1['oxoG'] = 'FAIL'
                                else:
                                    _dict_1['oxoG'] = 'PASS'
                            else:
                                _dict_1['oxoG'] = 'PASS'
                        else:
                            _dict_1['oxoG'] = 'NotDetermined'
                    elif (_dict_1.get('REF') == 'G') & (_dict_1.get('ALT') == 'T'):
                        if (int(_dict_1.get('tumor:alt_F2R1'))+int(_dict_1.get('tumor:alt_F1R2'))) > 0:
                            if (int(_dict_1.get('tumor:alt_F1R2'))/(int(_dict_1.get('tumor:alt_F1R2'))+int(_dict_1.get('tumor:alt_F2R1')))) >= _config.get('oxog details').get('foxog cutoff'):
                                if str(records[_dict_1.get('#CHROM')].seq)[int(_dict_1.get('POS'))-1:int(_dict_1.get('POS'))+2] in _oxog_dict.get('G>T'):
                                    _dict_1['oxoG'] = 'FAIL'
                                else:
                                    _dict_1['oxoG'] = 'PASS'
                            else:
                                _dict_1['oxoG'] = 'PASS'
                        else:
                            _dict_1['oxoG'] = 'NotDetermined'
                    else:
                        _dict_1['oxoG'] = 'NA'
            _res_dict[p] = _dict_1
            p = p + 1
        _res_df = pd.DataFrame(_res_dict).transpose()
        _mask_cols = _config.get("default fields constant")
        _rename_dict = {x: x + ':' + __var_caller for x in [l for l in _res_df.columns if l not in _mask_cols]}
        _res_df = _res_df.rename(columns = _rename_dict)
        return _res_df

if __name__ == '__main__':
    _master_dict = read_args()
    _config_file = _master_dict.get('config file')
    __sample = _master_dict.get('sample name')

    _config = {}
    with open(_config_file, 'r') as f:
        _config = json.load(f)

    if _config.get('oxog details').get('status') == 'report':
        records = SeqIO.to_dict(SeqIO.parse(open(_config.get('reference fasta')), 'fasta'))

    _primary_callers = [l for l in _master_dict.get('variant calls').keys() if _config.get('variant callers').get(l).get('set') == 'primary']
    _secondary_callers = [l for l in _master_dict.get('variant calls').keys() if _config.get('variant callers').get(l).get('set') == 'secondary']

    _columns_for_vcf_merge = _config.get("default fields constant")
    _res_df = pd.DataFrame()
    for _vc in _primary_callers:
        __input_vcf = _master_dict.get('variant calls').get(_vc)
        if len(_res_df) == 0:
            _res_df = getSomatic(__input_vcf, _vc)
        else:
            _res_df = _res_df.merge(getSomatic(__input_vcf, _vc), left_on = _columns_for_vcf_merge, right_on = _columns_for_vcf_merge, how = 'outer')
    for _vc in _secondary_callers:
        __input_vcf = _master_dict.get('variant calls').get(_vc)
        _res_df = _res_df.merge(getSomatic(__input_vcf, _vc), left_on = _columns_for_vcf_merge, right_on = _columns_for_vcf_merge, how = 'left')
    _res_df.insert(loc = 0, column = 'SAMPLE', value = __sample)
    _res_df.to_csv(_master_dict.get('output file'), sep = '\t', header = True, index = False, na_rep = 'NA')
