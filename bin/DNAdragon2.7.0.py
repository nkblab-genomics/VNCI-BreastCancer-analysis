# DNAdragon version 2.7.0 by Arnab Ghosh and Nidhan K Biswas, National institute of Biomedical Genomics
# The code is designed to run in SGE cluster, change "all.q" in line 245 with appropriate queue name
# See help by python DNAdragon2.7.0.py --help
# reading sample sheet
_req_py_modules = ['os', 'pandas', 'snakemake']

def chk_req_py_modules(_list_modules):
    _chk_var = True
    for _module in _list_modules:
        try:
            _dummy_var = __import__(_module)
        except:
            _chk_var = False
            print('can\'t import required python module -->',_module)
    return(_chk_var)

if chk_req_py_modules(_req_py_modules) == False:
    print('try installing required modules --> ', _req_py_modules)
    exit()

import os
import pandas as pd

def read_args():
    import argparse
    parser = argparse.ArgumentParser(description="one stop Exome and WGS analysis pipeline")
    parser.add_argument("--conf", help="tool config file", required=True)
    parser.add_argument("--info", help="sample sheet file", required=True)
    parser.add_argument("--resource", help="refseq, 1000Genome, dbsnp data", required=True)
    parser.add_argument("--workDir", help="working directory", required=True)
    parser.add_argument("--sgeLogDir", help="sge stdout/stderr dump directory", required=True)
    parser.add_argument("--maxJobPerNode", help="max job to be submitted per compute node", required=True)
    parser.add_argument("--keepAllFiles", help="whether to keep intermediate files (default: F)", choices=['T', 'F'], default='F')
    parser.add_argument("--dataType", help="Specify data type", choices=['Ex', 'Wg', 'Amp'], required=True)
    parser.add_argument("--gvcf", help="GVCF calling (default: F)", choices=['T', 'F'], default='F')
    parser.add_argument("--fixedSeqLen", help="analysis with fixed sequence length", default=0)
    args = parser.parse_args()
    global __tool_config
    global __samplesheet
    global __res_config
    global __work_dir
    global __sge_log
    global __master_node_max_job
    global __keep_all_files
    global __data_type
    global __gvcf_call
    global __fix_seq_len
    __tool_config = args.conf
    __samplesheet = args.info
    __res_config = args.resource
    __work_dir = args.workDir
    __sge_log = args.sgeLogDir
    __master_node_max_job = str(args.maxJobPerNode)
    __keep_all_files = args.keepAllFiles
    __data_type = args.dataType
    __gvcf_call = args.gvcf
    __fix_seq_len = args.fixedSeqLen

read_args()

def chk_samplesheet(samplesheet_file):
    supported_pl = ['illumina',]
    import pandas as pd
    import os
    _chk_var = True
    __samplesheet = pd.read_table(samplesheet_file, sep = '\t', header = 0, dtype=str)
    _header_elements = __samplesheet.columns.values.tolist()
    _index_elements = __samplesheet.index.values.tolist()
    arr_key = []
    _header_requirements = ['name', 'pl', 'fcl', 'lane', 'read1', 'read2']
    for header_elm in _header_requirements:
        if header_elm not in _header_elements:
            print(header_elm, 'not found in samplesheet header')
            _chk_var = False
    if _chk_var == True:
        for indx_elm in _index_elements:
            key_str = __samplesheet.loc[indx_elm, 'name'] + __samplesheet.loc[indx_elm, 'fcl'] + __samplesheet.loc[indx_elm, 'lane']
            if key_str not in arr_key:
                arr_key.append(key_str)
            else:
                print('duplicated samples present in samplesheet')
                _chk_var = False
            if __samplesheet.loc[indx_elm, 'pl'].lower() not in supported_pl:
                print('not supported ->', __samplesheet.loc[indx_elm, 'pl'])
                _chk_var = False
            if os.path.isfile(__samplesheet.loc[indx_elm, 'read1']) == True:
                if os.path.isfile(__samplesheet.loc[indx_elm, 'read2']) == False:
                    print('cant open ->',__samplesheet.loc[indx_elm, 'read2'])
                    _chk_var = False
            else:
                print('cant open ->',__samplesheet.loc[indx_elm, 'read1'])
                _chk_var = False
    return(_chk_var)

def chk_config_file(__tool_config):
    _chk_var = True
    import pandas as pd
    import os
    required_steps = ['bam_qc', 'bamfilter', 'merge', 'fastq_qual_filter', 'adpremove', 'fastq_qc', 'alignment', 'samsort', 'removedup', 'lirRealignerTargetCreator', 'lirIndelRealigner', 'bqsr_BaseRecalibrator', 'bqsr_PrintReads']
    if __gvcf_call == 'T':
        required_steps.append('gvcfcalling')
    required_fields = ['step', 'toolpath', 'threads', 'memory']
    __configfile = pd.read_table(__tool_config, sep = '\t', header = 0)
    _header_elements = __configfile.columns.values.tolist()
    _index_elements = __configfile.index.values.tolist()
    for elm in required_fields:
        if elm not in _header_elements:
            print('missing',elm,'from configuration file header')
            _chk_var = False
    if _chk_var == True:
        _step_elements = __configfile['step'].tolist()
        for elm in required_steps:
            if elm not in _step_elements:
                print('missing',elm,'from configuration file > \'step\'')
                _chk_var = False
    if _chk_var == True:
        for elm in _index_elements:
            if os.path.isfile(__configfile.loc[elm, 'toolpath']) == False:
                _chk_var = False
                print('can\'nt open ->', __configfile.loc[elm, 'toolpath'])
            if 'G' not in __configfile.loc[elm, 'memory']:
                _chk_var = False
                print('expected <int>G in memory ->', __configfile.loc[elm, 'step'])
            else:
                try:
                    _temp_var1 = int(__configfile.loc[elm, 'memory'].replace('G', ''))
                    _temp_var2 = float(__configfile.loc[elm, 'memory'].replace('G', ''))
                    if _temp_var1 != _temp_var2:
                        _chk_var = False
                        print('expected <int>G in memory ->', __configfile.loc[elm, 'step'])
                    elif _temp_var1 == 0:
                        _chk_var = False
                        print('expected at least 1G in memory ->', __configfile.loc[elm, 'step'])
                except:
                    _chk_var = False
                    print('expected <int>G in memory ->', __configfile.loc[elm, 'step'])
            try:
                _temp_var1 = int(__configfile.loc[elm, 'threads'])
                _temp_var2 = float(__configfile.loc[elm, 'threads'])
                if _temp_var1 != _temp_var2:
                    _chk_var = False
                    print('expected <int> in threads ->', __configfile.loc[elm, 'step'])
                elif _temp_var1 == 0:
                    _chk_var = False
                    print('expected at least 1 in threads ->', __configfile.loc[elm, 'step'])
            except:
                _chk_var = False
                print('expected <int> in threads ->', __configfile.loc[elm, 'step'])
    return(_chk_var)

def chk_resFile(resources_file, data_type):
    _chk_var = True
    required_fields = ['resources', 'path']
    if data_type == 'Wg':
        required_res = ['adapter_file', 'refseq_hg19_decoy', 'gb_1000G_indel', 'gb_mills_1000G_gold_indel', 'gb_dbsnp']
    else:
        required_res = ['adapter_file', 'refseq_hg19_decoy', 'gb_1000G_indel', 'gb_mills_1000G_gold_indel', 'gb_dbsnp', 'capture_bed']
    import pandas as pd
    import os
    __resfile = pd.read_table(resources_file, sep = '\t', header = 0)
    _header_elements = __resfile.columns.values.tolist()
    _index_elements = __resfile.index.values.tolist()
    for elm in required_fields:
        if elm not in _header_elements:
            _chk_var = False
            print('missing',elm,'from resources file')
    if _chk_var == True:
        _res_elements = __resfile['resources'].tolist()
        for elm in required_res:
            if elm not in _res_elements:
                _chk_var = False
                print('missing',elm,'from resources')
        for elm in _index_elements:
            if os.path.isfile(__resfile.loc[elm, 'path']) == False:
                _chk_var = False
                print('can\'nt open -->',__resfile.loc[elm, 'path'])
    return(_chk_var)

# exit if chk_samplesheet/chk_config_file/chk_resFile value is "False"
if chk_samplesheet(__samplesheet) == False:
    print('error found in samplesheet')
    exit()

if chk_config_file(__tool_config) == False:
    print('error found in config file')
    exit()

if chk_resFile(__res_config, __data_type) == False:
    print('error found in resources file')
    exit()

#exit()

my_db = pd.read_csv(__samplesheet, sep = '\t', dtype=str)
my_db['key'] = my_db['name'] + '_' + my_db['fcl'] + '_' + my_db['lane']
_sample_dict = my_db.set_index('key').transpose().to_dict()

_individuals = []
for keys in _sample_dict.keys():
    if _sample_dict.get(keys).get('name') not in _individuals:
        _individuals.append(_sample_dict.get(keys).get('name'))

_master_dict = {}
for indv in _individuals:
    scratch = []
    for keys in _sample_dict.keys():
        if _sample_dict.get(keys).get('name') == indv:
            scratch.append(keys)
    _master_dict[indv] = scratch

def write_master_workflow(_sample_array,_max_process,_work_dir,_sge_log_dir):
    _snake_file_name = _work_dir + '/Snakefile' 
    comp1 = "SAMPLES=["
    for l in _sample_array:
        comp1 = comp1 + "'" + l + "',"
    comp1 = comp1 + "]"

    with open(_snake_file_name, 'w') as snakef:
        snakef.write(comp1)
        snakef.write('\n')
        snakef.write('localrules: target, fire_job\n')
        snakef.write('rule target:\n')
        snakef.write('    input:\n')
        snakef.write('        expand("{sample}/{sample}_fq2bam_completed.flag", sample=SAMPLES)\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        echo "all done :)"\n')
        snakef.write('        \'\'\'\n')
        snakef.write('rule fire_job:\n')
        snakef.write('    input:\n')
        snakef.write('        snake_file="{sample}/Snakefile",\n')
        snakef.write('        run_dir=directory("{sample}")\n')
        snakef.write('    output:\n')
        snakef.write('        "{sample}/{sample}_fq2bam_completed.flag"\n')
        snakef.write('    run:\n')
        snakef.write('        import os\n')
        if __keep_all_files == 'T':
            snakef.write('        command = \'snakemake --rerun-incomplete --restart-times 3 --max-jobs-per-second 1 --notemp --latency-wait 3600 -s \' + str(input.snake_file) + \' -j ')
        else:
            snakef.write('        command = \'snakemake --rerun-incomplete --restart-times 3 --max-jobs-per-second 1 --latency-wait 3600 -s \' + str(input.snake_file) + \' -j ')
        snakef.write(_max_process)
        snakef.write(' -d \' + str(input.run_dir) + str(\' --cluster "qsub -pe smp {threads} -R y -l mem_free={params.qsub_mem} -e ')
        snakef.write(_sge_log_dir)
        snakef.write(' -o ')
        snakef.write(_sge_log_dir)
        snakef.write(' -q all.q"\')\n')
        snakef.write('        print(command)\n')
        snakef.write('        os.system(command)\n')
    snakef.close()

write_master_workflow(_individuals,__master_node_max_job,__work_dir,__sge_log)

def env_creation(_work_dir,_sample_array,_sample_dict):
    import os
    for l in _sample_array:
        _sample_dir_command = 'mkdir -p ' + _work_dir + '/' + l
        _sub_input_dir_command = 'mkdir -p ' + _work_dir + '/' + l + '/.input'
        os.system(_sample_dir_command)
        os.system(_sub_input_dir_command)
        for sub_sample in _master_dict.get(l):
            _input_r1_link_gen_command = 'ln -s -f ' + _sample_dict.get(sub_sample).get('read1') + ' ' + _work_dir + '/' + l + '/.input/' + sub_sample + '_R1.fastq.gz'
            _input_r2_link_gen_command = 'ln -s -f ' + _sample_dict.get(sub_sample).get('read2') + ' ' + _work_dir + '/' + l + '/.input/' + sub_sample + '_R2.fastq.gz'
            os.system(_input_r1_link_gen_command)
            os.system(_input_r2_link_gen_command)

env_creation(__work_dir,_individuals,_sample_dict)

def write_sub_workflow_sampleDef(_individual,_master_dict,_work_dir):
    _sample_dir = _work_dir + '/' + _individual
    _sample_snakefile = _work_dir + '/' + _individual + '/' + 'Snakefile'
    with open(_sample_snakefile, 'w') as snakef:
        comp1 = 'SAMPLES = ['
        for l in _master_dict.get(_individual):
            comp1 = comp1 + "'" + l + "',"
        comp1 = comp1 + "]"
        snakef.write(comp1)
        snakef.write('\n')
        snakef.write('import pandas as pd')
        snakef.write('\n')
        _temp_str = '__samplesheet = \'' + __samplesheet + '\''
        snakef.write(_temp_str)
        snakef.write('\n')
        snakef.write('my_db = pd.read_csv(__samplesheet, sep = \'\\t\', dtype=str)')
        snakef.write('\n')
        snakef.write('my_db[\'key\'] = my_db[\'name\'] + \'_\' + my_db[\'fcl\'] + \'_\' + my_db[\'lane\']')
        snakef.write('\n')
        snakef.write('myDict = my_db.set_index(\'key\').transpose().to_dict()')
        snakef.write('\n')
        _temp_str = 'sampleID = \'' + _individual + '\''
        snakef.write(_temp_str)
        snakef.write('\n')
        _temp_str = '__tool_config = \'' + __tool_config + '\''
        snakef.write(_temp_str)
        snakef.write('\n')
        snakef.write('tool_db = pd.read_csv(__tool_config, sep = \'\\t\', header = 0, index_col = 0)')
        snakef.write('\n')
        snakef.write('tool_db[\'qsub_memory\'] = ((tool_db[\'memory\'].str.replace(\'G\',\'\').astype(int) * 1024).astype(int) + 2048).astype(str) + \'M\'')
        snakef.write('\n')
        _temp_str = '__res_config = \'' + __res_config + '\''
        snakef.write(_temp_str)
        snakef.write('\n')
        snakef.write('res_db = pd.read_csv(__res_config, sep = \'\\t\', header = 0, index_col = 0)')
        snakef.write('\n')
        snakef.write('import gzip')
        snakef.write('\n')
        _temp_str = '_dummy_samp = \'' + __work_dir + '/' + _individual + '/.input/\''  + ' + SAMPLES[0] + ' + '\'_R1.fastq.gz\''
        snakef.write(_temp_str)
        snakef.write('\n')
        #_temp_str = '_temp_file_ = gzip.open(_dummy_samp, \'rb\')'
        #snakef.write(_temp_str)
        #snakef.write('\n')
        #snakef.write('_buffer_line_ = _temp_file_.readline()')
        #snakef.write('\n')
        #snakef.write('_buffer_line_ = _temp_file_.readline()')
        #snakef.write('\n')
        #snakef.write('_seq_length = len([l for l in _buffer_line_.decode("utf-8").replace(\'\\n\', \'\')])')
        #snakef.write('\n')
        #snakef.write('_temp_file_.close()')
        #snakef.write('\n')
    snakef.close()

def write_sub_workflow_main_body(_individual,_work_dir):
    _sample_dir = _work_dir + '/' + _individual
    _sample_snakefile = _work_dir + '/' + _individual + '/' + 'Snakefile'
    with open(_sample_snakefile, 'a') as snakef:
        snakef.write('\n')
        snakef.write('localrules: target\n')
        snakef.write('\n')
        snakef.write('rule target:\n')
        snakef.write('    input:\n')
        snakef.write('        ancient(str(sampleID) + \'_merge.bam\'),\n')
        snakef.write('        ancient(str(sampleID) + \'_mergeFilter.bam\'),\n')
        snakef.write('        ancient(expand("FastQC_REPORT_raw/{sample}_R1_fastqc.zip", sample=SAMPLES)),\n')
        snakef.write('        ancient(expand("FastQC_REPORT/{sample}_R1_filter_adprm_fastqc.zip", sample=SAMPLES)),\n')
        snakef.write('        ancient("Qualimap_REPORT/genome_results.txt"),\n')
        if __gvcf_call == 'T':
            snakef.write('        ancient(str(sampleID) + \'_HaplotypeCaller.g.vcf.gz\'),\n')
        snakef.write('    output:\n')
        snakef.write('        flagfile=str(sampleID) + \'_fq2bam_completed.flag\'\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        touch {output.flagfile}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule bam_qc:\n')
        snakef.write('    input:\n')
        snakef.write('        bam=ancient(str(sampleID) + \'_mergeFilter.bam\'),\n')
        snakef.write('        bai=ancient(str(sampleID) + \'_mergeFilter.bam\'),\n')
        if __data_type != 'Wg':
            snakef.write('        bed=res_db.loc[\'capture_bed\',\'path\']\n')
        snakef.write('    output:\n')
        snakef.write('        bam_qc_data="Qualimap_REPORT/genome_results.txt"\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'bam_qc\',\'toolpath\'],\n')
        snakef.write('        bam_qc_dir=\'Qualimap_REPORT\',\n')
        snakef.write('        qsub_mem=tool_db.loc[\'bam_qc\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'bam_qc\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'bam_qc\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        if __data_type != 'Wg':
            snakef.write('        {params.tool} bamqc -bam {input.bam} -outdir {params.bam_qc_dir} --java-mem-size={params.runtime_mem} -gff {input.bed} -nt {threads}\n')
        else:
            snakef.write('        {params.tool} bamqc -bam {input.bam} -outdir {params.bam_qc_dir} --java-mem-size={params.runtime_mem} -nt {threads}\n')
        snakef.write('        \'\'\'\n')
        if __gvcf_call == 'T':
            snakef.write('rule gvcfcalling:\n')
            snakef.write('    input:\n')
            snakef.write('        bam=ancient(str(sampleID) + \'_mergeFilter.bam\'),\n')
            snakef.write('        ref=res_db.loc[\'refseq_hg19_decoy\',\'path\'],\n')
            if __data_type != 'Wg':
                snakef.write('        bed=res_db.loc[\'capture_bed\',\'path\']\n')
            snakef.write('    output:\n')
            snakef.write('        gvcf=str(sampleID) + \'_HaplotypeCaller.g.vcf.gz\',\n')
            snakef.write('        tmpdir=temporary(directory("gvcf_tmpdir")),\n')
            snakef.write('    params:\n')
            snakef.write('        tool=tool_db.loc[\'gvcfcalling\',\'toolpath\'],\n')
            snakef.write('        qsub_mem=tool_db.loc[\'gvcfcalling\',\'qsub_memory\'],\n')
            snakef.write('        runtime_mem=tool_db.loc[\'gvcfcalling\',\'memory\']\n')
            snakef.write('    threads: int(tool_db.loc[\'bamfilter\',\'threads\'])\n')
            snakef.write('    shell:\n')
            snakef.write('        \'\'\'\n')
            snakef.write('        mkdir -p {output.tmpdir}\n')
            if __data_type != 'Wg':
                snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar {params.tool} HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.gvcf} -ERC GVCF --native-pair-hmm-threads {threads} --tmp-dir {output.tmpdir} -L {input.bed}\n')
            else:
                snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar {params.tool} HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.gvcf} -ERC GVCF --native-pair-hmm-threads {threads} --tmp-dir {output.tmpdir} \n')
            snakef.write('        \'\'\'\n')
        snakef.write('rule bamfilter:\n')
        snakef.write('    input:\n')
        snakef.write('        bam=ancient(str(sampleID) + \'_merge.bam\')\n')
        snakef.write('    output:\n')
        snakef.write('        bam=str(sampleID) + \'_mergeFilter.bam\'\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'bamfilter\',\'toolpath\'],\n')
        snakef.write('        chrom=\'chr{1..22} chrX chrY\',\n')
        snakef.write('        filter=\'-f 0x2 -F 0x100 -q 40\',\n')
        snakef.write('        qsub_mem=tool_db.loc[\'bamfilter\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'bamfilter\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'bamfilter\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        {params.tool} view {params.filter} {input.bam} {params.chrom} -@ {threads} -o {output.bam}\n')
        snakef.write('        {params.tool} index {output.bam}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('rule merge:\n')
        snakef.write('    input:\n')
        snakef.write('        ancient(expand("{sample}_mem_sorted_rmdup_realign_recal.bam", sample=SAMPLES)),\n')
        snakef.write('        ancient(expand("{sample}_mem_sorted_rmdup_realign_recal.bai", sample=SAMPLES))\n')
        snakef.write('    output:\n')
        snakef.write('        bam=str(sampleID) + \'_merge.bam\'\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'merge\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'merge\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'merge\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'merge\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        ls *_recal.bam > mergeBam.list\n')
        snakef.write('        echo -e \'java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar {params.tool} MergeSamFiles \\\\\' > mergeBam.sh\n')
        snakef.write('        while read listBamForMerge; do\n')
        snakef.write('            echo \'I=\'$listBamForMerge\' \\\\\' >> mergeBam.sh\n')
        snakef.write('        done < mergeBam.list\n')
        snakef.write('        echo \'O={output.bam} \\\\\' >> mergeBam.sh\n')
        snakef.write('        echo \'CREATE_INDEX=true \\\\\' >> mergeBam.sh\n')
        snakef.write('        echo \'ASSUME_SORTED=true \\\\\' >> mergeBam.sh\n')
        snakef.write('        echo \'VALIDATION_STRINGENCY=LENIENT\' >> mergeBam.sh\n')
        snakef.write('        chmod a+x mergeBam.sh\n')
        snakef.write('        ./mergeBam.sh\n')
        snakef.write('        rm -f mergeBam.list mergeBam.sh\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule fastq_qc_raw:\n')
        snakef.write('    input:\n')
        snakef.write('        read1=ancient(".input/{sample}_R1.fastq.gz"),\n')
        snakef.write('        read2=ancient(".input/{sample}_R2.fastq.gz")\n')
        snakef.write('    output:\n')
        snakef.write('        zip_file="FastQC_REPORT_raw/{sample}_R1_fastqc.zip",\n')
        snakef.write('        html_file="FastQC_REPORT_raw/{sample}_R1_fastqc.html",\n')
        snakef.write('        zip_file2="FastQC_REPORT_raw/{sample}_R2_fastqc.zip",\n')
        snakef.write('        html_file2="FastQC_REPORT_raw/{sample}_R2_fastqc.html",\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'fastq_qc\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'fastq_qc\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'fastq_qc\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'fastq_qc\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        mkdir -p FastQC_REPORT_raw/\n')
        snakef.write('        {params.tool} -o FastQC_REPORT_raw -f fastq {input.read1} {input.read2} -t {threads}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule fastq_qual_filter:\n')
        snakef.write('    input:\n')
        snakef.write('        read1=ancient(".input/{sample}_R1.fastq.gz"),\n')
        snakef.write('        read2=ancient(".input/{sample}_R2.fastq.gz")\n')
        snakef.write('    output:\n')
        snakef.write('        read1f=temporary("{sample}_R1_filter.fastq.gz"),\n')
        snakef.write('        read2f=temporary("{sample}_R2_filter.fastq.gz"),\n')
        snakef.write('        lowqual=temporary("{sample}_lowQ.fastq.gz")\n')
        snakef.write('    threads: int(tool_db.loc[\'fastq_qual_filter\',\'threads\'])\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'fastq_qual_filter\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'fastq_qual_filter\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'fastq_qual_filter\',\'memory\'],\n')
        snakef.write('        q_val=20,\n')
        snakef.write('        pct_hiqual_base=90,\n')
        snakef.write('        max_n_bases=5\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        {params.tool} {input.read1} {input.read2} {output.read1f} {output.read2f} {output.lowqual} {params.q_val} {params.pct_hiqual_base} {params.max_n_bases}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule adpremove:\n')
        snakef.write('    input:\n')
        snakef.write('        read1="{sample}_R1_filter.fastq.gz",\n')
        snakef.write('        read2="{sample}_R2_filter.fastq.gz",\n')
        snakef.write('        adapter=res_db.loc[\'adapter_file\',\'path\']\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'adpremove\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'adpremove\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'adpremove\',\'memory\'],\n')
        #snakef.write('        seq_length=_seq_length,\n')
        snakef.write('        ext=\'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15\'\n')
        snakef.write('    output:\n')
        snakef.write('        read1p=temporary("{sample}_R1_filter_adprm.fastq.gz"),\n')
        snakef.write('        read2p=temporary("{sample}_R2_filter_adprm.fastq.gz"),\n')
        snakef.write('        read1up=temporary("{sample}_R1_filter_unpair.fastq.gz"),\n')
        snakef.write('        read2up=temporary("{sample}_R2_filter_unpair.fastq.gz")\n')
        snakef.write('    threads: int(tool_db.loc[\'adpremove\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        #snakef.write('        SEQ_LEN=$(gunzip -c {input.read1} | head -2 | sed 1d | fold -w 1 | wc -l)\n')
        if __fix_seq_len == 0:
            snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar {params.tool} PE -threads {threads} -phred33 {input.read1} {input.read2} {output.read1p} {output.read1up} {output.read2p} {output.read2up} ILLUMINACLIP:{input.adapter}:2:30:10 {params.ext} MINLEN:"$(gunzip -c {input.read1} | head -2 | sed 1d | fold -w 1 | wc -l)"\n')
        else:
            snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar {params.tool} PE -threads {threads} -phred33 {input.read1} {input.read2} {output.read1p} {output.read1up} {output.read2p} {output.read2up} ILLUMINACLIP:{input.adapter}:2:30:10 {params.ext} MINLEN:')
            snakef.write(__fix_seq_len)
            snakef.write('\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule unzip:\n')
        snakef.write('    input:\n')
        snakef.write('        read1=ancient("{sample}_R1_filter_adprm.fastq.gz"),\n')
        snakef.write('        read2=ancient("{sample}_R2_filter_adprm.fastq.gz")\n')
        snakef.write('    output:\n')
        snakef.write('        read1=temporary("{sample}_R1_filter_adprm.fastq"),\n')
        snakef.write('        read2=temporary("{sample}_R2_filter_adprm.fastq")\n')
        snakef.write('    params:\n')
        snakef.write('        qsub_mem=\'500M\',\n')
        snakef.write('        runtime_mem=\'1G\'\n')
        snakef.write('    threads: 1\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        gunzip -c {input.read1} > {output.read1}\n')
        snakef.write('        gunzip -c {input.read2} > {output.read2}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule fastq_qc:\n')
        snakef.write('    input:\n')
        snakef.write('        read1=ancient("{sample}_R1_filter_adprm.fastq.gz"),\n')
        snakef.write('        read2=ancient("{sample}_R2_filter_adprm.fastq.gz")\n')
        snakef.write('    output:\n')
        snakef.write('        zip_file="FastQC_REPORT/{sample}_R1_filter_adprm_fastqc.zip",\n')
        snakef.write('        html_file="FastQC_REPORT/{sample}_R1_filter_adprm_fastqc.html"\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'fastq_qc\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'fastq_qc\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'fastq_qc\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'fastq_qc\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        mkdir -p FastQC_REPORT/\n')
        snakef.write('        {params.tool} -o FastQC_REPORT -f fastq {input.read1} {input.read2} -t {threads}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule alignment:\n')
        snakef.write('    input:\n')
        snakef.write('        read1=ancient("{sample}_R1_filter_adprm.fastq.gz"),\n')
        snakef.write('        read2=ancient("{sample}_R2_filter_adprm.fastq.gz"),\n')
        snakef.write('        ref=res_db.loc[\'refseq_hg19_decoy\',\'path\']\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'alignment\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'alignment\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'alignment\',\'memory\']\n')
        snakef.write('    output:\n')
        snakef.write('        temporary("{sample}_mem.sam")\n')
        snakef.write('    threads: int(tool_db.loc[\'alignment\',\'threads\'])\n')
        snakef.write('    run:\n')
        snakef.write('        import os\n')
        snakef.write('        id = str(input.read1).replace(\'_R1_filter_adprm.fastq.gz\', \'\')\n')
        snakef.write('        pl = myDict.get(id).get(\'pl\')\n')
        snakef.write('        pu = myDict.get(id).get(\'fcl\') + "." + myDict.get(id).get(\'lane\')\n')
        snakef.write('        lb = str(input.read1).replace(\'_R1_filter_adprm.fastq.gz\', \'\').split(\'_\')[0]\n')
        snakef.write('        sm = "Sample_" + lb\n')
        snakef.write('        rgtag = "-R \'@RG\\\\tID:" + id + "\\\\tPL:" + pl + "\\\\tPU:" + pu + "\\\\tLB:" + lb + "\\\\tSM:" + sm + "\'"\n')
        snakef.write('        command = str(params.tool) + " mem -t " + str(threads) + " -T 1 -M " + rgtag + " " + str(input.ref) + " \'<zcat " + str(input.read1) + "\' \'<zcat " + str(input.read2) + "\' > " + str(output)\n')
        snakef.write('        print(command)\n')
        snakef.write('        os.system(command)\n')
        snakef.write('\n')
        snakef.write('rule samsort:\n')
        snakef.write('    input:\n')
        snakef.write('        sam=ancient("{sample}_mem.sam")\n')
        snakef.write('    output:\n')
        snakef.write('        bam=temporary("{sample}_mem_sorted.bam"),\n')
        snakef.write('        bai=temporary("{sample}_mem_sorted.bai"),\n')
        snakef.write('        tmpdir=temporary(directory("{sample}_sort"))\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'samsort\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'samsort\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'samsort\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'samsort\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar {params.tool} SortSam INPUT={input.sam} OUTPUT={output.bam} SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR={output.tmpdir}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        if __data_type != 'Amp':
            snakef.write('rule removedup:\n')
            snakef.write('    input:\n')
            snakef.write('        bam=ancient("{sample}_mem_sorted.bam"),\n')
            snakef.write('        bai=ancient("{sample}_mem_sorted.bai")\n')
            snakef.write('    output:\n')
            snakef.write('        bam=temporary("{sample}_mem_sorted_rmdup.bam"),\n')
            snakef.write('        bai=temporary("{sample}_mem_sorted_rmdup.bai"),\n')
            snakef.write('        tmpdir=temporary(directory("{sample}_rmdup")),\n')
            snakef.write('        matrix_file=temporary("{sample}_mem_sorted_rmdup_matrix.txt")\n')
            snakef.write('    params:\n')
            snakef.write('        tool=tool_db.loc[\'removedup\',\'toolpath\'],\n')
            snakef.write('        qsub_mem=tool_db.loc[\'removedup\',\'qsub_memory\'],\n')
            snakef.write('        runtime_mem=tool_db.loc[\'removedup\',\'memory\']\n')
            snakef.write('    threads: int(tool_db.loc[\'removedup\',\'threads\'])\n')
            snakef.write('    shell:\n')
            snakef.write('        \'\'\'\n')
            snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar {params.tool} MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} REMOVE_DUPLICATES=true METRICS_FILE={output.matrix_file} VALIDATION_STRINGENCY=LENIENT TMP_DIR={output.tmpdir} CREATE_INDEX=true\n')
            snakef.write('        \'\'\'\n')
            snakef.write('\n')
        snakef.write('rule lirRealignerTargetCreator:\n')
        snakef.write('    input:\n')
        if __data_type == 'Amp':
            snakef.write('        bam=ancient("{sample}_mem_sorted.bam"),\n')
            snakef.write('        bai=ancient("{sample}_mem_sorted.bai"),\n')
        else:
            snakef.write('        bam=ancient("{sample}_mem_sorted_rmdup.bam"),\n')
            snakef.write('        bai=ancient("{sample}_mem_sorted_rmdup.bai"),\n')
        snakef.write('        ref=res_db.loc[\'refseq_hg19_decoy\',\'path\'],\n')
        snakef.write('        gb_1000G_indel=res_db.loc[\'gb_1000G_indel\',\'path\'],\n')
        snakef.write('        gb_mills_1000G_gold_indel=res_db.loc[\'gb_mills_1000G_gold_indel\',\'path\']\n')
        snakef.write('    output:\n')
        snakef.write('        interval_file=temporary("{sample}_mem_sorted_rmdup_forIndelRealigner.intervals"),\n')
        snakef.write('        tmpdir=temporary(directory("{sample}_lir"))\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'lirRealignerTargetCreator\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'lirRealignerTargetCreator\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'lirRealignerTargetCreator\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'lirRealignerTargetCreator\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar -Djava.io.tmpdir={output.tmpdir} {params.tool} -I {input.bam} -R {input.ref} -T RealignerTargetCreator -o {output.interval_file} -known {input.gb_1000G_indel} -known {input.gb_mills_1000G_gold_indel} -nt {threads}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule lirIndelRealigner:\n')
        snakef.write('    input:\n')
        if __data_type == 'Amp':
            snakef.write('        bam=ancient("{sample}_mem_sorted.bam"),\n')
            snakef.write('        bai=ancient("{sample}_mem_sorted.bai"),\n')
        else:
            snakef.write('        bam=ancient("{sample}_mem_sorted_rmdup.bam"),\n')
            snakef.write('        bai=ancient("{sample}_mem_sorted_rmdup.bai"),\n')
        snakef.write('        interval_file="{sample}_mem_sorted_rmdup_forIndelRealigner.intervals",\n')
        snakef.write('        ref=res_db.loc[\'refseq_hg19_decoy\',\'path\'],\n')
        snakef.write('        gb_1000G_indel=res_db.loc[\'gb_1000G_indel\',\'path\'],\n')
        snakef.write('        gb_mills_1000G_gold_indel=res_db.loc[\'gb_mills_1000G_gold_indel\',\'path\']\n')
        snakef.write('    output:\n')
        snakef.write('        bam=temporary("{sample}_mem_sorted_rmdup_realign.bam"),\n')
        snakef.write('        bai=temporary("{sample}_mem_sorted_rmdup_realign.bai"),\n')
        snakef.write('        tmpdir=temporary(directory("{sample}_lir"))\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'lirIndelRealigner\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'lirIndelRealigner\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'lirIndelRealigner\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'lirIndelRealigner\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -Djava.io.tmpdir={output.tmpdir} -jar {params.tool} -I {input.bam} -R {input.ref} -T IndelRealigner -targetIntervals {input.interval_file} -o {output.bam} -known {input.gb_1000G_indel} -known {input.gb_mills_1000G_gold_indel}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('rule bqsr_BaseRecalibrator:\n')
        snakef.write('    input:\n')
        snakef.write('        bam=ancient("{sample}_mem_sorted_rmdup_realign.bam"),\n')
        snakef.write('        bai=ancient("{sample}_mem_sorted_rmdup_realign.bai"),\n')
        snakef.write('        ref=res_db.loc[\'refseq_hg19_decoy\',\'path\'],\n')
        snakef.write('        gb_dbsnp=res_db.loc[\'gb_dbsnp\',\'path\']\n')
        snakef.write('    output:\n')
        snakef.write('        pre_recal_table=temporary("{sample}_mem_sorted_rmdup_realign_pre.recal.table"),\n')
        snakef.write('        tmpdir=temporary(directory("{sample}_bqsr"))\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'bqsr_BaseRecalibrator\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'bqsr_BaseRecalibrator\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'bqsr_BaseRecalibrator\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'bqsr_BaseRecalibrator\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -jar -Djava.io.tmpdir={output.tmpdir} {params.tool} -R {input.ref} -knownSites {input.gb_dbsnp} -I {input.bam} -T BaseRecalibrator -o {output.pre_recal_table}\n')
        snakef.write('        \'\'\'\n')
        snakef.write('\n')
        snakef.write('rule bqsr_PrintReads:\n')
        snakef.write('    input:\n')
        snakef.write('        bam=ancient("{sample}_mem_sorted_rmdup_realign.bam"),\n')
        snakef.write('        bai=ancient("{sample}_mem_sorted_rmdup_realign.bai"),\n')
        snakef.write('        pre_recal_table="{sample}_mem_sorted_rmdup_realign_pre.recal.table",\n')
        snakef.write('        ref=res_db.loc[\'refseq_hg19_decoy\',\'path\']\n')
        snakef.write('    output:\n')
        snakef.write('        bam=temporary("{sample}_mem_sorted_rmdup_realign_recal.bam"),\n')
        snakef.write('        bai=temporary("{sample}_mem_sorted_rmdup_realign_recal.bai"),\n')
        snakef.write('        tmpdir=temporary(directory("{sample}_bqsr"))\n')
        snakef.write('    params:\n')
        snakef.write('        tool=tool_db.loc[\'bqsr_PrintReads\',\'toolpath\'],\n')
        snakef.write('        qsub_mem=tool_db.loc[\'bqsr_PrintReads\',\'qsub_memory\'],\n')
        snakef.write('        runtime_mem=tool_db.loc[\'bqsr_PrintReads\',\'memory\']\n')
        snakef.write('    threads: int(tool_db.loc[\'bqsr_PrintReads\',\'threads\'])\n')
        snakef.write('    shell:\n')
        snakef.write('        \'\'\'\n')
        snakef.write('        java -Xms{params.runtime_mem} -Xmx{params.runtime_mem} -Djava.io.tmpdir={output.tmpdir} -jar {params.tool} -R {input.ref} -I {input.bam} -T PrintReads -BQSR {input.pre_recal_table} -o {output.bam}\n')
        snakef.write('        \'\'\'\n')
    snakef.close()

for indv in _individuals:
    write_sub_workflow_sampleDef(indv,_master_dict,__work_dir)
    write_sub_workflow_main_body(indv,__work_dir)

_final_run_command = 'snakemake -s ' + __work_dir + '/Snakefile' + ' -j 576 --restart-times 3 -d ' + __work_dir + '/'
os.system(_final_run_command)
