import os,sys

configfile: "conf/cami-low.yml"

out_dir = config['output_dir']
logs_dir = os.path.join(out_dir, 'logs')
basename = config['basename']

scaled = config.get('scaled', 1)
if not isinstance(scaled, list):
    scaled = [scaled]
sketch_scaled = min(scaled)


# find databases
database_dir = config.get('database_dir', '')
databases = config.get('databases', [])
dbs=[]
if databases:
    if not isinstance(databases, list):
        databases = [databases]
    for db in databases:
        # check file path
        if not os.path.isfile(db):
            db = os.path.join(database_dir, db)
            if not os.path.isfile(db):
                sys.stderr.write(f'Database {db} cannot be found in this dir or in database dir {database_dir}. Ignoring for now.')
                continue
        dbs.append(db)

if not dbs:
    sys.stderr.write('Please input databases to gather against! Use "databases: " in the config file.')
    sys.exit(-1)

databases = dbs

## get read_info
read_info = config.get('read_info')
if not isinstance(read_info, dict):
    sys.stderr.write(f'{read_info} must be a dictionary of trim_type: filepath')

rule all:
    input:
        #expand(os.path.join(out_dir, "sketch", '{basename}.{triminfo}.sig'), basename=basename, triminfo =read_info.keys()),
        expand(os.path.join(out_dir, "pfam-gather", "{basename}.{triminfo}.scaled{scaled}.pfam-gather.csv"), scaled =scaled, basename=basename, triminfo=read_info.keys())

# whole file sketch
rule sketch_input:
    input: lambda w: read_info[w.triminfo]
    output: os.path.join(out_dir, 'sketch', '{basename}.{triminfo}.sig')
    log: os.path.join(logs_dir, 'sketch', '{basename}.{triminfo}.log')
    benchmark: os.path.join(logs_dir, 'sketch', '{basename}.{triminfo}.benchmark')
    params:
        scaled=sketch_scaled,
        name= lambda w: f"{basename}.{w.triminfo}"
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch protein -p k=10,scaled={params.scaled},abund {input} --name {params.name} -o {output} 2> {log}
        """


# WHOLE FILE GATHER
rule file_gather:
    input:
        db=databases,
        reads_sig=os.path.join(out_dir, 'sketch', "{basename}.{triminfo}.sig"),
    output: os.path.join(out_dir, "pfam-gather", "{basename}.{triminfo}.scaled{scaled}.pfam-gather.csv"),
    log: os.path.join(logs_dir, 'pfam-gather', '{basename}.{triminfo}.scaled{scaled}.pfam-gather.log')
    benchmark: os.path.join(logs_dir, 'pfam-gather', '{basename}.{triminfo}.scaled{scaled}.pfam-gather.benchmark')
    shell:
        # sigh, I really want this to save a prefetch CSV, not the prefetch MATCHES
        """
        sourmash gather {input.reads_sig} {input.db} --scaled {wildcards.scaled} --protein --ksize 10 --threshold-bp 0 -o {output} 2> {log}
        """


# next: 
## READ BY READ GATHER
# check annotation against eggnog pfams CAMI_low_gs_read_annotations_with_eggnog.tsv.gz
# pfam key: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
