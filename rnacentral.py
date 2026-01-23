import duckdb
import hashlib
from pathlib import Path
import pandas as pd
import gzip

class QueryRnaCentral:
    def __init__(self, base_dir="current_release", db_path="rnacentral.duckdb"):
        self.base = Path(base_dir)
        self.con = duckdb.connect(db_path, read_only=True)
        # self._init_tables() COMMENT THIS IF THE DB EXISTS

    def _init_tables(self):
        # id mapping (tab-delimited; variable number of columns allowed)
        self.con.execute("""
            CREATE TABLE IF NOT EXISTS id_map AS
            SELECT
                column0        AS RNAcentral_ID,
                column1        AS external_db,
                column2        AS external_id,
                column3        AS taxid,
                column4        AS rna_type,
                column5        AS misc
            FROM read_csv_auto(?, delim='\t', header=False)
        """, [str(self.base / "id_mapping" / "id_mapping.tsv.gz")])

        # Rfam annotations (tab-delimited)
        self.con.execute("""
            CREATE TABLE IF NOT EXISTS rfam_map AS
            SELECT
                column0        AS RNAcentral_ID,
                column1        AS rfam_id,
                column2        AS score,
                column3        AS evalue,
                column4        AS seq_start,
                column5        AS seq_stop,
                column6        AS model_start,
                column7        AS model_stop,
                column8        AS rfam_description
            FROM read_csv_auto(?, delim='\t', header=False)
        """, [str(self.base / "rfam" / "rfam_annotations.tsv.gz")])

        # GO annotations: read with pandas to split RNAcentralID_TAXID -> RNAcentral_ID + taxid
        go_path = self.base / "go_annotations" / "rnacentral_rfam_annotations.tsv.gz"
        if go_path.exists():
            df_go = pd.read_csv(go_path, sep='\t', header=None, compression='gzip', names=['col0','go_id','rfam_id'], dtype=str)
            # split on last underscore to avoid splitting IDs that may contain underscores
            split_df = df_go['col0'].str.rsplit('_', n=1, expand=True)
            df_go['RNAcentral_ID'] = split_df[0]
            df_go['taxid'] = split_df[1]
            df_go = df_go[['RNAcentral_ID', 'taxid', 'go_id', 'rfam_id']]

            # register DataFrame in DuckDB and persist as table
            self.con.register('tmp_go_df', df_go)
            self.con.execute("CREATE TABLE IF NOT EXISTS go_map AS SELECT * FROM tmp_go_df")
            # optional: unregister to free name (duckdb will keep the table)
            self.con.unregister('tmp_go_df')

        # Sequences: parse FASTA.gz (headers like: >URS000149A9AF rRNA from 1 species)
        seq_path = self.base / "sequences" / "rnacentral_active.fasta.gz"
        if seq_path.exists():
            records = []
            with gzip.open(seq_path, 'rt') as fh:
                cur_id = None
                cur_type = None
                cur_misc = None
                seq_chunks = []
                for line in fh:
                    line = line.rstrip('\n')
                    if not line:
                        continue
                    if line.startswith('>'):
                        # store previous
                        if cur_id is not None:
                            records.append({
                                'RNAcentral_ID': cur_id,
                                'rna_type': cur_type,
                                'misc': cur_misc,
                                'sequence': ''.join(seq_chunks)
                            })
                        header = line[1:].strip()
                        parts = header.split()
                        cur_id = parts[0]
                        cur_type = parts[1] if len(parts) >= 2 else None
                        cur_misc = ' '.join(parts[2:]) if len(parts) >= 3 else None
                        seq_chunks = []
                    else:
                        seq_chunks.append(line.strip())
                # final record
                if cur_id is not None:
                    records.append({
                        'RNAcentral_ID': cur_id,
                        'rna_type': cur_type,
                        'misc': cur_misc,
                        'sequence': ''.join(seq_chunks)
                    })

            df_seq = pd.DataFrame(records, columns=['RNAcentral_ID', 'rna_type', 'misc', 'sequence'])
            self.con.register('tmp_seq_df', df_seq)
            self.con.execute("CREATE TABLE IF NOT EXISTS sequence_map AS SELECT * FROM tmp_seq_df")
            self.con.unregister('tmp_seq_df')

        # md5 mapping (tab-delimited)
        self.con.execute("""
            CREATE TABLE IF NOT EXISTS md5_map AS
            SELECT
                column0 AS RNAcentral_ID,
                column1 AS md5
            FROM read_csv_auto(?, delim='\t', header=False)
        """, [str(self.base / "md5" / "md5.tsv.gz")])

        # Indexes for fast lookups
        # id_map
        self.con.execute("CREATE INDEX IF NOT EXISTS idx_idmap_ext_tax ON id_map (external_id, taxid)")
        self.con.execute("CREATE INDEX IF NOT EXISTS idx_idmap_db_ext ON id_map (external_db, external_id)")
        self.con.execute("CREATE INDEX IF NOT EXISTS idx_idmap_urs ON id_map (RNAcentral_ID)")
        
        # md5_map
        self.con.execute("CREATE INDEX IF NOT EXISTS idx_md5 ON md5_map (md5)")
        self.con.execute("CREATE INDEX IF NOT EXISTS idx_md5_urs ON md5_map (RNAcentral_ID)")
        
        # go_map
        self.con.execute("CREATE INDEX IF NOT EXISTS idx_go_urs_tax ON go_map (RNAcentral_ID, taxid)")
        
        # rfam_map
        self.con.execute("CREATE INDEX IF NOT EXISTS idx_rfam_urs ON rfam_map (RNAcentral_ID)")
        
        # sequence_map
        self.con.execute("CREATE INDEX IF NOT EXISTS idx_seq_urs ON sequence_map (RNAcentral_ID)")


    def close(self):
        self.con.close()

    @staticmethod
    def seq_to_md5(seq):
        seq = seq.upper().replace("U", "T")
        return hashlib.md5(seq.encode()).hexdigest()

    
    # %%%%--------------------%%%% Queries %%%%--------------------%%%%

    def refseq_to_rnacentral(self, refseq_id, taxid):
        """
        Given RefSeq ID and taxid, return RNAcentral ID and RNA type.
        """
        return self.con.execute("""
            SELECT DISTINCT
                RNAcentral_ID,
                rna_type
            FROM id_map
            WHERE external_db = 'REFSEQ'
              AND external_id = ?
              AND taxid = ?
        """, [refseq_id, taxid]).df()

    def rnacentral_metadata(self, urs_id): 
        """ 
        Given RNAcentral ID, return GO terms and Rfam annotations (if available). 
        """ 
        return self.con.execute(""" 
            SELECT DISTINCT 
                g.go_id, 
                g.rfam_id AS go_rfam_id, 
                r.rfam_id, 
                r.rfam_description 
            FROM rfam_map r 
            LEFT JOIN go_map g 
              ON r.RNAcentral_ID = g.RNAcentral_ID 
            WHERE r.RNAcentral_ID = ? """, [urs_id]).df()
    
    '''
    function for taxID related GO terms only:
    def rnacentral_metadata(self, urs_id, taxid):
        """
        Given RNAcentral ID and taxid, return GO terms (filtered by taxid)
        and Rfam annotations (no taxid filtering, RFAM is universal).
        """
        return self.con.execute("""
            SELECT DISTINCT
                g.go_id,
                g.rfam_id AS go_rfam_id,
                g.taxid,
                r.rfam_id,
                r.rfam_description
            FROM rfam_map r
            LEFT JOIN go_map g
              ON r.RNAcentral_ID = g.RNAcentral_ID AND g.taxid = ?
            WHERE r.RNAcentral_ID = ?
        """, [taxid, urs_id]).df()'''
    
    
    def rnacentral_to_external_ids(self, urs_id, taxid):
        """
        Given RNAcentral ID, return all external database mappings.
        """
        return self.con.execute("""
            SELECT
                external_db,
                external_id,
                taxid,
                rna_type
            FROM id_map
            WHERE RNAcentral_ID = ?
              AND taxid = ?
        """, [urs_id, taxid]).df()
    
    
    def rnacentral_to_sequence(self, urs_id):
        """
        Given RNAcentral ID, return the RNA sequence.
        """
        return self.con.execute("""
            SELECT
                sequence
            FROM sequence_map
            WHERE RNAcentral_ID = ?
        """, [urs_id]).fetchone()

    
    #TODO: this function must to be tested and refactored, 
    # because sometimes RNA sequences are saved in RNAcentral with Ts instead of Us
    def sequence_to_rnacentral(self, sequence, taxid=None):
        """
        Given RNA sequence and optional taxid, return matching RNAcentral IDs.
        """
        md5 = self.seq_to_md5(sequence)
    
        if taxid is None:
            return self.con.execute("""
                SELECT DISTINCT
                    m.RNAcentral_ID
                FROM md5_map m
                WHERE m.md5 = ?
            """, [md5]).df()
    
        return self.con.execute("""
            SELECT DISTINCT
                m.RNAcentral_ID
            FROM md5_map m
            JOIN id_map i
              ON m.RNAcentral_ID = i.RNAcentral_ID
            WHERE m.md5 = ?
              AND i.taxid = ?
        """, [md5, taxid]).df()    


    def ensembl_to_rnacentral(self, ensembl_id, taxid=None):
        """
        Given Ensembl ID and optional taxid, return RNAcentral ID and RNA type.
        """
        if taxid is None:
            return self.con.execute("""
                SELECT DISTINCT
                    RNAcentral_ID,
                    rna_type,
                    taxid
                FROM id_map
                WHERE external_id = ?
                  AND external_db LIKE 'ENSEMBL%'
            """, [ensembl_id]).df()
    
        return self.con.execute("""
            SELECT DISTINCT
                RNAcentral_ID,
                rna_type,
                taxid
            FROM id_map
            WHERE external_id = ?
              AND taxid = ?
              AND external_db LIKE 'ENSEMBL%'
        """, [ensembl_id, taxid]).df()


