from rnacentral import QueryRnaCentral

rc = QueryRnaCentral()

# ------------------- Query RefSeq ID to RNAcentral -------------------
refseq_id = "NR_036832"  # example RefSeq
taxid = "58123"
print(f"\nRefSeq to RNAcentral (RefSeq={refseq_id}, taxid={taxid}):")
df_ref = rc.refseq_to_rnacentral(refseq_id, taxid)
print(df_ref)

# ------------------- Query RNAcentral metadata -------------------

urs_id = "URS000000000F"#"URS0000AB68F4"  # example RNAcentral ID 
print(f"\nRNAcentral metadata for URS={urs_id}:")
df_meta = rc.rnacentral_metadata(urs_id)
print(df_meta)

# ------------------- Query all external IDs for a URS and taxid -------------------
urs_id = "URS0000AB68F4"  # example RNAcentral ID 
taxid = "300269"
print(f"\nExternal IDs for RNAcentral ID={urs_id, taxid}:")
df_ext = rc.rnacentral_to_external_ids(urs_id, taxid)
print(df_ext)

# ------------------- Query Ensembl ID to RNAcentral -------------------
ensembl_id = "ENST00000618786"  # example
print(f"\nEnsembl to RNAcentral (Ensembl ID={ensembl_id}):")
df_ens = rc.ensembl_to_rnacentral(ensembl_id)
print(df_ens)

taxid= "9606"
ensembl_id = "ENST00000618786"  # example
print(f"\nEnsembl to RNAcentral (Ensembl ID={ensembl_id, taxid}):")
df_ens = rc.ensembl_to_rnacentral(ensembl_id, taxid)
print(df_ens)

# ------------------- Query sequence by RNAcentral ID -------------------
# URS is unique for a unique sequence!

urs_id = "URS00000478B7"
print(f"\nSequence for RNAcentral ID={urs_id}:")
seq_result = rc.rnacentral_to_sequence(urs_id)
print(seq_result)

# ------------------- List all distinct external databases -------------------
print("\nDistinct external databases in id_map:")
df_dbs = rc.con.execute("SELECT DISTINCT external_db FROM id_map ORDER BY external_db").df()
print(df_dbs)

# ------------------- Done -------------------
print('ALL SET AND RUNNING, CLOSING DB')
rc.close()
