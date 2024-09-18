import os
import sys
import warnings

import pandas as pd

gsm_url = (
	'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59086/matrix/'
	'GSE59086_series_matrix.txt.gz'
)

print('Downloading GSE59086 data from GEO...')

with warnings.catch_warnings(): # ignore mixed data type warnings -- we'll fix
	warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)
	gsm_df = pd.read_csv(gsm_url, sep='\t', index_col=0, skiprows=31).transpose()

gsm_df.insert(
	0,
	'TFI',
	gsm_df['!Sample_characteristics_ch1'].iloc[:, 0].map(
		lambda s: s.removeprefix('experiment: '))
)
gsm_df.drop(columns=[
	'!Sample_status', '!Sample_submission_date', '!Sample_last_update_date',
	'!Sample_type', '!Sample_channel_count', '!Sample_source_name_ch1',
	'!Sample_organism_ch1', '!Sample_characteristics_ch1',
	'!Sample_characteristics_ch1', '!Sample_characteristics_ch1',
	'!Sample_characteristics_ch1', '!Sample_treatment_protocol_ch1',
	'!Sample_growth_protocol_ch1', '!Sample_molecule_ch1',
	'!Sample_extract_protocol_ch1', '!Sample_label_ch1',
	'!Sample_label_protocol_ch1', '!Sample_taxid_ch1', '!Sample_hyb_protocol',
	'!Sample_scan_protocol', '!Sample_description', '!Sample_data_processing',
	'!Sample_platform_id', '!Sample_contact_name', '!Sample_contact_email',
	'!Sample_contact_laboratory', '!Sample_contact_institute',
	'!Sample_contact_address', '!Sample_contact_city', '!Sample_contact_state',
	'!Sample_contact_zip/postal_code', '!Sample_contact_country',
	'!Sample_supplementary_file', '!Sample_data_row_count', '!Sample_relation',
	'!series_matrix_table_begin', '!series_matrix_table_end',
], inplace=True)

print('Saving...')

reads_path = os.path.join(sys.path[0], 'out', 'GSE59086_reads.tsv')
gsm_df.drop(columns=['ID_REF', 'TFI']).to_csv(reads_path, sep='\t', index=False)

meta_path = os.path.join(sys.path[0], 'out', 'GSE59086_meta.tsv')
gsm_df['TFI'].to_csv(meta_path, sep='\t')

print('Done.')
