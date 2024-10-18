import pandas as pd

# all data, sampled data
file = './data/metadata_label_noscore.csv'
file1 = './data/2024-02-04_unmasked_SARS-COV-2_sampled_200.csv'

res_dir = './count_result/'

df = pd.read_csv(file)
# df1 = pd.read_csv(file1)

# all data
df = df[(df['label_0.99'] == 1)]

# single dimension
df_year = df.groupby(['Collection Year'])['Accession ID'].count()
df_year.to_excel(res_dir + 'metadata_0.99_year.xlsx')

df_host = df.groupby(['Host'])['Accession ID'].count()
df_host.to_excel(res_dir + 'metadata_0.99_host.xlsx')

df_continent = df.groupby(['Continent'])['Accession ID'].count()
df_continent.to_excel(res_dir + 'metadata_0.99_continent.xlsx')

# df_continent_country = df.groupby(['Continent', 'Country'])['Accession ID'].count()
# df_continent_country.to_excel(res_dir + 'metadata_0.99_continent_country.xlsx')

df_variant = df.groupby(['variant_label2'])['Accession ID'].count()
df_variant.to_excel(res_dir + 'metadata_0.99_variant.xlsx')

df_lineage = df.groupby(['variant_label2', 'Pango lineage'])['Accession ID'].count()
df_lineage.to_excel(res_dir + 'metadata_0.99_variant_lineage.xlsx')

# all data for human host
df = df[(df['label_0.99'] == 1) & (df['Host'] == 'Human')]

# single dimension
df_lineage = df.groupby(['Pango lineage'])['Accession ID'].count()
df_lineage.to_excel(res_dir + 'metadata_0.99_lineage_human.xlsx')

# two dimension
df_continent_variant = df.groupby(['Continent', 'variant_label2'])['Accession ID'].count()
df_continent_variant.to_excel(res_dir + 'metadata_0.99_continent_variant_human.xlsx')

df_ym_variant = df.groupby(['Year-Month', 'variant_label2'])['Accession ID'].count()
df_ym_variant.to_excel(res_dir + 'metadata_0.99_y_m_variant_human.xlsx')

# df_ym_continent = df.groupby(['Year-Month', 'Continent'])['Accession ID'].count()
# df_ym_continent.to_excel(res_dir + 'metadata_0.99_y_m_continent_human.xlsx')
#
# df_y_continent = df.groupby(['Collection Year', 'Continent'])['Accession ID'].count()
# df_y_continent.to_excel(res_dir + 'metadata_0.99_y_continent_human.xlsx')

# three dimension
df_y_continent_variant = df.groupby(['Collection Year', 'Continent', 'variant_label2'])['Accession ID'].count()
df_y_continent_variant.to_excel(res_dir + 'metadata_0.99_y_continent_variant_human.xlsx')

# df_continent_variant_lineage = df.groupby(['Continent', 'variant_label2', 'Pango lineage'])['Accession ID'].count()
# df_continent_variant_lineage.to_excel(res_dir + 'metadata_0.99_continent_variant_lineage.xlsx')
