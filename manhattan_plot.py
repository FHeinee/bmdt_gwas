import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_chromosome_results(chromosome: int, path: str, preprocessed: bool = False) -> pd.DataFrame:
    filename = path + f'/synthetic_v1_chr-{str(chromosome)}.csv'
    df = pd.read_csv(filename, header=None)
    df.columns = ['snp', 'pvalue', 'maf', 'hwe']if preprocessed else['snp', 'pvalue']

    df['chromosome'] = df['snp'].str.split(':').str[0].str[3:].astype(int)
    df['position'] = df['snp'].str.split(':').str[1].astype(int)
    df['log_pvalue'] = -np.log10(df['pvalue'])

    return df

def read_all_results(path: str, preprocessed:bool = False, maf_threshold:float = 0.05, hwe_threshold:float = 1e-6):
    combined_df = pd.DataFrame()

    # read all data for chromosomes 1 - 22 (end is exclusive)
    chromosomes = range(1, 23)

    for chromosome in chromosomes:
        print(f"Reading data for chromosome {str(chromosome)}")
        combined_df = pd.concat([
            combined_df,
            read_chromosome_results(chromosome,path, preprocessed)
        ], ignore_index=True)

    print(f"Read a total of {str(len(combined_df))}")

    if preprocessed:
        combined_df = combined_df[(combined_df['maf']>= maf_threshold) & (combined_df['hwe'] >= hwe_threshold)]
        print(f"After preprocessing step, {len(combined_df)} SNPs remain")


    return combined_df

def manhattan_plot(df: pd.DataFrame, preprocessed: bool =False):
    cumulative_offset = 0
    tick_positions = []
    peak_snps = []
    unique_chromosomes = sorted(df['chromosome'].unique())
    chromosome_colors = ['#1f77b4', '#ff7f0e']

    for chromosome in unique_chromosomes:
        chromosome_data = df[df['chromosome'] == chromosome &df['maf']>0.05] if preprocessed else df[df['chromosome'] == chromosome ]
        chromosome_data = chromosome_data.sort_values(by='position')
        df.loc[chromosome_data.index, 'cumulative_position'] = chromosome_data['position'] + cumulative_offset

        maximum_position = chromosome_data['position'].max()

        tick_position = (cumulative_offset + cumulative_offset + maximum_position) / 2
        tick_positions.append(tick_position)

        peak_snp_index = chromosome_data['log_pvalue'].idxmax()
        peak_snp = chromosome_data.loc[peak_snp_index]
        peak_snps.append({
            'snp': peak_snp['snp'],
            'log_pvalue': peak_snp['log_pvalue'],
            'position': peak_snp['position'],
            'cumulative_position': peak_snp['position'] + cumulative_offset,
        })

        cumulative_offset += maximum_position

    peak_snps_df = pd.DataFrame(peak_snps).sort_values(by='log_pvalue', ascending=False)
    first_3_peak_snps = peak_snps_df.head(3)

    plt.figure(figsize=(12,6))
    for chromosome_index, chromosome in enumerate(unique_chromosomes):
        chromosome_data = df[df['chromosome'] == chromosome]
        plt.scatter(
            chromosome_data['cumulative_position'], 
            chromosome_data['log_pvalue'],
            color=chromosome_colors[chromosome_index % len(chromosome_colors)], 
            alpha=0.6, 
            s=10, 
            label=f'Chr {chromosome}'
        )

    significance_value = -np.log10(5e-8)
    plt.axhline(significance_value, color='red', linestyle='--', label=f'Significance ({significance_value})')

    for index, row in first_3_peak_snps.iterrows():
        plt.text(row['cumulative_position'], row['log_pvalue'], f"{row['snp']}", fontsize=8)

    plt.xticks(ticks=tick_positions, labels=unique_chromosomes)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.ylim(0,60)
    plt.title('GWAS')
    plt.tight_layout()
    plt.show()

path = './data/preprocessed_full_AFR'
df = read_all_results(path = path,preprocessed=True)
manhattan_plot(df)