import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_chromosome_results(chromosome: int) -> pd.DataFrame:
    filename = f'./data/outputs/synthetic_small_v1_chr-{str(chromosome)}.csv'
    df = pd.read_csv(filename, header=None)
    df.columns = ['snp', 'pvalue']

    df['chromosome'] = df['snp'].str.split(':').str[0]
    df['position'] = df['snp'].str.split(':').str[1]
    df['position'] = df['position'].astype(int)
    df['log_pvalue'] = -np.log10(df['pvalue'])

    return df

def read_all_results():
    combined_df = pd.DataFrame()

    # read all data for chromosomes 1 - 22 (end is exclusive)
    chromosomes = range(1, 23)

    for chromosome in chromosomes:
        print(f"Reading data for chromosome {str(chromosome)}")
        combined_df = pd.concat([
            combined_df,
            read_chromosome_results(chromosome)
        ], ignore_index=True)

    print(f"Read a total of {str(len(combined_df))}")

    return combined_df


def manhattan_plot(df: pd.DataFrame):
    cumulative_offset = 0
    tick_positions = []
    unique_chromosomes = sorted(df['chromosome'].unique())
    chromosome_colors = ['#1f77b4', '#ff7f0e']

    for chromosome in unique_chromosomes:
        chromosome_data = df[df['chromosome'] == chromosome]
        chromosome_data = chromosome_data.sort_values(by='position')
        df.loc[chromosome_data.index, 'cumulative_position'] = chromosome_data['position'] + cumulative_offset

        cumulative_offset += chromosome_data['position'].max()
        #tick_positions.append((chromosome_data['cumulative_position'].iloc[0] + chromosome_data['cumulative_position'].iloc[-1]) / 2)

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

    # Add labels and title
    #plt.xticks(ticks=tick_positions, labels=unique_chromosomes)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.title('GWAS')
    plt.tight_layout()
    plt.show()
    
df = read_all_results()
manhattan_plot(df)