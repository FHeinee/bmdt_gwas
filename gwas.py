import pandas as pd 
import numpy as np
import math
from scipy.stats import linregress, chisquare

pheno_path = './synthetic_v1.pheno7'
sample_path = './synthetic_v1.sample'

def get_pheno_df():
    phenos = []
    with open(pheno_path, 'r') as file:
        for line in file:import pandas as pd 
import numpy as np
import math
from scipy.stats import linregress, chisquare

pheno_path = './synthetic_v1.pheno7'
sample_path = './synthetic_v1.sample'

def get_pheno_df():
    phenos = []
    with open(pheno_path, 'r') as file:
        for line in file:
            phenos.append(line.split())

    phenos_df = pd.DataFrame(phenos[1:], columns = phenos[0])
    phenos_df.set_index('Sample')

    return phenos_df

def get_fam_df(fam_path:str):
    fams = []
    fam_columns = ['fam_id','individual_id', 'father', 'mother', 'sex', 'pheno']
    with open(fam_path, 'r') as file:
        for line in file:
            fams.append(line.split())

    fam_df = pd.DataFrame(fams, columns=fam_columns)
    fam_df.set_index('individual_id')

    return fam_df

def get_ancestries_list():
    ancestries = []
    with open(sample_path, 'r') as file:
        for line in file:
            ancestries.append(line.strip())

    return ancestries

def get_bim_df(bim_path:str):
    bim_columns = ['chromosome', 'SNP', 'rel_pos', 'pos','Allele1','Allele2']
    bims = []
    with open(bim_path, 'r') as file:
        for line in file:
            bims.append(line.split())

    bim_df = pd.DataFrame(bims, columns = bim_columns)
    
    return bim_df


def decode_binary_string(byte_string) -> list[int]:
    # See https://www.cog-genomics.org/plink/1.9/formats#bed 
    results = []
    binary_rep = f"{byte_string:08b}"

    for i in reversed(range(1, 8, 2)):
        bit_pair = binary_rep[i] + binary_rep[i-1]  # Get the bit pair in order
        if bit_pair == "00":
            results.append(0) # homozygous for first allele
        elif bit_pair == "01":
            results.append(1) # heterozygous
        elif bit_pair == "11":
            results.append(2) # homozygous for second allele
        elif bit_pair == "10":
            results.append(3) # missing genotype (should never happen with our synthetic data)
        else:
            raise ValueError(f"Invalid bit pair: {bit_pair}")
    return results

def read_bed_file(bed_path: str, total_sample_count: int, sample_offset: int, sample_count: int):
    #iterating over bytes whith each byte storing the results of 4 samples
    if total_sample_count %4 !=0:
        raise RuntimeError("Number of samples must be a multiple of 4, reading "+bed_path)

    with open(bed_path, 'rb') as file:
        # Read the first 3 bytes to ignore the magic number
        (norm_bytes :=file.read(3))

        # 4 samples are stored in 1 byte, so we need to read total_sample_count/4 bytes
        bytes_per_snp = int(total_sample_count/4)
        
        # Read one entire snp at a time
        while (bytes_of_snp := file.read(bytes_per_snp)):
            samples = []

            current_byte = math.floor(sample_offset / 4)
            current_sample_index = current_byte * 4

            end_byte = math.ceil((sample_offset + sample_count) / 4)
                
            while current_byte < end_byte:
                byte = bytes_of_snp[current_byte]
                decoded_byte = decode_binary_string(byte)

                for sample in decoded_byte:
                    if current_sample_index >= sample_offset and current_sample_index < sample_offset + sample_count:
                        samples.append(sample)
                    current_sample_index += 1
                current_byte += 1

            yield samples

def calculate_minor_allele_frequency(snp: list[int]) -> float:
    total_allele_count = len(snp) * 2
    first_allele_count = snp.count(0) * 2 + snp.count(1)  # homozygous for first allele + heterozygous
    second_allele_count = snp.count(2) * 2 + snp.count(1)  # homozygous for second allele + heterozygous
    minor_allele_count = min(first_allele_count, second_allele_count)

    return minor_allele_count / total_allele_count

def calculate_hwe_pvalue(snp: list[int]) -> float:
    total_samples = len(snp)
    first_allele_homozygous_samples = snp.count(0)
    heterozygous_samples = snp.count(1)
    second_allele_homozygous_samples = snp.count(2)

    p = (2 * second_allele_homozygous_samples + heterozygous_samples) / (2 * total_samples)
    q = 1 - p

    expected_second_allele_homozygous_samples = p**2 * total_samples
    expected_heterozygous_samples = 2 * p * q * total_samples
    expected_first_allele_homozygous_samples = q**2 * total_samples

    observed = [second_allele_homozygous_samples, heterozygous_samples, first_allele_homozygous_samples]
    expected = [expected_second_allele_homozygous_samples, expected_heterozygous_samples, expected_first_allele_homozygous_samples]

    chi_squared, p_value = chisquare(f_obs=observed, f_exp=expected)

    if p_value < 1e-6:
        print(p_value)

    return p_value

def process_gwas(prefix: str, ancestry: str | None = None):
    bim_path = prefix + '.bim'
    bed_path = prefix + '.bed'
    fam_path = prefix + '.fam'

    fam_df = get_fam_df(fam_path)
    pheno_df = get_pheno_df()

    liability_df = pd.merge(
        fam_df,
        pheno_df,
        how = 'left',
        left_index=True,
        right_index=True
    )[list(fam_df.columns) + ['Phenotype(liability)']]

    if ancestry is not None:
        ancestries = get_ancestries_list()
        liability_df['ancestry'] = np.array(ancestries)
        liability = liability_df[liability_df['ancestry'] == ancestry]['Phenotype(liability)']

        # Note: This assumes that all samples of an ancestry are contiguous
        sample_offset = ancestries.index(ancestry)
        sample_count = ancestries.count(ancestry)
    else:
        liability = liability_df['Phenotype(liability)']
        sample_offset = 0
        sample_count = liability_df.shape[0]

    bim_df = get_bim_df(bim_path)
    snp_count = liability_df.shape[0]

    for snp_index, snp in enumerate(read_bed_file(bed_path, snp_count, sample_offset, sample_count)):
        liability =np.array(liability, dtype='float')

        snp_name = bim_df.iloc[snp_index]['SNP']
        if len(set(snp)) == 1:
            print(f'All samples are the same for SNP {snp_name}')
            continue
        slope, intercept, rvalue, pvalue, stderr = linregress(x = snp, y = liability)
        
        if snp_index % 1000 == 0:
            print(".", end='', flush=True)

        maf = calculate_minor_allele_frequency(snp)
        hwe_pvalue = calculate_hwe_pvalue(snp)

        with (open(f'outputs/{prefix}.csv', 'a')) as file:
            file.write(f'{snp_name},{pvalue},{maf},{hwe_pvalue}\n')
            phenos.append(line.split())

    phenos_df = pd.DataFrame(phenos[1:], columns = phenos[0])
    phenos_df.set_index('Sample')

    return phenos_df

def get_fam_df(fam_path:str):
    fams = []
    fam_columns = ['fam_id','individual_id', 'father', 'mother', 'sex', 'pheno']
    with open(fam_path, 'r') as file:
        for line in file:
            fams.append(line.split())

    fam_df = pd.DataFrame(fams, columns=fam_columns)
    fam_df.set_index('individual_id')

    return fam_df

def get_ancestries_list():
    ancestries = []
    with open(sample_path, 'r') as file:
        for line in file:
            ancestries.append(line.strip())

    return ancestries

def get_bim_df(bim_path:str):
    bim_columns = ['chromosome', 'SNP', 'rel_pos', 'pos','Allele1','Allele2']
    bims = []
    with open(bim_path, 'r') as file:
        for line in file:
            bims.append(line.split())

    bim_df = pd.DataFrame(bims, columns = bim_columns)
    
    return bim_df


def decode_binary_string(byte_string) -> list[int]:
    # See https://www.cog-genomics.org/plink/1.9/formats#bed 
    results = []
    binary_rep = f"{byte_string:08b}"

    for i in reversed(range(1, 8, 2)):
        bit_pair = binary_rep[i] + binary_rep[i-1]  # Get the bit pair in order
        if bit_pair == "00":
            results.append(0) # homozygous for first allele
        elif bit_pair == "01":
            results.append(1) # heterozygous
        elif bit_pair == "11":
            results.append(2) # homozygous for second allele
        elif bit_pair == "10":
            results.append(3) # missing genotype (should never happen with our synthetic data)
        else:
            raise ValueError(f"Invalid bit pair: {bit_pair}")
    return results

def read_bed_file(bed_path: str, total_sample_count: int, sample_offset: int, sample_count: int):
    #iterating over bytes whith each byte storing the results of 4 samples
    if total_sample_count %4 !=0:
        raise RuntimeError("Number of samples must be a multiple of 4, reading "+bed_path)

    with open(bed_path, 'rb') as file:
        # Read the first 3 bytes to ignore the magic number
        (norm_bytes :=file.read(3))

        # 4 samples are stored in 1 byte, so we need to read total_sample_count/4 bytes
        bytes_per_snp = int(total_sample_count/4)

       
        # Read one entire snp at a time
        while (bytes_of_snp := file.read(bytes_per_snp)):
            samples = []

            current_byte = math.floor(sample_offset / 4)
            current_sample_index = current_byte * 4

            end_byte = math.ceil((sample_offset + sample_count) / 4)
                
            while current_byte < end_byte:
                byte = bytes_of_snp[current_byte]
                decoded_byte = decode_binary_string(byte)

                for sample in decoded_byte:
                    if current_sample_index >= sample_offset and current_sample_index < sample_offset + sample_count:
                        samples.append(sample)
                    current_sample_index += 1
                current_byte += 1

            yield samples

def calculate_minor_allele_frequency(snp: list[int]) -> float:
    total_allele_count = len(snp) * 2
    first_allele_count = snp.count(0) * 2 + snp.count(1)  # homozygous for first allele + heterozygous
    second_allele_count = snp.count(2) * 2 + snp.count(1)  # homozygous for second allele + heterozygous
    minor_allele_count = min(first_allele_count, second_allele_count)

    return minor_allele_count / total_allele_count

def calculate_hwe_pvalue(snp: list[int]) -> float:
    total_samples = len(snp)
    first_allele_homozygous_samples = snp.count(0)
    heterozygous_samples = snp.count(1)
    second_allele_homozygous_samples = snp.count(2)

    p = (2 * second_allele_homozygous_samples + heterozygous_samples) / (2 * total_samples)
    q = 1 - p

    expected_second_allele_homozygous_samples = p**2 * total_samples
    expected_heterozygous_samples = 2 * p * q * total_samples
    expected_first_allele_homozygous_samples = q**2 * total_samples

    observed = [second_allele_homozygous_samples, heterozygous_samples, first_allele_homozygous_samples]
    expected = [expected_second_allele_homozygous_samples, expected_heterozygous_samples, expected_first_allele_homozygous_samples]
    chi_squared, p_value = chisquare(f_obs=observed, f_exp=expected)

    return p_value

def process_gwas(prefix: str, ancestry: str | None = None, saving_directory:str = "outputs/"):
    bim_path = prefix + '.bim'
    bed_path = prefix + '.bed'
    fam_path = prefix + '.fam'

    fam_df = get_fam_df(fam_path)
    pheno_df = get_pheno_df()

    liability_df = pd.merge(
        fam_df,
        pheno_df,
        how = 'left',
        left_index=True,
        right_index=True
    )[list(fam_df.columns) + ['Phenotype(liability)']]

    if ancestry is not None:
        ancestries = get_ancestries_list()
        liability_df['ancestry'] = np.array(ancestries)
        liability = liability_df[liability_df['ancestry'] == ancestry]['Phenotype(liability)']

        # Note: This assumes that all samples of an ancestry are contiguous
        sample_offset = ancestries.index(ancestry)
        sample_count = ancestries.count(ancestry)
    else:
        liability = liability_df['Phenotype(liability)']
        sample_offset = 0
        sample_count = liability_df.shape[0]

    bim_df = get_bim_df(bim_path)
    snp_count = liability_df.shape[0]

    for snp_index, snp in enumerate(read_bed_file(bed_path, snp_count, sample_offset, sample_count)):
        liability =np.array(liability, dtype='float')

        snp_name = bim_df.iloc[snp_index]['SNP']
        if len(set(snp)) == 1:
            print(f'All samples are the same for SNP {snp_name}')
            continue
        slope, intercept, rvalue, pvalue, stderr = linregress(x = snp, y = liability)
        
        if snp_index % 1000 == 0:
            print(".", end='', flush=True)

        maf = calculate_minor_allele_frequency(snp)
        try:
            hwe_pvalue = calculate_hwe_pvalue(snp)
        except:
            print(f"{snp_name} somehow no hwe")
            hwe_pvalue = -1.0

        with (open(saving_directory + f'{prefix}.csv', 'a')) as file:
            file.write(f'{snp_name},{pvalue},{maf},{hwe_pvalue}\n')