from gwas import process_gwas

for i in range(1,7):
    prefix = f'synthetic_small_v1_chr-{i}'
    print(f'\n\n========== starting {prefix} =====================\n\n')
    process_gwas(prefix)

    print(f'\n\n========== finished {prefix} =====================\n\n')