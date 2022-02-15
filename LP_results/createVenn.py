from venn import venn

def read_protein(file):
    protein_list=[]
    with open(file) as f:
        for line in f:
            s=line.strip().split('\t')[0]
            proteins=s.split(',')
            for p in proteins:
                if p.startswith('Rev_'):
                    continue
                else:
                    protein_list.append(p)
    return set(protein_list)

def read_protein_sipros(file):
    protein_list=[]
    with open(file) as f:
        for line_id,line in enumerate(f):
            if line_id<61:
                continue
            s = line.strip().split('\t')[0]
            s = s.replace('{','')
            s = s.replace('}', '')
            proteins = s.split(',')
            for p in proteins:
                if p.startswith('Rev_'):
                    continue
                else:
                    protein_list.append(p)
    return set(protein_list)

if __name__ == "__main__":
    sipros_set=read_protein_sipros('HumanGutFiber/test.pro.txt')
    LP_set=read_protein('HumanGutFiber/filteringProteins_LP_HumanGut.txt')
    pp_set=read_protein('HumanGutFiber/PP_human_gut_filteringProteins.txt')
    deeppep_set=read_protein('HumanGutFiber/DeepPep_HumanGut_filteringProteins.txt')
    metaLP_set=read_protein('HumanGutFiber/filteringProteinsHumanGut_metaLP.txt')

    benchmark_set={'SE-PI':sipros_set,'LP':LP_set,'PP':pp_set,'DP':deeppep_set,'metaLP':metaLP_set}

    second_frequency_dict = {}
    for key in pp_set:
        key = key.split('_')[0]
        if key in second_frequency_dict:
            second_frequency_dict[key] += 1
        else:
            second_frequency_dict[key] = 1

    metalp_frequency_dict = {}

    for key in metaLP_set:
        key = key.split('_')[0]
        if key in metalp_frequency_dict:
            metalp_frequency_dict[key] += 1
        else:
            metalp_frequency_dict[key] = 1

    sorted(second_frequency_dict.items(), key=lambda v: v[1])
    sorted(metalp_frequency_dict.items(), key=lambda v: v[1])


    venn(benchmark_set)
