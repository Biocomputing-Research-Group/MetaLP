
def pro2spec():
    proteins=dict()
    with open('human_gut_fiber_V2_Rev.fasta') as f:
        for line_id ,line in enumerate(f):
            if line_id%2==0:
                proID=line.strip().split( )[0][1:]
            else:
                protein=line.strip()
                proteins[proID]=protein

    pro2spec_dict=dict()

    for proID in proteins:
        if proteins[proID] in pro2spec_dict:
            pro2spec_dict[proteins[proID]].append(proID)
        else:
            pro2spec_dict[proteins[proID]]=[proID]

    proID2species_dict=dict()



    prefix2specie=dict()
    with open('protein_probability.txt') as f:
        for line_id,line in enumerate(f):
            if line_id==0:
                continue
            prefixes=line.strip().split('\t')[2].split(',')
            specie=line.strip().split('\t')[0]
            for prefix in prefixes:
                prefix2specie[prefix]=specie

    for value in pro2spec_dict.values():
        for proID in value:
            specie=[prefix2specie[x.split('_')[0]] for x in value if x.split('_')[0] in prefix2specie]
            if proID.split('_')[0] in prefix2specie:
                proID2species_dict[proID]=specie

    for protein in proID2species_dict:
        proID2species_dict[protein]=list(set(proID2species_dict[protein]))


    return proID2species_dict

def pro2specHuman():
    proteins=dict()
    with open('human_gut_fiber_V2_Rev.fasta') as f:
        for line_id ,line in enumerate(f):
            if line_id%2==0:
                proID=line.strip().split( )[0][1:]
            else:
                protein=line.strip()
                proteins[proID]=protein

    pro2spec_dict=dict()

    for proID in proteins:
        if proteins[proID] in pro2spec_dict:
            pro2spec_dict[proteins[proID]].append(proID)
        else:
            pro2spec_dict[proteins[proID]]=[proID]

    proID2species_dict=dict()


    prefix2specie=dict()
    with open('protein_probability_humanGut.txt') as f:
        for line_id,line in enumerate(f):
            if line_id==0:
                continue
            prefixes=line.strip().split('\t')[2].split(',')
            specie=line.strip().split('\t')[0]
            for prefix in prefixes:
                prefix2specie[prefix]=specie


    for value in pro2spec_dict.values():
        for proID in value:
            specie=[]
            for x in value:
                if x.split('_')[0] == 'Rev':
                    flag = '_'.join(x.split('_')[1:])
                    for prefix in prefix2specie:
                        if prefix in flag:
                            break
                        else:
                            prefix='Unknown'
                    specie.append(prefix2specie[prefix])

                else:
                    flag=x
                    for prefix in prefix2specie:
                        if prefix in flag:
                            break
                        else:
                            prefix = 'Unknown'
                    specie.append(prefix2specie[prefix])

            proID2species_dict[proID]=specie

    for protein in proID2species_dict:
        proID2species_dict[protein]=list(set(proID2species_dict[protein]))


    return proID2species_dict

def pro2specSoil():
    proteins=dict()
    with open('soil3_assembly_REV.fasta') as f:
        for line_id ,line in enumerate(f):
            if line_id%2==0:
                proID=line.strip().split( )[0][1:]
            else:
                protein=line.strip()
                proteins[proID]=protein

    pro2spec_dict=dict()

    for proID in proteins:
        if proteins[proID] in pro2spec_dict:
            pro2spec_dict[proteins[proID]].append(proID)
        else:
            pro2spec_dict[proteins[proID]]=[proID]

    proID2species_dict=dict()



    prefix2specie=dict()
    with open('soil3_to_bin.tsv') as f:
        for line_id,line in enumerate(f):
            if line_id==0:
                continue
            prefix=line.strip().split('\t')[0]
            specie=line.strip().split('\t')[1]
            prefix2specie[prefix]=specie
    for value in pro2spec_dict.values():
        for proID in value:
            specie=[]
            for x in value:
                if x.split('_')[0].startswith('Rev'):
                    if ('_'.join(x.split('_')[1:-1])) in prefix2specie:
                        specie.append(prefix2specie['_'.join(x.split('_')[1:-1])])
                    else:
                        specie.append('Unknown')
                else:
                    if ('_'.join(x.split('_')[0:-1])) in prefix2specie:
                        specie.append(prefix2specie['_'.join(x.split('_')[0:-1])])
                    else:
                        specie.append('Unknown')

            proID2species_dict[proID]=specie

    for protein in proID2species_dict:
        proID2species_dict[protein]=list(set(proID2species_dict[protein]))


    return proID2species_dict



if __name__=="__main__":
    #proID2species_dict=pro2spec()
    proID2species_dict = pro2specSoil()
    print('done')
