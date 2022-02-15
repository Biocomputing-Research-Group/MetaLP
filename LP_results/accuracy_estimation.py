def read_file(filename):
    peptide2protein_dict={}
    protein2peptide_dict={}
    pep_prob_dict={}
    with open(filename) as f:
        for line in f:
            line=line.strip().split('\t')
            peptide=line[0]
            protein=line[1]
            probability=line[2]
            if peptide not in peptide2protein_dict:
                peptide2protein_dict[peptide]=[protein]
            else:
                peptide2protein_dict[peptide].append(protein)
            if protein not in protein2peptide_dict:
                protein2peptide_dict[protein]=[peptide]
            else:
                protein2peptide_dict[protein].append(peptide)
            pep_prob_dict[peptide]=probability

    peptide2proteinGroup={}
    for key, value in protein2peptide_dict.items():
        if tuple(value) not in peptide2proteinGroup:
            peptide2proteinGroup[tuple(value)] = [key]
        else:
            peptide2proteinGroup[tuple(value)].append(key)

    proteinGrousoil3peptide={}
    for key, value in peptide2proteinGroup.items():
        proteinGrousoil3peptide[tuple(value)]=key


    return peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict


def PP_topk(fname,peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict):
    pp_list=[]
    '''
    with open(fname) as f:
        for line in f:
            s=line.strip().split('\t')
            proteins=s[0].split(',')
            score=float(s[1])
            if score<1:
                break
            else:
                pp_list.append([proteins,score])
    '''
    with open(fname) as f:
        for line_id,line in enumerate(f):
            s=line.strip().split('\t')
            proteins=s[0].split(',')
            score=float(s[1])
            if line_id>4499:
                break
            else:
                pp_list.append([proteins,score])
    target=[]
    decoy=[]
    target_degeneracy=0
    target_simple=0
    decoy_degeneracy=0
    decoy_simple=0
    for line in pp_list:
        TargetMatch=False
        for p in line[0]:
            if 'Rev_' not in p:
                TargetMatch=True
                break
        if TargetMatch:
            childpeptides=proteinGrousoil3peptide[tuple(line[0])]
            target.append(line[0])
            IF_degeneracy=False
            for pep in childpeptides:
                if len(peptide2protein_dict[pep])>1:
                    IF_degeneracy=True
                    break

            if IF_degeneracy:
                IF_count=False
                for pep in childpeptides:
                    if float(pep_prob_dict[pep])<0.9:
                        continue
                    else:
                        target_degeneracy+=1
                        break
            else:
                target_simple+=1

        else:
            childpeptides=proteinGrousoil3peptide[tuple(line[0])]
            decoy.append(line[0])
            IF_degeneracy=False
            for pep in childpeptides:
                if len(peptide2protein_dict[pep])>1:
                    IF_degeneracy=True
                    break

            if IF_degeneracy:
                IF_count=False
                for pep in childpeptides:
                    if float(pep_prob_dict[pep])<0.9:
                        continue
                    else:
                        decoy_degeneracy+=1
                        break
            else:
                decoy_simple+=1

    print('ProteinProphet result:')
    print('k of ProteinProphet is:'+str(len(pp_list)))
    print('# target degeneracy: '+str(target_degeneracy))
    print('# target simple: '+str(target_simple))
    print('# decoy degeneracy: '+str(decoy_degeneracy))
    print('# decoy simple: '+str(decoy_simple))
    print('++++++++++++++++++++++++++++++++++++++')
    return len(pp_list)

def benchmarking(fname,toolname,k,peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict):
    pro_list=[]
    with open(fname) as f:
        for line_id,line in enumerate(f):
            s=line.strip().split('\t')
            proteins=s[0].split(',')
            score=float(s[1])
            if line_id+1>k:
                break
            else:
                pro_list.append([proteins,score])
    target=[]
    decoy=[]
    target_degeneracy=0
    target_simple=0
    decoy_degeneracy=0
    decoy_simple=0
    for line in pro_list:
        TargetMatch=False
        for p in line[0]:
            if 'Rev_' not in p:
                TargetMatch=True
                break
        if TargetMatch:
            childpeptides=proteinGrousoil3peptide[tuple(line[0])]
            target.append(line[0])
            IF_degeneracy=False
            for pep in childpeptides:
                if len(peptide2protein_dict[pep])>1:
                    IF_degeneracy=True
                    break

            if IF_degeneracy:
                IF_count=False
                for pep in childpeptides:
                    if float(pep_prob_dict[pep])<0.9:
                        continue
                    else:
                        target_degeneracy+=1
                        break
            else:
                target_simple+=1

        else:
            childpeptides=proteinGrousoil3peptide[tuple(line[0])]
            decoy.append(line[0])
            IF_degeneracy=False
            for pep in childpeptides:
                if len(peptide2protein_dict[pep])>1:
                    IF_degeneracy=True
                    break

            if IF_degeneracy:
                IF_count=False
                for pep in childpeptides:
                    if float(pep_prob_dict[pep])<0.9:
                        continue
                    else:
                        decoy_degeneracy+=1
                        break
            else:
                decoy_simple+=1

    print(toolname+' result:')
    print('# target degeneracy: '+str(target_degeneracy))
    print('# target simple: '+str(target_simple))
    print('# decoy degeneracy: '+str(decoy_degeneracy))
    print('# decoy simple: '+str(decoy_simple))
    print('++++++++++++++++++++++++++++++++++++++')

def benchmarking_DeepPep(fname,toolname,k,peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict):
    pro_list=[]
    with open(fname) as f:
        for line_id,line in enumerate(f):
            s=line.strip().split('\t')
            proteins=s[0].split(',')
            newproteins=[]
            for p in proteins:
                newproteins.append(p.replace('DECOY0_','Rev_'))
            score=float(s[1])
            if line_id+1>k:
                break
            else:
                pro_list.append([newproteins,score])
    target=[]
    decoy=[]
    target_degeneracy=0
    target_simple=0
    decoy_degeneracy=0
    decoy_simple=0
    for line in pro_list:
        TargetMatch=False
        for p in line[0]:
            if 'Rev_' not in p:
                TargetMatch=True
                break
        if TargetMatch:
            if tuple(line[0]) in proteinGrousoil3peptide:
                childpeptides=proteinGrousoil3peptide[tuple(line[0])]
            else:
                continue
            target.append(line[0])
            IF_degeneracy=False
            for pep in childpeptides:
                if len(peptide2protein_dict[pep])>1:
                    IF_degeneracy=True
                    break

            if IF_degeneracy:
                IF_count=False
                for pep in childpeptides:
                    if float(pep_prob_dict[pep])<0.9:
                        continue
                    else:
                        target_degeneracy+=1
                        break
            else:
                target_simple+=1

        else:
            childpeptides=proteinGrousoil3peptide[tuple(line[0])]
            decoy.append(line[0])
            IF_degeneracy=False
            for pep in childpeptides:
                if len(peptide2protein_dict[pep])>1:
                    IF_degeneracy=True
                    break

            if IF_degeneracy:
                IF_count=False
                for pep in childpeptides:
                    if float(pep_prob_dict[pep])<0.9:
                        continue
                    else:
                        decoy_degeneracy+=1
                        break
            else:
                decoy_simple+=1

    print(toolname+' result:')
    print('# target degeneracy: '+str(target_degeneracy))
    print('# target simple: '+str(target_simple))
    print('# decoy degeneracy: '+str(decoy_degeneracy))
    print('# decoy simple: '+str(decoy_simple))
    print('++++++++++++++++++++++++++++++++++++++')

def benchmarking_sipros(fname1,fname2,toolname,k):
    pro_dict={}
    with open(fname2) as f:
        for line in f:
            s=line.strip().split('\t')
            peptide=s[0]
            protein=s[1]
            score=s[2]
            if protein not in pro_dict:
                pro_dict[protein]=[[peptide,score]]
            else:
                pro_dict[protein].append([peptide,score])

    pro_score_list=[]
    with open(fname1) as f:
        for line_id,line in enumerate(f):
            if line_id<61:
                continue
            s=line.strip().split('\t')[0]
            s=s.replace('{','')
            s=s.replace('}','')
            proteins=s.split(',')
            scores=[]
            for p in proteins:
                scores.append(pro_dict[p][0][1])

            pro_score_list.append([proteins,max(scores)])
        sorted_pro_score_list=sorted(pro_score_list,key=lambda x:x[1],reverse=True)

    pro_list=sorted_pro_score_list[0:k]
    target=[]
    decoy=[]
    target_degeneracy=0
    target_simple=0
    decoy_degeneracy=0
    decoy_simple=0
    for line in pro_list:
        TargetMatch=False
        for p in line[0]:
            if 'Rev_' not in p:
                TargetMatch=True
                break
        if TargetMatch:
            target.append(line[0])
            if len(line[0])>1:
                target_degeneracy+=1
            else:
                target_simple+=1
        else:
            decoy.append(line[0])
            if len(line[0])>1:
                decoy_degeneracy+=1
            else:
                decoy_simple+=1
    print(toolname+' result:')
    print('# target degeneracy: '+str(target_degeneracy))
    print('# target simple: '+str(target_simple))
    print('# decoy degeneracy: '+str(decoy_degeneracy))
    print('# decoy simple: '+str(decoy_simple))
    print('++++++++++++++++++++++++++++++++++++++')

if __name__ == "__main__":
    peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict=read_file('soil3.identification_1%.txt')
    k=PP_topk('soil3/PP_soil3_filteringProteins_Rev.txt',peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict)
    #benchmarking_sipros('soil3/test.pro.txt','soil3/soil3.identification_1%.txt','sipros',k)
    k=4500
    benchmarking('soil3/filteringProteins_LP_soil3_Rev.txt','LP',k,peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict)
    benchmarking('soil3/filteringProteins_LP2_soil3_Rev.txt', 'LP2', k,peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict)
    benchmarking_DeepPep('soil3/DeepPep_soil3_filteringProteins_Rev.txt', 'DeepPep', k,peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict)
    benchmarking('soil3/filteringProteinsSoil3_metaLP_Rev.txt', 'metaLP', k,peptide2protein_dict,protein2peptide_dict,peptide2proteinGroup,proteinGrousoil3peptide,pep_prob_dict)



