'''
@author SF
'''

import sys
import sipros_post_module
import time
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import math
import getopt

def read_file(filename):
    peptideID = {}
    proteinID = {}
    pep_prob_dict = {}
    peptide2protein_dict = {}
    protein2peptide_dict = {}
    with open(filename) as f:
        ID = 0
        proID = 0
        for line in f:
            line = line.strip().split('\t')
            peptide = line[0]
            protein = line[1]
            probability = line[2]
            if peptide not in peptideID:
                peptideID[peptide] = ID
                ID += 1
            if protein not in proteinID:
                proteinID[protein] = proID
                proID += 1
            if peptideID[peptide] not in peptide2protein_dict:
                peptide2protein_dict[peptideID[peptide]] = [proteinID[protein]]
            else:
                peptide2protein_dict[peptideID[peptide]].append(proteinID[protein])
            if proteinID[protein] not in protein2peptide_dict:
                protein2peptide_dict[proteinID[protein]] = [peptideID[peptide]]
            else:
                protein2peptide_dict[proteinID[protein]].append(peptideID[peptide])
            pep_prob_dict[peptideID[peptide]] = probability

    return peptideID, proteinID, pep_prob_dict, peptide2protein_dict, protein2peptide_dict

def read_pro2specie(filename):
    proID2species_dict={}
    with open(filename) as f:
        for line in f:
            s=line.strip().split('\t')
            proID=s[0]
            species=s[1].split(',')
            proID2species_dict[proID]=species

    return proID2species_dict


def RevMap(my_map):
    inv_map = {v: k for k, v in my_map.items()}
    return inv_map


def write_solution(model, peptideID, proteinID, pep_prob_dict, peptide2protein_dict, protein2peptide_dict,
                         proteinGroup2peptide, proteinGroupIDs, peptideIDPro_dict, proteinGroup_dict, VarIDX, VarIDX2,outputfile):
    start=time.time()
    pro_prob_dict = {}
    RevProteinID = RevMap(proteinID)
    paraV=[v.X for v in model.getVars()]
    for proteinGroupID in proteinGroup2peptide:
        sum_prob = []
        for peptide in proteinGroup2peptide[proteinGroupID]:
            for spec_prob in peptideIDPro_dict[peptide]:
                sum_prob.append((paraV[VarIDX[(peptide, proteinGroupID, spec_prob)]]) * spec_prob)
        pro_prob_dict[proteinGroupID] = sum(sum_prob)


    pro_prob_merge_dict = dict()
    for key in pro_prob_dict:
        proteinIDs = tuple(proteinGroupIDs[key])
        for proteinID in proteinIDs:
            if proteinID in pro_prob_merge_dict:
                if pro_prob_merge_dict[proteinID]>pro_prob_dict[key]:
                    continue
                else:
                    pro_prob_merge_dict[proteinID] = pro_prob_dict[key]
            else:
                pro_prob_merge_dict[proteinID] = pro_prob_dict[key]

    f = open(outputfile, 'w')
    for group in proteinGroup_dict.values():
        proteingn = []
        for protein in group:
            proteingn.append(RevProteinID[protein])
        f.write(','.join(proteingn))
        f.write('\t')
        f.write(str(pro_prob_merge_dict[protein]))
        f.write('\n')
    f.close()

def speciePro(filename, pro2otu, proteinID):
    proteinIDPro_dict = dict()
    prob_dict = dict()
    with open(filename) as f:
        for line_id, line in enumerate(f):
            if line_id == 0:
                continue
            prob = line.strip().split('\t')[1]
            specie = line.strip().split('\t')[0]
            prob_dict[specie] = float(prob)


    proID2species_dict = read_pro2specie(pro2otu)
    probability_base = min(prob_dict.values())
    for key in proteinID:
        if key.startswith('Rev_'):
            proID = key.replace('Rev_', '')
        else:
            proID = key
        # proID = key
        if proID in proID2species_dict:
            species = proID2species_dict[proID]
            speciesPro = [prob_dict[s] for s in species]
            proteinIDPro_dict[proteinID[key]] = speciesPro
        else:
            proteinIDPro_dict[proteinID[key]] = [probability_base]

    for sp in proteinIDPro_dict:
        prob = []
        for p in proteinIDPro_dict[sp]:
            prob.append(p/probability_base)
            proteinIDPro_dict[sp] = prob

    return proteinIDPro_dict


def speciePep(proteinIDPro_dict, peptide2protein_dict):
    peptideID2Pro = dict()
    revP = RevMap(proteinID)
    for peptide in peptide2protein_dict:
        proteins = peptide2protein_dict[peptide]
        peptideprotPro = []
        for p in proteins:
            peptideprotPro.extend(proteinIDPro_dict[p])
        peptideID2Pro[peptide] = list(set(peptideprotPro))

    return peptideID2Pro


divide = sipros_post_module.divide
FDR_parameter = 1.0


def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(TP)
    FDR_accept = True

    if FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))


def proteinFDR(filename, filteringfile, fdr_float):
    filtering_proteins = []
    data = []
    decoy = 0
    target = 0
    best_nums = [0, 0]
    cutoff_probability = 1000.0
    with open(filename) as f:
        for line in f:
            line = line.strip().split('\t')
            data.append([line[0], float(line[1])])
    sort_data = sorted(data, key=lambda data: data[1], reverse=True)
    for pro_prob in sort_data:
        pros = pro_prob[0].split(',')
        TargetMatch = False
        for p in pros:
            if 'Rev_' not in p:
                TargetMatch = True
                break
        if TargetMatch:
            target += 1
        else:
            decoy += 1

        (FDR_accept, FDR_value) = FDR_calculator(decoy, target)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (best_nums[0] + best_nums[1]) < (decoy + target):
                best_nums = [decoy, target]
                cutoff_probability = pro_prob[1]

    proteins = []
    for pro_prob in sort_data:
        if pro_prob[1] >= cutoff_probability:
            pros = pro_prob[0].split(',')
            TargetMatch = False
            for p in pros:
                if 'Rev_' in p:
                    TargetMatch = False
                else:
                    TargetMatch = True
                    break
            if TargetMatch:
               filtering_proteins.append(pro_prob)
               for p in pros:
                    if 'Rev_' not in p:
                        proteins.append(p)
    print('number of target protein groups within FDR ' + str(fdr_float) + ': ' + str(len(filtering_proteins)))
    print('number of target proteins within FDR ' + str(fdr_float) + ': ' + str(len(set(proteins))))
    with open(filteringfile, 'w') as f:
        for pro_prob in filtering_proteins:
            f.write(pro_prob[0])
            f.write('\t')
            f.write(str(pro_prob[1]))
            f.write('\n')


if __name__ == "__main__":
    argv=sys.argv[1:]
    try:
        opts,args=getopt.getopt(argv,"hi:f:o:g:d:p:")
    except:
        print("Error Option, using -h for help information.")
        sys.exit(1)
    if len(opts)==0:
        print("\n\nUsage:\n")
        print("-i\t Plain text file parsed from PeptideProphet result. The plain text contains peptide, peptide's parent protein and peptide's probability\n")
        print("-o\t Protein inference result\n")
        print("-g\t Protein inference result filtered by assigned FDR\n")
        print("-f\t Protein level fdr. Default: 0.01\n")
        print("-d\t tab delimited file which contains the mapping between protein identifier with the correspongding OTUs\n")
        print("-p\t The probability of OTUs\n")
        sys.exit(1)
    inputfile=""
    outputfile=""
    fdr=0.01
    pro2otu=""
    probability_file=""
    filteringfile=""
    for opt, arg in opts:
        if opt in ("-h"):
            print("\n\nUsage:\n")
            print("-i\t Plain text file parsed from PeptideProphet result. The plain text contains peptide, peptide's parent protein and peptide's probability\n")
            print("-o\t Protein inference result\n")
            print("-g\t Protein inference result filtered by assigned FDR\n")
            print("-f\t Protein level fdr. Default: 0.01\n")
            print("-d\t tab delimited file which contains the mapping between protein identifier with the correspongding OTUs\n")
            print("-p\t The probability of OTUs\n")
            sys.exit(1)
        elif opt in ("-i"):
            inputfile=arg
        elif opt in ("-o"):
            outputfile=arg
        elif opt in ("-g"):
            filteringfile=arg
        elif opt in ("-p"):
            probability_file=arg
        elif opt in ("-d"):
            pro2otu=arg
        else:
            fdr=float(arg)

    start = time.time()
    peptideID, proteinID, pep_prob_dict, peptide2protein_dict, protein2peptide_dict = read_file(inputfile)
    print("Number of protein:" + str(len(protein2peptide_dict)))
    print("Number of peptide:" + str(len(peptide2protein_dict)))
    proteinIDPro_dict = speciePro(probability_file, pro2otu, proteinID)
    peptideIDPro_dict = speciePep(proteinIDPro_dict,peptide2protein_dict)
    e = 0
    VarIDX = {}  # pij
    VarIDX2 = {}  # tj

    proteinGroup_dict = {}
    for key, value in protein2peptide_dict.items():
        if tuple(value) not in proteinGroup_dict:
            proteinGroup_dict[tuple(value)] = [key]
        else:
            proteinGroup_dict[tuple(value)].append(key)
    print('Number of protein groups:'+str(len(proteinGroup_dict)))
    proteinGroupIDs = {}
    proteinGroup2peptide = {}
    proteinGroupIDPro_dict={}
    groupIDs = 0
    for key, value in proteinGroup_dict.items():
        specie_probs = proteinIDPro_dict[value[0]]
        #proteinGroupIDs[groupIDs] = [value, specie_probs]
        proteinGroupIDs[groupIDs]=value
        proteinGroupIDPro_dict[groupIDs]=specie_probs
        #proteinGroup2peptide[groupIDs] = [list(key), specie_probs]
        proteinGroup2peptide[groupIDs] = list(key)
        groupIDs += 1
        for protein in value:
            if proteinIDPro_dict[protein] == specie_probs:
                continue
            else:
                #proteinGroupIDs[groupIDs] = [value, proteinIDPro_dict[protein]]
                proteinGroupIDs[groupIDs] = value
                proteinGroupIDPro_dict[groupIDs] = specie_probs
                #proteinGroup2peptide[groupIDs] = [list(key), proteinIDPro_dict[protein]]
                proteinGroup2peptide[groupIDs]=list(key)
                groupIDs += 1

    peptide2proteinGroup = {}
    for key, value in proteinGroup2peptide.items():
        for peptide in value:
            if peptide not in peptide2proteinGroup:
                peptide2proteinGroup[peptide] = [key]
            else:
                peptide2proteinGroup[peptide].append(key)

    model = gp.Model("LP")
    # Add variables: pij and tj
    p = []
    t = []
    index = 0
    for peptide in pep_prob_dict:
        for proteinGroup in peptide2proteinGroup[peptide]:
            for spe_prob in peptideIDPro_dict[peptide]:
                VarIDX[(peptide, proteinGroup, spe_prob)] = index
                p.append(model.addVar(0, 1, 0, vtype=GRB.CONTINUOUS))
                index += 1

    index = 0
    for proteinGroupID in proteinGroup2peptide:
        VarIDX2[proteinGroupID] = index
        t.append(model.addVar(0, GRB.INFINITY, 1.0, vtype=GRB.CONTINUOUS))
        index += 1

    model.update()
    # Add constrains
    for peptide in pep_prob_dict:
        Pepi = float(pep_prob_dict[peptide])
        Vs = []
        for proteinGroup in peptide2proteinGroup[peptide]:
            for spe_prob in peptideIDPro_dict[peptide]:
                Vs.append(p[VarIDX[(peptide, proteinGroup, spe_prob)]] * spe_prob)
        model.addConstr(gp.quicksum(Vs) <= Pepi + e)

    for peptide in pep_prob_dict:
        Pepi = float(pep_prob_dict[peptide])
        Vs = []
        for proteinGroup in peptide2proteinGroup[peptide]:
            for spe_prob in peptideIDPro_dict[peptide]:
                Vs.append(p[VarIDX[(peptide, proteinGroup, spe_prob)]] * spe_prob)
        model.addConstr(gp.quicksum(Vs) >= Pepi - e)

    for peptide in pep_prob_dict:
        for proteinGroup in peptide2proteinGroup[peptide]:
            model.addConstr(gp.quicksum(p[VarIDX[(peptide, proteinGroup, spec_prob)]] * spec_prob for spec_prob in
                                    proteinGroupIDPro_dict[proteinGroup]) <= t[VarIDX2[proteinGroup]])
    model.setObjective(gp.quicksum(t[j] for j in range(len(t))), GRB.MINIMIZE)

    model.optimize()
    print('calculation time:'+str(time.time()-start))
    start2=time.time()
    write_solution(model, peptideID, proteinID, pep_prob_dict, peptide2protein_dict, protein2peptide_dict,
                         proteinGroup2peptide, proteinGroupIDs,peptideIDPro_dict, proteinGroup_dict, VarIDX, VarIDX2,outputfile)
    print('solution time:'+str(time.time()-start2))
    print('epsilon:'+str(e))
    proteinFDR(outputfile, filteringfile,fdr)
    end = time.time()
    print('Completion Time:' + str(end - start))
