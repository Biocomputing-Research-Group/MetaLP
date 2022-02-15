'''
@author SF
'''

import sipros_post_module
import re
import sys
import getopt
class PSM:
    def __init__(self, filename, file, scan, ParentCharge, rank, MeasuredParentMass,CalculatedParentMass,Massdiff, rescore, PTM_score, IdentifiedPeptide,PSM_Label,
                 Proteins,Proteinname,ProteinCount):
        self.filename=filename
        self.file = file
        self.scan = scan
        self.ParentCharge = ParentCharge
        self.rank = rank
        self.MeasuredParentMass=MeasuredParentMass
        self.CalculatedParentMass=CalculatedParentMass
        self.Massdiff = Massdiff
        self.MassErrorPPM='NA'
        self.ScanType='HCD'
        self.SearchName='Deep learning'
        self.ScoringFunction='softmax'
        self.rescore = rescore
        self.DeltaZ='NA'
        self.DeltaP='NA'
        self.PTM_score = PTM_score
        self.IdentifiedPeptide = IdentifiedPeptide
        self.OriginalPeptide='NA'
        self.PSM_Label = PSM_Label
        self.Proteins = Proteins
        self.Proteinname=Proteinname
        self.ProteinCount=ProteinCount

class Peptide:
    def __init__(self):
        self.IdentifiedPeptide = ''
        self.ParentCharge = ''
        self.OriginalPeptide = ''
        self.ProteinNames = []
        self.Proteins = []
        self.ProteinCount = 0
        self.SpectralCount = 0
        self.BestScore = 0.0
        self.PSMs = []

    def add(self, oPsm):
        self.SpectralCount += 1
        if self.BestScore < oPsm.rescore:
            self.BestScore = oPsm.rescore
        self.PSMs.append('{0}_{1}_{2}_{3}'.format(oPsm.file, oPsm.scan, oPsm.ParentCharge, oPsm.rank))
        self.ScanType = 'HCD'
        self.SearchName = 'Deep learning'
        if oPsm.PSM_Label==True:
            self.TargetMatch='T'
        else:
            self.TargetMatch='F'

    def set(self, oPsm):
        self.IdentifiedPeptide = oPsm.IdentifiedPeptide
        self.ParentCharge = oPsm.ParentCharge
        self.OriginalPeptide = oPsm.OriginalPeptide
        self.ProteinNames = oPsm.Proteinname
        self.ProteinCount = oPsm.ProteinCount
        self.Proteins = oPsm.Proteins
        self.SpectralCount = 1
        self.BestScore = oPsm.rescore
        self.PSMs.append('{0}_{1}_{2}_{3}'.format(oPsm.file, oPsm.scan, oPsm.ParentCharge, oPsm.rank))
        self.ScanType = 'HCD'
        self.SearchName = 'Deep learning'
        if oPsm.PSM_Label==True:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'


def read_prophet(input_file, mix_version=False):
    PSMs=[]
    C_pattern = re.compile('C\[160\]')
    M_pattern = re.compile('M\[147\]')
    clean_pattern = re.compile('[">/]')
    scan_id = 0
    charge_id = ''
    original_pep = ''
    identified_pep = ''
    protein_l = []
    iProbability = 0.0
    ntt = 0
    nmc = 0

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("<spectrum_query "):
                count=0
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('spectrum='):
                        split_l_2 = one.split('=')
                        filename=split_l_2[-1].split('.')[0].replace('"','')+'.ms2'
                        prefix=split_l_2[-1].split('.')[0].split('_')
                        '''
                        if int(prefix[-2].replace('Run',''))==1:
                            file_id=int(prefix[-1])
                        else:
                            file_id=int(prefix[-1])+11
                        '''
                        file_id=prefix[-1].replace('soil','')
                    if one.startswith('start_scan='):
                        split_l_2 = one.split('=')
                        scan_id = int(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith('precursor_neutral_mass='):
                        split_l_2 = one.split('=')
                        MeasuredParentMass= float(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith('assumed_charge='):
                        split_l_2 = one.split('=')
                        charge_id = clean_pattern.sub('', split_l_2[-1])

                protein_l = []
                ntt = 2
                nmc = 0
            if line.startswith("<parameter name=\"ntt\""):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('value='):
                        split_l_2 = one.split('=')
                        ntt = int(clean_pattern.sub('', split_l_2[-1]))
            if line.startswith("<parameter name=\"nmc\""):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('value='):
                        split_l_2 = one.split('=')
                        nmc = int(clean_pattern.sub('', split_l_2[-1]))

            if line.startswith("<search_hit"):
                count+=1
                if count>1:
                    continue
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('peptide='):
                        split_l_2 = one.split('=')
                        original_pep = clean_pattern.sub('', split_l_2[-1])
                        identified_pep = original_pep
                        PTM_score=identified_pep.count('~')
                    if one.startswith("protein="):
                        split_l_2 = one.split('=')
                        protein_l.append(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith("calc_neutral_pep_mass="):
                        split_l_2 = one.split('=')
                        CalculatedParentMass=float(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith("massdiff="):
                        split_l_2 = one.split('=')
                        MassDiff=float(clean_pattern.sub('', split_l_2[-1]))

            if line.startswith("<modification_info modified_peptide"):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('modified_peptide='):
                        split_l_2 = one.split('=')
                        identified_pep = C_pattern.sub('C', (clean_pattern.sub('', split_l_2[-1])))
                        identified_pep = M_pattern.sub('M~', (clean_pattern.sub('', split_l_2[-1])))
            if line.startswith("<alternative_protein"):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('protein='):
                        split_l_2 = one.split('=')
                        tmp_str = clean_pattern.sub('', split_l_2[-1])
                        if tmp_str not in protein_l:
                            protein_l.append(tmp_str)
            if line.startswith("<peptideprophet_result "):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('probability='):
                        split_l_2 = one.split('=')
                        iProbability = float(clean_pattern.sub('', split_l_2[-1]))
            if line.startswith("</spectrum_query>"):
                PSM_Label = False
                for p in protein_l:
                    if 'Rev' not in p:
                        PSM_Label = True
                        break
                PSMs.append(PSM(filename,file_id,scan_id,charge_id,'NA',MeasuredParentMass, CalculatedParentMass, MassDiff,iProbability,PTM_score,identified_pep,PSM_Label,protein_l,','.join(protein_l),len(protein_l)))
                protein_l = []


    return PSMs


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

def re_rank(psm_list, consider_charge_bool=False):
    psm_new_list = []
    psm_dict = {}
    if consider_charge_bool:
        for oPsm in psm_list:
            sId = '{0}_{1}_{2}'.format(str(oPsm.file), str(oPsm.scan), str(oPsm.ParentCharge))
            if sId in psm_dict:
                if oPsm.rescore > psm_dict[sId].rescore:
                    psm_dict[sId] = oPsm
                elif oPsm.rescore == psm_dict[sId].rescore:
                    if abs(oPsm.Massdiff) < abs(psm_dict[sId].Massdiff):
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.Massdiff) == abs(psm_dict[sId].Massdiff):
                        # calculate PTM scores
                        if oPsm.PTM_score < psm_dict[sId].PTM_score:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTM_score == psm_dict[sId].PTM_score:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm
                            elif oPsm.IdentifiedPeptide.upper() == psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId].add_protein(oPsm.protein_list)

            else:
                psm_dict[sId] = oPsm
    else:
        for oPsm in psm_list:
            sId = '{0}_{1}'.format(str(oPsm.file), str(oPsm.scan))
            if sId in psm_dict:
                if oPsm.rescore > psm_dict[sId].rescore:
                    psm_dict[sId] = oPsm
                elif oPsm.rescore == psm_dict[sId].rescore:
                    if abs(oPsm.Massdiff) < abs(psm_dict[sId].Massdiff):
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.Massdiff) == abs(psm_dict[sId].Massdiff):
                        # calculate PTM scores
                        if oPsm.PTM_score < psm_dict[sId].PTM_score:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTM_score == psm_dict[sId].PTM_score:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm


            else:
                psm_dict[sId] = oPsm

    for _key, value in psm_dict.items():
        psm_new_list.append(value)

    return psm_new_list
def show_Fdr_Pep(psm_list, fdr_float, charge_left_given=-1, charge_right_given=-1):
    list_sorted = sorted(psm_list, key=lambda psm: (psm.rescore, psm.Massdiff, psm.PTM_score), reverse=True)

    peptide_set = set()
    decoy = 0
    target = 0
    best_nums = [0, 0]

    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        '''
        if charge_left_given != -1 and (
                oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        '''
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        #pep_str = oPsm.IdentifiedPeptide
        if pep_str not in peptide_set:
            if oPsm.PSM_Label:
                target += 1
                peptide_set.add(pep_str)
            elif not oPsm.PSM_Label:
                decoy += 1
                peptide_set.add(pep_str)
            else:
                sys.stderr.write('error 768.\n')

        (FDR_accept, FDR_value) = FDR_calculator(decoy, target)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (best_nums[0] + best_nums[1]) < (decoy + target):
                best_nums = [decoy, target]
                cutoff_probability = oPsm.rescore

    print('Number of Target Peptides at FDR '+str(fdr_float)+': '+str(best_nums[1]))
    print('Number of Decoy Peptides at FDR '+str(fdr_float)+': '+str(best_nums[0]))
    peptide = dict()
    for oPsm in list_sorted:
        '''
        if charge_left_given != -1 and (
                oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        '''
        pep_str=oPsm.IdentifiedPeptide+'_'+str(oPsm.ParentCharge)
        #pep_str = oPsm.IdentifiedPeptide
        if oPsm.rescore >= cutoff_probability:
            if pep_str in peptide:
                peptide[pep_str].add(oPsm)
            else:
                oPeptide=Peptide()
                oPeptide.set(oPsm)
                peptide[pep_str]=oPeptide

    # return set(psm_filtered_list)
    return peptide


if __name__ == "__main__":
    argv=sys.argv[1:]
    try:
        opts,args=getopt.getopt(argv,"hi:f:o:")
    except:
        print("Error Option, using -h for help information.")
        sys.exit(1)
    if len(opts)==0:
        print("\n\nUsage:\n")
        print("-i\t PeptideProphet output results, xml format")
        print("-o\t output file stored the parsing result from PeptideProphet, plain text format\n")
        print(
            "-f: Peptide level threshold, is used to keep high-quality peptide candidates for protein inference. Default: 0.01\n")
        sys.exit(1)
    inputfile=""
    outputfile=""
    fdr=0.01
    for opt, arg in opts:
        if opt in ("-h"):
            print("\n\nUsage:\n")
            print("-i\t PeptideProphet output results, xml format")
            print("-o\t output file stored the parsing result from PeptideProphet, plain text format\n")
            print("-f: Peptide level threshold, is used to keep high-quality peptide candidates for protein inference. Default: 0.01\n")
            sys.exit(1)
        elif opt in ("-i"):
            inputfile=arg
        elif opt in ("-o"):
            outputfile=arg
        else:
            fdr=float(arg)
    PSMs = read_prophet(inputfile)
    psm_list = sorted(PSMs, key=lambda psm: (psm.rescore, psm.Massdiff, psm.PTM_score), reverse=True)
    rank_list = re_rank(PSMs)
    filter_pep_list = show_Fdr_Pep(rank_list,fdr) 


    with open(outputfile,'w') as f:
        for pep in filter_pep_list:
            for protein in filter_pep_list[pep].Proteins:
                if filter_pep_list[pep].BestScore>=0.05:
                    f.write(pep+'\t'+protein+'\t'+str(filter_pep_list[pep].BestScore)+'\n')
