#!/usr/bin/python3
'''
    @release_date  : $release_date
    @version       : $release_version
    @author        : Sarath Chandra Dantu
    

     Copyright (C) 2020 Sarath Chandra Dantu & Alessandro Pandini

     This file is part of: Dynamics based Network cOmparisons in Python (DyNoPy)

     DyNoPy is a free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     DyNoPy is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with DyNoPy.  If not, see <http://www.gnu.org/licenses/>.

'''
comboDict={};
def aadict():
    aadict = {
        "A": 0,
        "R": 0,
        "N": 0,
        "D": 0,
        "C": 0,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": 0,
        "I": 0,
        "L": 0,
        "K": 0,
        "M": 0,
        "F": 0,
        "P": 0,
        "S": 0,
        "T": 0,
        "W": 0,
        "Y": 0,
        "V": 0,
        "-":0,
        "X":0
        }
    return aadict

'''
    to fix the dictionary bug
'''
def aalist():
    # not going to check for X
    aalist = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-"] 
    return aalist


def aacombos(aa):
    global comboDict;
    list_of_aa=aalist();
    for i in list_of_aa:
        if(i!="-"):
            comboDict[aa+i]=0;
def generate_combi_aa():
    global comboDict;
    comboDict={};
    list_of_aa=aalist();
    for aa in list_of_aa:
        if(aa!="-"):
            aacombos(aa)
    return comboDict
def funcaadict():
    aadict = {
        "G": [0],        
        "A": [0],
        "V": [0],
        "L": [0],
        "I": [0],
        "S": [0],
        "C": [0],
        "T": [0],
        "M": [0],
        "P": [0],
        "F": [0],
        "Y": [0],
        "W": [0],
        "H": [0],                
        "K": [0],
        "R": [0],        
        "N": [0],
        "Q": [0],
        "D": [0],
        "E": [0],
        "-": [0],
        "X": [0]
        }
    return aadict


def mutdict(key):
    mutdict = {
        0: ['i',1],
        1: ['j',1],
        2: ['both',2]
        }
    return mutdict[key]
def flagdict():
    flagdict = {
        "#pdbid"        : 0,
        "#nresidues"    : 0,
        "#npairs"       : 0
        }
    return flagdict

def aa2res(aa):
    aadict = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "C": "CYS",
        "Q": "GLN",
        "E": "GLU",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL", 
        }
    return aadict[aa]
def res2aa(residue):
    dict={"ALA":"A",
          "ARG":"R",
          "ASN":"N",
          "ASP":"D",
          "CYS":"C",
          "GLN":"Q",
          "GLU":"E",
          "GLY":"G",
          "HIS":"H",
          "HIP":"H",
          "HIE":"H",
          "ILE":"I",
          "LEU":"L",
          "LYS":"K",
          "MET":"M",
          "PHE":"F",
          "PRO":"P",
          "SER":"S",
          "THR":"T",
          "TYR":"Y",
          "TRP":"W",
          "VAL":"V",
          "DC":"DC",
          "sM":"M",
          "sN":"N",
          "sK":"K",
          "sA":"A",
          "sR":"R",
          "sL":"L",
          "sE":"E",
          "sI":"I",
          "sT":"T",
          "sP":"P",
          "sHd":"H",
          "sHe":"H",
          "sHp":"H",
          "sF":"F",
          "sS":"S",
          "sV":"V",
          "sQ":"Q",
          "sD":"D",
          "sY":"Y",
          "sG":"G",
          "sC":"C",
          "sW":"W",
          }
    if(residue in dict):
        return dict[residue]
    else:
        return "UNK"
def atom_mass(atom):
    atommass = {
        "H":1.008,
        "C":12.01,
        "N":14.0067,
        "O":16.001,
        "S":32.065,
        "P":30.97,
        "Na":22.989,
        "Mg":24.3050,
        "Cl":35.453,
        "Ca":40.078,
    }
    return atommass[atom]
def aafreq_from_literature(aa):
    '''
        kozlowski NAR 2016: doi:10.1093/nar/gkw978
    '''
    freq=0;
    aafreq = {
        "A": 10.06,
        "R": 5.88,
        "N": 3.58,
        "D": 5.59,
        "C": 0.94,
        "Q": 3.58,
        "E": 6.15,
        "G": 7.76,
        "H": 2.06,
        "I": 5.89,
        "L": 10.09,
        "K": 4.68,
        "M": 2.38,
        "F": 3.89,
        "P": 4.61,
        "S": 5.85,
        "T": 5.52,
        "W": 1.27,
        "Y": 2.94,
        "V": 7.27,
        "-":0 
        }
    if aa in aafreq:
        freq=aafreq[aa]/100.0;
    return freq
def get_ft_types_txt():
    '''
        https://web.expasy.org/docs/userman.html#FT_keys
    '''
    dict_ft_types={};
    
    dict_ft_types   =   {
            "SIGNAL"    :   [],
            "CA_BIND"   :   [],
            "ZN_FING"   :   [],
            "DNA_BIND"  :   [],
            "NP_BIND"   :   [],
            "ACT_SITE"  :   [],
            "METAL"     :   [],
            "BINDING"   :   [],
            "SITE"      :   [],
            "MOD_RES"   :   [],
            "CARBOHYD"  :   [],
            "DISULFID"  :   [],
            "CROSSLINK" :   [],
            "MUTAGEN"   :   [],
            }
    return dict_ft_types

def get_ft_code(ft_type,fmt="xml"):
    '''
        https://web.expasy.org/docs/userman.html#FT_keys

        0   :   FUNCTION
        1   :   SECONDARY STRUCTURE
        
    '''
    ft_code=9;
    dict_ft_types   =   {
            "ACTIVE_SITE"       :   0,
            "MUTAGENESIS_SITE"  :   0,
            "BINDING_SITE"      :   0,
            "NUCLEOTIDE_PHOSPHATE-BINDING_REGION" : 0,
            "SEQUENCE_CONFLICT" :   0,
            "REGION_OF_INTEREST":   0,
            "MODIFIED_RESIDUE"  :   0,
            "HELIX"             :   0,
            "STRAND"            :   0,
            "TURN"              :   0,
            "CHAIN"             :   2,
            }
    if(ft_type in dict_ft_types):
        ft_code =   dict_ft_types[ft_type]
    else:
        print('FT_CODE_NOT_FOUND FOR %s'%(ft_type))
        ft_code = 0
    return ft_code

def get_fn_code(ft_type,fmt="xml"):
    '''
        https://web.expasy.org/docs/userman.html#FT_keys

        0   :   FUNCTION
        1   :   SECONDARY STRUCTURE
        
    '''
    fn_code=-1;
    dict_ft_types   =   {
            "ACTIVE_SITE"       :   0,
            "MUTAGENESIS_SITE"  :   1,
            "BINDING_SITE"      :   2,
            "NUCLEOTIDE_PHOSPHATE-BINDING_REGION" : 3,
            "REGION_OF_INTEREST"    :   4,
            "SEQUENCE_CONFLICT" :   -1,
            "MODIFIED_RESIDUE"  :   -1,
            "HELIX"             :   "H",
            "STRAND"            :   "S",
            "TURN"              :   "T",
            }
    if(ft_type in dict_ft_types):
        fn_code =   dict_ft_types[ft_type]
    else:
        print('FN_CODE_NOT_FOUND FOR %s'%(ft_type))
        fn_code = -1
        
    return fn_code


def get_ss_code(ss_type="H"):
    ss_code=9;
    dict_ss =   {
            "HELIX"     :   "H",
            "STRAND"    :   "S",
            "TURN"      :   "T"
            }
    if(ss_code in dict_ss):
        ss_code=dict_ss[ss_type];
    else:
        print('SS_CODE_NOT_FOUND FOR %s'%(ss_type))
        exit()
    return ss_code

