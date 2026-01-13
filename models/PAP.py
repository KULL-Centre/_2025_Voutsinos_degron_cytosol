#!/usr/bin/env python
# (C) 2023-2025 by Kristoffer E. Joahnsson <kristoffer.johansson@bio.ku.dk>

import sys
import os
import argparse
import numpy as np
import time
# imports tensorfolw on request

__version__ = "0.2"

aa_one_nat = "ACDEFGHIKLMNPQRSTVWY"

def is_aa_one_nat(sequence, additional=""):
    for aa in sequence.upper():
        if not (aa in aa_one_nat or aa in additional.upper()):
            return(False)
    return(True)


def read_fasta(filename, comment_char="#;", extra_aa="", skip_nonnatural=False, verbose=0):
    """Flexible FASTA file reader without dependencies"""
    seq_list = []
    file_id = filename.split("/")[-1].split(".")[0]
    reading_fasta = False
    with open(filename,"r") as file_handle:
        for line in file_handle:
            words = line.split()
            if len(words) == 0:
                continue
            if words[0][0] in comment_char:
                continue
            # start reading new fasta record
            if words[0][0] == ">":
                if reading_fasta:
                    if is_aa_one_nat(seq, extra_aa) or not skip_nonnatural:
                        seq_list.append((seq_id,seq))
                    else:
                        if verbose > 0:
                            print("Skipping protein with non-standard amino acids: %s" % (seq_id), file=sys.stderr)
                seq_id = words[0][1:]
                seq = ""
                reading_fasta = True
            # continue reading fasta record
            elif reading_fasta:
                if len(words) > 1:
                    raise ValueError("Found FASTA line with more white space separated fields: '%s'" % (line.strip()))
                seq = seq + words[0]
            # read table-format (.seq) with single-column sequences or multi-column with sequences in second column
            else:
                if len(words) == 1:
                    seq = words[0]
                    if is_aa_one_nat(seq, extra_aa) or not skip_nonnatural:
                        seq_id = file_id+"%05d" % (len(seq_list))
                        seq_list.append((seq_id,seq))
                    else:
                        if verbose > 0:
                            print("Skipping line in non-FASTA single-column file with non-natural amino acid sequence in first column:", file=sys.stderr)
                            print("         %s" % (line.strip()), file=sys.stderr)
                elif len(words) > 1:
                    seq = words[1]
                    if is_aa_one_nat(seq, extra_aa) or not skip_nonnatural:
                        seq_list.append((words[0],seq))
                    else:
                        print("WARNING: Skipping line in non-FASTA multi-column file with non-natural amino acid sequence in second column:", file=sys.stderr)
                        print("         %s" % (line.strip()), file=sys.stderr)
        # append the last fasta record
        if reading_fasta:
            if is_aa_one_nat(seq, extra_aa) or not skip_nonnatural:
                seq_list.append((seq_id,seq))
            else:
                print("WARNING: Skipping protein with non-standard amino acids: %s" % (seq_id), file=sys.stderr)
    return(seq_list)


yeast17 = {}
yeast17["type"] = "glm"
yeast17["link"] = "sigmoid"
yeast17["size"] = 17
yeast17["intercept"] = -0.89102423
yeast17["comp"] = {}
yeast17["comp"]["C"] =   0.37721431
yeast17["comp"]["D"] =  -0.78986558
yeast17["comp"]["E"] =  -0.65124014
yeast17["comp"]["K"] =  -0.15518666
yeast17["comp"]["R"] =  -0.02030300
yeast17["comp"]["H"] =  -0.02110156
yeast17["comp"]["N"] =  -0.32782161
yeast17["comp"]["Q"] =  -0.17676485
yeast17["comp"]["A"] =   0.10844211
yeast17["comp"]["G"] =  -0.37594135
yeast17["comp"]["S"] =  -0.09627044
yeast17["comp"]["T"] =  -0.08533912
yeast17["comp"]["V"] =   0.43746326
yeast17["comp"]["M"] =   0.31182498
yeast17["comp"]["L"] =   0.53427787
yeast17["comp"]["I"] =   0.61465146
yeast17["comp"]["F"] =   0.52882600
yeast17["comp"]["Y"] =   0.45253658
yeast17["comp"]["W"] =   0.58693535
yeast17["comp"]["P"] =  -0.25880796

human25 = {}
human25["type"] = "glm"
human25["link"] = "sigmoid"
human25["size"] = 25
human25["intercept"] = 1.59844363
human25["comp"] = {}
human25["comp"]["C"] =  0.23501148 
human25["comp"]["D"] = -0.71078372 
human25["comp"]["E"] = -0.63879547 
human25["comp"]["K"] = -0.18827385 
human25["comp"]["R"] =  0.09256831  
human25["comp"]["H"] = -0.06207261  
human25["comp"]["N"] = -0.37443898  
human25["comp"]["Q"] = -0.30649578  
human25["comp"]["A"] = -0.14842730
human25["comp"]["G"] = -0.33494145
human25["comp"]["S"] = -0.29271972  
human25["comp"]["T"] = -0.17276016  
human25["comp"]["V"] =  0.24236702 
human25["comp"]["M"] =  0.23746524 
human25["comp"]["L"] =  0.44688517 
human25["comp"]["I"] =  0.50168095 
human25["comp"]["F"] =  0.62951898 
human25["comp"]["Y"] =  0.47749474 
human25["comp"]["W"] =  0.91842501  
human25["comp"]["P"] = -0.36074230  

human30 = {}
human30["type"] = "glm"
human30["link"] = "linear"
human30["size"] = 30
human30["intercept"] = 0.489921917
human30["comp"] = {}
human30["comp"]["E"] =  0.037084228
human30["comp"]["L"] = -0.043636559
human30["comp"]["F"] = -0.053407775
human30["comp"]["I"] = -0.051675763
human30["comp"]["D"] =  0.040739392
human30["comp"]["Y"] = -0.047630778
human30["comp"]["S"] =  0.009223921
human30["comp"]["W"] = -0.073252881
human30["comp"]["V"] = -0.030630140
human30["comp"]["P"] =  0.011450366
human30["comp"]["C"] = -0.032196916
human30["comp"]["R"] = -0.019694770
human30["comp"]["M"] = -0.030993800
human30["comp"]["Q"] =  0.011084338
human30["comp"]["G"] =  0.014319071
human30["comp"]["H"] = -0.009140511
human30["comp"]["N"] =  0.014679317
human30["comp"]["K"] =  0.0
human30["comp"]["A"] =  0.0
human30["comp"]["T"] =  0.0
# positions are one-based relative to CT, e.g. pos 1 is the CT/last/rightmost amino acid in sequence
human30["pssm"] = {}
human30["pssm"]["G1"] = -0.065740482
human30["pssm"]["A1"] = -0.052849574
human30["pssm"]["K1"] =  0.054402914
human30["pssm"]["R3"] = -0.049967679
human30["pssm"]["K2"] =  0.035361966
human30["pssm"]["A2"] = -0.032112630
human30["pssm"]["E1"] =  0.029693876
human30["pssm"]["L1"] =  0.026114745
human30["pssm"]["K5"] =  0.024547954
human30["pssm"]["A3"] = -0.024365990
human30["pssm"]["K4"] =  0.022683833
human30["pssm"]["V2"] = -0.027960791
human30["pssm"]["V3"] = -0.025205832
human30["pssm"]["R1"] = -0.018860502
human30["pssm"]["L5"] = -0.012282595
human30["pssm"]["D2"] =  0.011048007
human30["pssm"]["E2"] = -0.019685794
human30["pssm"]["G2"] = -0.019555017
human30["pssm"]["D3"] =  0.013405138
human30["pssm"]["T1"] = -0.021889475
human30["pair"] = {}
human30["pair"]["E2:E1"] = -0.230729238
human30["pair"]["G2:G1"] = -0.177261215
human30["pair"]["A2:A1"] = -0.088430778
human30["pair"]["G2:A1"] = -0.144060642
human30["pair"]["A4:A1"] = -0.073953252
human30["pair"]["R4:G1"] = -0.105240042
human30["pair"]["A2:G1"] = -0.108544265
human30["pair"]["R2:G1"] = -0.129695386
human30["pair"]["P2:G1"] = -0.121846803
human30["pair"]["A3:A1"] = -0.032161964
human30["pair"]["K2:G1"] = -0.144603640
human30["pair"]["R4:A1"] = -0.073407395
human30["pair"]["R3:E2"] = -0.051428389
human30["pair"]["R5:G2"] = -0.067951546
human30["pair"]["R3:P1"] = -0.077359054
human30["pair"]["K3:R1"] = -0.066072125
human30["triplet"] = {}
human30["triplet"]["P4:P2:A1"] = -0.283872981
human30["triplet"]["A4:A2:A1"] = -0.093776558

def calc_composition(seq, aa_param):
    """Calculate composition score of a sequence and correct for sequence length
    """
    if all([aa in aa_param["comp"].keys() for aa in seq]):
        score = np.sum([ aa_param["comp"][aa] for aa in seq]) * aa_param["size"]*1.0/len(seq)
    else:
        score = np.nan
    return(score)


def calc_glm(tile_list, model_param, verbose=0):
    """Calculate composition scores of a generalized linear model
    """
    score_list = []
    for tile in tile_list:
        comp = calc_composition(tile.seq, model_param)
        score_list.append( model["intercept"] + comp )
    score_list = np.array(score_list)

    if model_param["link"] == "linear":
        pass
    elif model_param["link"] == "sigmoid":
        score_list = 1/(1+np.exp(-score_list))
    else:
        raise ValueError("Unknown link function '%s'" % (model_param["link"]))
    return(score_list)


def calc_glm_ct(tile_list, model_param, include_ct=False, verbose=0):
    """Calculate composition and C-degron scores of a generalized linear model
    """
    # If model has CT PSSM terms, preprocess
    pssm_vec = None
    model_has_pssm = False
    if "pssm" in model_param.keys():
        pssm_terms = model_param["pssm"].keys()
        if len(pssm_terms) > 0:
            model_has_pssm = True
            pssm_aa_pos = [ (aa_pos[0],int(aa_pos[1:])) for aa_pos in pssm_terms]
            pssm_vec = np.array([ model_param["pssm"][term] for term in pssm_terms ])
    # If model has CT pair terms, preprocess
    pair_vec = None
    model_has_pair = False
    if "pair" in model_param.keys():
        pair_terms = model_param["pair"].keys()
        if len(pair_terms) > 0:
            model_has_pair = True
            pair_terms_split = [ term.split(':') for term in pair_terms ]
            pair_aa_pos = [ (aa_pos[0][0],aa_pos[1][0],int(aa_pos[0][1:]),int(aa_pos[1][1:])) for aa_pos in pair_terms_split ]
            pair_vec = np.array([ model_param["pair"][term] for term in pair_terms ])
    # If model has CT triplet terms, preprocess
    trip_vec = None
    model_has_trip = False
    if "triplet" in model_param.keys():
        trip_terms = model_param["triplet"].keys()
        if len(trip_terms) > 0:
            model_has_trip = True
            trip_terms_split = [ term.split(':') for term in trip_terms ]
            trip_aa_pos = [ (aa_pos[0][0],aa_pos[1][0],aa_pos[2][0],int(aa_pos[0][1:]),int(aa_pos[1][1:]),int(aa_pos[2][1:])) for aa_pos in trip_terms_split ]
            trip_vec = np.array([ model_param["triplet"][term] for term in trip_terms ])
    # evaluate sequences
    score_list = []
    ct_list = []
    for tile in tile_list:
        # if is_aa_one_nat(tile.seq):
        comp = calc_composition(tile.seq, model_param)
        pssm = 0.0
        pair = 0.0
        trip = 0.0
        if tile.is_ct():
            if model_has_pssm:
                pssm = np.sum(pssm_vec * np.array([ tile.seq[-pos]==aa for (aa,pos) in pssm_aa_pos ]))
            if model_has_pair:
                pair = np.sum(pair_vec * np.array([ tile.seq[-pos1]==aa1 and tile.seq[-pos2]==aa2 for (aa1,aa2,pos1,pos2) in pair_aa_pos ]))
            if model_has_trip:
                trip = np.sum(trip_vec * np.array([ tile.seq[-pos1]==aa1 and tile.seq[-pos2]==aa2 and tile.seq[-pos3]==aa3 for (aa1,aa2,aa3,pos1,pos2,pos3) in trip_aa_pos ]))
        score_list.append( model["intercept"] + comp + pssm + pair + trip)
        ct_list.append( pssm + pair + trip )
            
    score_list = np.array(score_list)
    ct_list = np.array(ct_list)
    if model_param["link"] == "linear":
        pass
    elif model_param["link"] == "sigmoid":
        score_list = 1/(1+np.exp(-score_list))
    else:
        raise ValueError("Unknown link function '%s'" % (model_param["link"]))
    return((score_list,ct_list))    


cnn2w1 = {}
cnn2w1['type'] = 'keras'
cnn2w1['arch'] = 'cnn2w'
cnn2w1['size'] = 30
cnn2w1['cnn_res'] = 1
cnn2w1['nf_int'] = 25
cnn2w1['ct_res'] = 5
cnn2w1['nf_ct'] = 25
cnn2w1['nn'] = 50
cnn2w1['act'] = 'elu'
cnn2w1['weight_file'] = './pap_weights/cnn2w1.keras'

cnn2w7 = {}
cnn2w7['type'] = 'keras'
cnn2w7['arch'] = 'cnn2w'
cnn2w7['size'] = 30
cnn2w7['cnn_res'] = 7
cnn2w7['nf_int'] = 25
cnn2w7['ct_res'] = 5
cnn2w7['nf_ct'] = 25
cnn2w7['nn'] = 50
cnn2w7['act'] = 'elu'
cnn2w7['weight_file'] = './pap_weights/cnn2w7.keras'

class KerasModel:
    """Abstract base class of Kreas models
    """
    aa2i = {'A':0, 'C':1, 'D':2, 'E':3, 'F':4, 'G':5, 'H':6, 'I':7, 'K':8, 'L':9,
            'M':10, 'N':11, 'P':12, 'Q':13, 'R':14, 'S':15, 'T':16, 'V':17, 'W':18, 'Y':19}
    def str2mat(self, s):
        n = len(s)
        v = np.zeros((n,len(self.aa2i)))
        for (i,aa) in enumerate(s):
            v[i,self.aa2i[aa] ] = 1
        return v
    def seqs2mat(self, seq_list):
        seqs_mat = np.array([self.str2mat(s) for s in seq_list])
        return(seqs_mat)
    def seqs2mat_len_sort(self, seq_list):
        ret_dict = {}
        for seq in seq_list:
            if not len(seq) in ret_dict.leys():
                ret_dict[len(seq)] = []
            ret_dict[len(seq)].append(self.str2mat(seq))
        for key in ret_dict.keys():
            ret_dict[key] = np.array(ret_dict[key])
        return(ret_dict)
    
class KerasModelCNN2W(KerasModel):
    """The PAP two-way convolutional neural network (CNN2w) model implemented in Keras
    """
    def __init__(self, model_param):
        from tensorflow import keras
        self.ct_res = model_param['ct_res']
        self.int_res = model_param['cnn_res']
        
        input_int  = keras.layers.Input(shape=(None, len(self.aa2i)), name='input_int', dtype="float32")
        internal   = keras.layers.Conv1D(model_param['nf_int'], model_param['cnn_res'], strides=1, padding='valid', name='int_cnn',
                                         activation=model_param['act'], trainable=False)(input_int)
        internal   = keras.layers.GlobalAveragePooling1D(data_format='channels_last', keepdims=False)(internal)
        internal   = keras.layers.Dense(model_param['nn'], name='int2', activation=model_param['act'], trainable=False)(internal)
        internal   = keras.layers.Dense(model_param['nn'], name='int3', activation=model_param['act'], trainable=False)(internal)
        output_int = keras.layers.Dense(1, name='latent_int', activation=None, trainable=False)(internal)
        
        input_ct   = keras.layers.Input(shape=(self.ct_res*len(self.aa2i),), name='input_ct', dtype="float32")
        cterm      = keras.layers.Dense(model_param['nf_ct'], name='ct1', activation=model_param['act'], trainable=False)(input_ct)
        cterm      = keras.layers.Dense(model_param['nn'], name='ct2', activation=model_param['act'], trainable=False)(cterm)
        cterm      = keras.layers.Dense(model_param['nn'], name='ct3', activation=model_param['act'], trainable=False)(cterm)
        output_ct  = keras.layers.Dense(1, name='latent_ct', activation=None, trainable=False)(cterm)

        # build the full model in order to load weights
        output_full = keras.layers.Add()([output_int,output_ct])
        model_full  = keras.models.Model(inputs=[input_int,input_ct], outputs=output_full)
        if not os.path.isfile(model_param['weight_file']):
            raise ValueError("Cannot find weights file '%s'" % (model_param['weight_file']))
        # model_full.load_weights(model_param['weight_file'], skip_mismatch=False, by_name=True, options=None)
        model_full.load_weights(model_param['weight_file'])
        
        # build models for each channel to be used later
        self.model_int = keras.models.Model(inputs=input_int, outputs=output_int)
        self.model_ct = keras.models.Model(inputs=input_ct, outputs=output_ct)

    def predict(self, tile_list, include_ct, verbose=0):
        tiles_nres = [ len(tile) for tile in tile_list ]
        tiles_is_ct = [ tile.is_ct() for tile in tile_list ]
        
        if len(set(tiles_nres)) == 1:
            # if all tiles have same length they may be evaluated in one call
            nres = tiles_nres[0]
            if nres < self.int_res:
                raise ValueError("Cannot evaluate tile of length %d with CNN layer that considers %d positions" % (prev_nres,self.int_res))
            seqs_mat = self.seqs2mat([ tile.seq for tile in tile_list ])
            scores_int = self.model_int.predict(seqs_mat, verbose=0).flatten()
            if include_ct and any(tiles_is_ct):
                seqs_mat_ct = seqs_mat[:,(nres-self.ct_res):nres,:]
                seqs_oh_ct = seqs_mat_ct.reshape(*seqs_mat_ct.shape[:-2], -1)
                scores_ct = self.model_ct.predict(seqs_oh_ct, verbose=0).flatten()
                scores_ct = scores_ct * np.array(tiles_is_ct)
            else:
                scores_ct = np.zeros(len(tile_list))
        else:
            # if tiles have different lengths, collect consequtive tiles of same lengths
            scores_int = []
            scores_ct = []
            tiles_same_len = []
            prev_nres = None
            calc_ct = False
            # signal for last evaluation
            tiles_nres.append(0)
            # count number of evaluations to give hint
            n_eval = 0; n_eval_warn = 1000
            for it in range(len(tile_list)+1):
                if tiles_nres[it] != prev_nres and not prev_nres is None:
                    n_eval += 1
                    # if tile has length different from the previous, evaluate
                    if prev_nres < self.int_res:
                        raise ValueError("Cannot evaluate tile of length %d with CNN layer that considers %d positions" % (prev_nres,self.int_res))
                    if verbose > 1:
                        print("    evaluate %d tiles of length %d" % (len(tiles_same_len),prev_nres), file=sys.stderr)
                    if verbose > 0 and n_eval % n_eval_warn == 0:
                        print("    doing keras evaluation number %d, consider ordering input tiles by length to improve speed" % (n_eval), file=sys.stderr)
                    seqs_mat = self.seqs2mat(tiles_same_len)
                    scores_int_same_len = self.model_int.predict(seqs_mat, verbose=0).flatten()
                    scores_int.extend(list(scores_int_same_len))
                    if calc_ct and include_ct:
                        if prev_nres < self.ct_res:
                            raise ValueError("Cannot evaluate CT term of tile length %d with model that has %d CT positions" % (prev_nres,self.ct_res))
                        seqs_mat_ct = seqs_mat[:,(prev_nres-self.ct_res):prev_nres,:]
                        seqs_oh_ct = seqs_mat_ct.reshape(*seqs_mat_ct.shape[:-2], -1)
                        scores_ct_same_len = self.model_ct.predict(seqs_oh_ct, verbose=0).flatten()
                        scores_ct.extend(list(scores_ct_same_len))
                    else:
                        scores_ct.extend([0.0]*len(tiles_same_len))
                    if it < len(tile_list):
                        # reset collection with current tile
                        tiles_same_len = [ tile_list[it].seq ]
                        prev_nres = tiles_nres[it]
                        calc_ct = tiles_is_ct[it]
                else:
                    # if tile has same length as the previous, continue collecting for faster evaluation
                    tiles_same_len.append( tile_list[it].seq )
                    prev_nres = tiles_nres[it]
                    calc_ct = calc_ct or tiles_is_ct[it]

            assert len(scores_int) == len(tile_list)
            assert len(scores_ct) == len(tile_list)
            scores_int = np.array(scores_int)
            scores_ct = np.array(scores_ct) * np.array(tiles_is_ct)
        scores_full = scores_int+scores_ct
        return((scores_full,scores_ct))

    
# Keep the model object as a global model to avoid loading it many times
keras_model = None
def calc_keras(tile_list, model_param, include_ct, verbose=0):
    global keras_model
    if not 'tensorflow' in sys.modules:
        if verbose > 0:
            print("Requested Keras-based model - importing tensorflow", file=sys.stderr)
        t1 = time.time()
        import tensorflow
        if verbose > 0:
            print("    done importing version %s, took %.0f seconds" % (tensorflow.__version__, time.time()-t1), file=sys.stderr)

    if keras_model is None:
        if model_param['arch'] == 'cnn2w':
            keras_model = KerasModelCNN2W(model_param)
        else:
            raise ValueError("Unknown keras architecture '%s'" % (model_param['arch']))
        
    (score_list,ct_list) = keras_model.predict(tile_list, include_ct, verbose=verbose)
    return (score_list,ct_list)


def calc_score(tile_list, model_param, include_ct=True, verbose=0):
    """Calculate scores using the appropriate model type
    """
    # call a score function based on 'type' key of model
    if model_param['type']=='glm' and include_ct is False:
        # simple and perhaps slightly faster special case of calc_glm_ct
        score_list = calc_glm(tile_list, model_param, verbose=verbose)
        ct_list = np.zeros(len(tile_list))
    elif model_param['type']=='glm':
        # generalised linear model
        (score_list,ct_list) = calc_glm_ct(tile_list, model_param, include_ct, verbose=verbose)
    elif model_param['type']=='keras':
        # neural network model implemented in Keras with tensorflow back-end
        (score_list,ct_list) = calc_keras(tile_list, model_param, include_ct, verbose=verbose)
    else:
        raise ValueError("Unknown model type '%s'" % (model_param['type']))
    return (score_list,ct_list)


class Tile:
    """Class that describes a peptide tile as a segment of a protein with a single (center) position representing the tile
    """
    def __init__(self, seq, begin=0, center_index=None, parent_id=None, parent_nres=None):
        self.seq = seq
        self.begin = begin
        self.end = begin+len(seq)
        if center_index is None:
            self.center_index = int( np.ceil( len(seq)/2.0 ) )
        else:
            self.center_index = center_index
        self.parent_id = parent_id
        self.parent_nres = parent_nres
        self.wt = ['']*len(self.seq)
        if not self.check():
            raise ValueError("Tile is inconsistent, check sequence length ans indices (begin,end and center)")
    def __len__(self):
        return(len(self.seq))
    def __eq__(self, other):
        is_eq = self.seq==other.seq and self.begin==other.begin
        if not self.parent_id is None:
            is_eq = is_eq and self.parent_id==other.parent_id
        return(is_eq)
    def __str__(self):
        s = "Tile %s, residue %d-%d (1-based first-last) with center at %s%d" % (self.seq,self.begin+1,self.end,self.seq[self.center_index],
                                                                                 self.center_index+1)
        n_mut = sum([wt!='' for wt in self.wt])
        if n_mut > 0:
            s = s+". Tile has %d mutation(s)" % (n_mut)
        if not self.parent_id is None:
            s = s+". Parent is '%s'" % (self.parent_id)
        return(s)
    def __repr__(self):
        return(self.__str__())
    def copy(self):
        return(Tile(self.seq, self.begin, center_index=self.center_index, parent_id=self.parent_id, parent_nres=self.parent_nres))
    def check(self):
        return(len(self.seq)==self.end-self.begin and self.center_index>=0 and self.center_index<len(self.seq))
    def is_nt(self):
        return(self.begin==0)    
    def is_ct(self):
        if self.parent_nres is None:
            return(None)
        else:
            return(self.end==self.parent_nres)
    def mutate(self,resi,new_aa):
        """Returns a copy of a tile with the given mutation"""
        if resi < 0 or resi >= len(self.seq):
            raise ValueError("Tile mutation residue %d out of range" % (resi))
        if new_aa == self.seq[resi] or not new_aa in aa_one_nat:
            raise ValueError("Cannot mutate tile %s%d to %s" % (self.seq[resi],resi,new_aa))
        if self.wt[resi] != '':
            raise ValueError("Tile position %d already mutated (%s -> %s)" % (resi, self.wt[resi], self.seq[resi]))
        mut_tile = self.copy()
        mut_tile.wt[resi] = mut_tile.seq[resi]
        mut_tile.seq = mut_tile.seq[:resi]+new_aa+mut_tile.seq[(resi+1):]
        return(mut_tile)
    def is_mutated(self):
        return(any([wt!='' for wt in self.wt]))


def tile_sequence(sequence, tile_size, terminal_tiles=False, seq_id=None):
    ## Example: protein length 8 and tiles length 5 or 6. 'X' is the central position,
    ## in case of even sized tiles, the central amino acid is right of the center
    ## tile   indices                length  left-of-cent  rigth-of-cent
    ##   CT               |5-6-X|    3       2             0
    ##   CT             |4-5-X-7|    4       2             1
    ##   Full         |3-4-X-6-7|    5       2             2
    ##   Full       |2-3-X-5-6|      5       2             2
    ##   Full     |1-2-X-4-5|        5       2             2
    ##   Full   |0-1-X-3-4|          5       2             2 
    ##   NT     |0-X-2-3|            4       1             2
    ##   NT     |X-1-2|              3       0             2
    ## Protein   0 1 2 3 4 5 6 7     8
    ##   NT     |X-1-2|              3       0             2
    ##   NT     |0-X-2-3|            4       1             2
    ##   NT     |0-1-X-3-4|          5       3             2
    ##   Full   |0-1-2-X-4-5|        6       3             2
    ##   Full     |1-2-3-X-5-6|      6       3             2
    ##   Full       |2-3-4-X-6-7|    6       3             2
    ##   CT           |3-4-5-X-7|    5       3             1
    ##   CT             |4-5-6-X|    4       3             0
    n_res = len(sequence)
    
    # in case of even sized tiles, the central amino acid is right of the center
    tile_left_of_center = int( np.floor( tile_size/2.0 ) )
    tile_right_of_center = int( np.ceil( tile_size/2.0 ) ) -1

    # tile indices in full sequence
    tiles = []
        
    # short NT tiles if requested
    if terminal_tiles:
        # include center aa and everthing to the right of the center
        for end_index in range(1+tile_right_of_center, tile_size):
            tiles.append( Tile(sequence[0:end_index], 0, center_index=end_index-tile_right_of_center-1,
                               parent_nres=len(sequence), parent_id=seq_id) )
        
    # full tiles
    for begin_index in range(n_res-tile_size+1):
        tiles.append( Tile(sequence[begin_index:(begin_index+tile_size)], begin_index, center_index=tile_left_of_center,
                           parent_nres=len(sequence), parent_id=seq_id) )

    # short CT tiles if requested
    if terminal_tiles:
        # from where full tiles ended and until a tile is only the central position and positions left of this
        for begin_index in range(n_res-tile_size+1, n_res-tile_left_of_center):
            tiles.append( Tile(sequence[begin_index:n_res], begin_index, center_index=tile_left_of_center,
                               parent_nres=len(sequence), parent_id=seq_id) )
    return(tiles)


def saturation_mut_expand(wt_tiles):
    """Construct a saturation mutagenesis library where only the central position of tiles are mutated
    """
    tiles = []
            
    # make new lists that are expanded with variants
    for wt_tile in wt_tiles:
        # always append wt first in tile_list - variants in order of
        tiles.append(wt_tile.copy())
        wt_aa = wt_tile.seq[wt_tile.center_index]

        # mutants to test
        aa_muts = aa_one_nat.replace(wt_aa,'')     
        
        # append mutants - avoid wt here
        for aa_mut in aa_muts:
            tiles.append(wt_tile.mutate(wt_tile.center_index, aa_mut))
            
    return(tiles)


def guess_uniprot_from_filename(filename):
    """Robust method thar attempts to guess the uniprot entry from a prism-like file name
    """
    if filename.count('_') >=2:
        # assume file name is prism_method_uniprot[-_]*.txt and that uniprot does not contain - or _
        words = filename.split("_",2)
        if words[0] == 'prism':
            uniprot = words[2].split('_',1)[0].split('-',1)[0]
            return(uniprot)
    elif '_' in filename:
        # assume file name is uniprot_some_name
        words = filename.split("_",1)
        return(words[0])
    return("unknown")


if __name__ == "__main__":
    # argument parser
    arg_parser = argparse.ArgumentParser(description="Predict degradation/abundance scores of proteins and peptides assuming these are solvent/PQC exposed in the cell")
    # keyword arguments
    models_str = "Model to use\n"+"cnn2w1: 2-way convolutional neural network\n"+"    Default tile size 30\n"+"    Predicts abundance scores with zero most unstable and un-tagged GFP ~0.6\n"
    # human30: Linear regression model
    #     Default tile size 30
    #     Predicts abundance scores with zero most unstable and un-tagged GFP ~0.6
    # human25: Logistic regression model trained on the 
    #     Default tile size 25
    #     Predicts degron scores with zero stable and one most unstable
    # yeast17: QCDPred trained on yeast-based screen
    #     Default tile size 17
    #     Predicts degron scores with zero stable and one most unstable
    arg_parser.add_argument("-m", "--model", choices=["yeast17", "human25", "human30", "cnn2w1", "cnn2w7"], default="cnn2w1",
                            help="Model to use")
    # when using -s 0 with keras-based models, it is highly recommended to order input sequences by length (keras uses arrays for input)
    arg_parser.add_argument("-s", "--tile-size", type=int, default=None,
                            help="Number of amino acids per tile. Set to zero to keep sequence lengths as inputed. Default depends on the model used")
    arg_parser.add_argument("-c","--ignore-ct-term", default=False, action='store_true',
                            help="By default, C-degron terms (if available in model) are applied to the last tile of each sequence")
    arg_parser.add_argument("-t", "--term-method", choices=["full-tiles", "full-coverage"], default="full-tiles",
                            help="Method to handle terminals. Either considers tiles of full length only (default) or cover the entire sequence "+
                            "by evaluating tiles down to half length. If length is even the central amino acid is the one to the right of the center")
    arg_parser.add_argument('--saturation-mut', dest='saturation', default=False, action='store_true',
                            help="Perform satuation mutagenesis of the central amino acid in each tile. The output will typically have 20 lines per tile"+
                            "with the central amino acid mutated. If tile length is even the central amino acid is the one to the right of the center")
    arg_parser.add_argument('--skip-nonnatural', default=False, action='store_true',
                            help="Ignore tiles that contains other then the 20 standard amino acids")
    arg_parser.add_argument("--output-format", choices=["std", "prism"], default="std",
                            help="Output to stdout or write prism files in directory 'prism_pap'")
    arg_parser.add_argument('-v', '--verbose', action='count', default=0)
    # positional argument
    arg_parser.add_argument("seq_input", nargs="+", metavar="SEQ",
                            help="Sequence(s) in text or FASTA format. More inputs should by whitespace separated")
    args = arg_parser.parse_args()

    # set model parameters
    model = cnn2w1
    if args.model == "cnn2w7":
        model = cnn2w7
    elif args.model == "human30":
        model = human30
    elif args.model == "human25":
        model = human25
    elif args.model == "yeast17":
        model = yeast17

    # set tile size
    tile_size = None
    if args.tile_size is None:
        tile_size = model['size']
    else:
        tile_size = args.tile_size
        
    # if writing prism files is requested, make a output directory if not already present
    prism_dir = None
    if args.output_format == "prism":
        prism_dir = "prism_pap"
        if not args.saturation:
            print("WARNING: Cannot write prism formatted output without variants (use option --saturation)", file=sys.stderr)
        elif not os.path.isdir(prism_dir):
            os.mkdir(prism_dir)
        
    # extra_aa="XU" 
    extra_aa=""

    # read all sequences into memory
    input_list = []
    for seq_input in args.seq_input:
        if is_aa_one_nat(seq_input.upper(), extra_aa):
            input_list.append(("seq%05d" % (len(input_list)), seq_input.upper()))
        elif seq_input[-4:].lower() in [".seq",".fas"] or seq_input[-6:].lower() == ".fasta":
            if not os.path.isfile(seq_input):
                print("ERROR: Cannot find file %s" % (seq_input), file=sys.stderr)
                sys.exit(2)
            input_list.extend(read_fasta(seq_input, extra_aa=extra_aa, skip_nonnatural=args.skip_nonnatural, verbose=args.verbose))
        else:
            print("ERROR: Argument %s seems to be neither a protein sequence nor a FASTA file (.seq, .fas or .fasta)" % (seq_input), file=sys.stderr)
            sys.exit(2)
    if args.verbose > 0:
        print("Loaded %d input sequences" % (len(input_list)), file=sys.stderr)

    # check if names of input sequences are unique
    name_list = np.array([ elem[0] for elem in input_list ])
    name_uniq_list = np.unique(name_list, return_counts=True)
    if len(name_uniq_list[0]) != len(name_list):
        if args.output_format == "prism":
            # for prism-format output, sequence names need to be unique
            n_changed_names = 0
            for (name,count) in zip(name_uniq_list[0],name_uniq_list[1]):
                if count > 1:
                    indices = np.where(name_list==name)[0]
                    for (ic,iname) in enumerate(indices):
                        assert input_list[iname][0] == name
                        n_changed_names += 1
                        new_name = "%s_%02d" % (name,ic)
                        input_list[iname] = (new_name,input_list[iname][1])
                        if args.verbose > 1:
                            print("Changing name of input sequence  %s  ->  %s" % (name,new_name))
            if args.verbose > 0:
                print("Changed %d protein names to make unique file names" % (n_changed_names), file=sys.stderr)
            # update name_list
            name_list = np.array([ elem[0] for elem in input_list ])
        else:
            # for other output, give a warning
            print("WARNING: Names of input sequences are not unique", file=sys.stderr)    

    # process one input sequence at a time and collect into batches of same-sized tiles
    batch_size = 1e5
    batch_tiles = []
    for si in range(len(input_list)):
        (aaseq_name,aaseq) = input_list[si]

        # always skip short sequences
        n_res = len(aaseq)
        if n_res < tile_size:
            if args.verbose > 0:
                print("Skip sequence '%s' length %d < tile size %d" % (aaseq_name,n_res,tile_size), file=sys.stderr)
            continue
        
        # make tiles from sequence
        if tile_size > 0:
            tiles = tile_sequence(aaseq, tile_size, terminal_tiles=args.term_method=="full-coverage", seq_id=aaseq_name)
        else:
            tiles = [ Tile(aaseq, parent_nres=len(aaseq), parent_id=aaseq_name) ]

        # if saturation mutagenesis is requested, expand tile list to include mutants
        if args.saturation:
            tiles = saturation_mut_expand(tiles)
            
        batch_tiles.extend(tiles)
        
        # Evaluate and output if batch is full
        if len(batch_tiles) >= batch_size or si+1 >= len(input_list):
            # Calculate scores of all tiles
            time_begin = time.time()
            (tile_degron_score,tile_ct_score) = calc_score(batch_tiles, model, include_ct=(not args.ignore_ct_term), verbose=args.verbose)
            if args.verbose > 0:
                timing = time.time() - time_begin
                print("Evaluated %d tiles in %.2f seconds, %.4f s/tile" % (len(batch_tiles), timing, timing/len(batch_tiles)), file=sys.stderr)
        
            # Output
            if not args.saturation:
                if args.output_format == "std":
                    for (tile,score,ct_score) in zip(batch_tiles,tile_degron_score,tile_ct_score):
                        print("%18s  %30s  %8.5f  %8.5f  %4d  %4d  %1s  %4d" % (tile.parent_id, tile.seq, ct_score, score, tile.begin+1, tile.end,
                                                                                tile.seq[tile.center_index], tile.begin+tile.center_index+1))
                elif args.output_format == "prism":
                    raise NotImplementedError("Cannot do prism files of residues yet (only variant data files using --saturation-mut)")
            else:
                # also output WT AA and score difference for sasturation setting
                wt_score = None
                if args.output_format == "std":
                    for (tile,score,ct_score) in zip(batch_tiles,tile_degron_score,tile_ct_score):
                        wt_aa = None
                        if not tile.is_mutated():
                            # assume WT is always first in a sequence of mutants
                            wt_score = score
                            wt_aa = tile.seq[tile.center_index]
                        else:
                            wt_aa = tile.wt[tile.center_index]
                            assert tile.seq[tile.center_index] != wt_aa
                        print("%18s  %30s  %8.5f  %8.5f  %4d  %4d  %1s  %4d  %1s  %8.5f" % (tile.parent_id, tile.seq, ct_score, score, tile.begin+1, tile.end,
                                                                                            tile.seq[tile.center_index], tile.begin+tile.center_index+1,
                                                                                            wt_aa, score-wt_score))
                elif args.output_format == "prism":
                    # make output for each protein in batch
                    batch_names = np.array([tile.parent_id for tile in batch_tiles])
                    for prot_name in np.unique(batch_names):
                        # indices in batch lists that belongs to protein
                        prot_indices = np.where(batch_names==prot_name)[0]
                        # find protein name among input sequences
                        ii = np.where(np.array([elem[0] for elem in input_list])==prot_name)[0]
                        assert len(ii)==1
                        prot_seq = input_list[ii[0]][1]
                        prism_filename = "%s/prism_pap_%s.txt" % (prism_dir,prot_name)
                        with open(prism_filename, "wt") as prismfile:
                            prismfile.write("# --------------------\n")
                            prismfile.write("# version: 1\n")
                            prismfile.write("# protein:\n")
                            prismfile.write("#     name: %s\n" % (prot_name))
                            uniprot = guess_uniprot_from_filename(prot_name)
                            prismfile.write("#     uniprot: %s\n" % (uniprot))
                            prismfile.write("#     sequence: %s\n" % (prot_seq))
                            prismfile.write("# pap:\n")
                            prismfile.write("#     model: %s\n" % (args.model))
                            prismfile.write("#     tile_size: %d\n" % (tile_size))
                            prismfile.write("#     version: %s\n" % (__version__))
                            prismfile.write("# columns:\n")
                            prismfile.write("#     resi_first: First residue index of tile (one-based)\n")
                            prismfile.write("#     resi_last: Last residue index included in tile (one-based)\n")
                            prismfile.write("#     score_var_ct: Latent C-degron term as applied\n")
                            prismfile.write("#     score_var: Score of variant tile being a degron\n")
                            prismfile.write("#     delta_score: Score difference of variant - wt of tile being a degron\n")
                            prismfile.write("# --------------------\n")
                            prismfile.write("# \n")
                            prismfile.write("%7s %10s %10s %12s %12s %12s\n" % ("variant","resi_first","resi_last","score_var_ct","score_var","delta_score"))
                            for i in prot_indices:
                                cent_resi = batch_tiles[i].center_index
                                if not batch_tiles[i].is_mutated():
                                    # assume WT is always first in a sequence of mutants
                                    wt_score = tile_degron_score[i]
                                    variant = "%s%d=" % (batch_tiles[i].seq[cent_resi], batch_tiles[i].begin+cent_resi+1)
                                else:
                                    variant = "%s%d%s" % (batch_tiles[i].wt[cent_resi], batch_tiles[i].begin+cent_resi+1, batch_tiles[i].seq[cent_resi])
                                # if tile_degron_score[i] is np.nan or wt_score is np.nan:
                                #     continue
                                prismfile.write("%-7s %10d %10d %12.5f %12.5f %12.5f\n" % (variant, batch_tiles[i].begin+1, batch_tiles[i].end, tile_ct_score[i],
                                                                                           tile_degron_score[i], tile_degron_score[i]-wt_score))
                    if args.verbose > 0:
                        print("Wrote %d prism-formatted files in directory %s" % (len(np.unique(batch_names)),prism_dir), file=sys.stderr)
            # reset to make new batch
            batch_tiles = []

