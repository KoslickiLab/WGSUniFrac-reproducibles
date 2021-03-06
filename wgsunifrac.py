import logging
from collections import defaultdict
from ete3 import NCBITaxa
from ete3 import TreeNode
import copy
import numpy as np
import dendropy
from biom import Table
import pandas as pd
import random
import seaborn as sns
from scipy.stats import halfnorm
from skbio import DistanceMatrix #to install: pip install scikit-bio
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
import re
import itertools as it
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('TkAgg') #has to be enabled for plots to be shown
from skbio.stats.ordination import pcoa
ncbi = NCBITaxa()
import os

#classes
class Prediction:
    def __init__(self):
        pass

    @property
    def rank(self):
        return self.__rank

    @property
    def taxid(self):
        return self.__taxid

    @property
    def percentage(self):
        return self.__percentage

    @property
    def taxpath(self):
        return self.__taxpath

    @property
    def taxpathsn(self):
        return self.__taxpathsn

    @rank.setter
    def rank(self, rank):
        self.__rank = rank

    @taxid.setter
    def taxid(self, taxid):
        self.__taxid = taxid

    @percentage.setter
    def percentage(self, percentage):
        self.__percentage = percentage

    @taxpath.setter
    def taxpath(self, taxpath):
        self.__taxpath = taxpath

    @taxpathsn.setter
    def taxpathsn(self, taxpathsn):
        self.__taxpathsn = taxpathsn

    def get_dict(self):
        return self.__dict__

    def get_pretty_dict(self):
        return {property.split("_")[3]: value for property, value in self.__dict__.items()}

    def get_metadata(self):
        return {'rank': self.__rank, 'taxpath': self.__taxpath, 'taxpathsn': self.__taxpathsn}


#profiling
class Profile(object):
    def __init__(self, sample_metadata=None, profile=None, branch_length_fun=lambda x: 1 / x):
        self.sample_metadata = sample_metadata
        self.profile = profile
        self._data = dict()
        # Stick in the root node just to make sure everything is consistent
        self._data["-1"] = dict()
        self._data["-1"]["rank"] = None
        self._data["-1"]["tax_path"] = list()
        self._data["-1"]["tax_path_sn"] = list()
        self._data["-1"]["abundance"] = 0
        self._data["-1"]["descendants"] = list()
        self._header = list()
        self._tax_id_pos = None
        self._rank_pos = None
        self._tax_path_pos = None
        self._tax_path_sn_pos = None
        self._abundance_pos = None
        self._eps = .0000000000000001  # This is to act like zero, ignore any lines with abundance below this quantity
        self._all_keys = ["-1"]
        self._merged_flag = False
        self.root_len = 1  # the length you want between the "root" of "-1" and the superkingdom level (eg. Bacteria)
        self.branch_len_func = branch_length_fun  # Given a node n at depth d in the tree, branch_len_func(d)
        # is how long you want the branch length between n and ancestor(n) to be
        self._data["-1"]["branch_length"] = self.root_len
        self.parse_file()  # TODO: this sets all the branch lengths to 1 currentlclass Profile(object):


    def parse_file(self):
        _data = self._data
        _all_keys = self._all_keys
        _header = self._header
        for k, v in self.sample_metadata.items():
            _header.append('{}:{}'.format(k, v))

        # populate all the correct keys
        for prediction in self.profile:
            _all_keys.append(prediction.taxid.strip())

        # crawl over all profiles tax_path and create the ancestors and descendants list
        for prediction in self.profile:
            tax_id = prediction.taxid.strip()
            tax_path = prediction.taxpath.strip().split("|")  # this will be a list, join up late
            if tax_id not in _data:
                _data[tax_id] = dict()
            else:
                raise Exception(f"Improperly formatted profile: row starting with {tax_id} shows up more than once")
            _data[tax_id]["tax_path"] = tax_path

            # populate abundance
            _data[tax_id]["abundance"] = prediction.percentage

            # populate tax path sn
            if not (prediction.taxpathsn is None):  # might not be present
                _data[tax_id]["tax_path_sn"] = prediction.taxpathsn.strip().split(
                    "|")  # this will be a list, join up later

            # populate the rank
            _data[tax_id]["rank"] = prediction.rank.strip()

            # populate the branch length
            _data[tax_id]["branch_length"] = self.tax_path_to_branch_len(tax_path, self.branch_len_func, self.root_len)

            # Find the ancestors
            if len(tax_path) <= 1:  # note, due to the format, we will never run into the case tax_path == []
                _data[tax_id]["ancestor"] = "-1"  # no ancestor, it's a root
            else:  # go from the bottom up, looking for an ancestor that is an acceptable key
                ancestor = "-1"  # this is the default
                tax_path_rev = tax_path[::-1]
                for potential_ancestor in tax_path_rev:
                    if potential_ancestor != tax_id and potential_ancestor in _all_keys:
                        ancestor = potential_ancestor
                        break  # you found the ancestor, so can quit looking
                _data[tax_id]["ancestor"] = ancestor

            # Create a placeholder descendant key initialized to [], just so each tax_id has a descendant key associated to it
            if "descendants" not in _data[tax_id]:  # if this tax_id doesn't have a descendant list,
                _data[tax_id]["descendants"] = list()  # initialize to empty list

        self._add_descendants()
        self._delete_missing()  # make sure there aren't any missing internal nodes

    def _add_descendants(self):
        """
        Idea here is to look at all the ancestors of each key, and make the key the descendant of that ancestor
        Returns
        -------
        None: modifies Profile in place
        """
        _data = self._data
        _all_keys = self._all_keys
        for prediction in self.profile:
            tax_id = prediction.taxid.strip()  # the tax ID we are looking at
            ancestor = _data[tax_id]['ancestor']  # the tax ID's ancestor
            if tax_id not in _data[ancestor]['descendants']:
                _data[ancestor]['descendants'].append(
                    tax_id)  # so make the tax ID we're looking at the descendant of the ancestor

    def _delete_missing(self):
        """
        Deletes from the descendants all those taxids that aren't keys in the profile (i.e. there is no line that starts with that taxID)
        Returns
        -------
        none: modifies Profile in place
        """
        for key in self._data:
            clean_descendants = []
            for descendant in self._data[key]["descendants"]:
                if descendant in self._all_keys:  # if it's one of the taxids that the line starts with, add it
                    clean_descendants.append(descendant)
                else:
                    pass  # don't include the taxids that aren't actually in the final tax tree
            self._data[key]["descendants"] = clean_descendants
        return

    def write_file(self, out_file_name=None):
        if out_file_name is None:
            raise Exception
        _data = self._data
        keys = _data.keys()
        # This will be annoying to keep things in order...
        # Let's iterate on the length of the tax_path since we know that will be in there
        tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
        fid = open(out_file_name, 'w')
        # Write the header
        for head in self._header:
            fid.write("%s\n" % head)

        # Loop over length of tax_path and write data
        # always make the output tax_id, rank, tax_path, tax_path_sn, abundance in that order
        for path_length in range(1, tax_path_lengths + 1):
            for key in keys:
                if len(_data[key]["tax_path"]) == path_length and _data[key]["abundance"] > self._eps:
                    line_data = _data[key]
                    fid.write("%s\t" % key)
                    if self._rank_pos is not None:
                        fid.write("%s\t" % line_data["rank"])
                    fid.write("%s\t" % "|".join(line_data["tax_path"]))
                    if self._tax_path_sn_pos is not None:
                        fid.write("%s\t" % "|".join(line_data["tax_path_sn"]))
                    fid.write("%f\n" % line_data["abundance"])
        fid.close()
        return

    def threshold(self, threshold=None):
        if threshold is None:
            raise Exception
        _data = self._data
        keys = _data.keys()
        for key in keys:
            if _data[key]["abundance"] < threshold:
                _data[key]["abundance"] = 0
        return

    def _subtract_down(self):
        # helper function to push all the weights up by subtracting
        # NOTE: when subtracting, need to start at root and go down
        # NOTE: when adding, need to start at leaves and go up
        _data = self._data
        keys = _data.keys()
        # This will be annoying to keep things in order...
        # Let's iterate on the length of the tax_path since we know that will be in there
        tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
        for path_length in range(1, tax_path_lengths):  # eg tax_path_lengths = 5, use 1,2,3,4 since we stop at leaves
            for key in keys:
                if len(_data[key]["tax_path"]) == path_length:
                    descendants = _data[key]["descendants"]  # get all descendants
                    for descendant in descendants:
                        _data[key]["abundance"] -= _data[descendant]["abundance"]  # subtract the descendants abundance

    def _add_up(self):
        # helper function to push all the weights up by subtracting
        # NOTE: when subtracting, need to start at root and go down
        # NOTE: when adding, need to start at leaves and go up
        _data = self._data
        keys = _data.keys()
        # This will be annoying to keep things in order...
        # Let's iterate on the length of the tax_path since we know that will be in there
        tax_path_lengths = max([len(_data[key]["tax_path"]) for key in keys])
        for path_length in range(tax_path_lengths, 1,
                                 -1):  # eg tax_path_lengths = 5, use 5,4,3,2, since we stop at roots
            for key in keys:
                if len(_data[key]["tax_path"]) == path_length:
                    ancestor = _data[key]["ancestor"]
                    if ancestor in _data:  # don't do anything if this is a/the root node
                        _data[ancestor]["abundance"] += _data[key]["abundance"]  # add the descendants abundance

    def normalize(self):
        # Need to really push it up while subtracting, then normalize, then push up wile adding
        # self._push_up(operation="subtract")
        self._subtract_down()
        _data = self._data
        keys = _data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                _data[key]["abundance"] /= total_abundance
                _data[key]["abundance"] *= 100  # make back into a percentage
        # self._push_up(operation="add")
        self._add_up()
        return

    def merge(self, other):
        # Warning: not checking for taxonomic consistency
        if not isinstance(other, Profile):
            print("Only works with other Profiles")
            raise Exception
        if self._merged_flag is False:
            self._header.insert(0, "# This is a merged file, ignore files in headers below")
            self._merged_flag = True
        _data = self._data
        _other_data = other._data
        other_keys = _other_data.keys()
        for key in other_keys:
            if key in _data:
                _data[key]["abundance"] += _other_data[key]["abundance"]  # if already in there, add abundances
            else:
                _data[key] = copy.copy(_other_data[key])  # otherwise use the whole thing

    @staticmethod
    def tax_path_to_branch_len(tax_path, func, root_len=1):
        """
        This function modifies the branch lengths based on the input tax_path.
        intent is: ["2", "", "123", "456"] would result in a branch length of func(4)
        Parameters
        ----------
        tax_path : a list of strings (tax ID's)
        func : a function whose argument is the depth in the tree of a tax ID, and whose output is the branch length
               from the tax ID to its ancestor.
        root_len : how long you want the root of the tree "-1" to be to the descendants (eg. "-1" -> "Bacteria")
        Returns
        -------
        float
        """
        # eg. "-1" -> "Bacteria" should have a branch length of root_len
        if not tax_path:
            return root_len
        else:
            depth_in_tree = len(tax_path)  # this takes into account that the tax_path doesn't include the root of "-1"
            return func(depth_in_tree)

    def make_unifrac_input_and_normalize(self, other):
        if not isinstance(other, Profile):
            raise Exception
        _data = self._data
        _other_data = other._data

        _data_keys = _data.keys()
        tax_path_lengths1 = max([len(_data[key]["tax_path"]) for key in _data_keys])
        _other_data_keys = _other_data.keys()
        tax_path_lengths2 = max([len(_other_data[key]["tax_path"]) for key in _other_data_keys])
        tax_path_lengths = max(tax_path_lengths1, tax_path_lengths2)
        all_keys = set(_data_keys)
        all_keys.update(_other_data_keys)  # all the taxID's in the union of self and other profile
        nodes_in_order = []
        for path_length in range(tax_path_lengths, 0, -1):
            for key in all_keys:
                if key in _data:
                    if len(_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
                elif key in _other_data:
                    if len(_other_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
        # Make the graph
        # Put the root at the very end
        if '-1' in nodes_in_order:
            nodes_in_order.pop(nodes_in_order.index('-1'))
            nodes_in_order.append('-1')
        else:
            nodes_in_order.append('-1')
        Tint = dict()
        lint = dict()
        for key in nodes_in_order:
            if key in _data:
                if "ancestor" in _data[key]:  # If ancestor is not in there, then it's an ancestor
                    ancestor = _data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _data[key]["branch_length"]
            elif key in _other_data:
                if "ancestor" in _other_data[key]:
                    ancestor = _other_data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _other_data[key]["branch_length"]
        nodes_to_index = dict(
            zip(nodes_in_order, range(len(nodes_in_order))))  # maps '45202.15' -> 0 (i.e taxID to integer index)

        # Now need to change over to the integer-based indexing
        Tint2 = dict()
        lint2 = dict()
        nodes_in_order2 = []
        for key in nodes_in_order:
            if key in Tint:
                ancestor = Tint[key]
                Tint2[nodes_to_index[key]] = nodes_to_index[ancestor]
                if (key, ancestor) in lint:
                    lint2[nodes_to_index[key], nodes_to_index[ancestor]] = lint[key, ancestor]
            nodes_in_order2.append(nodes_to_index[key])

        # Next make the probability distributions
        # Would be nice if I could find a non-destructive way to subtract up and normalize

        # Do it for P
        self._subtract_down()
        keys = _data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                _data[key]["abundance"] /= total_abundance  # Should be a fraction, summing to 1
        P = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _data:
                P[key_ind] = _data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        for key in keys:
            if total_abundance > 0:
                _data[key]["abundance"] *= 100
        self._add_up()

        # Next do for Q
        other._subtract_down()
        keys = _other_data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _other_data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                _other_data[key]["abundance"] /= total_abundance  # should be a fraction, summing to 1
        Q = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _other_data:
                Q[key_ind] = _other_data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        for key in keys:
            if total_abundance > 0:
                _other_data[key]["abundance"] *= 100
        other._add_up()

        return Tint2, lint2, nodes_in_order2, nodes_to_index, P, Q

    def make_unifrac_input_no_normalize(self, other):
        if not isinstance(other, Profile):
            raise Exception
        _data = self._data
        _other_data = other._data

        _data_keys = _data.keys()
        tax_path_lengths1 = max([len(_data[key]["tax_path"]) for key in _data_keys])
        _other_data_keys = _other_data.keys()
        tax_path_lengths2 = max([len(_other_data[key]["tax_path"]) for key in _other_data_keys])
        tax_path_lengths = max(tax_path_lengths1, tax_path_lengths2)
        all_keys = set(_data_keys)
        all_keys.update(_other_data_keys)  # all the taxID's in the union of self and other profile
        nodes_in_order = []
        for path_length in range(tax_path_lengths, 0, -1):
            for key in all_keys:
                if key in _data:
                    if len(_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
                elif key in _other_data:
                    if len(_other_data[key]["tax_path"]) == path_length:
                        if key not in nodes_in_order:
                            nodes_in_order.append(key)
        # Make the graph
        # Put the root at the very end
        if '-1' in nodes_in_order:
            nodes_in_order.pop(nodes_in_order.index('-1'))
            nodes_in_order.append('-1')
        else:
            nodes_in_order.append('-1')
        Tint = dict()
        lint = dict()
        for key in nodes_in_order:
            if key in _data:
                if "ancestor" in _data[key]:  # If ancestor is not in there, then it's an ancestor
                    ancestor = _data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _data[key]["branch_length"]
            elif key in _other_data:
                if "ancestor" in _other_data[key]:
                    ancestor = _other_data[key]["ancestor"]
                    Tint[key] = ancestor
                    lint[key, ancestor] = _other_data[key]["branch_length"]
        nodes_to_index = dict(
            zip(nodes_in_order, range(len(nodes_in_order))))  # maps '45202.15' -> 0 (i.e taxID to integer index)

        # Now need to change over to the integer-based indexing
        Tint2 = dict()
        lint2 = dict()
        nodes_in_order2 = []
        for key in nodes_in_order:
            if key in Tint:
                ancestor = Tint[key]
                Tint2[nodes_to_index[key]] = nodes_to_index[ancestor]
                if (key, ancestor) in lint:
                    lint2[nodes_to_index[key], nodes_to_index[ancestor]] = lint[key, ancestor]
            nodes_in_order2.append(nodes_to_index[key])

        # Next make the probability distributions
        # Would be nice if I could find a non-destructive way to subtract up and normalize

        # Do it for P
        self._subtract_down()
        keys = _data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                # _data[key]["abundance"] /= total_abundance  # Should be a fraction, summing to 1
                pass
        P = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _data:
                P[key_ind] = _data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        # for key in keys:
        #    if total_abundance > 0:
        #        _data[key]["abundance"] *= 100
        self._add_up()

        # Next do for Q
        other._subtract_down()
        keys = _other_data.keys()
        total_abundance = 0
        for key in keys:
            total_abundance += _other_data[key]["abundance"]
        # print(total_abundance)
        for key in keys:
            if total_abundance > 0:
                # _other_data[key]["abundance"] /= total_abundance  # should be a fraction, summing to 1
                pass
        Q = np.zeros(len(nodes_in_order))
        for key_ind in range(len(nodes_in_order)):
            key = nodes_in_order[key_ind]
            if key in _other_data:
                Q[key_ind] = _other_data[key]["abundance"]

        # Make back into percentages and add the mass back up (effectively normalizing the vector)
        # for key in keys:
        #    if total_abundance > 0:
        #        _other_data[key]["abundance"] *= 100
        other._add_up()

        return Tint2, lint2, nodes_in_order2, nodes_to_index, P / 100., Q / 100.

class Node:
    def __init__(self, name, abundance=0, tax=''):
        self.name = name
        self.tax = tax
        self.abundance = abundance

# loading file
def open_profile_from_tsv(file_path, normalize):
    header = {}
    column_name_to_index = {}
    profile = []
    samples_list = []
    predictions_dict = {}
    reading_data = False
    got_column_indices = False

    with open(file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0 or line.startswith("#"):
                continue
            line = line.rstrip('\n')

            # parse header with column indices
            if line.startswith("@@"):
                for index, column_name in enumerate(line[2:].split('\t')):
                    column_name_to_index[column_name] = index
                index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn = get_column_indices(column_name_to_index)
                got_column_indices = True
                reading_data = False
                continue

            # parse header with metadata
            if line.startswith("@"):
                # if last line contained sample data and new header starts, store profile for sample
                if reading_data:
                    if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
                        if len(profile) > 0:
                            samples_list.append((header['SAMPLEID'], header, profile))
                            profile = []
                            predictions_dict = {}
                    else:
                        logging.getLogger('opal').critical(
                            "Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.\n".format(
                                file_path))
                        raise RuntimeError
                    header = {}
                reading_data = False
                got_column_indices = False
                key, value = line[1:].split(':', 1)
                header[key.upper()] = value.strip()
                continue

            if not got_column_indices:
                logging.getLogger('opal').critical(
                    "Header line starting with @@ in file {} is missing or at wrong position.\n".format(file_path))
                raise RuntimeError

            reading_data = True
            row_data = line.split('\t')

            taxid = row_data[index_taxid]
            # if there is already a prediction for taxon, only sum abundance
            if taxid in predictions_dict:
                prediction = predictions_dict[taxid]
                prediction.percentage += float(row_data[index_percentage])
            else:
                if float(row_data[index_percentage]) == .0:
                    continue
                prediction = Prediction()
                predictions_dict[taxid] = prediction
                prediction.taxid = row_data[index_taxid]
                prediction.rank = row_data[index_rank]
                prediction.percentage = float(row_data[index_percentage])
                prediction.taxpath = row_data[index_taxpath]
                if isinstance(index_taxpathsn, int):
                    prediction.taxpathsn = row_data[index_taxpathsn]
                else:
                    prediction.taxpathsn = None
                profile.append(prediction)

    # store profile for last sample
    if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
        if reading_data and len(profile) > 0:
            samples_list.append((header['SAMPLEID'], header, profile))
    else:
        logging.getLogger('opal').critical(
            "Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.\n".format(
                file_path))
        raise RuntimeError

    if normalize:
        normalize_samples(samples_list)

    return samples_list

def get_column_indices(column_name_to_index):
    if "TAXID" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("TAXID"))
        raise RuntimeError
    if "RANK" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("RANK"))
        raise RuntimeError
    if "PERCENTAGE" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("PERCENTAGE"))
        raise RuntimeError
    if "TAXPATH" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("TAXPATH"))
        raise RuntimeError
    index_taxid = column_name_to_index["TAXID"]
    index_rank = column_name_to_index["RANK"]
    index_percentage = column_name_to_index["PERCENTAGE"]
    index_taxpath = column_name_to_index["TAXPATH"]
    if "TAXPATHSN" in column_name_to_index:
        index_taxpathsn = column_name_to_index["TAXPATHSN"]
    else:
        index_taxpathsn = None
    return index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn


def normalize_samples(samples_list):
    for sample in samples_list:
        sample_id, sample_metadata, profile = sample
        sum_per_rank = defaultdict(float)
        for prediction in profile:
            sum_per_rank[prediction.rank] += prediction.percentage
        for prediction in profile:
            if prediction.percentage > 0:
                prediction.percentage = (prediction.percentage / sum_per_rank[prediction.rank]) * 100.0

#main functions
def EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q):
    '''
    (Z, diffab) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the weighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
    F[(i,j)] == num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node
    nodes_in_order[i] to the node nodes_in_order[j].
    '''
    num_nodes = len(nodes_in_order)
    Z = 0
    diffab = dict()
    partial_sums = P - Q
    for i in range(num_nodes - 1):
        val = partial_sums[i]
        partial_sums[Tint[i]] += val
        if val != 0:
            diffab[(i, Tint[i])] = lint[i, Tint[i]] * val  # Captures diffab
        Z += lint[i, Tint[i]] * abs(val)
    return (Z, diffab)

def pairwise_unifrac(dir, plot_title="plot", alpha=-1, show=False):
    '''
    Computes pairwise unifrac distance among profiles in a given directory
    :param dir: a directory containing profile files
    :return: a matrix of pairwise distances
    '''
    cur_dir = os.getcwd()
    file_lst = os.listdir(dir) #list files in the directory
    #print(file_lst)
    os.chdir(dir)
    if '.DS_Store' in file_lst:
        file_lst.remove('.DS_Store')
    sample_lst = [os.path.splitext(profile)[0].split('.')[0] for profile in file_lst] #e.g.env1sam10
   #create metadata
    metadata = dict()
    for name in sample_lst:
        env = name[3]
        metadata[name] = {'environment': env}
    # enumerate sample_lst, for filling matrix
    id_dict = dict()
    for i, id in enumerate(file_lst):
        id_dict[id] = i
    #initialize matrix
    dim = len(file_lst)
    dist_matrix = np.zeros(shape=(dim, dim))
    for pair in it.combinations(file_lst, 2):
        id_1,id_2 = pair[0], pair[1]
        i,j = id_dict[id_1], id_dict[id_2]
        profile_list1 = open_profile_from_tsv(id_1, False)
        profile_list2 = open_profile_from_tsv(id_2, False)
        name1, metadata1, profile1 = profile_list1[0]
        name2, metadata2, profile2 = profile_list2[0]
        profile1 = Profile(sample_metadata=metadata1, profile=profile1, branch_length_fun=lambda x: x**alpha)
        profile2 = Profile(sample_metadata=metadata2, profile=profile2, branch_length_fun=lambda x: x**alpha)
        #(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_no_normalize(profile2)
        (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
        (weighted, _) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
        dist_matrix[i][j] = dist_matrix[j][i] = weighted
    os.chdir(cur_dir)
    label = list(map(lambda x: x[3], sample_lst))
    # print(label)
    #score = silhouette_score(dist_matrix, label, metric="precomputed")
    #print("silhouette score is: %d" %score)
    #view_df = dm.to_data_frame()
    #view_df.to_csv('./data/test_df.txt', sep="\t")
    if show is True:
        df = pd.DataFrame.from_dict(metadata, orient='index')
        dm = DistanceMatrix(dist_matrix, sample_lst)
        dist_pc = pcoa(dm)
        dist_pc.plot(df=df, column="environment", cmap="Set1", title=plot_title)
        print("showing plot")
        plt.show()
    #plt.savefig('data/wgs_pcoa.png')
    return sample_lst, dist_matrix, metadata


#exp 1
def run_one(dist_dict, tax_dict, num_org, num_sample, Range, dissimilarity, run, out_dir):
    '''
    Create one run of experiment for exp1
    '''
    #create directory
    cur_dir = os.getcwd()
    dir_name = out_dir + "/env2r" + str(Range) + "d" + str(dissimilarity) + "-" + str(run)
    if os.path.exists(dir_name):
        print("directory exists")
        return
    os.mkdir(dir_name)
    os.chdir(dir_name)
    data = create_data_simple(num_org=num_org, num_sample=num_sample, sample_range=Range,
                              similarity=dissimilarity, distance_dict=dist_dict, tax_dict=tax_dict)
    updated_data = create_biom_table(dir_name, data, 'otu_table.tsv', True)
    os.mkdir("profiles")
    for key, value in list(updated_data.items()):
        filename = "{}{}".format(key, '.profile')
        create_profile(value, 'profiles', filename)
    os.chdir(cur_dir)

def create_data_simple(num_org, num_sample, sample_range, distance_dict, tax_dict, similarity=-1):
    '''
    :param num_org: number of organisms per sample
    :param num_sample: number of samples per environment
    :param sample_range: the spread of sample
    :param similarity: how close are the two nodes chosen
    :param tax_dict: otu:taxid dict
    :return:
    '''
    rank_list = (["superkingdom", "phylum", "class", "order", "family", "genus", "species"])
    node1 = random.choice(list(distance_dict.keys()))
    node2 = distance_dict[node1][similarity]
    #for test only
    #node2 = random.choice(list(distance_dict.keys()))
    print(node1)
    print(node2)
    data_dict = dict()
    env1_nodes = distance_dict[node1][:sample_range - 1]
    env2_nodes = distance_dict[node2][:sample_range - 1]
    remove_list = set()
    #update the 2 lists above to contain Node object instead
    for i, node in enumerate(env1_nodes):
        #create Nodes, update tax
        new_node = Node(name=node, tax=tax_dict[node])
        if int(new_node.tax) != ncbi.get_lineage(new_node.tax)[-1]:
            new_node.tax = ncbi.get_lineage(new_node.tax)[-1]
        env1_nodes[i] = new_node
        if check_rank(new_node.tax) is False:
            remove_list.add(new_node)
            print("remove ", new_node.tax)
    for i, node in enumerate(env2_nodes):
        new_node = Node(name=node, tax=tax_dict[node])
        if int(new_node.tax) != ncbi.get_lineage(new_node.tax)[-1]:
            new_node.tax = ncbi.get_lineage(new_node.tax)[-1]
        env2_nodes[i] = new_node
        if check_rank(new_node.tax) is False:
            print("remove ", new_node.tax)
            remove_list.add(new_node)
    #create sample
    if len(remove_list) > 0:
        env1_nodes = list(filter(lambda i:i not in remove_list, env1_nodes))
        env2_ndoes = list(filter(lambda i:i not in remove_list, env2_nodes))
    print(len(env1_nodes), " env1 nodes")
    print(len(env2_nodes), " env2 nodes")

    if num_org >= sample_range:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = copy.deepcopy(env1_nodes)
            value2 = copy.deepcopy(env2_nodes)
            random.shuffle(value1)
            random.shuffle(value2)
            data_dict[env1_key] = value1
            data_dict[env2_key] = value2
    else:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = copy.deepcopy(random.sample(env1_nodes, num_org))
            value2 = copy.deepcopy(random.sample(env2_nodes, num_org))
            data_dict[env1_key] = value1
            data_dict[env2_key] = value2
    return data_dict

def create_biom_table(table_id, data, filename, normalize=False):
    '''
    to be called after obtaining data by calling create_data
    :param sample_metadata:
    :param table_id:
    :param data: dictionary in the form of sample_id:list of Nodes
    :return:
    '''
    otus = []
    sample_id = []  # column index
    for key, value in list(data.items()):
        sample_id.append(key)
        value_name = list(map(lambda x: x.name, value))
        otus = otus + value_name
    otu_ids = list(set(otus))  # row index unique otus
    print('total {} otus'.format(len(otu_ids)))
    df = pd.DataFrame(columns=sample_id, index=otu_ids)
    for key, value in list(data.items()): #key = sample id, value = list of Nodes
        for x, node in enumerate(value, 1):
            ab = 100. / (1.5 ** x)
            ab = ab + halfnorm.rvs()
            df.at[node.name, key] = ab
    df = df.fillna(.0)
    #print(df)
    table = Table(data=df.to_numpy(), observation_ids=otu_ids, sample_ids=sample_id, observation_metadata=None,
                  sample_metadata=None, table_id=table_id)
    normed = table.norm(axis='sample', inplace=False)
    for key, value in list(data.items()):
        for node in value:
            node.abundance = normed.get_value_by_ids(node.name, key) * 100
    with open(filename, "w") as f:
        if normalize:
            normed.to_tsv(direct_io=f)
        else:
            table.to_tsv(direct_io=f)
    return data

def run_one(dist_dict, tax_dict, num_org, num_sample, Range, dissimilarity, run, out_dir):
    #create directory
    cur_dir = os.getcwd()
    dir_name = out_dir + "/env2r" + str(Range) + "d" + str(dissimilarity) + "-" + str(run)
    if os.path.exists(dir_name):
        print("directory exists")
        return
    os.mkdir(dir_name)
    os.chdir(dir_name)
    data = create_data_simple(num_org=num_org, num_sample=num_sample, sample_range=Range,
                              similarity=dissimilarity, distance_dict=dist_dict, tax_dict=tax_dict)
    updated_data = create_biom_table(dir_name, data, 'otu_table.tsv', True)
    os.mkdir("profiles")
    for key, value in list(updated_data.items()):
        filename = "{}{}".format(key, '.profile')
        create_profile(value, 'profiles', filename)
    os.chdir(cur_dir)

def create_profile(node_list, outdir, filename):
    '''

    :param node_list: list of Node objects
    :param outdir:
    :param filename:
    :return:
    '''
    rank_list = (["superkingdom", "phylum", "class", "order", "family", "genus", "species"])
    ab_dict=dict()
    for node in node_list:
        lin = ncbi.get_lineage(node.tax) #lineage
        lin_dict = ncbi.get_rank(lin) #id:rank
        #check if all ranks are present
        for rank in rank_list:
            if rank not in lin_dict.values():
                continue
        lin_dict_reverse = {y: x for x, y in lin_dict.items()} #rank:id
        if node.tax != lin[-1]: #in case the id is obsolete and translated to another
            node.tax = lin[-1]
        if lin_dict[node.tax] != "species":
            node.tax = lin_dict_reverse["species"]
        if node.tax not in ab_dict:
            ab_dict[node.tax] = .0
        ab_dict[node.tax]+=node.abundance
    unique_nodes = list(map(lambda x: Node(name=x, tax=x, abundance=ab_dict[x]), list(ab_dict.keys())))
    id_list = list(ab_dict.keys())
    outfile = outdir + '/' + filename
    f = open(outfile, "w+")
    f.write("# Taxonomic Profiling Output\n"
            "@SampleID:SAMPLEID\n"
            "@Version:0.9.1\n"
            "@Ranks:superkingdom|phylum|class|order|family|genus|species\n"
            "@TaxonomyID:ncbi-taxonomy_DATE\n"
            "@@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE\n")
    rank_dict = dict()
    for rank in rank_list:
        rank_dict[rank] = []
    rank_list.reverse()
    #update abundance information
    #print("here are the nodes")
    #print(len(unique_nodes))
    #print(list(map(lambda x: x.tax, unique_nodes)))
    #print(list(map(lambda x:x.abundance, unique_nodes)))
    for id in id_list:
        lin_list = ncbi.get_lineage(id)  # get lineage
        lin_dict = ncbi.get_rank(lin_list)  # create dict id:rank
        lin_dict_reverse = {y: x for x, y in lin_dict.items()}  # reverse dict rank:id
        cur_node = _get_node_from_taxid(id, unique_nodes)
        cur_abund = cur_node.abundance
        rank_dict["species"].append(cur_node)
        for rank in rank_list[1:]:
            cur_taxid = lin_dict_reverse[rank]
            cur_node = _get_node_from_taxid(cur_taxid, unique_nodes)
            if cur_node is None:
                cur_node = Node(name=rank+str(cur_taxid), tax=cur_taxid)
                unique_nodes.append(cur_node)
            cur_node.abundance = cur_node.abundance+cur_abund
            if cur_node not in rank_dict[rank]:
                rank_dict[rank].append(cur_node)
            #cur_abund = cur_node.abundance
    rank_list.reverse()
    #for k, v in rank_dict.items():
    #    print(k)
    #    print(len(list(v)))
    #print out
    for rank in rank_list:
        rank_pos = rank_list.index(rank)
        for node in rank_dict[rank]:
            lin_list = ncbi.get_lineage(node.tax)
            lin_dict = ncbi.get_rank(lin_list)
            lin_dict_reverse = {y: x for x, y in lin_dict.items()}  # reverse dict rank:id
            superkingdom = lin_dict_reverse["superkingdom"]
            taxpath = str(superkingdom)
            namepath = ncbi.get_taxid_translator([superkingdom])[superkingdom]
            for r in rank_list[1:rank_pos+1]:
                taxid = lin_dict_reverse[r]
                taxpath = taxpath + "|" + str(taxid)
                name = ncbi.get_taxid_translator([taxid])[taxid]
                cur_node = _get_node_from_taxid(taxid, unique_nodes)
                namepath = namepath + "|" + name
            f.writelines([str(node.tax), "\t", rank, "\t", taxpath, "\t", namepath, "\t", str(cur_node.abundance)])
            f.write("\n")
    f.close()
    return

def get_dataframe(dir, alpha, save_as):
    col_names = ["range", "dissimilarity", "silhouette", "Calinski-Harabasz", "Davies-Bouldin", "data_type", "sample_id"]
    file_lst = os.listdir(dir)
    cur_dir = os.getcwd()
    save_file_name = cur_dir + '/' + save_as
    print("Dataframe will be saved as: ", save_file_name)
    open(save_file_name, 'w+')
    os.chdir(dir)
    sil_score_16s = []
    sil_score_wgs = []
    calinski_16s = []
    calinski_wgs = []
    davies_16s = []
    davies_wgs = []
    Range = []
    similarity = []
    for file in file_lst:
        os.chdir(file) #individual run
        #get 16s score
        rg = int(re.findall("env2r(.*)d", file)[0])
        Range.append(rg)
        sim = int(re.findall("d(.*)-", file)[0])
        if sim == -1:
            sim = 35461
        similarity.append(sim)
        dist_matrix_16s = pd.read_table("distance-matrix.tsv", index_col=0)
        label_16s = list(map(lambda x: x[3], dist_matrix_16s))
        score_16s = silhouette_score(dist_matrix_16s, label_16s, metric="precomputed")
        sil_score_16s.append(score_16s)
        calinski_16s.append(calinski_harabasz_score(dist_matrix_16s, label_16s))
        davies_16s.append(davies_bouldin_score(dist_matrix_16s, label_16s))
        #get wgs score
        sample_lst_wgs, dist_matrix_wgs, metadata = pairwise_unifrac('profiles', alpha=alpha, show=False)
        label_wgs = list(map(lambda x: x[3], sample_lst_wgs)) #which environment
        score_wgs = silhouette_score(dist_matrix_wgs, label_wgs, metric="precomputed")
        sil_score_wgs.append(score_wgs)
        calinski_wgs.append(calinski_harabasz_score(dist_matrix_wgs, label_wgs))
        davies_wgs.append(davies_bouldin_score(dist_matrix_wgs, label_wgs))
        os.chdir('..')
    os.chdir(cur_dir)
    df_16s = pd.DataFrame(columns=col_names, index=range(len(file_lst)))
    df_16s["data_type"] = "16s"
    df_16s["sample_id"] = file_lst
    df_16s["range"] = Range
    df_16s["dissimilarity"] = similarity
    df_16s["silhouette"] = sil_score_16s
    df_16s["Calinski-Harabasz"] = calinski_16s
    df_16s["Davies-Bouldin"] = davies_16s
    df_wgs = pd.DataFrame(columns=col_names, index=range(len(file_lst)))
    df_wgs["data_type"] = "wgs"
    df_wgs["sample_id"] = file_lst
    df_wgs["range"] = Range
    df_wgs["dissimilarity"] = similarity
    df_wgs["silhouette"] = sil_score_wgs
    df_wgs["Calinski-Harabasz"] = calinski_wgs
    df_wgs["Davies-Bouldin"] = davies_wgs
    df_combined = pd.concat([df_16s, df_wgs])
    df_combined.to_csv(save_as, sep="\t")
    return df_combined

#GTDB experiments
def parse_taxonomy_file(tax_file, generate_all_tax=False, outfile=None):
    '''

    :param tax_file:
    :param generate_all_tax: if set to True, generate a file with unique taxonomy names to be used for taxonkit
    :return: a nested dictionary of taxonomy for each organism
    '''
    if generate_all_tax==True:
        tax_set = set()
    tax_dict = dict()
    with open(tax_file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            id_tax = line.split('\t')
            tax_path = id_tax[1].split(';')
            item_dict = dict()
            for rank_name in tax_path:
                rank = rank_name.split('__')[0]
                name = rank_name.split('__')[1]
                if generate_all_tax == True:
                    tax_set.add(name)
                if rank == 'd':
                    item_dict['superkingdom'] = name
                elif rank == 'p':
                    item_dict['phylum'] = name
                elif rank == 'c':
                    item_dict['class'] = name
                elif rank == 'o':
                    item_dict['order'] = name
                elif rank == 'f':
                    item_dict['family'] = name
                elif rank == 'g':
                    item_dict['genus'] = name
                elif rank == 's':
                    item_dict['species'] = name
                else:
                    print("unknown rank: ", rank)
            tax_dict[id_tax[0]] = item_dict
    if generate_all_tax is True:
        with open(outfile, 'w+') as f:
            for item in tax_set:
                f.write(item)
                f.write('\n')
    return tax_dict

def get_name_taxid_dict(mapping_file):
    '''

    :param mapping_file: bac120_all_taxons_taxids_valid.txt
    :return:
    '''
    tax_dict=dict()
    with open(mapping_file, 'r') as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            name = line[0]
            taxid = line[1]
            tax_dict[name] = taxid
    return tax_dict

def get_GTDBid_taxid_dict(taxonomy_dict, name_taxid_dict):
    GTDBid_taxid_dict = dict()
    for GTDBid in taxonomy_dict:
        name = taxonomy_dict[GTDBid]['species']
        if name in name_taxid_dict:
            taxid = name_taxid_dict[name]
            GTDBid_taxid_dict[GTDBid] = taxid
    return GTDBid_taxid_dict

def create_GTDB_data(distance_dict, similarity, sample_range, num_org, num_sample):
    node1 = random.choice(list(distance_dict.keys()))
    node2 = distance_dict[node1][similarity]
    print(node1)
    print(node2)
    data_dict = dict()
    print(len(distance_dict.keys()))
    env1_nodes = distance_dict[node1][:sample_range]
    env2_nodes = distance_dict[node2][:sample_range]
    for i, node in enumerate(env1_nodes):
        #create Nodes, update tax
        new_node = Node(name=node)
        env1_nodes[i] = new_node
    for i, node in enumerate(env2_nodes):
        new_node = Node(name=node)
        env2_nodes[i] = new_node
    # create sample
    if num_org >= sample_range:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = copy.deepcopy(env1_nodes)
            value2 = copy.deepcopy(env2_nodes)
            random.shuffle(value1)
            random.shuffle(value2)
            data_dict[env1_key] = value1
            data_dict[env2_key] = value2
    else:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = copy.deepcopy(random.sample(env1_nodes, num_org))
            value2 = copy.deepcopy(random.sample(env2_nodes, num_org))
            data_dict[env1_key] = value1
            data_dict[env2_key] = value2
    return data_dict

def create_GTDB_biom_table(table_id, data, filename, normalize=False):
    otus = []
    sample_id = []  # column index
    for key, value in list(data.items()):
        sample_id.append(key)
        value_name = list(map(lambda x: x.name, value))
        otus = otus + value_name
    otu_ids = list(set(otus))  # row index unique otus
    print('total {} otus'.format(len(otu_ids)))
    df = pd.DataFrame(columns=sample_id, index=otu_ids)
    for key, value in list(data.items()):  # key = sample id, value = list of Nodes
        for x, node in enumerate(value, 1):
            ab = 100. / (1.5 ** x)
            ab = ab + halfnorm.rvs()
            df.at[node.name, key] = ab
    df = df.fillna(.0)
    print(df)
    table = Table(data=df.to_numpy(), observation_ids=otu_ids, sample_ids=sample_id, observation_metadata=None,
                  sample_metadata=None, table_id=table_id)
    normed = table.norm(axis='sample', inplace=False)
    for key, value in list(data.items()):
        for node in value:
            node.abundance = normed.get_value_by_ids(node.name, key) * 100
    with open(filename, "w") as f:
        if normalize:
            normed.to_tsv(direct_io=f)
        else:
            table.to_tsv(direct_io=f)
    return data

def create_GTDB_profile(node_list, outdir, filename, tax_dict, name_tax_dict, GTDBid_taxid_dict):
    '''
    node_lst: a list of Node object with abundance. No taxid yet. Only species
    tax_dict: a nested dict {GTDB id: {rank : scientific name}}
    name_tax_dict: {scientific name : taxid}
    GTDBid_taxid_dict: {GTDBid : taxid}
    '''
    taxpath_dict = dict()
    ab_dict = dict() # { taxid : abundance }, all ranks
    rank_list = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    rank_dict = dict() #{ rank : set of taxids belonging to this rank}
    for rank in rank_list:
        rank_dict[rank] = set()
    for node in node_list:
        node.tax = GTDBid_taxid_dict[node.name]
        #print(node.tax)
        #ab_dict[node.tax] = node.abundance #added all the species level abundances
        #print(ab_dict)
        lineage = tax_dict[node.name] #{rank : scientific name}
        #debug-----
        #tax_lin = dict()
        #for rank in lineage:
        #    tax_lin[rank] = name_tax_dict[lineage[rank]]
        #print(tax_lin)
        #debug-----
        rank_dict['species'].add(node.tax)
        for rank in lineage:
            taxid = name_tax_dict[lineage[rank]]
            rank_dict[rank].add(taxid)
            if taxid not in ab_dict:
                ab_dict[taxid] = node.abundance #if new, abundance is the abundance of the species
            else: #if exists, add up
                ab_dict[taxid] += node.abundance
    #print(ab_dict)
    #tax path for each node
    for node in node_list: #species
        lineage = tax_dict[node.name]	#{rank : scientific name}
        taxpath = lineage["superkingdom"]
        for rank in rank_list[1:]:
            taxpath = taxpath + "|" + lineage[rank]
        taxpath_dict[node.tax] = taxpath
        #add other paths in this lineage, starting from phylum
        for rank in rank_list: #from phylum
            taxpath = lineage["superkingdom"]
            name = lineage[rank]
            taxid = name_tax_dict[name]
            rank_pos = rank_list.index(rank)
            for other_rank in rank_list[1:rank_pos+1]:
                taxpath = taxpath + "|" + lineage[other_rank]
            taxpath_dict[taxid] = taxpath
    #print(taxpath_dict)
    #print out
    # print
    outfile = outdir + '/' + filename
    f = open(outfile, "w+")
    f.write("# Taxonomic Profiling Output\n"
            "@SampleID:SAMPLEID\n"
            "@Version:0.9.1\n"
            "@Ranks:superkingdom|phylum|class|order|family|genus|species\n"
            "@TaxonomyID:ncbi-taxonomy_DATE\n"
            "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
    for rank in rank_list:
        for taxid in rank_dict[rank]:
            namepath = taxpath_dict[taxid]
            namepath_lst = namepath.split('|')
            taxpath_lst = list(map(lambda x: name_tax_dict[x], namepath_lst))
            taxpath = "|".join(taxpath_lst)
            f.writelines([str(taxid), "\t", rank, "\t", taxpath, "\t", namepath, "\t", str(ab_dict[taxid])])
            f.write("\n")
    f.close()
    return

def node_selection(distance_dict, dissimilarity, sample_range, num_org, num_sample):
    '''

    :param distance_dict: obtained from second distance matrix
    :param dissimilarity:
    :param sample_range:
    :param num_org:
    :param num_sample:
    :return:
    '''
    node1 = random.choice(list(distance_dict.keys()))
    node2 = distance_dict[node1][dissimilarity]
    print(node1)
    print(node2)
    data_dict = dict()
    print(len(distance_dict.keys()))
    env1_nodes = distance_dict[node1][:sample_range]
    env2_nodes = distance_dict[node2][:sample_range]
    for i, node in enumerate(env1_nodes):
        # create Nodes
        new_node = Node(name=node)
        env1_nodes[i] = new_node
    for i, node in enumerate(env2_nodes):
        new_node = Node(name=node)
        env2_nodes[i] = new_node
    # create sample
    if num_org >= sample_range:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = copy.deepcopy(env1_nodes)
            value2 = copy.deepcopy(env2_nodes)
            random.shuffle(value1)
            random.shuffle(value2)
            data_dict[env1_key] = value1
            data_dict[env2_key] = value2
    else:
        for i in range(num_sample):
            env1_key = "{}{}".format('env1sam', i)
            env2_key = "{}{}".format('env2sam', i)
            value1 = copy.deepcopy(random.sample(env1_nodes, num_org))
            value2 = copy.deepcopy(random.sample(env2_nodes, num_org))
            data_dict[env1_key] = value1
            data_dict[env2_key] = value2
    # add abundance
    for key, value in list(data_dict.items()):  # key = sample id, value = list of Nodes
        total = 0.
        for x, node in enumerate(value, 1):
            ab = 100. / (1.5 ** x)
            ab = ab + halfnorm.rvs()
            total += ab
            node.abundance = ab
        for node in value: #normalize
            node.abundance = node.abundance/total * 100.
    return data_dict

def create_ncbi_GTDB_profile(node_list, outdir, filename, name_tax_dict, GTDBid_taxid_dict):
    '''
    Create profile based on NCBI lineage
    :param node_list: a list of Node objects with abundances. No taxid. Only species. Obtained from node_selection
    :param outdir:
    :param filename:
    :param tax_dict:
    :param name_tax_dict:
    :param GTDBid_taxid_dict:
    :return:
    '''

    ab_dict = dict() # {taxid:abundance}
    rank_list = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    rank_dict = dict()
    taxpath_dict = dict()
    tax_name_dict = {y:x for x, y in name_tax_dict.items()} #taxid:sci name
    for rank in rank_list:
        rank_dict[rank] = set()
    for node in node_list:
        node.tax = GTDBid_taxid_dict[node.name]
        #print(node.tax)
        lin = ncbi.get_lineage(node.tax) #lineage
        lin_dict = ncbi.get_rank(lin) #id:rank
        lin_dict_reverse = {y:x for x, y in lin_dict.items()} #rank:id
        if node.tax != lin[-1]:
            node.tax = lin[-1]
        if lin_dict[node.tax] != "species":
            node.tax = lin_dict_reverse["species"]
        rank_dict['species'].add(node.tax)
        for rank in rank_list:
            taxid = lin_dict_reverse[rank]
            rank_dict[rank].add(taxid)
            if taxid not in ab_dict:
                ab_dict[taxid] = node.abundance
            else:
                ab_dict[taxid] += node.abundance
        taxpath = str(lin_dict_reverse["superkingdom"])
        for rank in rank_list[1:]:
            taxpath = taxpath + "|" + str(lin_dict_reverse[rank])
        taxpath_dict[node.tax] = taxpath
        for rank in rank_list:
            taxpath = str(lin_dict_reverse["superkingdom"])
            taxid = lin_dict_reverse[rank]
            rank_pos = rank_list.index(rank)
            for other_rank in rank_list[1:rank_pos+1]:
                taxpath = taxpath + '|' + str(lin_dict_reverse[other_rank])
            taxpath_dict[taxid] = taxpath

    outfile = outdir + '/' + filename
    f = open(outfile, "w+")
    f.write("# Taxonomic Profiling Output\n"
            "@SampleID:SAMPLEID\n"
            "@Version:0.9.1\n"
            "@Ranks:superkingdom|phylum|class|order|family|genus|species\n"
            "@TaxonomyID:ncbi-taxonomy_DATE\n"
            "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
    #print(rank_dict)
    for rank in rank_list:
        for taxid in rank_dict[rank]:
            taxpath = taxpath_dict[taxid]
            taxpath_lst = taxpath.split('|')
            namepath_lst = list(map(lambda x: ncbi.get_taxid_translator([taxid])[taxid], taxpath_lst))
            namepath = "|".join(namepath_lst)
            f.writelines([str(taxid), "\t", rank, "\t", taxpath, "\t", namepath, "\t", str(ab_dict[taxid])])
            f.write("\n")
    f.close()
    return

def get_GTDB_dataframe(dir, alpha, save_as):
    col_names = ["range", "dissimilarity", "silhouette", "Calinski-Harabasz", "Davies-Bouldin", "data_type", "sample_id"]
    cur_dir = os.getcwd()
    #dir_lst = os.listdir(dir)[1:]
    dir_lst = os.listdir(dir)
    os.chdir(dir)
    sil_score_16s = []
    sil_score_wgs = []
    calinski_16s = []
    calinski_wgs = []
    davies_16s = []
    davies_wgs = []
    Range = []
    dissimilarity = []
    for exp_dir in dir_lst:
        os.chdir(exp_dir)  # individual run
        # get 16s score
        rg = int(re.findall("r(.*)d", exp_dir)[0])
        Range.append(rg)
        sim = int(re.findall("d(.*)-", exp_dir)[0])
        dissimilarity.append(sim)
        dist_matrix_16s = pd.read_table("distance-matrix.tsv", index_col=0)
        label_16s = list(map(lambda x: x[3], dist_matrix_16s))
        score_16s = silhouette_score(dist_matrix_16s, label_16s, metric="precomputed")
        print(score_16s)
        sil_score_16s.append(score_16s)
        calinski_16s.append(calinski_harabasz_score(dist_matrix_16s, label_16s))
        davies_16s.append(davies_bouldin_score(dist_matrix_16s, label_16s))
        # get wgs score
        sample_lst_wgs, dist_matrix_wgs, metadata = pairwise_unifrac('profiles', alpha=alpha, show=False)
        label_wgs = list(map(lambda x: x[3], sample_lst_wgs))
        score_wgs = silhouette_score(dist_matrix_wgs, label_wgs, metric="precomputed")
        sil_score_wgs.append(score_wgs)
        calinski_wgs.append(calinski_harabasz_score(dist_matrix_wgs, label_wgs))
        davies_wgs.append(davies_bouldin_score(dist_matrix_wgs, label_wgs))
        os.chdir('..')
    df_16s = pd.DataFrame(columns=col_names, index=range(len(dir_lst)))
    df_16s["data_type"] = "16s"
    df_16s["sample_id"] = dir_lst
    df_16s["range"] = Range
    df_16s["dissimilarity"] = dissimilarity
    df_16s["silhouette"] = sil_score_16s
    df_16s["Calinski-Harabasz"] = calinski_16s
    df_16s["Davies-Bouldin"] = davies_16s
    df_wgs = pd.DataFrame(columns=col_names, index=range(len(dir_lst)))
    df_wgs["data_type"] = "wgs"
    df_wgs["sample_id"] = dir_lst
    df_wgs["range"] = Range
    df_wgs["dissimilarity"] = dissimilarity
    df_wgs["silhouette"] = sil_score_wgs
    df_wgs["Calinski-Harabasz"] = calinski_wgs
    df_wgs["Davies-Bouldin"] = davies_wgs
    df_combined = pd.concat([df_16s, df_wgs])
    print(df_combined)
    os.chdir(cur_dir)
    df_combined.to_csv(save_as, sep="\t")
    return df_combined

def get_GTDB_dataframe_2(dir, alpha, save_as):
    col_names = ["range", "dissimilarity", "silhouette", "Calinski-Harabasz", "Davies-Bouldin", "data_type",
                 "sample_id"]
    cur_dir = os.getcwd()
    dir_lst = os.listdir(dir)
    os.chdir(dir)
    sil_score_gtdb = []
    sil_score_ncbi = []
    calinski_gtdb = []
    calinski_ncbi = []
    davies_gtdb = []
    davies_ncbi = []
    Range = []
    dissimilarity = []
    for exp_dir in dir_lst:
        os.chdir(exp_dir)  # individual run
        rg = int(re.findall("r(.*)d", exp_dir)[0])
        Range.append(rg)
        sim = int(re.findall("d(.*)-", exp_dir)[0])
        dissimilarity.append(sim)
        # get gtdb score
        sample_lst_gtdb, dist_matrix_gtdb, metadata = pairwise_unifrac('GTDB_profiles', alpha=alpha, show=False)
        label_gtdb = list(map(lambda x: x[3], sample_lst_gtdb))
        score_gtdb = silhouette_score(dist_matrix_gtdb, label_gtdb, metric="precomputed")
        sil_score_gtdb.append(score_gtdb)
        calinski_gtdb.append(calinski_harabasz_score(dist_matrix_gtdb, label_gtdb))
        davies_gtdb.append(davies_bouldin_score(dist_matrix_gtdb, label_gtdb))
        #get ncbi score
        sample_lst_ncbi, dist_matrix_ncbi, metadata = pairwise_unifrac('NCBI_profiles', alpha=alpha, show=False)
        label_ncbi = list(map(lambda x: x[3], sample_lst_ncbi))
        score_ncbi = silhouette_score(dist_matrix_ncbi, label_ncbi, metric="precomputed")
        sil_score_ncbi.append(score_ncbi)
        calinski_ncbi.append(calinski_harabasz_score(dist_matrix_ncbi, label_ncbi))
        davies_ncbi.append(davies_bouldin_score(dist_matrix_ncbi, label_ncbi))
        os.chdir('..')
    df_gtdb = pd.DataFrame(columns=col_names, index=range(len(dir_lst)))
    df_gtdb["data_type"] = "GTDB"
    df_gtdb["sample_id"] = dir_lst
    df_gtdb["range"] = Range
    df_gtdb["dissimilarity"] = dissimilarity
    df_gtdb["silhouette"] = sil_score_gtdb
    df_gtdb["Calinski-Harabasz"] = calinski_gtdb
    df_gtdb["Davies-Bouldin"] = davies_gtdb
    df_ncbi = pd.DataFrame(columns=col_names, index=range(len(dir_lst)))
    df_ncbi["data_type"] = "NCBI"
    df_ncbi["sample_id"] = dir_lst
    df_ncbi["range"] = Range
    df_ncbi["dissimilarity"] = dissimilarity
    df_ncbi["silhouette"] = sil_score_ncbi
    df_ncbi["Calinski-Harabasz"] = calinski_ncbi
    df_ncbi["Davies-Bouldin"] = davies_ncbi
    df_combined = pd.concat([df_gtdb, df_ncbi])
    print(df_combined)
    os.chdir(cur_dir)
    df_combined.to_csv(save_as, sep="\t")
    return df_combined

#simulated
def get_grinder_abundances_for_both(sample_num, org_num, out_dir, env_num, rnge, dist):
    distance_dict = get_dist_dict('data/grinder_distance_matrix_match_primer.txt')
    otu_acc_dict = get_dict_from_file("data/mapping_file2.txt", 0, 2)
    
    #if not os.path.exists(out_dir):
    #    os.mkdir(out_dir)
    cur_dir = os.getcwd()
    os.chdir(out_dir)
    if not os.path.exists('16s_abundance_files'):
        os.mkdir('16s_abundance_files')
    if not os.path.exists('wgs_abundance_files'):
        os.mkdir('wgs_abundance_files')
    # choose env nodes
    env_nodes = []
    for i in range(env_num):
        print(env_nodes)
        if len(env_nodes) == 0:
            node1 = random.choice(list(distance_dict.keys()))  # env 1
            env_nodes.append(node1)
        else:
            prev_node = env_nodes[-1]
            #print('prev_node is', prev_node)
            cur_node = distance_dict[prev_node][dist]
            #print('current node is', cur_node)
            neighbor = 1
            while cur_node in env_nodes:
                # if cur_node already in env_nodes, pick one close enough
                cur_node = distance_dict[cur_node][neighbor]
                neighbor += 1
            env_nodes.append(cur_node)
    # 16s abundance files
    os.chdir('16s_abundance_files')
    for env in range(env_num):
        for sampl in range(sample_num):
            file_name = 'env' + str(env) + '-sample' + str(sampl) + '.txt'
            print(file_name)
            sample_nodes = random.sample(distance_dict[env_nodes[env]][:rnge], org_num)
            get_abundance_file(sample_nodes, file_name)
    wgs_dir = '../wgs_abundance_files/'
    for file in os.listdir():
        wgs_file_name = wgs_dir + 'wgs-' + file
        with open(file,'r') as f:
            with open(wgs_file_name, 'w+') as g:
                for line in f.readlines():
                    (otu, abundance) = line.strip().split('\t')
                    g.writelines([str(otu_acc_dict[otu]), '\t', abundance, '\n'])
    os.chdir(cur_dir)

def get_abundance_file(node_list, file_name, abund_fun="exp", factor=1.5):
    '''
    writes out a file with first column being otu and second being the relative abundance
    :param node_list: a list of otus
    :param abund_fun: abundance function, can be uniform or exponential
    :return:
    '''
    if abund_fun == "uniform":
        with open(file_name, 'w') as f:
            for otu in node_list:
                f.writelines([str(otu), '\t', str(1./len(node_list)), '\n'])
    else:
        abundances = []
        for i in range(len(node_list)):
            abundances.append(100. / (factor ** (i+1)) + halfnorm.rvs())
        normed = list(map(lambda x: round(float(x)/sum(abundances),6), abundances))
        with open(file_name, 'w') as f:
            for i, otu in enumerate(node_list):
                f.writelines([str(otu), '\t', str(normed[i]), '\n'])

#OGU experiments
def get_grinder_abundance_for_ogu(sample_num, org_num, out_dir, env_num, rnge, distance_dict, dissim):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    os.chdir(out_dir)
    if dissim == 3000:
        dissim = -1
    if not os.path.exists('abundance_files'):
        os.mkdir('abundance_files')
    # choose env nodes
    env_nodes = []
    for i in range(env_num):
        print(env_nodes)
        if len(env_nodes) == 0:
            node1 = random.choice(list(distance_dict.keys()))  # env 1
            env_nodes.append(node1)
        else:
            prev_node = env_nodes[-1]
            #print('prev_node is', prev_node)
            cur_node = distance_dict[prev_node][dissim]
            #print('current node is', cur_node)
            neighbor = 1
            while cur_node in env_nodes:
                # if cur_node already in env_nodes, pick one close enough
                cur_node = distance_dict[cur_node][neighbor]
                neighbor += 1
            env_nodes.append(cur_node)
    
    os.chdir('abundance_files')
    for env in range(env_num):
        for sampl in range(sample_num):
            file_name = 'env' + str(env) + '-sample' + str(sampl) + '.txt'
            print(file_name)
            sample_nodes = random.sample(distance_dict[env_nodes[env]][:rnge], org_num)
            get_abundance_file(sample_nodes, file_name)
    

def get_ogu_vs_wgsunifrac_df(dir, save):
    '''
    Produce a file with the following columns: range, dissimilarity, silhouette, data_type, sample_id
    param: dir containing all experiment files, which contains profiles and ogu matrix.
    '''
    col_names = ["range", "dissimilarity", "setup", "Silhouette", "method", "sample_id"]
    cur_dir = os.getcwd()
    exp_lst = os.listdir(dir)
    if '.DS_Store' in exp_lst:
        exp_lst.remove('.DS_Store')
    #exp_lst = list(set(map(lambda x:re.findall("(.*)-", x)[0], exp_lst))) #remove -1, -2 etc. suffix
    os.chdir(dir)
    sil_score_ogu = []
    sil_score_wgs = []
    Range = []
    dissimilarity = []
    setup_list = []
    for exp in exp_lst:
        rg = int(re.findall("r(.*)d", exp)[0])
        Range.append(rg)
        dissim = int(re.findall("d(.*)-", exp)[0])
        dissimilarity.append(dissim)
        setup = str(re.findall("(.*)-", exp)[0])
        setup_list.append(setup)
        os.chdir(exp)  # individual run
        #get ogu score
        dist_matrix_ogu = pd.read_table("exported/distance-matrix.tsv", index_col=0)
        label_ogu = list(map(lambda x: x[3], dist_matrix_ogu)) #which environment
        score_ogu = silhouette_score(dist_matrix_ogu, label_ogu, metric="precomputed")
        sil_score_ogu.append(score_ogu)
        #get wgs score
        sample_lst_wgs, dist_matrix_wgs, metadata = pairwise_unifrac('profiles', alpha=-1, show=False)
        label_wgs = list(map(lambda x: x[3], sample_lst_wgs))
        score_wgs = silhouette_score(dist_matrix_wgs, label_wgs, metric="precomputed")
        sil_score_wgs.append(score_wgs)
        os.chdir('..')
        #ave_sil_ogu = np.mean(this_sil_score_ogu)
        #ave_sil_wgs = np.mean(this_sil_score_wgs)
        print(exp)
    
    df_ogu = pd.DataFrame(columns=col_names, index=range(len(exp_lst)))
    df_ogu["method"] = "OGU"
    df_ogu["sample_id"] = exp_lst
    df_ogu["range"] = Range
    df_ogu["dissimilarity"] = dissimilarity
    df_ogu["Silhouette"] = sil_score_ogu
    df_ogu["setup"] = setup_list
    df_wgs = pd.DataFrame(columns=col_names, index=range(len(exp_lst)))
    df_wgs["method"] = "WGSUniFrac"
    df_wgs["sample_id"] = exp_lst
    df_wgs["range"] = Range
    df_wgs["dissimilarity"] = dissimilarity
    df_wgs["Silhouette"] = sil_score_wgs
    df_wgs["setup"] = setup_list
    df_combined = pd.concat([df_ogu, df_wgs])
    print(df_combined)
    os.chdir(cur_dir)
    df_combined.to_csv(save, sep="\t")
    return df_combined
    return

def get_ogu_vs_wgsunifrac_plot(dataframe_file, x, save):
    '''
    :param: dataframe_file. A dataframe file with the following columns: range, dissimilarity, silhouette, data_type, sample_id
    :param: file_name: file name to be saved as
    :param: x-asix: range or dissimilarity
    '''
    df = pd.read_table(dataframe_file, index_col=0)
    #sns.set_theme(style="ticks", palette="pastel")
    sns.lineplot(x=x, y="Silhouette", hue="method", data=df)
    #sns.lineplot(x=x, y="Silhouette", hue="method", data=df, err_style="bars", ci="sd")
    plt.savefig(save)

#Real data
def get_metadata_from_real_data_partial(meta_file, profile_dir, by):
    file_lst = os.listdir(profile_dir)
    sample_lst = list(map(lambda x: x.split('.')[0], file_lst))
    metadata=dict()
    if by == "bodysites":
        by_col = "HMgDB_sample_site_1"
    elif by == "study":
        by_col = "BioprojectID"
    df = pd.read_csv(meta_file, usecols=["library_id", by_col])
    id_sites_dict=df.set_index('library_id').to_dict()[by_col]
    for id in sample_lst:
        metadata[id] = {'environment': id_sites_dict[id]}
    return metadata

def just_pairwise_unifrac(dir, alpha, save_as):
    '''
    :param dir: directory containing the .profile files
    :param alpha factor for branch length function  x**alpha
    :return: a dataframe of pairwise distance matrix
    '''
    if save_as is None:
        save_as = "pairwise_WGSUniFrac_matrix.csv"
    cur_dir = os.getcwd()
    file_lst = os.listdir(dir)  # list files in the directory
    # print(file_lst)
    os.chdir(dir)
    if '.DS_Store' in file_lst:
        file_lst.remove('.DS_Store')
    sample_lst = [os.path.splitext(profile)[0].split('.')[0] for profile in file_lst] #e.g.env1sam10. i.e.filenames without extension
    #print(sample_lst)
    # enumerate sample_lst, for filling matrix
    id_dict = dict()
    for i, id in enumerate(file_lst):
        id_dict[id] = i
    # initialize matrix
    dim = len(file_lst)
    dist_matrix = np.zeros(shape=(dim, dim))
    count=0
    for pair in it.combinations(file_lst, 2): #all pairwise combinations
        #to keep the running less boring
        count+=1
        if count % 100 == 0:
            print(count)
        id_1, id_2 = pair[0], pair[1]
        i, j = id_dict[id_1], id_dict[id_2]
        profile_list1 = open_profile_from_tsv(id_1, False)
        profile_list2 = open_profile_from_tsv(id_2, False)
        name1, metadata1, profile1 = profile_list1[0]
        name2, metadata2, profile2 = profile_list2[0]
        profile1 = Profile(sample_metadata=metadata1, profile=profile1, branch_length_fun=lambda x: x ** alpha)
        profile2 = Profile(sample_metadata=metadata2, profile=profile2, branch_length_fun=lambda x: x ** alpha)
        # (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_no_normalize(profile2)
        (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = profile1.make_unifrac_input_and_normalize(profile2)
        (weighted, _) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
        dist_matrix[i][j] = dist_matrix[j][i] = weighted
    os.chdir(cur_dir)
    pd.DataFrame(data=dist_matrix, index=sample_lst, columns=sample_lst).to_csv(save_as, sep="\t")
    return dist_matrix, sample_lst

### helper functions
def get_dist_dict(file):
    '''
    Creates distance_dict to be used for generation of raw data.
    :param file: a text file containing ranked distance information of the nodes. For each row, the first position is the node and the subsequent entries are all other nodes ranked with respect to the distance from the first node in ascending order.
    :return: a dist_dict in the form of node:list of nodes
    '''
    node_dict = dict()
    count=0
    with open(file,'r') as f:
        for line in f.readlines():
            line = line.strip()
            nodes = line.split('\t')
            node_dict[nodes[0]] = nodes[1:]
            count+=1
            #print(len(node_dict))
    return node_dict

def get_dict_from_file(file, key_col, val_col):
    '''
    A generic function that parses a file to obtain a dictionary, with entries in the key_col being the keys and entries in val_col being the values.
    Also removes the '>' symbol in front of NCBI IDs
    :param file: a text file from which the dict is to be generated
    :param key_col: key column
    :param val_col: value column
    :return: a dict with entries in the key_col being the keys and entries in val_col being the values.
    '''
    _dict = dict()
    with open(file,'r') as f:
        for line in f.readlines():
            line = line.strip()
            item = line.split('\t')
            #print(item)
            if len(item) > val_col:
                key = item[key_col]
                value = item[val_col]
            if value[0] == ">":
                value = value[1:]
            _dict[key] = value
    return _dict

def _get_node_from_taxid(taxid, node_list):
    for node in node_list:
        if int(node.tax) == int(taxid):
            return node

def check_rank(id):
    '''
    Check if all ranks are present
    :param id: a taxid to be checked
    :return: nothing if any of the ranks is missing. otherwise return id
    '''
    rank_list = (["superkingdom", "phylum", "class", "order", "family", "genus", "species"])
    try:
        lineage = ncbi.get_lineage(id)  # get lineage dict
    except ValueError:
        return False
    ranks = ncbi.get_rank(lineage).values()
    for r in rank_list:
        if r not in ranks:
            print("rank %s not present" % r)
            return False
    return True

def get_plot_from_file(file, x, y, palette, save):
    df = pd.read_table(file, index_col=0)
    print(df)
    sns.set_theme(style="ticks", palette="pastel")
    sns.boxplot(x=x, y=y, hue="data_type", data=df, palette=palette)
    plt.savefig(save)

def get_pcoa(dist_matrix, sample_lst, metadata, plot_title, save):
    df = pd.DataFrame.from_dict(metadata, orient='index')
    print(len(metadata.keys()))
    dm = DistanceMatrix(dist_matrix, sample_lst)
    print(dm)
    dist_pc = pcoa(dm)
    dist_pc.plot(df=df, column="environment", cmap="Set1", title=plot_title, axis_labels=('PC1', 'PC2', 'PC3'))
    #label = list(map(lambda x: metadata[x], sample_lst))
    #print(label)
    #score = silhouette_score(dist_matrix, label, metric="precomputed")
    #print("silhouette score is: %d" % score)
    plt.savefig(save)
