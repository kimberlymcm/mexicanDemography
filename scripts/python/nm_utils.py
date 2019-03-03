""" Utility functions and constants for project """

import numpy

import dadi


DATA_DICT_FN = "/home/kfm/kfm_projects/NA/NA_data/getIntrons/remaskingNAH/NA_CHB_exomes_20150826.data_dict"

# Maps population to size used in dadi
POP_DICT = {"TAR": 26, "MYA": 14, "TRQ": 16, "HUI": 24}

NSPEC = 6

PTS_L = [40, 50, 60] # Grid size we search over



def compare_to_self(fss):
    """ Compares the real data to the real data """

    # What is the point of this?
    ll_model = map(lambda i: dadi.Inference.ll_multinom(fss[i], fss[i]), range(NSPEC))
    print "data to data: by marginal", ll_model
    print "all marginals", numpy.sum(ll_model)


def compare_to_base(model, fss):
    """ Get results of base model """

    # Sum of the likelihoods of the models based on the startparms.
    # Compares model based on starting params to the real data
    ll_low_model = numpy.sum(map(lambda i: dadi.Inference.ll_multinom(model[i], fss[i]), range(NSPEC)))
    print "Coarse grid model Model log-likelihood: ", ll_low_model
    # The optimal value of theta given the model.
    thetas = map(lambda i: dadi.Inference.optimal_sfs_scaling(model[i], fss[i]), range(NSPEC))
    print "Initial thetas are: " + str(thetas)
