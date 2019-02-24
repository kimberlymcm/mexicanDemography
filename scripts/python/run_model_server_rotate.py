"""
   For each possible branching topology for the four native american populations,
   optimizes demographic parameters
"""


import sys
import argparse

import numpy
from numpy import array

import dadiHacks
import marginal_optimization
import dadi
import nat_model_SG
from nm_utils import POP_DICT, DATA_DICT_FN, NSPEC, PTS_L, compare_to_self, compare_to_base



pops = [["TAR", "HUI", "TRQ", "MYA"],
        ["TRQ", "MYA", "TAR", "HUI"],
        ["TRQ", "HUI", "TAR", "MYA"],
        ["TRQ", "TAR", "HUI", "MYA"],
        ["MYA", "TAR", "HUI", "TRQ"]]

popSizes = [[26, 24, 16, 14],
            [16, 14, 26, 24],
            [16, 24, 26, 14],
            [16, 26, 24, 14],
            [14, 26, 24, 16]]


def create_fss(dd, k, divergence):
    """ Create true frequency spectra """

    fux_table = "fux_table_" + str(divergence) + ".dat"
    fss = []
    fs = dadi.Spectrum.from_data_dict_corrected(dd, [pop_list[0], pop_list[3]], [POP_SIZES[k][0], POP_SIZES[k][3]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, [pop_list[0], pop_list[2]], [POP_SIZES[k][0], POP_SIZES[k][2]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, [pop_list[1], pop_list[3]], [POP_SIZES[k][1], POP_SIZES[k][3]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, [pop_list[1], pop_list[2]], [POP_SIZES[k][1], POP_SIZES[k][2]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, [pop_list[0], pop_list[1]], [POP_SIZES[k][0], POP_SIZES[k][1]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, [pop_list[3], pop_list[2]], [POP_SIZES[k][3], POP_SIZES[k][2]], fux_table)
    fss.append(fs.copy())
    return fss



def main(args):

    # TODO (kmcmanus): Break this function into smaller pieces


for k in range(5):


    dd = dadi.Misc.make_data_dict(DATA_DICT_FN)

    fss = create_fss(dd, k, args.divergence)

    nss = [fs.sample_sizes for fs in fss]
    NSPEC = 6
    #PTS_L = [30,40,50] #kk
    PTS_L = [40, 50, 60]
        print "PTS_L: ",PTS_L

    #Prepare the custom model
    func = nat_model_SG.Kimberly_method_fewerSizes

    #              [Tb,   Nb,   T_s1, T_TarHui, T_TrqMya, N_Tar, N_Hui, N_Trq, N_Mya]
    start_params = [0.09, 0.25, 0.01, 0.01,     0.01,     0.2,   0.2,   0.2,   0.2]
    fixed_params = [0.09, None, None, None,     None,     None,  None,  None,  None]
    upper_bound =  [1,    0.5,  0.02, 0.02,     0.02,     0.5,   0.5,   0.5,   0.5]
    lower_bound =  [0,    0.1,  1e-7, 1e-7,     0.0003,   0.05,  0.05,  0.05,  0.05]

    # Make the extrapolating version of our demographic model function.
    func_ex = dadiHacks.make_extrap_func_multispec(func, extrap_log=True)
    print "running model"
    model = func_ex(start_params, nss, PTS_L) # Calculate the model AFS
    print "number of reads in the data" # perform a sanity check
    print [fs.sum() for fs in fss]

    compare_to_self(fss)

    compare_to_base(model, fss)

    ratios = [fs.data.sum() / len(dd) for fs in fss]

    outFile = "20190218_4pop_Tb0.09_rotate_" + str(pops[k][0]) + "_" + str(pops[k][1]) + "_" + str(pops[k][2]) + "_" + str(pops[k][3]) + ".txt"
    f = open(outFile, "a")

    for i in range(4):
        p0 = dadi.Misc.perturb_params(start_params, lower_bound=lower_bound, upper_bound=upper_bound, fold=2)
        popts=[]
        lls=[]
        lthetas=[]
        for j in range(2):    
            popt = marginal_optimization.optimize_log_fmin(p0, fss, func_ex, PTS_L,
                                                           fixed_params=fixed_params,
                                                           lower_bound=lower_bound,
                                                           upper_bound=upper_bound,
                                                           nmarginals=NSPEC)
            popt = marginal_optimization.optimize_log_fmin(popt, fss, func_ex, PTS_L,
                                                           fixed_params=fixed_params,
                                                           lower_bound=lower_bound,
                                                           upper_bound=upper_bound,
                                                           nmarginals=NSPEC)
            model = func_ex(popt, nss, PTS_L)
            thetas=map(lambda i: dadi.Inference.optimal_sfs_scaling(model[i], fss[i]),range(NSPEC))
            lthetas.append(thetas)
            ll_model=numpy.sum(map(lambda i: dadi.Inference.ll_multinom(model[i], fss[i]), range(NSPEC)))
            lls.append(ll_model)
            popts.append(popt)
        best=lls.index(max(lls))    
        #converting parameters
        ll_model = func_ex(popts[best], nss, PTS_L)
        thetas=map(lambda i: dadi.Inference.optimal_sfs_scaling(model[i], fss[i]),range(NSPEC))
        #print "thetas are", thetas 
        effectTheta=(numpy.array(lthetas[best])/numpy.array(ratios)).mean()

        end = nat_model_SG.Kimberly_convert_params_fewerSizes(popts[best], theta=effectTheta, gen_time=29.0, L=8889201, mu=1.25e-8)
        result = str(lls[best]) + "\t" + str(effectTheta) + "\t" + "\t".join(map(str,end)) + "\n"
        print result
        f.write(result)
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Optimize demographic parameters for rotating topologies")
    parser.add_option("--out")
    parser.add_option("--inFile")
    parser.add_option("--divergence", default=0.0121, type=float)
    args = parser.parse_args()
    main(args)