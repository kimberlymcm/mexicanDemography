"""
   For each possible branching topology for the four native american populations,
   optimizes demographic parameters
"""


import argparse

import numpy

import dadiHacks
import marginal_optimization
import dadi
import nat_model_SG
from nm_utils import DATA_DICT_FN, NSPEC, PTS_L, compare_to_self, compare_to_base


POPS = [['TAR', 'HUI', 'TRQ', 'MYA'],
        ['TAR', 'TRQ', 'MYA', 'HUI'],
        ['TAR', 'MYA', 'HUI', 'TRQ'],
        ['HUI', 'TAR', 'TRQ', 'MYA'],
        ['HUI', 'TRQ', 'MYA', 'TAR'],
        ['HUI', 'MYA', 'TAR', 'TRQ'],
        ['TRQ', 'TAR', 'MYA', 'HUI'],
        ['TRQ', 'MYA', 'HUI', 'TAR'],
        ['TRQ', 'HUI', 'TAR', 'MYA'],
        ['MYA', 'TRQ', 'TAR', 'HUI'],
        ['MYA', 'TAR', 'HUI', 'TRQ'],
        ['MYA', 'HUI', 'TRQ', 'TAR']]

POP_SIZES = [[26, 24, 16, 14],
             [26, 16, 14, 24],
             [26, 14, 24, 16],
             [24, 26, 16, 14],
             [24, 16, 14, 26],
             [24, 14, 26, 16],
             [16, 26, 14, 24],
             [16, 14, 24, 26],
             [16, 24, 26, 14],
             [14, 16, 26, 24],
             [14, 26, 24, 16],
             [14, 24, 16, 26]]


OUTFILE_PREFIX = "20190218_branching_Tb0.09_"



def create_fss(dd, k, divergence, pop_list):
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
    k = 0
    for pop_list in POPS:

        dd = dadi.Misc.make_data_dict(DATA_DICT_FN)

        fss = create_fss(dd, k, args.divergence, pop_list)

        nss = [fs.sample_sizes for fs in fss]

        #Prepare the custom model
        func = nat_model_SG.Kimberly_method_branching

        #              (Tb,   Nb,  T_s1,  T_s2, T_s3,   N_Tar, N_Hui, N_Trq, N_Mya)
        # Starting params from the bootstrapping code
        start_params = [0.09, 0.25, 0.01, 0.01, 0.01,   0.2,   0.2,   0.2,   0.2]
        fixed_params = [0.09, None, None, None, None,   None,  None,  None,  None]
        upper_bound =  [1,    0.5,  0.02, 0.02, 0.02,   0.5,   0.5,   0.5,   0.5]
        lower_bound =  [0,    0.1,  1e-7, 1e-7, 0.0003, 0.05,  0.05,  0.05,  0.05]


        # Make the extrapolating version of our demographic model function.
        func_ex = dadiHacks.make_extrap_func_multispec(func, extrap_log=True)
        print "running model"
        model = func_ex(start_params, nss, PTS_L) # Calculate the model AFS
        print "number of reads in the data" # perform a sanity check
        print [fs.sum() for fs in fss]

        compare_to_self(fss)

        compare_to_base(model, fss)

        ratios = [fs.data.sum() / len(dd) for fs in fss]

        outFile = "20190218_branching_Tb0.09_" + str(pop_list[0]) + "_" + str(pop_list[1]) + "_" + str(pop_list[2]) + "_" + str(pop_list[3]) + ".txt"
        f = open(outFile, 'a')

        for i in range(4):
            p0 = dadi.Misc.perturb_params(start_params, lower_bound=lower_bound, upper_bound=upper_bound, fold=2)
            popts = []
            lls = []
            lthetas = []
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
                thetas = map(lambda i: dadi.Inference.optimal_sfs_scaling(model[i], fss[i]), range(NSPEC))
                lthetas.append(thetas)
                ll_model = numpy.sum(map(lambda i: dadi.Inference.ll_multinom(model[i], fss[i]), range(NSPEC)))
                lls.append(ll_model)
                popts.append(popt)

            best = lls.index(max(lls))    
            ll_model = func_ex(popts[best], nss, PTS_L)
            thetas = map(lambda i: dadi.Inference.optimal_sfs_scaling(model[i], fss[i]), range(NSPEC))
            #print "thetas are", thetas
            effectTheta = (numpy.array(lthetas[best]) / numpy.array(ratios)).mean()
            #print "effect theta", effectTheta
            #print "params" + " ".join(map(str, popt)) + "\n"
            end = nat_model_SG.Kimberly_convert_params_fewerSizes(popts[best], theta=effectTheta, gen_time=29.0, L=8889201, mu=1.25e-8)
            result = str(lls[best]) + "\t" + str(effectTheta) + "\t" + "\t".join(map(str, end)) + "\n"
            print result
            f.write(result)
        f.close()
        k = k + 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Optimize demographic parameters for branching topologies")
    parser.add_argument("--out")
    parser.add_argument("--inFile")
    parser.add_argument("--divergence", default=0.0121, type=float)
    args = parser.parse_args()
    main(args)