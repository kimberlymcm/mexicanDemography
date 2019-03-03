"""
    Run the best performing 4 population dadi model for the bootstrapping results
"""

import sys
sys.path.append('/home/kfm/kfm_projects/NA/mexicanDemography/scripts/python')

import argparse

import numpy

import dadi_hacks
import nat_model_SG
import marginal_optimization
import dadi

from nm_utils import NSPEC, PTS_L, POP_DICT, compare_to_self, compare_to_base


def create_fss(dd, divergence):
    """ Create true frequency spectra """

    fux_table = "/home/kfm/kfm_projects/NA/NA_data/getIntrons/dadi_other/fux_table_" + str(divergence) + ".dat"
    fss = []
    fs = dadi.Spectrum.from_data_dict_corrected(dd, ["TAR", "MYA"], [POP_DICT["TAR"], POP_DICT["MYA"]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, ["TAR", "TRQ"], [POP_DICT["TAR"], POP_DICT["TRQ"]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, ["HUI", "MYA"], [POP_DICT["HUI"], POP_DICT["MYA"]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, ["HUI", "TRQ"], [POP_DICT["HUI"], POP_DICT["TRQ"]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, ["TAR", "HUI"], [POP_DICT["TAR"], POP_DICT["HUI"]], fux_table)
    fss.append(fs.copy())
    fs = dadi.Spectrum.from_data_dict_corrected(dd, ["MYA", "TRQ"], [POP_DICT["MYA"], POP_DICT["TRQ"]], fux_table)
    fss.append(fs.copy())
    return fss


def output_results(outFile, result):

    print result
    f = open(outFile, "a")
    f.write(result)
    f.close()


def main(args):

    # TODO (kmcmanus): Break this function into smaller pieces
    outFile = "../../results/bootstrapping/results_0.09Tb/20190224_bootstrap." + str(args.out)
    inFile = "../../data/bootstrapping/input_bootstrap/Bootstrap_" + str(args.inFile) + ".txt"

    dd = dadi.Misc.make_data_dict(inFile)
    fss = create_fss(dd, args.divergence)

    nss = [fs.sample_sizes for fs in fss]

    # Prepare the custom model
    func = nat_model_SG.Kimberly_method_fewerSizes

    #             (Tb,   Nb,   T_s1, T_TarHui, T_TrqMya, N_Tar, N_Hui, N_Trq, N_Mya)
    startparams = [0.09, 0.25, 0.01, 0.01,     0.01,     0.2,   0.2,   0.2,   0.2]
    fixedparams = [0.09, None, None, None,     None,     None,  None,  None,  None]
    upper_bound = [1,    0.5,  0.02, 0.02,     0.02,     0.5,   0.5,   0.5,   0.5]
    lower_bound = [0,    0.1,  1e-7, 1e-7,     0.0003,   0.05,  0.05,  0.05,  0.05]

    # Make the extrapolating version of our demographic model function.
    func_ex = dadi_hacks.make_extrap_func_multispec(func, extrap_log=True)
    print "running model"
    model = func_ex(startparams, nss, PTS_L) # Calculate the model AFS
    print "number of reads in the data" # Perform a sanity check
    print [fs.sum() for fs in fss]

    compare_to_self(fss)

    compare_to_base(model, fss)

    ratios = [fs.data.sum() / len(dd) for fs in fss]

    popts = []
    lls = []
    lthetas=[]
    for i in range(4):
        print(i)

        # Slightly perturb params
        p0 = dadi.Misc.perturb_params(startparams, lower_bound=lower_bound, upper_bound=upper_bound, fold=2)
        for j in range(2):
            # Jan 2019: Looks like I run it twice just to optimize the params a little more.	
            popt = marginal_optimization.optimize_log_fmin(p0, fss, func_ex, PTS_L,
                                                           fixed_params=fixedparams,
                                                           lower_bound=lower_bound,
                                                           upper_bound=upper_bound,
                                                           nmarginals=NSPEC)
            popt = marginal_optimization.optimize_log_fmin(popt, fss, func_ex, PTS_L,
                                                           fixed_params=fixedparams,
                                                           lower_bound=lower_bound,
                                                           upper_bound=upper_bound,
                                                           nmarginals=NSPEC)
            model = func_ex(popt, nss, PTS_L)
            thetas = map(lambda i: dadi.Inference.optimal_sfs_scaling(model[i], fss[i]), range(NSPEC))
            lthetas.append(thetas)
            ll_model = numpy.sum(map(lambda i: dadi.Inference.ll_multinom(model[i], fss[i]), range(NSPEC)))
            lls.append(ll_model)
            popts.append(popt)

    # Jan 2019: It seems like the lines below should not have been indented, but I guess this way still works.
    best = lls.index(max(lls))
    ll_model = func_ex(popts[best], nss, PTS_L) # Converting parameters
    thetas = map(lambda i: dadi.Inference.optimal_sfs_scaling(model[i], fss[i]), range(NSPEC))
    effectTheta = (numpy.array(lthetas[best]) / numpy.array(ratios)).mean()
    # print "effect theta", effectTheta
    # print "params" + " ".join(map(str, popt)) + "\n"
    # end = nat_model_SG.Kimberly_convert_params_fewerSizes(popts[best],theta=effectTheta, L=8880201, mu=1.25e-8)
    result = str(lls[best]) + "\t" + str(effectTheta) + "\t" + "\t".join(map(str, popts[best])) + "\n"

    output_results(outFile, result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create fux table")
    parser.add_argument("--out", default="1")
    parser.add_argument("--inFile", default="1")
    parser.add_argument("--divergence", default=0.0121, type=float)
    args = parser.parse_args()
    main(args)

