"""
    Additional code for dadi.
    Initial script from Simon Gravel.
    I made very few changes.
"""


import os
import logging
import functools

import numpy
import dadi

logger = logging.getLogger("Inference")


def make_extrap_func_multispec(func, extrap_x_l=None, extrap_log=False, fail_mag=10):
    """
    Generate a version of func that extrapolates to infinitely many gridpoints.
	Modification from the initial dadi file to allow multiple marginal to be returned by func

    func: A function that returns a single scalar or array and whose last
        non-keyword argument is 'pts': the number of default_grid points to use
        in calculation.  
    extrap_x_l: An explict list of x values to use for extrapolation. If not 
        provided, the extrapolation routine will look for '.extrap_x'
        attributes on the results of func. The method Spectrum.from_phi will
        add an extrap_x attribute to resulting Spectra, equal to the x-value
        of the first non-zero grid point. An explicit list is useful if you
        want to override this behavior for testing.
    fail_mag:  Simon Gravel noted that there can be numerical instabilities in
        extrapolation when working with large spectra that have very small
        entires (of order 1e-24). To avoid these instabilities, we ignore the 
        extrapolation values (and use the input result with the smallest x) 
        if the extrapolation is more than fail_mag orders of magnitude away
        from the smallest x input result.

    Returns a new function whose last argument is a list of numbers of grid
    points and that returns a result extrapolated to infinitely many grid
    points.
    """
    x_l_from_results = (extrap_x_l is None)

    def extrap_func(*args, **kwargs):
        #print args
        #print kwargs
        # Separate pts (or pts_l) from arguments
        if 'pts' not in kwargs:
            other_args, pts_l = args[:-1], args[-1]
        else:
            other_args = args
            pts_l = kwargs['pts']
            del kwargs['pts']

        if 'no_extrap' in kwargs:
            no_extrap = True
            del kwargs['no_extrap']
        else:
            no_extrap = False

        if numpy.isscalar(pts_l):
            pts_l = [pts_l]

        # Create a sub-function that fixes all other arguments and only
        # takes in pts.
        partial_func = functools.partial(func, *other_args, **kwargs)
        

        ##
        ## The commented-out code implements distribution of fs calculations
        ## among multiple processors. Unfortunately, it doesn't work with
        ## iPython, because pickling functions is very fragile in the
        ## interactive interpreter. If we had some sort of data structure for
        ## defining models, then this would be much easier to implement. Note
        ## that there might still be issues on Windows systems, because of its
        ## poor-man's version of fork().
        ##
        #import cPickle, multiprocessing
        #try:
        #    # Test whether sub-function is picklable. If not, pool.map will
        #    # hang.
        #    cPickle.dumps(partial_func)
        #    pool = multiprocessing.Pool()
        #    result_l = pool.map(partial_func, pts_l)
        #    pool.close()
        #except (cPickle.PicklingError, TypeError):
        #    print('Function passed to extrap func must be picklable for '
        #          'multi-processor executation.')
        #    import sys
        #    sys.exit()
	
        allresults = map(partial_func, pts_l)
        # number of sfss in the results
        nres=len(allresults[0])
        if no_extrap:
            return allresults

        if x_l_from_results:
            try:
            	# infer the grid properties from the first SFS. That may be a bit dangerous.
                x_l = [r[0].extrap_x for r in allresults]
            except AttributeError:
                raise ValueError("Extrapolation function error: No explicit extrapolation x_l provided, and results do not have 'extrap_x' attributes. If this is an FS extrapolation, check your from_phi method.")
        else:
            x_l = extrap_x_l
	allouts = []
	# loop over the number of spectra in the model, and perform extrapolation for each
        for nspec in range(nres):
            # within this loop, everything is as in the original make_extrap_log
            result_l = [result[nspec] for result in allresults] 
            if extrap_log:
                result_l = [numpy.log(result[nspec]) for result in allresults] 
                
       	    # Extrapolate
            if len(pts_l) == 1:
                ex_result = result_l[0]
            elif len(pts_l) == 2:
                ex_result = dadi.Numerics.linear_extrap(result_l, x_l)
            elif len(pts_l) == 3:
                ex_result = dadi.Numerics.quadratic_extrap(result_l, x_l)
            elif len(pts_l) == 4:
                ex_result = dadi.Numerics.cubic_extrap(result_l, x_l)
            elif len(pts_l) == 5:
                ex_result = dadi.Numerics.quartic_extrap(result_l, x_l)
            elif len(pts_l) == 6:
                ex_result = dadi.Numerics.quintic_extrap(result_l, x_l)
            else:
                raise ValueError('Number of calculations to use for extrapolation '
                                 'must be between 1 and 6')

            if extrap_log:
                ex_result = numpy.exp(ex_result)
    
            # Simon Gravel noted that there can be numerical instabilities in
            # extrapolation when working with large spectra that have very small
            # entires (of order 1e-24).
            # To avoid these instabilities, we ignore the extrapolation values
            # if it is too different from the input values.
            if len(pts_l) > 1:
                # Assume the best input value comes from the smallest grid.
                best_result = result_l[numpy.argmin(x_l)]
                if extrap_log:
                     best_result = numpy.exp(best_result)

                # The extrapolation is deemed to have failed if it results in a
                # value more than fail_mag orders of magnitude away from the 
                # best input value.
                extrap_failed = abs(numpy.log10(ex_result/best_result)) > fail_mag
                if numpy.any(extrap_failed):
                    logger.warn('Extrapolation may have failed. Check resulting '
                            'frequency spectrum for unexpected results.')

                # For entries that fail, use the "best" input result.
                ex_result[extrap_failed] = best_result[extrap_failed]
            allouts.append(ex_result)
             
        return allouts

    extrap_func.func_name = func.func_name
    extrap_func.func_doc = func.func_doc

    return extrap_func


def make_extrap_log_func(func, extrap_x_l=None):
    """
    Generate a version of func that extrapolates to infinitely many gridpoints.

    Note that extrapolation here is done on the *log* of the function result,
    so this will fail if any returned values are < 0. It does seem to be better
    behaved for SFS calculation.

    func: A function whose last argument is the number of Numerics.default_grid 
          points to use in calculation and that returns a single scalar or 
          array.
    extrap_x_l: An explict list of x values to use for extrapolation. If not 
         provided, the extrapolation routine will look for '.extrap_x'
         attributes on the results of func. The method Spectrum.from_phi will
         add an extrap_x attribute to resulting Spectra, equal to the x-value
         of the first non-zero grid point. An explicit list is useful if you
         want to override this behavior for testing.

    Returns a new function whose last argument is a list of numbers of grid
    points and that returns a result extrapolated to infinitely many grid
    points.
    """
    return make_extrap_func(func, extrap_x_l=extrap_x_l, extrap_log=True)