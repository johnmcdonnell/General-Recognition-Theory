__author__ = "John McDonnell"
__copyright__ = "Copyright (C) 2011 John McDonnell"

# Questionable: angle conversions.
# Questionable: alpha negloglike

import os
import numpy as np
import numpy.linalg as npla
import scipy.stats as stat
import scipy.optimize as optimize

# Globals
SINIT = 10 # initial noise param for the Fisher fitting method.
ZLIMIT = 7  # z-scores past this value will be truncated.

##################################################
## Helper functions
##################################################
def angle2xy(anglevec):
    """
    For a vector of angles, finds their x/y location on a unit circle.
    """
    return np.cos(anglevec), np.sin(anglevec)

def xy2angle(xyvec):
    """
    For a unit vector of x/y coordinates, finds their theta (in polar coords).
    """
    if xyvec[1] > 0:
        return np.arccos( xyvec[0] )
    else:
        return 2*np.pi - np.arccos( xyvec[0] )

def norm_old_params(oldparams):
    """
    WARNING: May not work in the 3d case.
    """
    # TODO: Test if this is right.
    return oldparams / npla.norm( oldparams[1:-1], 2 )

def old2new_2dparams(oldparams):
    """
    WARNING: May not work in the 3d case.
    """
    # TODO: Test if this is right.
    return oldparams / npla.norm( oldparams[1:-1], 2 )

def findindices(data):
    respset = list( set( data[:,-1] ) )
    assert len( respset ) == 2
    respset.sort()
    
    aindices = data[:,-1] == respset[0]
    bindices = data[:,-1] == respset[1]
    return aindices, bindices

def truncate_zscores( zscores, z_limit=7 ):
    """
    Takes a list of z_scores and removes the extreme values *in place*.
    """
    for i in xrange(len(zscores)):
        if zscores[i] < -z_limit:
            zscores[i] = -z_limit
        elif zscores[i] > z_limit:
            zscores[i] = z_limit

##################################################
## Statistical functions
##################################################
def aic( negloglike, r ):
    """
    The famouse Akaike information criterion.
    r is the number of parameters in the model.
    """
    return 2 * ( negloglike + r )

def bic( negloglike, r, n ):
    """
    The Bayesian information criterion.
    r is the number of parameters in the model.
    n is the number of observations.
    """
    return 2 * negloglike + r * np.log(n)

def lindecisbnd( k, ma, mb ):
    # seems to be fine for 1d and 2d cases
    if k.size > 1:
        assert np.rank(k) == np.size( k, 1 ), "Covariance matrix not full rank."
        kinv = npla.inv( k )
    else:
        kinv = 1/k
    udiff = mb - ma
    a = np.inner(udiff, kinv)
    normconstant = npla.norm( a, 2 ) # Matlab norm is the 2-norm!
    anorm = a / normconstant
    bnorm = -.5 * np.inner(a, (mb+ma)) / normconstant
    return np.hstack((anorm, bnorm))

def fisherdiscrim( subjdata ):
    # seems to be fine for 1d and 2d cases.
    subjdata = np.array( subjdata )
    
    aindices, bindices = findindices( subjdata )
    
    alength = len( aindices )
    blength = len( bindices )
    if alength == 0 or blength == 0:
        midix = len(subjdata) // 2
        aindices = range(midix)
        bindices = range(midix, len(subjdata))
        alength = len( aindices )
        blength = len( bindices )
    meana = np.mean( subjdata[aindices, 0:-1], 0 )
    meanb = np.mean( subjdata[bindices, 0:-1], 0 )
    ka = np.cov( subjdata[aindices, 0:-1].T )
    kb = np.cov( subjdata[bindices, 0:-1].T )
    kpooled = ((alength-1)*ka + (blength-1)*kb) / (alength+blength-2)
    
    return lindecisbnd( kpooled, meana, meanb )

##################################################
## GLC specific functions
##################################################
def var_from_b( b, data ):
    """
    Takes the b variable and the data and returns the variance as dictated by
    Equation 15 in Ashby (1992).
    """
    return np.inner( np.inner( b, np.cov( data[:,:-1].T ) ), b )

def negloglike_glc( params, data, z_limit=7 ):
    """
    Computes the negative log likelihood for a linear discriminant with the
    given parameters for the given dataset.
    
    Datapoints should be in the form [dim1, dim2, ..., label], where label is
    either 0 or 1.
    
    Params should be in the form [ sigma, b0, b1, ..., c0 ]
    """
    mu_hx = np.inner( params[1:-1], data[:,:-1] ) + params[-1]
    zscores = mu_hx / np.sqrt( params[0] )
    
    # Truncate extreme z-scores (this function is in place)
    truncate_zscores( zscores, z_limit )
    
    aindices, bindices = findindices( data )
    
    log_a_likes = np.log( 1-stat.norm.cdf( zscores[aindices] ) )
    log_b_likes = np.log( stat.norm.cdf( zscores[bindices] ) )
    
    negloglike = -(np.sum(log_a_likes)+np.sum(log_b_likes))
    return negloglike

def negloglike_reduced( params, data, z_limit=7, a=None ):
    params = np.array( params )
    dims = len(data[0])-1
    
    if dims==2 and len(params) == 3:
        # Convert param[1] from angle to a1/a2 notation in 2d case.
        params = np.hstack([params[0], angle2xy( params[1] ), params[-1]])
    elif dims==1 and a:
        # a param fixed in the 1d case
        params = np.hstack([ params[0], a, params[1] ])
    
    if dims==3:
        raise Exception, "More than 2 dimensions not yet supported."
    #sigma = np.sqrt( var_from_b( z_coefs[:-1], data[:,:-1] ) )
    return negloglike_glc( params, data, z_limit )

def get_params_from_alpha( alpha, data ):
    labels = data[:,-1]
    stims = data[:,:-1]
    
    aindices, bindices = findindices( data )
    
    # Finding sample statistics
    # NOTE: this will break with sample size <= 2
    acov = np.cov( stims[ aindices ].T )
    bcov = np.cov( stims[ bindices ].T )
    amu = np.mean( stims[ aindices ], 0 )
    bmu = np.mean( stims[ bindices ], 0 )
    
    # Solving for 'V' in Fukunaga (1990) Equation 4.45.
    # Declaring the params of 'solve' to be matrices makes it work even in the
    # 1d case.
    b = np.array( np.linalg.solve( np.matrix(alpha * acov + (1-alpha) * bcov),
                                   np.matrix(bmu - amu).T
                                 )
                ).ravel()
    
    # Solving for sigma_i, Fukunaga (1990) Equation 4.20
    sigha = np.inner(np.inner(b, acov), b)
    sighb = np.inner(np.inner(b, bcov), b)
    
    # Solving for v_0, Fukunaga (1990) Equation 4.46
    c0 = -(alpha * sigha * np.inner( b, bmu) + \
         (1-alpha) * sighb * np.inner(b, amu)) /  \
         (alpha*sigha + (1-alpha)*sighb)
    
    # To plot each intermediate fit line during optimization, uncomment this:
    #plotline( b[0], b[1], c0 )
    return b, c0

def negloglike_alpha( alpha, data, z_limit=7 ):
    """
    Single-parameter minimization suggested by Ashby (1992), explained in more
    detail in Fukunaga(1990).
    """
    b, c0 = get_params_from_alpha( alpha, data )
    variance_hx = np.inner( np.inner( b.T, np.cov(data[:,:-1].T) ), b )
    
    return negloglike_glc( np.hstack([ variance_hx, b, c0 ]), data, z_limit )

def fit_GLC_alpha( subjdata ):
    xopt, fval, ierr, numfunc = optimize.fminbound(negloglike_alpha
                                                   , 0, 1
                                                   , args=[subjdata]
                                                   , full_output=True )
    b, c0 = get_params_from_alpha( xopt, subjdata )
    sigmasq = var_from_b( b, subjdata )
    alpha_params = np.hstack([ sigmasq, b, c0 ])
    return xopt, fval, alpha_params 

def fit_GLC_allparams( subjdata, inparams, z_limit=ZLIMIT ):
    """
    b (inparams[2:-1]) is general case. (Alfonso-Reese method requires unit vector)
    """
    #thesedata=subjdata
    dims = len(subjdata[0]) - 1
    
    if dims == 2:
        x0 = inparams
        #bounds = [(.00001,np.inf), (-np.inf,np.inf),(-np.inf,np.inf)]
        optargs = (subjdata, z_limit)
    elif dims==1:
        if len( inparams ) != 3:
            raise Exception, "1d fit needs 3 starting params."
        x0 = inparams
        likfun = negloglike_glc
        #bounds = [(.00001,np.inf), (-np.inf,np.inf)]
        print "init: ", inparams
        print "running with x0 = ", x0
        
        optargs = ( subjdata, z_limit )
    
    else:
        raise Exception, "More than 2 dims not yet supported."
    
    xopt, fopt, iter, im, sm = optimize.fmin( func=negloglike_glc
                                            , x0 = x0 
                                            , args = optargs
                                            , full_output=True
                                            )
    
    if dims==1 and angleparam:
        xopt =  np.hstack([xopt[0], angleparam, xopt[1]])
    if dims==2:
        xopt = np.hstack([ xopt[0], angle2xy(xopt[1]), xopt[2] ])
    return xopt, fopt

def fit_GLC_grid( subjdata, z_limit=ZLIMIT ):
    #thesedata=subjdata
    dims = len(subjdata[0])-1
    
    bounds = [(.00001,50)] + [(-25,25) for _ in range( dims+1 ) ]
    optargs = ( subjdata, z_limit, None, True)
    
    xopt, fopt, grid, Jout = optimize.brute(func=negloglike_reduced 
                                              , ranges = bounds
                                              , args = optargs
                                              , Ns = 5
                                              , full_output=True
                                             )
    
    from scipy.ndimage.filters import minimum_filter
    from scipy.ndimage.morphology import generate_binary_structure
    
    neighborhood = generate_binary_structure( dims+2, dims+2 )
    local_mins = minimum_filter( Jout, footprint=neighborhood ) == Jout
    min_coords = np.array([ g[local_mins] for g in grid ]).T
    xoptglobal = xopt
    foptglobal = fopt
    for coords in min_coords:
        xopt, fopt, iter, im, sm = optimize.fmin(func=negloglike_reduced 
                                                  , x0 = coords
                                                  , args = optargs
                                                  , full_output=True
                                                 )
        if fopt < foptglobal:
            xoptglobal = xopt
            foptglobal = fopt
    
    return xoptglobal, foptglobal

##################################################
## Generic modelfit class, to be inherited.
##################################################
class modelfit:
    def __init__( self, subjdata, r=None ):
        self.data = np.array( subjdata )
        self.n = len( self.data )
        self.r = r
        self.negloglike = 9999999
        self.aic = None
        self.bic = None
    
    def calcnegloglike( self ):
        """
        Overload to build new model classes.
        """
        pass
    
    def calcaic( self ):
        self.aic = aic( self.negloglike, self.r )
        return self.aic 
    
    def calcbic( self ):
        self.bic = bic( self.negloglike, self.r, self.n )
        return self.bic 
    
    def getnegloglike( self ):
        return self.negloglike
    
    def getaic( self ):
        if not self.aic:
            return self.calcaic()
        else:
            return self.aic
    
    def getbic( self ):
        if not self.bic:
            return self.calcbic()
        else:
            return self.bic

##################################################
## Classes for the different models.
##################################################
class nullmodel( modelfit ):
    def __init__( self, subjdata ):
        modelfit.__init__( self, subjdata, 0 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
    
    def calcnegloglike ( self ):
        self.negloglike = -np.log( .5 ) * self.n
        self.params = {}

class interceptmodel( modelfit ):
    def __init__( self, subjdata ):
        modelfit.__init__( self, subjdata, 1 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
    
    def calcnegloglike( self ):
        resps = self.data[:,-1]
        assert len(set(resps)) <= 2, "There can only be two."
        numa = list(resps).count(resps[0])
        numb = self.n - numa
        proba = float(numa) / self.n
        self.negloglike = -(np.log(proba) * numa + np.log(1-proba) * numb)
        self.params = dict( proba=proba)

class onedfit( modelfit ):
    """
    Fits the data with a 1-dimensional GLC model, initially using the scalar
    fitting method suggested by Ashby (1992) then using that fit as a starting
    point for a Nelder-Meade descent.
    
    This method seems to be strictly better than the other methods, although it
    is not guaranteed to find a global minimum.
    """
    def __init__( self, subjdata, r=None ):
        modelfit.__init__( self, subjdata, 2 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
    
    def calcnegloglike( self ):
        logliks = []
        params = []
        for dim in range( len( self.data[0] ) -1 ):
            # Remove one of the dimensions from the dataset.
            thesedata = self.data[:, [dim, -1]]
            
            # Optimize using the one-parameter optimization method.
            alpha, foptalpha, alpha_params = fit_GLC_alpha( self.data )
            
            negglc = negloglike_glc( alpha_params, self.data )
            assert negloglike_glc( alpha_params, self.data ) == foptalpha
            
            # Minimize using hill climbing from that point:
            print "Alpha params: ", alpha_params
            xopt, fopt = fit_GLC_allparams( self.data, alpha_params  )
            print "Sigma^2 inaccuracy: ", xopt[0] - var_from_b( xopt[1:-1], self.data )
            print "Alpha negloglike: ", negloglike_glc( alpha_params, self.data )
            assert fopt < foptalpha, "%f, %f" % (fopt, foptalpha)
            
            logliks.append( fopt )
            params.append( xopt )
            #logliks.append( foptalpha )
            #params.append( alpha_params )
        self.dim = np.argmin( logliks )
        self.negloglike = logliks[self.dim]
        self.params = params[self.dim]
        self.results = dict(
            logliks = logliks, ps=params
        )
        thisxy = np.zeros(2)
        thisxy[self.dim] = self.params[1]
        self.fullparams = np.hstack([ self.params[0], thisxy, self.params[-1] ])

class twodfit( modelfit ):
    """
    Fits the data with a 2-dimensional GLC model, initially using the scalar
    fitting method suggested by Ashby (1992) then using that fit as a starting
    point for a Nelder-Meade descent.
    
    This method seems to be strictly better than the other methods, although it
    is not guaranteed to find a global minimum.
    """
    def __init__( self, subjdata, r=None ):
        modelfit.__init__( self, subjdata, 3 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
        self.dim=-999
    
    def calcnegloglike( self ):
        # Fit using the Ashby method on one parameter (alpha)
        alpha, foptalpha, alpha_params = fit_GLC_alpha( self.data )
        
        # Minimize using hill climbing from that point:
        xopt, fopt = fit_GLC_allparams( self.data, alpha_params  )
        
        # Record the results.
        self.negloglike = fopt
        self.params = xopt
        self.paramsconv = np.hstack([ xopt[0], xy2angle(xopt[1:3]), xopt[3] ])

class onedfit_fisher( modelfit ):
    """
    Fits the data with a 1-dimensional GLC model, using an optimization on 2
    paramters with the Fisher discriminant as the starting point for
    optimization.
    """
    def __init__( self, subjdata, r=None ):
        self.sinit = SINIT
        modelfit.__init__( self, subjdata, 2 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
    
    def calcnegloglike( self ):
        logliks = []
        params = []
        for dim in range( len( self.data[0] ) -1 ):
            thesedata = self.data[:, [dim, -1]]
            fisher_coeffs = fisherdiscrim(thesedata)
            #print fisher_coeffs 
            raw_params = [self.sinit] + list( fisher_coeffs )
            start_params = norm_old_params(raw_params);
            #print "Searching for fit"
            xopt, fopt = fit_GLC_allparams( thesedata, start_params )
            logliks.append( fopt )
            params.append( xopt )
        self.dim = np.argmin( logliks )
        self.negloglike = logliks[self.dim]
        self.params = params[self.dim]
        self.results = dict(
            logliks = logliks, ps=params
        )
        thisxy = np.zeros(2)
        thisxy[self.dim] = self.params[1]
        self.fullparams = np.hstack([ self.params[0], thisxy, self.params[-1] ])

class twodfit_fisher( modelfit ):
    """
    Fits the data with a 2-dimensional GLC model, using an optimization on 3
    paramters with the Fisher discriminant as the starting point for
    optimization.
    
    This method does not seem to consistently find the best fit.
    """
    def __init__( self, subjdata, r=None ):
        self.sinit = SINIT
        modelfit.__init__( self, subjdata, 3 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
        self.dim=-999
    
    def calcnegloglike( self ):
        fisher_coeffs = fisherdiscrim(self.data)
        raw_params = [self.sinit] + list( fisher_coeffs )
        start_params = norm_old_params(raw_params);
        #print "Searching for fit"
        print "Fisher params: ", start_params
        xopt, fopt = fit_GLC_allparams( self.data, start_params )
        self.negloglike = fopt
        self.params = xopt
        self.paramsconv = np.hstack([ xopt[0], xy2angle(xopt[1:3]), xopt[3] ])

class onedfit_alpha( modelfit ):
    """
    Fits the data with a 1-dimensional GLC model using the one-parameter
    fitting method suggested by Ashby (1992).
    
    This method does not seem to consistently find the best fit, although
    qualitatively it seems to be an improvement on the Fisher method (for
    example, if there is a 1D solution almost as good as the 2D solution, it
    tends to find it.
    """
    def __init__( self, subjdata, r=None ):
        modelfit.__init__( self, subjdata, 2 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
    
    def calcnegloglike( self ):
        logliks = []
        params = []
        for dim in range( len( self.data[0] ) -1 ):
            thesedata = self.data[:, [dim, -1]]
            alpha, fopt = fit_GLC_alpha( thesedata )
            xopt = get_params_from_alpha( alpha, thesedata )
            xopt = np.hstack([ var_from_b( xopt[0], thesedata ), xopt[0]])
            logliks.append( fopt )
            params.append( xopt )
        self.dim = np.argmin( logliks )
        self.negloglike = logliks[ self.dim ]
        self.params = params[ self.dim ]
        self.results = dict(
            logliks = logliks, 
            ps=params
        )
        thisxy = np.zeros(2)
        thisxy[self.dim] = self.params[1]
        self.fullparams = np.hstack([ self.params[0], thisxy, self.params[-1] ])

class twodfit_alpha( modelfit ):
    """
    Fits the data with a 2-dimensional GRT model.
    
    This method does not seem to consistently find the best fit, although
    qualitatively it seems to be an improvement on the Fisher method (for
    example, if there is a 1D solution almost as good as the 2D solution, it
    tends to find it.
    """
    def __init__( self, subjdata, r=None ):
        modelfit.__init__( self, subjdata, 3 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
        self.dim=-999
    
    def calcnegloglike( self ):
        #print "Searching for fit"
        alpha, fopt = fit_GLC_alpha( self.data )
        xopt = get_params_from_alpha( alpha, self.data )
        self.negloglike = fopt
        self.params = np.hstack(xopt )
        #self.paramsconv = np.hstack([ xopt[0], xy2angle(xopt[1:3]), xopt[3] ])

# Plotting functions
def plotdata( data ):
    import pylab as pl
    pl.close()
    ch1 = np.array([ datum for datum in data if datum[-1]==0 ])
    ch2 = np.array([ datum for datum in data if datum[-1]==1 ])
    pl.plot( ch1[:,0], ch1[:,1], 'xb')
    pl.plot( ch2[:,0], ch2[:,1], '+g')

def plotline( a, b, intercept, **plotkwargs ):
    import pylab as pl
    ax = pl.axis()
    if b == 0:
        xs = [-intercept/a, -intercept/a]
        ys = ax[2:4]
    else:
        xs = ax[:2]
        ys = [ (-a*x - intercept)/b for x in xs ]
    pl.plot( xs, ys, **plotkwargs )

# Test code.
toydata = np.array([[ 0.1, 0, 1], [0, .9, 1], [.9, 0, 0], [1.1, 1, 0], [1, 1, 1]])
toy1 = toydata[:,[0,2]]

def plotalpha( fn ):
    import pylab as pl
    data = np.loadtxt( os.popen("sed 's/ch//g' %s |  awk NF==10 | awk '{if \
                                ($4==1) { print $0} }'" % fn), usecols=[4,5,9] )
    data[:,-1] -= 1
    alphas = np.linspace( 0, 1, 100 )
    negloglikes = [ negloglike_alpha( alpha, data ) for alpha in alphas ]
    pl.plot( alphas, negloglikes )
    pl.title( fn )
    pl.grid()

def plotalpha_all(fns):
    import pylab as pl
    for i in range( len( fns ) ):
        pl.subplot( 5, 5, i )
        plotalpha( fns[i] )


def compare_fits( fns ):
    models = [onedfit, twodfit, onedfit_alpha, twodfit_alpha, onedfit_fisher, twodfit_fisher ]
    modelnames = ["onedfit", "twodfit", "onedfit_alpha", "twodfit_alpha", "onedfit_fisher", "twodfit_fisher" ]
    results = {}
    for fn in fns:
        data = np.loadtxt( os.popen("sed 's/ch//g' %s |  awk NF==10 | \
                                    awk '{if ($4==1) { print $0} }'" % fn),
                          usecols=[4,5,9] )
        data[:,-1] -= 1
        results[fn] = [ m( data ) for m in models ]
    for fn in fns:
        print fn
        print [ fit.negloglike for fit in results[fn] ]
        for modelname, fit in zip( modelnames, results[fn] ):
            print modelname
            print "params: ", fit.params
            if 'twod' in modelname:
                print "2d var should be: ", var_from_b( fit.params[1:3], data )
            print "negloglike: ", fit.negloglike
            print "bic: ", fit.bic
            print "dim: ", fit.dim
        

def test(fn, twod=True, dim=0):
    data = np.loadtxt( os.popen("sed 's/ch//g' %s |  awk NF==10 | awk '{if \
                                ($4==1) { print $0} }'" % fn), usecols=[4,5,9] )
    data[:,-1] -= 1
    
    if twod:
        plotdata( data )
        alpha, foptalpha = fit_GLC_alpha( data )
        print "Final alpha: ", alpha
        print "alphanegloglik: ", foptalpha
        alpha_params = np.hstack( get_params_from_alpha( alpha, data ) ) 
        sigma_init = var_from_b( alpha_params[:-1], data )
        alpha_params = np.hstack( [[sigma_init], alpha_params] )
        print "Params from alpha optimization: ", alpha_params
        outputs = optimize.fmin(  negloglike_reduced
                                , x0 = alpha_params
                                , args=[ insert_ones_col( data ) ]
                                , full_output=True 
                               )
        xmin=outputs[0]
        fopt=outputs[1]
        print "Params: ", xmin
        print "negloglike: ", fopt
        print "alphaerror: ", foptalpha - fopt
        #thisfit = twodfit( data )
        
        #params = thisfit.params
        #b, c0 = get_params_from_alpha( alpha, data )
        import pylab as pl
        plotline( *alpha_params[1:], linestyle='--', color='k', label="Result of alpha opt" )
        plotline( *xmin[1:], linestyle='-', color='r', label="Final result." )
        pl.legend(loc=0)
        pl.title( fn )
    else:
        #plotdata( data )
        data = data[:,[dim,2]]
        alpha, foptalpha = fit_GLC_alpha( data )
        print "Final alpha: ", alpha
        print foptalpha
        #thisfit = twodfit( data )
        alpha_params = np.hstack( get_params_from_alpha( alpha, data ) ) 
        sigma_init = var_from_b( alpha_params[:-1], data )
        alpha_params = np.hstack( [[sigma_init], alpha_params] )
        print "Params from alpha optimization: ", alpha_params
        outputs = optimize.fmin(  negloglike_reduced
                                , x0 = alpha_params
                                , args=[ insert_ones_col( data ) ]
                                , full_output=True 
                               )
        xmin=outputs[0]
        fopt=outputs[1]
        print "Params: ", xmin
        print "negloglike: ", fopt
        print "alphaerror: ", foptalpha - fopt
        #thisfit = twodfit( data )
        
        #params = thisfit.params
        b, c0 = get_params_from_alpha( alpha, data )
        #plotline( b[0], b[1], c0 )
    #return thisfit


fns = [ 'data/'+fn for fn in os.listdir('data') if fn[-4:]==".dat" ]

