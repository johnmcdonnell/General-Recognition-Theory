__author__ = "John McDonnell"
__copyright__ = "Copyright (C) 2011 John McDonnell"

import os
import numpy as np
import numpy.linalg as npla
import scipy.stats as stat
import scipy.optimize as optimize
from copy import copy

# Globals
sinit = 10 # initial noise param
ZLIMIT = 7

# Helper functions
def angle2xy(anglevec):
    """
    For a vector of angles, finds their x/y location on a unit circle.
    """
    return np.cos(anglevec), np.sin(anglevec)

def xy2angle(xyvec):
    """
    For a vector of x/y coordinates, finds their theta (in polar coords).
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

def insert_ones_col(data):
    return np.column_stack( [ data[:,:-1], np.ones(len(data)), data[:,-1] ] )

# Statistical functions
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
    aindices = [ i for i in range(len(subjdata)) if subjdata[i,-1] != 0 ]
    bindices = [ i for i in range(len(subjdata)) if subjdata[i,-1] == 0 ]
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

def negloglike_glc( params, data, z_limit=7, a=None ):
    # seems to be fine for 1d and 2d cases.
    #print "Data to negloglike:", data
    #print "Params to negloglike:", params
    params = np.array( params )
    #print "Params into negloglike: ", params
    dims = len(data[0])-2
    #if any([ np.isnan(x) for x in params ]):
    #    raise Exception, "Wtf, param was nan!?!?"
    #print data
    
    if dims==2:
        #if any([ x>.0001 for x in params ]):
        #    return np.inf
        # Convert param[1] from angle to a1/a2 notation in 2d case.
        xy = angle2xy( params[1] ) # just converting angle back to a1/a2
        z_coefs = np.hstack([xy, params[2]])
        #print z_coefs
    elif dims==1:
    #    # a param fixed in the 1d case
        assert a, "a is required."
        if params[0] < .001:
            return np.inf
        params = np.hstack([ params[0], a, params[1] ])
        z_coefs = copy( params[1:] )
    else:
        raise Exception, "More than 2 dimensions not yet supported."
    
    # Normalize coefficients to units of SD:
    z_coefs /= params[0]
    #print "final zcoefs", z_coefs
    
    # Compute z-scores for data
    zscores = np.inner( np.matrix(data[:,:-1]), z_coefs )
    #print "zscores", zscores
    
    # Truncate extreme z-scores
    zscoreold = copy(zscores)
    zscores[zscores < (- z_limit)] = -z_limit
    zscores[zscores > (z_limit)] = z_limit
    #print "final zscores", zscores
    
    labels = data[:,-1]
    assert all( [label in [0,1] for label in labels ] )
    aindices = labels == 1
    bindices = labels == 0
    
    log_a_probs = np.log( 1-stat.norm.cdf( zscores[aindices] ) )
    if any([not np.isfinite(x) for x in log_a_probs ]):
        #print "caught div by zero."
        return np.inf
    
    log_b_probs = np.log( stat.norm.cdf( zscores[bindices] ) )
    if any([not np.isfinite(x) for x in log_b_probs ]):
        #print "caught div by zero."
        return np.inf
    
    negloglike = -(np.sum(log_a_probs)+np.sum(log_b_probs))
    #print "negloglike Params: ", params
    #print "calculated negloglik: ", negloglike 
    return negloglike

def getparams( alpha, data ):
    labels = data[:,-1]
    stims = data[:,:-1]
    assert all( [label in [0,1] for label in labels ] )
    aindices = labels == 1
    bindices = labels == 0
    acov = np.cov( stims[ aindices ].T )
    bcov = np.cov( stims[ bindices ].T )
    amu = np.mean( stims[ aindices ], 0 )
    bmu = np.mean( stims[ bindices ], 0 )
    # Declaring the params of 'solve' to be matrices makes it work even in the
    # 1d case.
    b = np.linalg.solve( np.matrix(alpha * acov + (1-alpha) * bcov), 
             np.matrix(bmu - amu ) )
    bsigAb = np.inner(np.inner(b, acov), b)
    bsigBb = np.inner(np.inner(b, bcov), b)
    c0 = -(alpha * bsigAb * np.inner( b, bmu) + \
         (1-alpha) * bsigBb * np.inner(b, amu)) /  \
         (alpha*bsigAb + (1-alpha)*bsigBb)
    
    #plotline( b[0], b[1], c0 )
    return b, c0

def negloglike_alpha( alpha, data, z_limit=7 ):
    b, c0 = getparams( alpha, data )
    
    labels = data[:,-1]
    stims = data[:,:-1]
    aindices = labels == 1
    bindices = labels == 0
    mu = np.inner( b, stims ) + c0  # NOTE: I'm very uncertain about this.
    stdeviation = np.inner( np.inner( b, np.cov(stims.T) ), b )
    zscores = mu / np.sqrt( stdeviation ) # TODO: subtracting 1 here does a lot, why? -1 # hmmm
    #zscores = np.array([ (stim-mu) / np.sqrt(stdeviation) for stim in stims ])
    for i in xrange(len(zscores) ):
        if zscores[i] > z_limit:
            zscores[i] = z_limit
        elif zscores[i] < -z_limit:
            zscores[i] = -z_limit
    
    log_a_probs = np.log( 1-stat.norm.cdf( zscores[aindices] ) )
    log_b_probs = np.log( stat.norm.cdf( zscores[bindices] ) )
    print "Loga, then logb"
    print zscores[aindices]
    print zscores[bindices]
    
    negloglike = -(np.sum(log_a_probs)+np.sum(log_b_probs))
    return negloglike

def fit_GLC_alpha( subjdata ):
    xopt, fval, ierr, numfunc = optimize.fminbound(negloglike_alpha
                                                   , 0, 1
                                                   , args=[subjdata]
                                                   , full_output=True )
    return xopt, fval


def fit_GLC( subjdata, inparams, z_limit=ZLIMIT ):
    #thesedata=subjdata
    thesedata = insert_ones_col( subjdata )
    dims = len(inparams) - 2
    
    print inparams
    if dims == 2:
        x0 = [inparams[0], xy2angle(inparams[1:3]), inparams[-1]]
        #bounds = [(.00001,np.inf), (-np.inf,np.inf),(-np.inf,np.inf)]
        optargs = (thesedata, z_limit)
    elif dims==1:
        angleparam = inparams[1]
        #x0 = (inparams[0], inparams[2])
        #bounds = [(.00001,np.inf), (-np.inf,np.inf)]
        print "init: ", inparams
        a = inparams[1]
        x0 = inparams[[0,2]]
        print "running with x0 = ", x0
        
        optargs = ( thesedata, z_limit, a)
    
    else:
        raise Exception, "More than 2 dims not yet supported."
    
    xopt, fopt, iter, im, sm = optimize.fmin(func=negloglike_glc 
                                                   , x0 = x0 
                                                   , args = optargs
                                                   , full_output=True
                                                  )
    
    
    if dims==1:
        xopt =  np.hstack([xopt[0], a, xopt[1]])
    if dims==2:
        xopt = np.hstack([ xopt[0], angle2xy(xopt[1]), xopt[2] ])
    return xopt, fopt

# Generic modelfit class, to be inherited.
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

# Classes for the different models.
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
    Fits the data with a 1-dimensional GRT model.
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
            fisher_coeffs = fisherdiscrim(thesedata)
            #print fisher_coeffs 
            raw_params = [sinit] + list( fisher_coeffs )
            start_params = norm_old_params(raw_params);
            #print "Searching for fit"
            xopt, fopt = fit_GLC( thesedata, start_params )
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

class twodfit( modelfit ):
    """
    Fits the data with a 2-dimensional GRT model.
    """
    def __init__( self, subjdata, r=None ):
        modelfit.__init__( self, subjdata, 3 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
    
    def calcnegloglike( self ):
        logliks = []
        params = []
        fisher_coeffs = fisherdiscrim(self.data)
        print "Fisher: ", fisher_coeffs 
        raw_params = [sinit] + list( fisher_coeffs )
        print "raw_params", raw_params
        start_params = norm_old_params(raw_params);
        print "start_params", start_params
        #print "Searching for fit"
        xopt, fopt = fit_GLC( self.data, start_params )
        self.negloglike = fopt
        self.params = xopt
        self.paramsconv = np.hstack([ xopt[0], xy2angle(xopt[1:3]), xopt[3] ])

class onedfit_alpha( modelfit ):
    """
    Fits the data with a 1-dimensional GRT model.
    TODO: fill me in!
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
            xopt = getparams( alpha, thesedata )
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
    """
    def __init__( self, subjdata, r=None ):
        modelfit.__init__( self, subjdata, 3 )
        self.calcnegloglike()
        self.calcaic()
        self.calcbic()
    
    def calcnegloglike( self ):
        logliks = []
        params = []
        #print "Searching for fit"
        alpha, fopt = fit_GLC_alpha( self.data )
        xopt = getparams( alpha, self.data )
        self.negloglike = fopt
        self.params = xopt
        #self.paramsconv = np.hstack([ xopt[0], xy2angle(xopt[1:3]), xopt[3] ])

# Plotting functions
def plotdata( data ):
    import pylab as pl
    pl.close()
    ch1 = np.array([ datum for datum in data if datum[-1]==0 ])
    ch2 = np.array([ datum for datum in data if datum[-1]==1 ])
    pl.plot( ch1[:,0], ch1[:,1], 'xb')
    pl.plot( ch2[:,0], ch2[:,1], '+g')

def plotline( a, b, intercept ):
    import pylab as pl
    ax = pl.axis()
    if b == 0:
        xs = [-intercept/a, -intercept/a]
        ys = ax[2:4]
    else:
        xs = ax[:2]
        ys = [ (-a*x - intercept)/b for x in xs ]
    pl.plot( xs, ys )


 
# Test code.
toydata = np.array([[ 0.1, 0, 1], [0, .9, 1], [.9, 0, 0], [1.1, 1, 0], [1, 1, 1]])
toyconv = insert_ones_col( toydata )
toy1 = toyconv[:,[0,2,3]]
def test(fn, twod=True):
    data = np.loadtxt( os.popen("sed 's/ch//g' %s |  awk NF==10 | awk '{if \
                                ($4==1) { print $0} }'" % fn), usecols=[4,5,9] )
    data[:,-1] -= 1
    
    twod = False
    if twod:
        plotdata( data )
        alpha, fopt = fit_GLC_alpha( data )
        print alpha
        print fopt
        #thisfit = twodfit( data )
        
       #params = thisfit.params
        b, c0 = getparams( alpha, data )
        plotline( b[0], b[1], c0 )
    else:
        #plotdata( data )
        data = data[:,[0,2]]
        alpha, fopt = fit_GLC_alpha( data )
        print "Final alpha: ", alpha
        print fopt
        #thisfit = twodfit( data )
        
        #params = thisfit.params
        b, c0 = getparams( alpha, data )
        #plotline( b[0], b[1], c0 )
    #return thisfit


