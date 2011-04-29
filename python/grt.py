__author__ = "John McDonnell"
__copyright__ = "Copyright (C) 2011 John McDonnell"


import numpy as np
import numpy.linalg as npla
import scipy.stats as stat
import scipy.optimize as optimize

# Globals
sinit = 10 # Noise param of some kind?

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
        return np.acos( xyvec[0] )
    else:
        return 2*np.pi - np.acos( xyvec[0] )

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
    if k.size > 1:
        assert np.rank(k) == np.size( k, 1 ), "Covariance matrix not full rank."
    if k.size > 1:
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

def negloglike_glc( params, data, z_limit ):
    if params[0] < .001 or params[0] > 500:
        return 999999
    
    # Convert param[1] from angle to a1/a2 notation in 2d case.
    if len(data[0])==4:
        xy = angle2xy( params[1] ) # just converting angle back to a1/a2
        z_coefs = np.hstack([xy, params[2]])
    else:
        z_coefs = params[1:]
    
    # Normalize coefficients to units of SD:
    z_coefs /= params[0]
    
    # Compute z-scores for data
    zscores = np.inner( np.matrix(data[:,:-1]), z_coefs )
    
    # Truncate extreme z-scores
    for i, score in enumerate( xrange(len(zscores)) ):
        if score < -z_limit:
            zscores[i] = -z_limit
        if score > z_limit:
            zscores[i] = z_limit
    
    aindices = [ i for i in range(len(data)) if data[i,-1] == 1 ]
    bindices = [ i for i in range(len(data)) if data[i,-1] != 1 ]
    
    log_a_probs = np.log( 1-stat.norm.cdf( zscores[aindices] ) )
    log_b_probs = np.log( stat.norm.cdf( zscores[bindices] ) )
    
    negloglike = -(np.sum(log_a_probs)+np.sum(log_b_probs))
    return negloglike

def fit_GLC( subjdata, inparams, z_limit=7 ):
    thesedata = np.column_stack( [ subjdata[:,:-1], np.ones(len(subjdata)),
                                 subjdata[:,-1] ] )
    
    xopt, fopt, iter, funcalls, warnflag = optimize.fmin( negloglike_glc, inparams, (thesedata, z_limit), xtol=.001, full_output=True)
    return xopt, fopt

def norm_old_params(oldparams):
    """
    WARNING: May not work in the 3d case.
    """
    # TODO: Test if this is right.
    return oldparams / npla.norm( oldparams[1:-1], 2 )

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

class nullmodel( modelfit ):
    def __init__( self, subjdata ):
        modelfit.__init__( self, subjdata, 0 )
        self.calcnegloglike()
    
    def calcnegloglike ( self ):
        self.negloglike = -np.log( .5 ) * self.n
        self.params = {}

class interceptmodel( modelfit ):
    def __init__( self, subjdata ):
        modelfit.__init__( self, subjdata, 1 )
        self.calcnegloglike()
    
    def calcnegloglike( self ):
        numa = np.sum( self.data[:,-1] )
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
        self.negloglike = np.min( logliks )
        self.params = dict(
            dim = self.dim,
            params = params[self.dim]
            ,logliks = logliks, ps=params
        )

def fitsubj( subjdata ):
    """
    Trial format: [x0,x1,...,y]
    
    fittype is as follows:
     0 : 0d with fixed p=0.5
     1 : 0d with floating p
     2 : 1d: x axis
     3 : 1d: y axis
     4 : 2d
    """
    pass

toydata = np.array([[ 0.1, 0, 1], [0, .9, 1], [.9, 0, 0], [1.1, 1, 0]])
def test():
    thisfit = onedfit( toydata )
    return thisfit


