
import numpy as np
import scipy.stats as stat
import scipy.optimize as optimize

# Globals
sinit = 10 # Noise param of some kind?

def aic( negloglike, r ):
    return 2 * ( negloglike + r )

def bic( negloglike, r, n ):
    return 2 * negloglike + r * np.log(n)

def lindecisbnd( k, ma, mb ):
    assert np.rank(k) == np.size( k, 1 ), "Covariance matrix not full rank."
    kinv = np.inv( k )
    udiff = mb - ma
    a = udiff  * kinv # TODO: Make sure this works in the 2d case.
    normconstant = np.norm( a )
    a /= normconstant
    b = -.5 * a * (mb+ma ) / normconstant
    return a,b

def fisherdiscrim1d( subjdata, dim ):
    aindices = [ i for i in range(len(subjdata)) if subjdata[i,-1] == 0 ]
    bindices = [ i for i in range(len(subjdata)) if subjdata[i,-1] != 1 ]
    alength = len( aindices )
    blength = len( bindices )
    if alength == 0 or blength == 0:
        midix = len(subjdata) // 2
        aindices = range(midix)
        bindices = range(midix, len(subjdata))
        alength = len( aindices )
        blength = len( bindices )
    meana = np.mean( subjdata[aindices, dim] )
    meanb = np.mean( subjdata[bindices, dim] )
    ka = np.cov( subjdata[aindices, dim] )
    kb = np.cov( subjdata[bindices, dim] )
    kpooled = ((alength-1)*ka + (blength-1)*kb) / (alength+blength-2)
    
    return lindecisbnd( kpooled, meana, meanb )

def negloglike_glc( params, data, z_limit ):
    # Not worked on yet...
    if params[0] < .001 or params[0] > 500:
        return 999999
    z_coefs = params[1:] / params[0]
    zscores = np.multiply( np.matrix(data[:,:-1]).T, z_coefs ).T
    for i, score in enumerate( xrange(len(zscores)) ):
        if score < -z_limit:
            zscores[i] = -z_limit
        if score > z_limit:
            zscores[i] = z_limit
    
    aindices = [ i for i in range(len(data)) if data[i,-1] == 0 ]
    bindices = [ i for i in range(len(data)) if data[i,-1] != 1 ]
    
    log_a_probs = np.log( 1 - stat.norm.cdf( zscores[aindices] ) )
    log_b_probs = np.log( stat.norm.cdf( zscores[bindices] ) )
    
    return -(np.sum(log_a_probs)+np.sum(log_b_probs))

def negloglike_1dglc( params, data, z_limit ):
    if params[0] < .001 or params[0] > 500:
        return 999999
    z_coefs = params[1:] / params[0]
    zscores = np.multiply( np.matrix(data[:,:-1]).T, z_coefs ).T
    for i, score in enumerate( xrange(len(zscores)) ):
        if score < -z_limit:
            zscores[i] = -z_limit
        if score > z_limit:
            zscores[i] = z_limit
    
    aindices = [ i for i in range(len(data)) if data[i,-1] == 0 ]
    bindices = [ i for i in range(len(data)) if data[i,-1] != 1 ]
    
    log_a_probs = np.log( 1 - stat.norm.cdf( zscores[aindices] ) )
    log_b_probs = np.log( stat.norm.cdf( zscores[bindices] ) )
    
    return -(np.sum(log_a_probs)+np.sum(log_b_probs))


def fit_1dGLC( subjdata, dim, inparams, z_limit ):
    thisdata = np.column_stack( [ subjdata[:,dim], np.ones(len(subjdata)),
                                 subjdata[:,-1] ] )
    
    x0 = inparams
    
    xopt, fopt, iter, funcalls, warnflag = optimize.fmin( negloglike_1dglc, x0, (thisdata, z_limit), xtol=.001, full_output=True)
    return xopt, fopt


def norm_old_params(oldparams):
    """
    WARNING: May not work in the 3d case.
    """
    # TODO: Test if this is right.
    return oldparams / np.norm( oldparams[1:-1] )

class modelfit:
    def __init__( self, subjdata, negloglikefun, r ):
        self.data = self.subjdata
        self.n = len( self.data )
        self.r = r
        self.negloglike = None
        self.aic = None
        self.bic = None
        self.negloglikefun = negloglikefun
    
    def calcnegloglike( self ):
        """
        Overload to build new model classes.
        """
        pass
    
    def calcaic( self ):
        if not self.negloglike:
            self.calcnegloglike()
        self.aic = aic( self.negloglike, self.r )
        return self.aic 
    
    def calcbic( self ):
        if not self.negloglike:
            self.calcnegloglike()
        self.bic = bic( self.negloglike, self.r, self.n )
        return self.bic 
    
    def getnegloglike( self ):
        if self.negloglike:
            return self.negloglike
        else:
            return self.calcnegloglike()
    
    def getaic( self ):
        if self.negloglike:
            return self.negloglike
        else:
            return self.calcaic()
    
    def getbic( self ):
        if self.negloglike:
            return self.negloglike
        else:
            return self.calcbic()

class nullmodel( modelfit ):
    def __init__( self, subjdata ):
        modelfit.__init__( self, subjdata, 0 )
    
    def calcnegloglike ( self ):
        if self.negloglike:
            return self.negloglike
        else:
            self.negloglike = -np.log( .5 ) * self.n
            self.params = [0]
            return self.negloglike

class interceptmodel( modelfit ):
    def __init__( self, subjdata ):
        modelfit.__init__( self, subjdata, 1 )
    
    def calcnegloglike( self ):
        if self.negloglike:
            return self.negloglike
        else:
            numa = np.sum( self.subjdata[:,-1] )
            numb = self.n - numa
            proba = float(numa) / self.n
            self.negloglike = -(np.log(proba) * numa + np.log(1-proba) * numb)
            self.params = [ proba ]
            return self.negloglike

class onedfit( modelfit ):
    def __init__( self, subjdata ):
        modelfit.__init__( self, subjdata, 2 )
    
    def fit_1dCLG( self ):
        pass
    
    def calcnegloglike( self ):
        if self.negloglike:
            return self.negloglike
        else:
            for dim in range( len( self.subjdata[0] ) -1 ):
                fisher_coeffs = fisherdiscrim1d(self.subjdata, dim)
                raw_params = [sinit] + list( fisher_coeffs )
                start_params = norm_old_1dparams(raw_params);
                #print "Searching for fit"


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
