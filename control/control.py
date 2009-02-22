import numpy as npy
import slicot

def lyap(A, C):
    """Solves the continous-time Lyapunov equation
    
    $ A^T X + X A + C = 0 $
    
    Parameters
    ----------

    A : array_like, shape (N, N), real
       A matrix

    C : array_like, shape (N, N), symmetric, real
       C matrix

    Returns
    -------

    x : ndarray, shape (M,M)
       Solution matrix X

    Examples
    --------

    >>> lyap(npy.array([[1., 2.], [3., 4.]]), npy.array([[5., 6.], [6., 8.]]))
    array([[-0.1, -0.8],
           [-0.8, -0.6]])

    Raises
    ------
    LinAlgError

       
    """
    
    a,u,c,scale,sep,ferr,wr,wi,info = slicot.sb03md(dico='C', 
                                                    job='X', 
                                                    fact='N', 
                                                    trana='N', 
                                                    a=A, 
                                                    c=-C)

    if info != 0:
        raise npy.linalg.LinAlgError()

    return c

def gram(sys, gramtype):
    """
    Computes the controllability and observability grammians of an lti system

    
    Parameters
    ----------

    sys : signal.ltisys.lti object
        lti system in state-space form

    gramtype: string
        'c' or 'o' to compute either controllability or observability grammians

    Example
    -------
    
    >>> import scipy.signal.ltisys as ltisys
    >>> s = ltisys.lti(npy.mat([[-1,0],[0,-1]]), npy.mat([[1],[1]]), \
        npy.mat([[1, -1]]), 0)
    >>> gram(s, 'c')
    array([[ 0.5,  0.5],
           [ 0.5,  0.5]])
    >>> gram(s, 'o')
    array([[ 0.5, -0.5],
           [-0.5,  0.5]])
    
    """
    if gramtype == 'c':
        return lyap(sys.A.T, npy.dot(sys.B, sys.B.T))
    elif gramtype == 'o':
        return lyap(sys.A, npy.dot(sys.C.T, sys.C))
