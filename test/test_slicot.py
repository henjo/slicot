import numpy as npy
import slicot

from nose.tools import *
from numpy.testing import assert_array_almost_equal

def test_sb03md():
    A = npy.array([[1., 2.], [3., 4.]])
    C = -npy.array([[5., 6.], [6., 8.]])

    a,u,c,scale,sep,ferr,wr,wi,info = slicot.sb03md(dico='C', job='X', fact='N', trana='N', a=A, c=C)
    
    assert_equal(info, 0)
    
    ## Test if eigenvalues of A are correct
    Alambda,Aeigvec = npy.linalg.eig(A)
    assert_array_almost_equal(wr, Alambda)
    
    ## Test result
    assert_array_almost_equal(c, npy.array([[-0.1, -0.8],[-0.8, -0.6]]))
    
    

