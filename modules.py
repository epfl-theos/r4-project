#!/usr/bin/env python
# coding: utf-8

# In[3]:


import aiida
from aiida import orm, load_profile
from aiida.orm import StructureData
import sys, os
from aiida.orm.utils import load_entity, load_code, load_computer, load_group, load_node
from collections import Counter
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from aiida.orm import QueryBuilder

from collections import Counter
import pandas as pd
from pandas.plotting import scatter_matrix
import seaborn as sns
import pandas as pd
from pandas.plotting import scatter_matrix
import ase
from ase.io import read, write
from ase.atoms import Atoms
from ase.calculators import calculator
from ase.calculators.calculator import Calculator, kpts2ndarray
from ase import Atoms
from operator import itemgetter
from ase.calculators.abinit import Abinit

# Librascal
from rascal.representations import SphericalInvariants as SOAP

sys.path.append('./kernel-tutorials/')
# Local Utilities for Notebook
from utilities.general import FPS, center_matrix, normalize_matrix, load_variables
from utilities.kernels import linear_kernel, gaussian_kernel, center_kernel
from utilities.general import load_variables, get_stats
from utilities.plotting import (
    plot_base, 
    plot_projection,
    plot_regression,
    plot_simple,
    get_cmaps,
    table_from_dict,
    check_mirrors,
)
from utilities.general import FPS, center_matrix, normalize_matrix, load_variables
from utilities.kernels import linear_kernel, gaussian_kernel, center_kernel
cmaps = get_cmaps()
plt.style.use("./kernel-tutorials/utilities/kernel_pcovr.mplstyle")
dbl_fig = (2 * plt.rcParams["figure.figsize"][0], plt.rcParams["figure.figsize"][1])

#skcosmo

from skcosmo.preprocessing import SparseKernelCenterer as SKC
from skcosmo.decomposition import PCovR
import random

from sklearn.datasets import load_iris
from sklearn.linear_model import Ridge
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.preprocessing import StandardScaler
from chemiscope import write_input
import json
from cur import cur_decomposition
from numpy.linalg import pinv as inv

# import pymatgen
import freud
import aiida
import ase
from aiida import orm, load_profile 
import matplotlib.pyplot as plt
from aiida.orm import QueryBuilder
from collections import Counter
import pandas as pd
from pandas.plotting import scatter_matrix
import seaborn as sns
import ase
from math import pi
import numpy as np
from ase.atoms import Atoms
from ase.calculators import calculator
from ase.calculators.calculator import Calculator, kpts2ndarray
from ase import Atoms
from operator import itemgetter 
from sklearn.preprocessing import LabelEncoder
from sklearn import preprocessing 
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from mpl_toolkits.mplot3d import Axes3D
import random
from aiida.orm import ArrayData
from aiida.orm import TrajectoryData
from aiida.orm import StructureData
from ase import Atoms
import matplotlib
from matplotlib import cm
from matplotlib.colors import Colormap as cmap
import matplotlib.cm
import numpy as np
import rowan
# import plato.draw.fresnel
from scipy.spatial import Voronoi, voronoi_plot_2d, ConvexHull, convex_hull_plot_2d, SphericalVoronoi
import landlab
from landlab import VoronoiDelaunayGrid
from landlab.grid.voronoi  import simple_poly_area
from scipy.spatial import Delaunay
from landlab import RasterModelGrid
from shapely.geometry import LineString
from shapely.ops import polygonize, unary_union


# In[4]:


"""CUR decomposition:matrix approximation of a matrix A in matrices C, U, and R such that:
    C is made from columns of A, 
    R is made from rows of A, 
    the product CUR closely approximates A."""

def decision(probability):
    return random.random() < probability
def colselect(A, k , row = False , eps =1):
    c = ( k*np.log(k))/(eps*eps)
    m, n = A.shape[0], A.shape[1]
    u, s, vh = np.linalg.svd(A,full_matrices= False )
    vh = vh[:k,:]
    probs = (1/ k)*(vh**2).sum(axis =0)
    probs = [min(1,c*p) for p in probs]
    idxs = [decision(p) for p in probs]
    cols= A[:,idxs]
    included_idx= [i for i in range(n) if idxs[i]]
    if row :    
        return cols.T, included_idx
    return cols, included_idx

def cur_decompose (A,k,e=1, return_idx = False ) :
    m, n = A.shape[ 0 ] , A.shape[ 1 ]
    if k>min (m, n) :
        return [ ] , [ ] , [ ]
    C, included_cols = colselect (A, k , False , eps=e )
    R, included_rows = colselect(A.T, k , True , eps=e )
    U = inv(C) @ A @ inv(R)
    if return_ixd:
        return C,U,R, included_cols, included_rows
    return C,U,R

def  give_cur_vals(A, k, N=10):
    c,u,r,cols,rows = cur_decompose(A,k,return_idx=True)
    err=give_error(A, c@u@r)
    for i in range (N):
        ctmp, utmp, rtmp ,c_tmp, r_tmp=cur_decompose(A, k, return_idx=True)
        err_temp=give_error(A, ctmp@utmp@rtmp)
        if err_temp<err:
            err=err_temp
            c=ctmp
            u=utmp
            r=rtmp
            cols=c_tmp
            rows=r_tmp
    return c,u,r,err,cols, rows
def give_cur_results(A, upto=10):
    erros=[]
    ks=[i for i in range(1, upto+1)]
    for k in ks:
        a,b,c,err,rows, cols= give_cur_vals(A,k)
        errors.append(err)
    return errors 

def plot_cur_error(A, upto =10):
    errs= give_cur_results (A, upto)
    x = [ i for i in range ( 1 , upto +1)]
    plt.plot(x , errs,'r-')
    plt.show() 


# In[5]:


def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)


# In[ ]:




