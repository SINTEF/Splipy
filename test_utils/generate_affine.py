import numpy as np
import subprocess
import datetime
import os
from general import *

def do_translate(f, dim, pardim):
    obj_name   = ['', 'crv', 'surf', 'vol']
    amount = np.random.rand(dim) * 10
    amount = np.floor(amount*10)/10
    f.write('        %s2 += np.'%(obj_name[pardim]) + repr(amount) + '\n')
    return amount

def expect_translate(amount):
    result = ''
    for i in range(len(amount)):
        result += '        pt2[...,%d] -= %f\n' % (i,amount[i])
    return result

def do_scale(f, dim, pardim):
    obj_name   = ['', 'crv', 'surf', 'vol']
    amount = np.random.rand(dim) * 10 + 2
    amount = np.floor(amount*10)/10
    f.write('        %s2 *= np.'%(obj_name[pardim]) + repr(amount) + '\n')
    return amount

def expect_scale(amount):
    result = ''
    for i in range(len(amount)):
        result += '        pt2[...,%d] /= %f\n' % (i,amount[i])
    return result




# dump large sets of control points
np.set_printoptions(threshold=np.inf)
# fetch random state
state = np.random.get_state() # might consider dumping this for reproducibility
# open file and write headers
f = open('affine_test.py', 'w')
f.write('# --- Automatic generated test file  ---\n')
f.write('# Generator    : ' + os.path.basename(__file__) + '\n')
f.write('# Date         : ' + str(datetime.date.today()) + '\n')
f.write('# Git revision : ' + subprocess.check_output(['git', 'rev-parse', 'HEAD']))
f.write("""
import numpy as np
from GeoMod import Volume, Surface, Curve, BSplineBasis
from math import sqrt
import unittest


class TestAffine(unittest.TestCase):
""")

evaluate  = [None, evaluate_curve, evaluate_surface, evaluate_volume]
do        = [do_translate,     do_scale    ]
expect    = [expect_translate, expect_scale]
name      = ['translate',      'scale'     ]
for baseP in [2,5]:
    for dim in [2,3]:
        for rational in [True, False]:
            for periodic in [-1,0,1]:
                for pardim in range(1,4):
                    for j in range(len(do)):
                        p  = np.array([baseP]*pardim) + np.random.randint(0,3,pardim)
                        if periodic >= p[0]-1:
                            continue
                        if dim < 3 and pardim > 2:
                            continue
                        n  = p + 3
                        n += np.random.randint(-2,3, pardim)
                        cp   = gen_controlpoints(n, dim, rational, periodic)
                        knot = [gen_knot(n[i], p[i], (i==0)*(periodic+1)-1) for i in range(pardim)]
                        f.write('    def test_%s_'%name[j] + get_name(n, p, dim, rational, periodic) + '(self):\n')
                        if cp.shape[0] > 30:
                            cp_str = repr(cp).replace('\n', '')
                        else:
                            cp_str = repr(cp)
                        f.write('        controlpoints = np.'+cp_str+'\n')
                        write_basis(f, p, knot, periodic)
                        write_object_creation(f, rational, pardim)
                        if name[j] is not 'scale':
                            new_dim = min(dim+np.random.randint(0,2),3) # shake it up by sometime requesting 3D-manipulation on 2D geometries
                        else:
                            new_dim = dim
                        params = do[j](f, new_dim, pardim)
                        f.write(evaluate[pardim]())
                        f.write(expect[j](params))
                        if new_dim != dim:
                            f.write('        allZero           = pt2\n')
                            f.write('        allZero[...,:-1] -= pt \n')
                            f.write('        self.assertAlmostEqual(np.max(allZero), 0.0)\n\n')
                        else:
                            f.write('        self.assertAlmostEqual(np.max(pt-pt2), 0.0)\n\n')

f.write("""
if __name__ == '__main__':
    unittest.main()
""")
print('affine_test.py written')
