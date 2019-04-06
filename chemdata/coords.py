import math
import numpy as np
from numpy.linalg import norm

import ase
from itertools import combinations

import logging


DEBUG = False

def cartesian_to_zmatrix(positions, zmatrix_dict = None,
    initial_num = 0, indices = None):
    def get_zmat_data(zmatrix_dict, keywords):
        return zmatrix_dict[keywords] if zmatrix_dict is not None \
            and keywords in zmatrix_dict else []
    shown_length  = get_zmat_data(zmatrix_dict, 'shown_length')
    shown_angle   = get_zmat_data(zmatrix_dict, 'shown_angle')
    shown_dihedral= get_zmat_data(zmatrix_dict, 'shown_dihedral')
    same_length   = get_zmat_data(zmatrix_dict, 'same_length')
    same_angle    = get_zmat_data(zmatrix_dict, 'same_angle')
    same_dihedral = get_zmat_data(zmatrix_dict, 'same_dihedral')
    natoms = len(positions)
    #shown_length = []
    #shown_angle = []
    #shown_dihedral = []
    if indices is None:
        indices = np.array(range(len(positions)))
    shown_length.sort()
    #for a0, a1, a2, a3 in in_shown_dihedral:
    #   i0, i1, i2, i3 = indices
    #for a0, a1 in shown_length:
    #   if indices[a0] > indices[a1]:
    #       indices[a0], indices[a1] = indices[a1], indices[a0]
    zmatrix = np.array([[[-1, -1], [-1, -1], [-1, -1]]]*natoms).tolist()
    same_bond_variables = [''] * len(same_length)
    variables = {}
    if isinstance(positions, ase.Atoms):
        atoms = positions.copy()
    else:
        atoms = ase.Atoms(symbols = 'X'*natoms, positions = positions)[indices]
    for ai in range(len(atoms)):
        if ai == 0:
            continue
        elif ai == 1:
            zmatrix[ai][0] = [0, atoms.get_distance(0, 1) ]
            continue
        for a0, a1 in shown_length:
            a0, a1 = indices[a0], indices[a1]
            if DEBUG: print(a0, a1)
            if ai == a1:
                alpha = 'R_'+str(a0+initial_num)+'_'+str(a1+initial_num)
                write_variable = True
                for same_length, index in zip(same_length, \
                        range(len(same_length))):
                    # print((a0, a1), same_length)
                    if (a0, a1) in same_length:
                        # print("UES")
                        if same_bond_variables[index] == '':
                            same_bond_variables[index] = alpha
                            if DEBUG: print(index, same_bond_variables)
                        else:
                            alpha = same_bond_variables[index]
                            write_variable = False
                        break
                zmatrix[ai][0] = [a0, alpha ]
                if write_variable:
                    variables[alpha] = [(a0, a1), atoms.get_distance(a0, a1)]
                break

        a0 = -1
        a1 = -1
        a2 = -1
        a0 = zmatrix[ai][0][0]
        if a0 == -1:
            a0 = 0
            dist = atoms.get_distance(ai, a0)
            if DEBUG: print('dist:', ai, a0, dist)
            zmatrix[ai][0] = [a0, dist]

        a1 = zmatrix[ai][1][0]
        if  a1 == -1:
            for a1 in range(0, ai):
                if not a1 in [a0]:
                    break
            if a1 == -1:
                raise ValueError('a1 is still -1')
            angle = atoms.get_angle(ai, a0, a1)
            if DEBUG: print('angle:', ai, a0, a1, angle)
            zmatrix[ai][1] = [a1, angle]
        a2 = zmatrix[ai][2][0]
        if ai >= 3 and a2 == -1:
            for a2 in range(0, ai):
                if not a1 in [a0, a1]:
                    break
            if a2 == -1:
                raise ValueError('a2 is still -1')
            dihedral = atoms.get_dihedral(ai, a0, a1, a2)
            if DEBUG: print('dihedral:', dihedral)
            zmatrix[ai][2] = [a2, dihedral]
    if initial_num != 0:
        for zmat in zmatrix:
            for zmat_x in zmat:
                if zmat_x[0] != -1:
                    zmat_x[0] += initial_num
    if DEBUG: print(zmatrix, variables, indices)
    return zmatrix, variables, indices





def cos(theta, arc = False):
    factor = 1 if arc else math.pi/180.0
    return math.cos(theta * factor)


def sin(theta, arc = False):
    factor = 1 if arc else math.pi/180.0
    return math.sin(theta * factor)

def acos(result, arc = False):
    factor = 1 if arc else 180.0/math.pi
    return math.acos(result) * factor


def asin(result, arc = False):
    factor = 1 if arc else 180.0/math.pi
    return math.asin(result) * factor




def cartesian_to_spherical(pos_o, pos_s):
    pos_o = np.array(pos_o)
    pos_s = np.array(pos_s)
    if DEBUG: print('cartesian to spherical:', pos_o, pos_s)
    v_os = pos_s - pos_o
    if norm(v_os) < 0.01:
        return (0, 0, 0)
    x, y, z = v_os
    length = np.linalg.norm(v_os)
    theta = acos(z/length)
    xy_length = math.sqrt(x*x+y*y)
    if DEBUG: print('xy_length', theta, xy_length)
    if xy_length <  0.05:
        phi_x = 0.0
        phi_y = 0.0
    else:
        phi_x = acos(x/xy_length)
        phi_y = asin(y/xy_length)
    if y>=0: phi =  phi_x
    else:    phi = -phi_x
    return (length, theta, phi)


def spherical_to_cartesian(pos_o, length, space_angle, space_angle0 = (0, 0)):
    theta , phi  = space_angle
    theta0, phi0 = space_angle0
    if DEBUG: print('sperical to cartesian:', theta, phi)
    pos_site = np.array(pos_o) + length * \
        np.array([sin(theta+theta0) * cos(phi+phi0), \
                     sin(theta+theta0) * sin(phi+phi0), \
                     cos(theta+theta0)])
    return pos_site


def rotate_site_angle(site_angle, theta, phi):
    for site_angle_i in site_angle:
        theta_i, phi_i = site_angle_i
        site_angle_i = [theta_i+theta, phi_i+phi]
    return site_angle


def transform_ijk(pos_o, pos_i, vector, direction = 'to_ijk'):
    assert(direction in ['to_ijk', 'to_xyz'])
    pos_o = np.array(pos_o)
    pos_i = np.array(pos_i)
    if direction == 'to_ijk':
        return np.dot(pos_i - pos_o, np.asarray(np.mat(vector).T))
    elif direction == 'to_xyz':
        return np.dot(pos_i, np.array(vector)) + pos_o


def get_cartesian_ijk(pos_o, pos_z, pos_x = None):
    pos_o = np.array(pos_o)
    pos_z = np.array(pos_z)
    if not pos_x is None:
        pos_x = np.array(pos_x)
    v_z = pos_z - pos_o
    if norm(v_z) < 0.01:
        # print('get_cartesian_ijk norm too small, maybe linear')
        return np.eye(3)
    k = v_z / np.linalg.norm(v_z)
    if pos_x is None:
        pos_x = np.array([np.random.random(), np.random.random(), np.random.random()])
        #pos_x = np.array([100, 80, 60])
    v_ox = pos_x - pos_o
    v_x = v_ox - np.dot(v_ox, k) * k
    if ((v_x==np.array([0,0,0])).all()):
        pos_x = np.array([np.random.random(), np.random.random(), np.random.random()])
        v_ox = pos_x - pos_o
        v_x = v_ox - np.dot(v_ox, k) * k
    i = v_x / np.linalg.norm(v_x)
    j = np.cross(k, i)
    return (i, j, k)


class LG_transform_obj(object):
    def __init__(self, posO, posZ = None, posX=None, vector = None):
        import copy
        if posZ is None and vector is None:
            raise ValueError('You need to give posZ or vector')
        self.posO = copy.copy(posO)
        if vector is not None:
            self.vector = vector
        else:
            self.vector = get_cartesian_ijk(posO, posZ, posX)
    def get_local(self, pos_global):
        return transform_ijk(self.posO, pos_global, self.vector)
    def get_global(self, pos_local):
        return transform_ijk(self.posO, pos_local, self.vector, 'to_xyz')

def get_cartesian_ijk_with_sphere(pos0, pos1, pos2, angle1, angle2):
    vec1  = pos1 - pos0
    vec2  = pos2 - pos0
    cvec1 = vec1/norm(vec1)
    cvec2 = vec2/norm(vec2)
    cvec3 = np.cross(cvec1, cvec2)
    svec1 = spherical_to_cartesian([0,0,0],1,angle1)
    svec2 = spherical_to_cartesian([0,0,0],1,angle2)
    svec3 = np.cross(svec1, svec2)
    angle_diff = abs(vector_angle(svec1, svec2)-vector_angle(cvec1, cvec2))
    # print('angle_diff:', angle_diff)
    sx, sy, sz = np.linalg.solve([svec1, svec2, svec3], np.eye(3))
    cx, cy, cz = np.dot([sx,sy,sz], [cvec1, cvec2, cvec3])
    return get_cartesian_ijk([0,0,0],cz,cx), angle_diff



def vector_angle(a, b):
    return acos(np.dot(a, b)/(norm(a)*norm(b)))



def input_standard_pos_transform(inp_pos, std_pos, t_vals,
        std_to_inp=True, is_coord = False, DEBUG=False):
    t_vals  = np.array(t_vals).copy()
    std_O   = np.array(std_pos)[-1].copy()
    inp_O   = np.array(inp_pos)[-1].copy()
    std_pos = np.array(std_pos).copy() - std_O
    inp_pos = np.array(inp_pos).copy() - inp_O
    natoms= len(inp_pos)
    if not is_coord:
        inp_O = std_O = np.array([0,0,0])

    R_mat = None
    # return std_pos, inp_pos
    for selection in combinations(range(natoms-1), 3):
        selection = list(selection)
        std_m = std_pos[selection]
        inp_m = inp_pos[selection]
        if np.linalg.det(std_m) > 0.01 and np.linalg.det(inp_m) > 0.01:
            # std_m * R_mat = inp_m
            # R_mat = std_m^-1 * inp_m
            R_mat = np.dot(np.linalg.inv(std_m) , inp_m)
            if DEBUG:
                print('selections:', selection)
                # print(std_m, np.linalg.det(std_m))
                # print(inp_m, np.linalg.det(inp_m))
            break
    if R_mat is None:
        # dimision is less than 3
        for selection in combinations(range(natoms-1), 2):
            std_v0 = std_pos[selection[0]]
            std_v1 = std_pos[selection[1]]
            std_v2 = np.cross(std_v0, std_v1)
            std_m  = np.array([std_v0, std_v1, std_v2])
            inp_v0 = inp_pos[selection[0]]
            inp_v1 = inp_pos[selection[1]]
            inp_v2 = np.cross(inp_v0, inp_v1)
            inp_m  = np.array([inp_v0, inp_v1, inp_v2])
            if np.linalg.det(std_m) > 0.01:
                R_mat = np.dot(np.linalg.inv(std_m) , inp_m)
                if DEBUG:
                    print('selections:', selection)
                break
    if R_mat is None:
        # 2 atoms
        std_v = std_pos[0]
        inp_v = inp_pos[0]
        R = np.cross(std_v, inp_v)
        R = R / np.linalg.norm(R)
        if DEBUG:
            print('stdv, inpv:', std_v, inp_v, '\nR:', R)
        if std_to_inp:
            return np.cross(R, t_vals-std_O)+inp_O
        else:
            return np.cross(t_vals-inp_O, R)+std_O
    else:
        # testification
        if DEBUG:
            assert((np.dot(std_pos, R_mat)-inp_pos < 0.001).all())
            print('test complete')
        if std_to_inp:
            return np.dot(t_vals - std_O, R_mat) + inp_O
        else:
            return np.dot(t_vals - inp_O, np.linalg.inv(R_mat)) + std_O



def compute_distance_matrix(X):
    dists = np.sum(np.square(X), axis = 1) \
        + np.transpose([np.sum(np.square(X), axis = 1)]) -2*np.dot(X, X.T)
    dists = np.sqrt(abs(dists))
    np.fill_diagonal(dists, 0)
    return dists


def dist_change_mat(X, positions):
    X = X.copy()
    positions = positions.copy()
    dists0 = compute_distance_matrix(positions)
    dists1 = compute_distance_matrix(X+positions)
    return dists1 - dists0

def freq_dist_change_mat(XX, positions):
    XX = XX.copy()
    dists0 = compute_distance_matrix(positions)
    return np.array([compute_distance_matrix(x) for x in XX+positions]) - dists0


if __name__ == '__main__':
    np.set_printoptions(precision=3, suppress=True, linewidth=100)


    inp_pos = np.array([
       [ 0.  ,  0.  ,  0.  ],
       [ 0.  ,  0.  , -1.11],
       [ 0.92, -0.61, -0.11],
       [-0.99, -0.49, -0.11],
       [ 0.07,  1.1 , -0.11],
       [ 0.  ,  0.  ,  1.76]])

    std_pos = np.array([
       [ 0.7579,  0.    ,  0.    ],
       [ 1.8679,  0.    ,  0.    ],
       [ 0.8679, -0.8569, -0.6959],
       [ 0.8679,  1.0326, -0.3922],
       [ 0.8679, -0.1758,  1.0881],
       [-1.0021,  0.    ,  0.    ]])

    inp_vec = inp_pos[[0,1]]
    std_vec = std_pos[[0,1]]

    print('inp_vec', inp_vec)
    print('std_vec', std_vec)
    out_inp_vec = input_standard_pos_transform(inp_pos, std_pos, std_vec,
        is_coord=True, DEBUG =True)
    print('out_vec', out_inp_vec)
    print('\n'*5)


    inp_pos = np.array([
       [-1.4951,  0.7264,  0.    ],
       [-1.5688, -0.6965,  0.    ],
       [-2.6259, -0.0039,  0.    ],
       [-0.868 ,  1.8174,  0.    ]])
    std_pos = np.array([
       [-0.0018,  0.475 ,  0.    ],
       [ 0.6433,  1.7454,  0.    ],
       [-0.6183,  1.6716,  0.    ],
       [-0.0018, -0.7834,  0.    ]])
    inp_vec = inp_pos[1]
    std_vec = std_pos[1]

    print('inp_vec', inp_vec)
    print('std_vec', std_vec)
    out_inp_vec = input_standard_pos_transform(inp_pos, std_pos, std_vec,
        is_coord=True, DEBUG =True)
    print('out_vec', out_inp_vec)
    print('\n'*5)


    # CO
    inp_pos = np.array([[ 0.493 ,  0.    ,  0.    ],
       [-0.6573,  0.    ,  0.    ]])

    std_pos = np.array([[ 0.    ,  0.    ,  0.493 ],
       [ 0.    ,  0.    , -0.6573]])


    inp_vec = inp_pos[[0,1]]
    std_vec = std_pos[[0,1]]

    print('inp_vec', inp_vec)
    print('std_vec', std_vec)
    out_inp_vec = input_standard_pos_transform(inp_pos, std_pos, std_vec,
        is_coord=True, DEBUG =True)
    print('out_vec', out_inp_vec)
    print('\n'*5)


def Rotation_matrix(k, theta, radians = False):
    """  使用罗德里格旋转公式 (Rodrigues' rotation formula )
    k is the unit vector of rotation axis;
    v is the rotated vector;
    Rotation Vector:
      R = Ecos(theta) + (1-cos(theta))*k*k^T + sin(theta)[[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]];

    Reference:
    https://baike.baidu.com/item/%E7%BD%97%E5%BE%B7%E9%87%8C%E6%A0%BC%E6%97%8B%E8%BD%AC%E5%85%AC%E5%BC%8F/18878562?fr=aladdin
    """
    k /= norm(k)
    if not radians:
        theta = math.radians(theta)
    kx = k[0]; ky = k[1]; kz = k[2];
    # k_c = k[np.newaxis].T # which now is column vector
    k_outer = np.outer(k, k)
    R = np.identity(3)*math.cos(theta) + (1-math.cos(theta))*k_outer + \
        math.sin(theta)*np.matrix([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    return np.array(R)


