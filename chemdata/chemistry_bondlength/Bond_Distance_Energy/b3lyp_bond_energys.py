import os 
import sys

# sys.path.append('/home/sky/christina')
import numpy as np 
import mse
import gaussian
import chem_db 
import configparser
import itertools 
import pandas as pd 


import pdb 


def test_atom_dist(a0, a1, order=1): 
	# pdb.set_trace() 
	print('run', (a0, a1, order)) 
	config = configparser.ConfigParser()
	config.read('../config') 
	atoms = mse.Mol([a0,a1]) 
	atoms.filename = atoms.formula+'.com' 
	#atoms.config = config 
	bond = chem_db.get_bond(a0,a1)  
	bond.order = order 
	bond_dict = {'bonds':{(0,1):bond}, 'half_bonds':{}, 
		'bond_order_matrix':np.array([[0,1],[1,0]]), 
		'bond_matrix':np.array([[0,1],[1,0]])} 
	atoms._bond_dict = bond_dict 
	dist  = 0.3 
	calc  = gaussian.Gaussian() 
	atoms.calc = calc 
	calc.method = 'b3lyp/genecp' 
	calc.method = 'b3lyp/6-311G(d,p)' 
	calc.method_i = -1 
	calc.run_type = 'sp' 
	calc.run_in_situ = True 
	calc.freq_analysis = False 
	calc.cpu_mem = (8, '10GB') 
	dist_list = [] 
	energy_list = [] 
	# print(atoms.bond_dict) 
	while dist < 5: 
		atoms[1].position = [dist, 0, 0] 
		atoms.calc.filename = atoms.formula+'_'+'%2.1f' %(dist)+'.com'  
		atoms.start_calc() 
		try: 
			energy = atoms.get_potential_energy() 
			dist_list.append(dist) 
			energy_list.append(energy) 
		except: 
			pass 
		if dist < 2: 
			dist += 0.1
		else: 
			dist += 0.2 
	dist_dict = {'dist':dist_list, 'energy': energy_list} 
	df = pd.DataFrame(dist_dict) 
	df.to_csv(atoms.formula+'_'+str(order), index=False) 
	# return result 

def all_combinations(run_list = []): 
	for a0, a1 in itertools.combinations_with_replacement(range(1, 21), 2): 
		if run_list != [] and not (a0, a1) in run_list: 
			continue 
		dirname = '%02d-%02d' %(a0, a1)
		if not os.path.exists(dirname): 
			os.makedirs(dirname) 
		if not os.path.exists(dirname+'/Done'): 
			os.chdir(dirname) 
			# try: 
			test_atom_dist(a0, a1) 
			#except: 
			#	print((a0, a1), 'error') 
			os.system('touch Done') 
			os.chdir('..') 


if __name__ == '__main__': 
	all_combinations() 





