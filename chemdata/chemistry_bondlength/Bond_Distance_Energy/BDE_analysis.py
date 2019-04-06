import os 
import pandas as pd 
from ase import Atoms 
from scipy import interpolate 
from scipy.optimize import fsolve 
import logging 
import json 
import pdb 

logging.basicConfig(level=logging.ERROR) 
logging.basicConfig(filename = 'error.log') 
console = logging.StreamHandler() 
console.setLevel(logging.INFO) 


if __name__ != '__main__': 
	exit() 


max_E = 3.0 
bond_dist_dict = {"max_E": max_E} 
if not os.path.exists('csv'): 
	os.makedirs('csv') 

for a0 in range(1, 100): 
	for a1 in range(1, 100): 
		a = Atoms([a0, a1]) 
		s0, s1 = a[0].symbol, a[1].symbol 
		dirname = '%02d-%02d' %(a0, a1) 
		if os.path.exists(dirname): 
			logging.info('Info: '+dirname+' '+a.get_chemical_formula()) 
			csv = dirname +'/'+ a.get_chemical_formula() +'_1'  
			csv_new = 'csv/'+a.get_chemical_formula() 
			if not os.path.exists(csv): 
				logging.info('dirname: '+dirname+' cvs :'+csv+' not exists') 
				continue 
			df = pd.read_csv(csv) 
			logging.info(df) 
			if len(df) == 0: 
				logging.error(dirname+' dataframe is empty') 
				# os.system('rm -f '+dirname+'/Done') 
				continue 
			try: 
				df['energy'] -= [min(df['energy'])]*len(df['energy']) 
			except: 
				logging.error('dirname: '+dirname+' '+a.get_chemical_formula()+' ERROR') 
				continue 
			try:
				df_eq_dist = df[df['energy'] == 0]['dist']
				df_eq_dist = df['dist'][df_eq_dist.index[0]]
				eq_dist = float(df_eq_dist)  
			except: 
				logging.error(eq_dist) 
				exit() 
				
			df.to_csv(csv_new, float_format='%.2f', index=False) 
			# qudratic spline  
			qud_spl = interpolate.splrep(df['dist'], df['energy'], k=2) 
			def f(x): 
				return float(interpolate.splev(x, qud_spl)[0]) - max_E 
			solution0 = fsolve(f, eq_dist-0.5)[0] 
			solution1 = fsolve(f, eq_dist+0.5)[0] 
			if (solution0 >= solution1): 
				logging.critical(dirname+' min>max') 
			df_valid = df[solution0 <= df['dist']][df['dist'] <= solution1] 
			# logging.info('\n'+str(df_valid)) 
			logging.info('valid df set:'+str(len(df_valid))) 
			logging.info('min: %.2f; max: %.2f' %(solution0, solution1)) 
			logging.info('') 
			if not str(s0) in bond_dist_dict: 
				bond_dist_dict[str(s0)] = {} 
			bond_dist_dict[str(s0)][str(s1)] = {
					'dist': round(eq_dist,2), 
					'min' : round(solution0,2), 
					'max' : round(solution1,2), 
					'max_over_precent': round((solution1-eq_dist)/eq_dist*100,0)} 
json.dump(bond_dist_dict, open('bond_dict.'+str(max_E)+'.json', 'w'), indent=4) 


