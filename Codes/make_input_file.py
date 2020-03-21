import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--pbin', nargs='*', type=float, default=[])
parser.add_argument('--rbin', nargs='*', type=float, default=[])
args = parser.parse_args()

bins = len(args.rbin) - 1

f = open('input_ExoOccurrence.txt','w+')
f.write('path_to_obs = FGK_planets.dat # \n')

b = 0
f.write('param_to_fit =')
while b < bins:
	b += 1
	f.write(' r' + str(b))
f.write(' # \n')
b = 0
f.write('param_to_sim =')
while b < bins:
	b += 1
	f.write(' r' + str(b))
f.write(' # \n')
b = 0
f.write('prior_func =')
while b < bins:
	b += 1
	f.write(' flat_prior')
f.write(' # \n')

b = 0
while b < bins:
	b += 1
	f.write('r' + str(b) + '_prior_par_name = pmin pmax # \n')
	pmin, pmax = args.pbin[0], args.pbin[1]
	rmin, rmax = args.rbin[b-1], args.rbin[b]
	fmax = 2*np.log2(pmax/pmin)*np.log2(rmax/rmin)
	f.write('r' + str(b) + '_prior_par_val = 0.0 ' + str(fmax) + ' # \n')
	f.write('r' + str(b) + '_lim = 0.0 ' + str(fmax) + ' # \n')
	f.write('r' + str(b) + ' = 0.05 # \n')

f.write('M = 200 # \n')
f.write('Mini = 500 # \n')
f.write('delta = 0.1 # \n')
f.write('qthreshold = 0.75 # \n')
f.write('file_root = exo_sims_ # \n')
f.write('screen = 0 # \n')
f.write('ncores = 1 # \n')
f.write('split_output = 1 # \n')
f.write('simulation_func = simulate_catalogue # \n')
f.write('distance_func = distance # \n')
f.write('dist_dim = 1 # \n')
f.close()
