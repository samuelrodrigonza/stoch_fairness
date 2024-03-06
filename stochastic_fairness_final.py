"""
bilevel tDINDP_model Python script

Author: Samuel Rodriguez-Gonzalez, PhD Candidate, The University of Oklahoma
Last modified: 2023-03-30


"""

from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import random
import openpyxl as opxl
import networkx as nx
import pandas as pd
from math import log
import winsound
import time
import ast
import tdINDP_parameters_build

N, A, T, S, L, R, K, N_star_k, N_k, N_prime_k, A_k, A_prime_k, L_k, v, h, p, M_plus, M_minus, alpha, beta, gamma, g, f, q, c, u, b = tdINDP_parameters_build.load_tdINDP_parameters()

# new_M_minus = sum(c[(i,j,k,l)] for k in K for l in L_k[k] for (i,j) in A_k[k])
# print(f'\nnew_M_minus: {new_M_minus}\n')

T = [0]

u = tdINDP_parameters_build.find_capacities()
# for key in M_minus.keys():
# 	M_minus[key] = min(1.5e6,M_minus[key])
# for key in M_plus.keys():
# 	M_plus[key] = min(1.5e6,M_minus[key]) 

N_k_d = {k:[i for i in N_k[k] if b[(i,k,k,0)] < 0] for k in K}
# print(f'N_k_d: {N_k_d}\n')
N_k_s = {k:[i for i in N_k[k] if b[(i,k,k,0)] > 0]  for k in K}
# print(f'N_k_s: {N_k_s}\n')
N_k_nod = {k:[i for i in N_k[k] if b[(i,k,k,0)] == 0] for k in K}
N_k_d_nod = {k:list(set(N_k_d[k]+N_k_nod[k])) for k in K}

# print('{')
# for k in K:
# 	print(f'{k}:')
# 	print('{')
# 	for i in N_k_d[k]:
# 		print(f'{i}:,')
# 	print('},')
# print('}')
# aaaa

SVI = {
'Gas':
{
'Gas_node_3':[0.208,0.4504],
'Gas_node_4':[0.208,0.4504],
'Gas_node_5':[0.4505,0.6909],
'Gas_node_6':[0.208,0.4504],
'Gas_node_7':[0,0.2079],
'Gas_node_8':[0,0.2079],
'Gas_node_9':[0,0.9993],
'Gas_node_10':[0,0.2079],
'Gas_node_11':[0,0.2079],
'Gas_node_12':[0,0.2079],
'Gas_node_13':[0.4505,0.6909],
'Gas_node_14':[0.4505,0.6909],
'Gas_node_15':[0,0.2079]
},
'Power':
{
'Power_node_9':[0.691,0.8594],
'Power_node_10':[0.8595,0.9993],
'Power_node_11':[0,0.2079],
'Power_node_12':[0.691,0.8594],
'Power_node_13':[0,0.2079],
'Power_node_14':[0,0.9993],
'Power_node_15':[0.8595,0.9993],
'Power_node_16':[0.691,0.8594],
'Power_node_17':[0,0.9993],
'Power_node_18':[0.8595,0.9993],
'Power_node_19':[0,0.9993],
'Power_node_20':[0.208,0.4504],
'Power_node_21':[0.8595,0.9993],
'Power_node_22':[0.691,0.8594],
'Power_node_23':[0.8595,0.9993],
'Power_node_24':[0.691,0.8594],
'Power_node_25':[0,0.4504],
'Power_node_26':[0.691,0.8594],
'Power_node_27':[0.8595,0.9993],
'Power_node_28':[0.691,0.9993],
'Power_node_29':[0.691,0.8594],
'Power_node_30':[0.8595,0.9993],
'Power_node_31':[0,0.6909],
'Power_node_32':[0.208,0.6909],
'Power_node_33':[0,0.6909],
'Power_node_34':[0,0.9993],
'Power_node_35':[0.691,0.8594],
'Power_node_36':[0,0.9993],
'Power_node_37':[0,0.2079],
'Power_node_38':[0.8595,0.9993],
'Power_node_39':[0.208,0.4504],
'Power_node_40':[0.8595,0.9993],
'Power_node_41':[0.208,0.4504],
'Power_node_42':[0.208,0.4504],
'Power_node_43':[0.208,0.4504],
'Power_node_44':[0.208,0.4504],
'Power_node_45':[0.8595,0.9993],
'Power_node_46':[0,0.2079],
'Power_node_47':[0,0.2079],
'Power_node_48':[0.208,0.4504],
'Power_node_49':[0,0.4504],
'Power_node_50':[0,0.2079]
},
'Water':
{
'Water_node_15':[0.4505,0.6909],
'Water_node_16':[0,0.2079],
'Water_node_17':[0.8595,0.9993],
'Water_node_18':[0,0.9993],
'Water_node_19':[0.8595,0.9993],
'Water_node_20':[0,0.9993],
'Water_node_21':[0.8595,0.9993],
'Water_node_22':[0.208,0.9993],
'Water_node_23':[0,0.2079],
'Water_node_24':[0.4505,0.6909],
'Water_node_25':[0,0.2079],
'Water_node_26':[0.691,0.8594],
'Water_node_27':[0,0.4504],
'Water_node_28':[0.8595,0.9993],
'Water_node_29':[0.208,0.4504],
'Water_node_30':[0.8595,0.9993],
'Water_node_31':[0.208,0.4504],
'Water_node_32':[0.8595,0.9993],
'Water_node_33':[0.691,0.9993],
'Water_node_34':[0,0.4504],
'Water_node_35':[0,0.2079],
'Water_node_36':[0,0.9993],
'Water_node_37':[0,0.2079],
'Water_node_38':[0.208,0.4504],
'Water_node_39':[0,0.2079],
'Water_node_40':[0,0.2079],
'Water_node_41':[0,0.2079],
'Water_node_42':[0,0.2079],
'Water_node_43':[0.4505,0.9993],
'Water_node_44':[0,0.2079],
'Water_node_45':[0.4505,0.6909],
'Water_node_46':[0,0.2079],
'Water_node_47':[0.4505,0.6909],
'Water_node_48':[0,0.2079]
}
}

# print(print(f'N_k_d_nod: {N_k_d_nod}\n'))
# for c in model.getGenConstrs():
# 	print(f'\nGenconstrName: {c.GenconstrName}')
# 	print(f'FuncPieceError: {c.FuncPieceError}')
# 	print(f'FuncPieceLength: {c.FuncPieceLength}')
# 	print(f'FuncPieceRatio: {c.FuncPieceRatio}')
# 	print(f'FuncPieces: {c.FuncPieces}')
# 	print(f'GenConstrType: {c.GenConstrType}')
# 	if model.status == 3:
# 		print(f'IISGenConstr: {c.IISGenConstr}')
# 		print(f'IISGenConstrForce: {c.IISGenConstrForce}')

# for var in objectives[1].keys(): 
# 	if objectives[1][var] != 0:
# 		print(f'{var}: {objectives[1][var]}')
# 		print()

N_prime_k, A_prime_k = tdINDP_parameters_build.disrput(0)

# print(f'K = {K}')
# print(f'L_k = {L_k}')
# print(f'T = {T}')
print(f'N_k_s = {N_k_s}')
print(f'N_k_nod = {N_k_nod}')
# print(f'N_k_d = {N_k_d}')
# print(f'A_k = {A_k}')

b_temp = {}
for (i,k,l,_) in b.keys():
	b_temp[(i,k,l)] = b[(i,k,l,0)]
b = b_temp
u_temp = {}
for (i,j,k,_) in u.keys():
	u_temp[(i,j,k)] = u[(i,j,k,0)]
u = u_temp

scarce = 0.25

for k in K:
	thy_sum = 0
	for i in N_k_s[k]:
		b[(i,k,k)] = b[(i,k,k)]*scarce

demand = {k:{i:abs(b[(i,k,l)]) for l in L_k[k] for i in N_k_d[k]} for k in K} 

# print(demand)

# print(A_k)
# print(gamma)
# aaaaa
# print(f'\nscarce = {scarce}')

random.seed(42)

same_probs = True
# same_probs = False

# utility_type = 'linear'
utility_type = 'power'

num_pieces = int(1)
pieces_length = 1
num_focus = 3

# rerun_x_start = True

# file_path = f"utiltiy_{utility_type}_scarce_{scarce}_x_start.txt"
# # Open the file and read its content
# if os.path.exists(file_path):
# 	with open(file_path, 'r') as file:
# 		content = file.read()
# 		# Safely evaluate the literal expression using ast.literal_eval
# 		x_start = ast.literal_eval(content)
# 	num_scenarios = 25
# else:
# 	num_scenarios = 1
# 	x_start = {}

num_scenarios = 150
OMEGA = [omega for omega in range(num_scenarios)]


if same_probs:
	#Same probabilities
	rho = {omega:1/len(OMEGA) for omega in OMEGA}
else:
	#Different probabilities
	rho = {omega:random.randint(1, 10) for omega in OMEGA}
	rho_sum = sum(rho[omega] for omega in OMEGA)
	rho = {omega:rho[omega]/rho_sum for omega in OMEGA}

rho_sum = sum(rho[omega] for omega in OMEGA)

if len(OMEGA) > 1:
	if SVI:
		print(f'SVI: {True}')
		alpha = {(i,k,l,omega): random.uniform(SVI[k][i][0], SVI[k][i][1]) for k in K for l in L_k[k] for i in N_k_d[k] for omega in OMEGA}
	else:
		alpha = {(i,k,l,omega): random.uniform(0, 1) for k in K for l in L_k[k] for i in N_k_d[k] for omega in OMEGA}
else:
	alpha = {(i,k,l,omega): 1 for k in K for l in L_k[k] for i in N_k_d[k] for omega in OMEGA}
print()


def EV_and_EEVS(OMEGA = OMEGA, N=N, A=A, T=T, S=S, L=L, R=R, K=K, N_star_k=N_star_k, N_k=N_k, N_prime_k=N_prime_k, A_k=A_k, A_prime_k=A_prime_k, L_k=L_k, v=v, h=h, p=p, M_plus=M_plus, M_minus=M_minus, alpha=alpha, beta=beta, gamma=gamma, g=g, f=f, q=q, c=c, u=u, b=b):
	
	# N_k_d = {k:[i for i in N_k[k] if b[(i,k,k,0)] < 0] for k in K}
	# # print(f'N_k_d: {N_k_d}\n')
	# N_k_s = {k:[i for i in N_k[k] if b[(i,k,k,0)] > 0]  for k in K}
	# # print(f'N_k_s: {N_k_s}\n')
	# N_k_nod = {k:[i for i in N_k[k] if b[(i,k,k,0)] == 0] for k in K}
	# N_k_d_nod = {k:list(set(N_k_d[k]+N_k_nod[k])) for k in K}

	# cake_l_k = {k:{l:sum(b[(i,k,l,0)] for i in N_k[k] if b[(i,k,l,0)] > 0) for l in L_k[k]} for k in K}
	# capacity = {(i,k,l,0):sum(u[(j,i,k,0)] for j in N_k[k] if (j,i) in A_k[k]) for k in K for l in L_k[k] for i in N_k_d[k]}
	# reachability_dict = {k:{l:{} for l in L_k[k]} for k in K}


	'''
	
	Expected value
	
	'''
	
	print('Expected value')
	model = Model('stochastic_tdINDP')
	
	'''
	Gurobi parameters

	'''
	# model.setParam('InfUnbdInfo',1)
	model.setParam('OutputFlag', 0)
	model.setParam('FuncPieceRatio', -1)
	model.setParam('FuncPieces',num_pieces)
	model.setParam('FuncPieceLength',pieces_length)
	model.setParam('DualReductions',0)
	model.setParam('NumericFocus',num_focus)
	# model.setParam(GRB.Param.NonConvex, 2)
	# model.setParam('MIPGap',0.01)
	# model.setParam('OptimalityTol',0.01)
	# model.setParam('IntegralityFocus',1) 
	
	'''
	first stage variables
	'''
	
	x_star = {(i,j,k,l): model.addVar(vtype = GRB.CONTINUOUS, name = 'x_star_'+str((i,j,k,l))) for k in K for l in L_k[k] for (i,j) in A_k[k] if i in N_k_s[k]}
	# w = {(i,k): model.addVar(vtype = GRB.BINARY, name = 'w_'+str((i,k))) for k in K for i in N_k[k]}
	# y = {(i,j,k): model.addVar(vtype = GRB.BINARY, name = 'y_'+str((i,j,k))) for k in K for (i,j) in A_k[k]}
	
	model.update()
	
	'''
	second stage variables
	'''

	proportional_fairness = model.addVar(vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'proportional_fairness')
	phi = {(i,k,l):model.addVar(vtype = GRB.CONTINUOUS, lb=0, ub = abs(b[(i,k,l)]),name="phi_"+str((i,k,l))) for k in K for l in L_k[k] for i in N_k_d[k]}
	log_phi = {(i,k,l):model.addVar(vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=math.log(abs(b[(i,k,l)])), name="log_phi_"+str((i,k,l))) for k in K for l in L_k[k] for i in N_k_d[k]}
	delta_plus = {(i,k,l): model.addVar(vtype = GRB.CONTINUOUS, name = 'delta_plus_'+str((i,k,l))) for k in K for l in L_k[k] for i in N_k_s[k]}
	delta_minus = {(i,k,l): model.addVar(vtype = GRB.CONTINUOUS, name = 'delta_minus_'+str((i,k,l))) for k in K for l in L_k[k] for i in N_k_d[k]}
	x = {(i,j,k,l): model.addVar(vtype = GRB.CONTINUOUS, name = 'x_'+str((i,j,k,l))) for k in K for l in L_k[k] for (i,j) in A_k[k] if i in N_k_d_nod[k]}

	model.update()
	
	'''
	Objective functions:
	
	'''
	if utility_type == 'linear':
		# Linear
		model.addConstr(proportional_fairness == quicksum(log_phi[(i,k,l)] + math.log(0.5) for k in K for l in L_k[k] for i in N_k_d[k]), name = 'proportional_fairness_constraint')
	else:
		# Power
		model.addConstr(proportional_fairness == quicksum(0.5*log_phi[(i,k,l)] for k in K for l in L_k[k] for i in N_k_d[k]), name = 'proportional_fairness_constraint')
	
	model.setObjective(proportional_fairness, GRB.MAXIMIZE)
	model.update()
	
	'''

	Constraints

	'''
	
	
	for k in K:
		for l in L_k[k]:
			
			for i in N_k_s[k]:
				model.addConstr(quicksum(x_star[(i,j,k,l)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == b[(i,k,l)] - delta_plus[(i,k,l)], name = 'supply_flow_balance_'+str((i,k,l)))
			
			
	for k in K:
		for (i,j) in A_k[k]:
			# model.addConstr(y[(i,j,k)] == y[(j,i,k)], name = 'bidir_1_'+str((i,j,k)))
			if i in N_k_s[k]:
				# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*y[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
				# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(i,k)], name = 'node_i_functionality_'+str((i,j,k)))
				# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(j,k)], name = 'node_j_functionality_'+str((i,j,k)))

				model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
				model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_i_functionality_'+str((i,j,k)))
				model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_j_functionality_'+str((i,j,k)))


	for k in K:
		for k_tilde in K:
			for j in N_k[k_tilde]:
				if (i,j,k,k_tilde) in gamma.keys():
					model.addConstr(quicksum(w[(i,k)]*gamma[(i,j,k,k_tilde)] for i in N_k[k]) <= w[(j,k_tilde)], name = 'interdependent_functionality_'+str((k,k_tilde,j)))

	

	for k in K:
		for (i,j) in A_k[k]:
			if i in N_k_d_nod[k]:
				# model.addConstr(quicksum(x[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*y[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
				# model.addConstr(quicksum(x[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(i,k)], name = 'node_i_functionality_'+str((i,j,k)))
				# model.addConstr(quicksum(x[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(j,k)], name = 'node_j_functionality_'+str((i,j,k)))
				model.addConstr(quicksum(x[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
				model.addConstr(quicksum(x[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_i_functionality_'+str((i,j,k)))
				model.addConstr(quicksum(x[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_j_functionality_'+str((i,j,k)))
					
	for k in K:
		for l in L_k[k]:
			for i in N_k_d[k]:
				model.addConstr(quicksum(x[(i,j,k,l)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == b[(i,k,l)] + delta_minus[(i,k,l)], name = 'demand_flow_balance_'+str((i,k,l)))
				model.addConstr(phi[(i,k,l)] == abs(b[(i,k,l)]) - delta_minus[(i,k,l)], name = 'capture_phi_'+str((i,k,l)))
				model.addGenConstrLog(phi[(i,k,l)], log_phi[(i,k,l)], name = 'log_phi_constr_'+str((i,k,l)))
					
			for i in N_k_nod[k]:
				model.addConstr(quicksum(x[(i,j,k,l)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == 0, name = 'transshipment_flow_balance_'+str((i,k,l)))
	
	model.update()
	
	model.optimize()

	x_star_EV = {var.varName:var.x for var in model.getVars()}
	OF_EV = model.objVal

	del model

	''' 

	Expectation of the expected value

	'''

	x_star_EEVS = {}
	OF_omega = []
	for omega in OMEGA:
		print(f'Expectation of the expected value, scenario {omega+1} out of {len(OMEGA)}')

		model = Model('stochastic_tdINDP')
		
		'''
		Gurobi parameters

		'''
		# model.setParam('InfUnbdInfo',1)
		model.setParam('OutputFlag', 0)
		model.setParam('FuncPieceRatio', -1)
		model.setParam('FuncPieces',num_pieces)
		model.setParam('FuncPieceLength',pieces_length)
		model.setParam('DualReductions',0)
		model.setParam('NumericFocus',num_focus)
		# model.setParam(GRB.Param.NonConvex, 2)
		# model.setParam('MIPGap',0.01)
		# model.setParam('OptimalityTol',0.01)
		# model.setParam('IntegralityFocus',1) 
		
		'''
		first stage variables
		'''
		
		x_star = {(i,j,k,l): model.addVar(vtype = GRB.CONTINUOUS, name = 'x_star_'+str((i,j,k,l))) for k in K for l in L_k[k] for (i,j) in A_k[k] if i in N_k_s[k]}
		# w = {(i,k): model.addVar(vtype = GRB.BINARY, name = 'w_'+str((i,k))) for k in K for i in N_k[k]}
		# y = {(i,j,k): model.addVar(vtype = GRB.BINARY, name = 'y_'+str((i,j,k))) for k in K for (i,j) in A_k[k]}
		
		model.update()
		
		'''
		second stage variables
		'''

		proportional_fairness = model.addVar(vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'proportional_fairness_'+str(omega))
		phi = {(i,k,l,omega):model.addVar(vtype = GRB.CONTINUOUS, lb=0, ub = abs(b[(i,k,l)]),name="phi_"+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k]}
		log_phi = {(i,k,l,omega):model.addVar(vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=math.log(abs(b[(i,k,l)])), name="log_phi_"+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k]}
		delta_plus = {(i,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'delta_plus_'+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_s[k]}
		delta_minus = {(i,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'delta_minus_'+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k]}
		x = {(i,j,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'x_'+str((i,j,k,l,omega))) for k in K for l in L_k[k] for (i,j) in A_k[k] if i in N_k_d_nod[k]}

		model.update()
		
		'''
		Objective functions:
		
		'''
		if utility_type == 'linear':
			# Linear
			model.addConstr(proportional_fairness == quicksum(log_phi[(i,k,l,omega)] + math.log(alpha[(i,k,l,omega)]) for k in K for l in L_k[k] for i in N_k_d[k]), name = 'proportional_fairness_constraint')
		else:
			# Power
			model.addConstr(proportional_fairness == quicksum(alpha[(i,k,l,omega)]*log_phi[(i,k,l,omega)] for k in K for l in L_k[k] for i in N_k_d[k]), name = 'proportional_fairness_constraint')
		
		model.setObjective(proportional_fairness, GRB.MAXIMIZE)
		model.update()
		
		'''

		Constraints

		'''
		
		for k in K:
			for l in L_k[k]:
				
				for i in N_k_s[k]:
					model.addConstr(quicksum(x_star[(i,j,k,l)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == b[(i,k,l)] - delta_plus[(i,k,l,omega)], name = 'supply_flow_balance_'+str((i,k,l,omega)))
				
				
		for k in K:
			for (i,j) in A_k[k]:
				# model.addConstr(y[(i,j,k)] == y[(j,i,k)], name = 'bidir_1_'+str((i,j,k)))
				if i in N_k_s[k]:
					# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*y[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
					# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(i,k)], name = 'node_i_functionality_'+str((i,j,k)))
					# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(j,k)], name = 'node_j_functionality_'+str((i,j,k)))

					model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
					model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_i_functionality_'+str((i,j,k)))
					model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_j_functionality_'+str((i,j,k)))


		for k in K:
			for k_tilde in K:
				for j in N_k[k_tilde]:
					if (i,j,k,k_tilde) in gamma.keys():
						model.addConstr(quicksum(w[(i,k)]*gamma[(i,j,k,k_tilde)] for i in N_k[k]) <= w[(j,k_tilde)], name = 'interdependent_functionality_'+str((k,k_tilde,j)))

		

		for k in K:
			for (i,j) in A_k[k]:
				if i in N_k_d_nod[k]:
					# model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*y[(i,j,k)], name = 'arc_functionality_'+str((i,j,k,omega)))
					# model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*w[(i,k)], name = 'node_i_functionality_'+str((i,j,k,omega)))
					# model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*w[(j,k)], name = 'node_j_functionality_'+str((i,j,k,omega)))
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_i_functionality_'+str((i,j,k)))
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_j_functionality_'+str((i,j,k)))
						
		for k in K:
			for l in L_k[k]:
				for i in N_k_d[k]:
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == b[(i,k,l)] + delta_minus[(i,k,l,omega)], name = 'demand_flow_balance_'+str((i,k,l,omega)))
					model.addConstr(phi[(i,k,l,omega)] == abs(b[(i,k,l)]) - delta_minus[(i,k,l,omega)], name = 'capture_phi_'+str((i,k,l,omega)))
					model.addGenConstrLog(phi[(i,k,l,omega)], log_phi[(i,k,l,omega)], name = 'log_phi_constr_'+str((i,k,l,omega)))
						
				for i in N_k_nod[k]:
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == 0, name = 'transshipment_flow_balance_'+str((i,k,l,omega)))
		
		for k in K:
			for (i,j) in A_k[k]:
				for l in L_k[k]:
					if i in N_k_s[k]:
						model.addConstr(model.getVarByName('x_star_'+str((i,j,k,l))) == x_star_EV['x_star_'+str((i,j,k,l))], name = 'fix_first_stage_x_star_'+str((i,j,k,l)))
							
		model.update()

		model.optimize()

		x_star_EEVS[omega] = {var.varName:var.x for var in model.getVars()}
		OF_omega.append(model.objVal) 

		del model
	
	OF_EEVS = sum(OF for OF in OF_omega)/len(OMEGA)

	return x_star_EV, OF_EV, x_star_EEVS, OF_EEVS

def WS(x_start = {}, OMEGA = OMEGA, N=N, A=A, T=T, S=S, L=L, R=R, K=K, N_star_k=N_star_k, N_k=N_k, N_prime_k=N_prime_k, A_k=A_k, A_prime_k=A_prime_k, L_k=L_k, v=v, h=h, p=p, M_plus=M_plus, M_minus=M_minus, alpha=alpha, beta=beta, gamma=gamma, g=g, f=f, q=q, c=c, u=u, b=b):

	x_star_WS = {}
	OF_omega = []
	for omega in OMEGA:
		print(f'Wait and see, scenario {omega+1} out of {len(OMEGA)}')

		model = Model('stochastic_tdINDP')
		
		'''
		Gurobi parameters

		'''
		# model.setParam('InfUnbdInfo',1)
		model.setParam('OutputFlag', 0)
		model.setParam('FuncPieceRatio', -1)
		model.setParam('FuncPieces',num_pieces)
		model.setParam('FuncPieceLength',pieces_length)
		model.setParam('DualReductions',0)
		model.setParam('NumericFocus',num_focus)
		# model.setParam(GRB.Param.NonConvex, 2)
		# model.setParam('MIPGap',0.01)
		# model.setParam('OptimalityTol',0.01)
		# model.setParam('IntegralityFocus',1) 
		
		'''
		first stage variables
		'''
		
		x_star = {(i,j,k,l): model.addVar(vtype = GRB.CONTINUOUS, name = 'x_star_'+str((i,j,k,l))) for k in K for l in L_k[k] for (i,j) in A_k[k] if i in N_k_s[k]}
		# w = {(i,k): model.addVar(vtype = GRB.BINARY, name = 'w_'+str((i,k))) for k in K for i in N_k[k]}
		# y = {(i,j,k): model.addVar(vtype = GRB.BINARY, name = 'y_'+str((i,j,k))) for k in K for (i,j) in A_k[k]}
		
		model.update()
		
		'''
		second stage variables
		'''

		proportional_fairness = model.addVar(vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'proportional_fairness_'+str(omega))
		phi = {(i,k,l,omega):model.addVar(vtype = GRB.CONTINUOUS, lb=0, ub = abs(b[(i,k,l)]),name="phi_"+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k]}
		log_phi = {(i,k,l,omega):model.addVar(vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=math.log(abs(b[(i,k,l)])), name="log_phi_"+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k]}
		delta_plus = {(i,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'delta_plus_'+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_s[k]}
		delta_minus = {(i,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'delta_minus_'+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k]}
		x = {(i,j,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'x_'+str((i,j,k,l,omega))) for k in K for l in L_k[k] for (i,j) in A_k[k] if i in N_k_d_nod[k]}

		model.update()
		
		'''
		Objective functions:
		
		'''
		if utility_type == 'linear':
			# Linear
			model.addConstr(proportional_fairness == quicksum(log_phi[(i,k,l,omega)] + math.log(alpha[(i,k,l,omega)]) for k in K for l in L_k[k] for i in N_k_d[k]), name = 'proportional_fairness_constraint')
		else:
			# Power
			model.addConstr(proportional_fairness == quicksum(alpha[(i,k,l,omega)]*log_phi[(i,k,l,omega)] for k in K for l in L_k[k] for i in N_k_d[k]), name = 'proportional_fairness_constraint')
		
		model.setObjective(proportional_fairness, GRB.MAXIMIZE)
		model.update()
		
		'''

		Constraints

		'''
		
		for k in K:
			for l in L_k[k]:
				for i in N_k_s[k]:
					model.addConstr(quicksum(x_star[(i,j,k,l)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == b[(i,k,l)] - delta_plus[(i,k,l,omega)], name = 'supply_flow_balance_'+str((i,k,l,omega)))
				
				
		for k in K:
			for (i,j) in A_k[k]:
				# model.addConstr(y[(i,j,k)] == y[(j,i,k)], name = 'bidir_1_'+str((i,j,k)))
				if i in N_k_s[k]:
					# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*y[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
					# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(i,k)], name = 'node_i_functionality_'+str((i,j,k)))
					# model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(j,k)], name = 'node_j_functionality_'+str((i,j,k)))

					model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
					model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_i_functionality_'+str((i,j,k)))
					model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_j_functionality_'+str((i,j,k)))


		for k in K:
			for k_tilde in K:
				for j in N_k[k_tilde]:
					if (i,j,k,k_tilde) in gamma.keys():
						model.addConstr(quicksum(w[(i,k)]*gamma[(i,j,k,k_tilde)] for i in N_k[k]) <= w[(j,k_tilde)], name = 'interdependent_functionality_'+str((k,k_tilde,j)))

		

		for k in K:
			for (i,j) in A_k[k]:
				if i in N_k_d_nod[k]:
					# model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*y[(i,j,k)], name = 'arc_functionality_'+str((i,j,k,omega)))
					# model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*w[(i,k)], name = 'node_i_functionality_'+str((i,j,k,omega)))
					# model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*w[(j,k)], name = 'node_j_functionality_'+str((i,j,k,omega)))
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_i_functionality_'+str((i,j,k)))
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)], name = 'node_j_functionality_'+str((i,j,k)))
						
		for k in K:
			for l in L_k[k]:
				for i in N_k_d[k]:
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == b[(i,k,l)] + delta_minus[(i,k,l,omega)], name = 'demand_flow_balance_'+str((i,k,l,omega)))
					model.addConstr(phi[(i,k,l,omega)] == abs(b[(i,k,l)]) - delta_minus[(i,k,l,omega)], name = 'capture_phi_'+str((i,k,l,omega)))
					model.addGenConstrLog(phi[(i,k,l,omega)], log_phi[(i,k,l,omega)], name = 'log_phi_constr_'+str((i,k,l,omega)))
						
				for i in N_k_nod[k]:
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == 0, name = 'transshipment_flow_balance_'+str((i,k,l,omega)))
		
		if x_start:
			for k in K:
				for l in L_k[k]:
					for i in N_k_d[k]:
						model.getVarByName('phi_'+str((i,k,l,omega))).Start = x_star_EV['phi_'+str((i,k,l))]
						model.update()
							
		model.update()

		model.optimize()

		x_star_WS[omega] = {var.varName:var.x for var in model.getVars()}
		OF_omega.append(model.objVal) 

		del model
	
	OF_WS = sum(OF for OF in OF_omega)/len(OMEGA)

	return x_star_WS, OF_WS

def DE(x_start = {}, OMEGA = OMEGA, N=N, A=A, T=T, S=S, L=L, R=R, K=K, N_star_k=N_star_k, N_k=N_k, N_prime_k=N_prime_k, A_k=A_k, A_prime_k=A_prime_k, L_k=L_k, v=v, h=h, p=p, M_plus=M_plus, M_minus=M_minus, alpha=alpha, beta=beta, gamma=gamma, g=g, f=f, q=q, c=c, u=u, b=b):
	
	# N_k_d = {k:[i for i in N_k[k] if b[(i,k,k,0)] < 0] for k in K}
	# # print(f'N_k_d: {N_k_d}\n')
	# N_k_s = {k:[i for i in N_k[k] if b[(i,k,k,0)] > 0]  for k in K}
	# # print(f'N_k_s: {N_k_s}\n')
	# N_k_nod = {k:[i for i in N_k[k] if b[(i,k,k,0)] == 0] for k in K}
	# N_k_d_nod = {k:list(set(N_k_d[k]+N_k_nod[k])) for k in K}

	# cake_l_k = {k:{l:sum(b[(i,k,l,0)] for i in N_k[k] if b[(i,k,l,0)] > 0) for l in L_k[k]} for k in K}
	# capacity = {(i,k,l,0):sum(u[(j,i,k,0)] for j in N_k[k] if (j,i) in A_k[k]) for k in K for l in L_k[k] for i in N_k_d[k]}
	# reachability_dict = {k:{l:{} for l in L_k[k]} for k in K}

	print(f'Deterministic equivalent')

	model = Model('stochastic_tdINDP')
	
	'''
	Gurobi parameters

	'''
	# model.setParam('InfUnbdInfo',1)
	# model.setParam('OutputFlag', 0)
	model.setParam('FuncPieceRatio', -1)
	model.setParam('FuncPieces',num_pieces)
	model.setParam('FuncPieceLength',pieces_length)
	model.setParam('DualReductions',0)
	model.setParam('NumericFocus',num_focus)
	# model.setParam(GRB.Param.NonConvex, 2)
	# model.setParam('MIPGap',0.01)
	# model.setParam('OptimalityTol',0.01)
	# model.setParam('IntegralityFocus',1) 
	
	'''
	first stage variables
	'''
	
	x_star = {(i,j,k,l): model.addVar(vtype = GRB.CONTINUOUS, name = 'x_star_'+str((i,j,k,l))) for k in K for l in L_k[k] for (i,j) in A_k[k] if i in N_k_s[k]}
	w = {(i,k): model.addVar(vtype = GRB.BINARY, name = 'w_'+str((i,k))) for k in K for i in N_k[k]}
	y = {(i,j,k): model.addVar(vtype = GRB.BINARY, name = 'y_'+str((i,j,k))) for k in K for (i,j) in A_k[k]}
	
	model.update()
	
	'''
	second stage variables
	'''

	proportional_fairness = {omega:model.addVar(vtype = GRB.CONTINUOUS, lb = -GRB.INFINITY, name = 'proportional_fairness_'+str(omega)) for omega in OMEGA}
	phi = {(i,k,l,omega):model.addVar(vtype = GRB.CONTINUOUS, lb=0, ub = abs(b[(i,k,l)]),name="phi_"+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k] for omega in OMEGA}
	log_phi = {(i,k,l,omega):model.addVar(vtype=GRB.CONTINUOUS,lb=-GRB.INFINITY,ub=math.log(abs(b[(i,k,l)])), name="log_phi_"+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k] for omega in OMEGA}
	delta_plus = {(i,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'delta_plus_'+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_s[k] for omega in OMEGA}
	delta_minus = {(i,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'delta_minus_'+str((i,k,l,omega))) for k in K for l in L_k[k] for i in N_k_d[k] for omega in OMEGA}
	x = {(i,j,k,l,omega): model.addVar(vtype = GRB.CONTINUOUS, name = 'x_'+str((i,j,k,l,omega))) for k in K for l in L_k[k] for (i,j) in A_k[k] if i in N_k_d_nod[k] for omega in OMEGA}

	model.update()
	
	'''
	Objective functions:
	
	'''
	if utility_type == 'linear':
		# Linear
		for omega in OMEGA:
			model.addConstr(proportional_fairness[omega] == quicksum(log_phi[(i,k,l,omega)] + math.log(alpha[(i,k,l,omega)]) for k in K for l in L_k[k] for i in N_k_d[k]), name = 'proportional_fairness_constraint')
	else:
		# Power
		for omega in OMEGA:
			model.addConstr(proportional_fairness[omega] == quicksum(alpha[(i,k,l,omega)]*log_phi[(i,k,l,omega)] for k in K for l in L_k[k] for i in N_k_d[k]), name = 'proportional_fairness_constraint')
	
	model.setObjective(quicksum(rho[omega]*proportional_fairness[omega] for omega in OMEGA), GRB.MAXIMIZE)
	model.update()
	
	'''

	Constraints

	'''
	
	for k in K:
		for l in L_k[k]:
			
			for i in N_k_s[k]:
				for omega in OMEGA:
					model.addConstr(quicksum(x_star[(i,j,k,l)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == b[(i,k,l)] - delta_plus[(i,k,l,omega)], name = 'supply_flow_balance_'+str((i,k,l,omega)))
			
			
	for k in K:
		for (i,j) in A_k[k]:
			# model.addConstr(y[(i,j,k)] == y[(j,i,k)], name = 'bidir_1_'+str((i,j,k)))
			if i in N_k_s[k]:
				model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*y[(i,j,k)], name = 'arc_functionality_'+str((i,j,k)))
				model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(i,k)], name = 'node_i_functionality_'+str((i,j,k)))
				model.addConstr(quicksum(x_star[(i,j,k,l)] for l in L_k[k]) <= u[(i,j,k)]*w[(j,k)], name = 'node_j_functionality_'+str((i,j,k)))


	for k in K:
		for k_tilde in K:
			for j in N_k[k_tilde]:
				if (i,j,k,k_tilde) in gamma.keys():
					model.addConstr(quicksum(w[(i,k)]*gamma[(i,j,k,k_tilde)] for i in N_k[k]) <= w[(j,k_tilde)], name = 'interdependent_functionality_'+str((k,k_tilde,j)))

	

	for k in K:
		for (i,j) in A_k[k]:
			for omega in OMEGA:
				if i in N_k_d_nod[k]:
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*y[(i,j,k)], name = 'arc_functionality_'+str((i,j,k,omega)))
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*w[(i,k)], name = 'node_i_functionality_'+str((i,j,k,omega)))
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for l in L_k[k]) <= u[(i,j,k)]*w[(j,k)], name = 'node_j_functionality_'+str((i,j,k,omega)))
					
	for k in K:
		for l in L_k[k]:
			for i in N_k_d[k]:
				for omega in OMEGA:
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == b[(i,k,l)] + delta_minus[(i,k,l,omega)], name = 'demand_flow_balance_'+str((i,k,l,omega)))
					model.addConstr(phi[(i,k,l,omega)] == abs(b[(i,k,l)]) - delta_minus[(i,k,l,omega)], name = 'capture_phi_'+str((i,k,l,omega)))
					model.addGenConstrLog(phi[(i,k,l,omega)], log_phi[(i,k,l,omega)], name = 'log_phi_constr_'+str((i,k,l,omega)))
					
			for i in N_k_nod[k]:
				for omega in OMEGA:
					model.addConstr(quicksum(x[(i,j,k,l,omega)] for j in N_k[k] if (i,j) in A_k[k])-(quicksum(x[(j,i,k,l,omega)] for j in N_k_d_nod[k] if (j,i) in A_k[k]) + quicksum(x_star[(j,i,k,l)] for j in N_k_s[k] if (j,i) in A_k[k])) == 0, name = 'transshipment_flow_balance_'+str((i,k,l,omega)))
	
	if x_start:
		# model.addConstr(proportional_fairness <= x_start['proportional_fairness'])
		for k in K:
			for l in L_k[k]:
				for i in N_k_d[k]:
					for omega in OMEGA:
						model.getVarByName('phi_'+str((i,k,l,omega))).Start = x_star_EV['phi_'+str((i,k,l))]
						model.update()
	model.update()

	model.params.LogFile=f'LOG_FILE_scarce_{scarce}_utility_type_{utility_type}_num_scenarios_{len(OMEGA)}_equal_probabilities_{same_probs}_nat_log_linear_pieces_length_{pieces_length}'
	model.optimize()

	print(f'\nscarce: {scarce*100}, Utility type: {utility_type}, Num. scenarios: {len(OMEGA)}, Equal probabilities?: {same_probs}, nat. log. linear pieces length: {pieces_length}')
	# print()
	# print(f'Scenario probabilities:\n{rho}')
	# print()
	# print(f'sum of probabilities: {rho_sum}')

	proportional_fairness_values = [model.getVarByName('proportional_fairness_'+str(omega)).x for omega in OMEGA]
		
	# Calculate mean and variance
	mean_value = np.mean(proportional_fairness_values)
	variance_value = np.var(proportional_fairness_values)
	std_value = np.std(proportional_fairness_values)

	# print()
	# print(f'mean proportional_fairness value: {mean_value}, proportional_fairness variance: {variance_value}, proportional_fairness std: {std_value}\n')
	# with open("monitoring_amount_"+str(PI[0])+"_horizon_size_"+str(len_T)+".txt", "w") as file:
	# 	file.write(f'\n-----\nHorizon size: {len_T}\n-----\n')

	for k in K:
		# print(f'\n---------------Network: {k}-----------------')
		for l in L_k[k]:
			for i in N_k_d[k]:
				# print()
				# print(f'node: {i}\n')
				allocation = []
				for omega in OMEGA:
					# print(f'allocation value in scenario {omega}: {model.getVarByName("phi_"+str((i,k,l,omega))).x}')
					allocation.append(model.getVarByName("phi_"+str((i,k,l,omega))).x)
				# Calculate mean and variance
				mean_value = np.mean(allocation)
				variance_value = np.var(allocation)
				std_value = np.std(allocation)
				# print()
				# print(f'mean allocation value: {mean_value}, allocation variance: {variance_value}, allocation std: {std_value}')


	results = {}

	results['proportional_fairness'] = [model.getVarByName('proportional_fairness_'+str(omega)).x for omega in OMEGA]

	for k in K:
		for l in L_k[k]:
			for i in N_k_d[k]:
				results['alpha_'+str((i,k,l))] = [alpha[(i,k,l,omega)] for omega in OMEGA]
				results['phi_'+str((i,k,l))] = [model.getVarByName('phi_'+str((i,k,l,omega))).x for omega in OMEGA]
				results['utility_of_phi_'+str((i,k,l))] = [max((model.getVarByName('phi_'+str((i,k,l,omega))).x**alpha[(i,k,l,omega)]).real,0) for omega in OMEGA]
				results['utility_of_demand_'+str((i,k,l))] = [abs(b[(i,k,l)])**alpha[(i,k,l,omega)] for omega in OMEGA]

	for k in K:
		for l in L_k[k]:
			for i in N_k_d[k]:
				results['log_phi_'+str((i,k,l))] = [model.getVarByName('log_phi_'+str((i,k,l,omega))).x for omega in OMEGA]

	for k in K:
		for l in L_k[k]:
			for i in N_k_s[k]:
				results['delta_plus_'+str((i,k,l))] = [model.getVarByName('delta_plus_'+str((i,k,l,omega))).x for omega in OMEGA]

	for k in K:
		for l in L_k[k]:
			for i in N_k_d[k]:
				results['delta_minus_'+str((i,k,l))] = [model.getVarByName('delta_minus_'+str((i,k,l,omega))).x for omega in OMEGA]

	for k in K:
		for l in L_k[k]:
			for i in A_k[k]:
				if i in N_k_d_nod[k]:
					results['x_'+str((i,j,k,l))] = [model.getVarByName('x_'+str((i,j,k,l,omega))).x for omega in OMEGA]			

	for k in K:
		for l in L_k[k]:
			for i in A_k[k]:
				if i in N_k_s[k]:
					results['x_star_'+str((i,j,k,l))] = [model.getVarByName('x_star_'+str((i,j,k,l))).x for omega in OMEGA]

	# for k in K:
	# 	for l in L_k[k]:
	# 		for i in N_k[k]:
	# 			results['w_'+str((i,k))] = [model.getVarByName('w_'+str((i,k))).x for omega in OMEGA]

	# for k in K:
	# 	for (i,j) in A_k[k]:
	# 		results['y_'+str((i,j,k))] = [model.getVarByName('y_'+str((i,j,k))).x for omega in OMEGA]

	results = pd.DataFrame(results, columns = results.keys())

	results = results.transpose()

	def compute_statistics(row):
	    return pd.Series({
	        'mean': np.mean(row),
	        'std': np.std(row),
	        'min': np.min(row),
	        'max': np.max(row),
	        'median': np.median(row),
	        'mode': row.mode().iloc[0],
	        'q1': np.percentile(row, 25),
	        'q2': np.percentile(row, 50),
	        'q3': np.percentile(row, 75)
	    })

	# Apply the function to each row
	summary_stats = results.apply(compute_statistics, axis=1)

	# Concatenate the original DataFrame with the new statistics
	results = pd.concat([results, summary_stats], axis=1)

	results.to_csv(f'scarce_{scarce}_utility_type_{utility_type}_num_scenarios_{len(OMEGA)}_equal_probabilities_{same_probs}_nat_log_linear_pieces_length_{pieces_length}.csv')

	if len(OMEGA) == 1:
		x_start = {(i,k,l):model.getVarByName('phi_'+str((i,k,l,omega))).x for k in K for l in L_k[k] for i in N_k_d[k]}
		with open(f"utiltiy_{utility_type}_scarce_{scarce}_x_start.txt", "w") as file:
			file.write(f'{x_start}')

	x_star_DE = {var.varName:var.x for var in model.getVars()}
	OF_DE = model.objVal
	
	return x_star_DE, OF_DE



x_star_EV, OF_EV, x_star_EEVS, OF_EEVS = EV_and_EEVS()
x_star_WS, OF_WS = WS()
x_star_DE, OF_DE = DE(x_start = {})

# if scarce != 0.75:
# 	x_star_DE, OF_DE = DE(x_start = x_star_EV)
# else:
# 	x_star_DE, OF_DE = DE(x_start = {})

print()
print(f'OF_EEVS: {OF_EEVS}')
print()
print(f'OF_WS: {OF_WS}')
print()
print(f'OF_DE: {OF_DE}')
print()
print(f'Expectation of the Expected Value Solution (EVVS): {OF_EEVS}')
print()
print(f'Expected Value of Perfect Information (EVPI): {OF_WS-OF_DE}')
print()
print(f'Value of the Stochastic Solution (VSS): {OF_DE-OF_EEVS}')
print()

with open(f"scarce_{scarce}_utility_type_{utility_type}_num_scenarios_{len(OMEGA)}_equal_probabilities_{same_probs}_nat_log_linear_pieces_length_{pieces_length}.txt", "w") as file:
	file.write(f"scarce_{scarce}_utility_type_{utility_type}_num_scenarios_{len(OMEGA)}_equal_probabilities_{same_probs}_nat_log_linear_pieces_length_{pieces_length}")
	file.write(f'\nOF_EEVS: {OF_EEVS}')
	file.write('\n')
	file.write(f'OF_WS: {OF_WS}')
	file.write('\n')
	file.write(f'OF_DE: {OF_DE}')
	file.write('\n')
	file.write(f'Expectation of the Expected Value Solution (EVVS): {OF_EEVS}')
	file.write('\n')
	file.write(f'Expected Value of Perfect Information (EVPI): {OF_WS-OF_DE}')
	file.write('\n')
	file.write(f'Value of the Stochastic Solution (VSS): {OF_DE-OF_EEVS}')
	
duration = 500# milliseconds
freq = 1000  # Hz
winsound.Beep(freq, duration)
sys.exit()
