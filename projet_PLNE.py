#Projet MOGPL
#Binh Thanh LUONG - 3504859
#Melissa YAYA - 3673113

import numpy as np
from gurobipy import *
import time
import matplotlib.pyplot as plt
import sys
import os.path
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def read_file(filename):
	"""
	Lit le fichier passe en parametre et retourne les sequences correspondantes
	"""
	file = open (filename, "r")
	lignes = []
	colonnes = []
	l = True
	for line in file.read().split("\n"):
		if l :
			if not(line):
				lignes.append([])
				continue
			if "#" in line:
				l = False
				continue
			lignes.append(list(map(int, line.split(" "))))
		else:
			if not(line):
				colonnes.append([])
				continue
			colonnes.append(list(map(int, line.split(" "))))
	file.close()
	del colonnes[len(colonnes)-1]
	return np.array(lignes),np.array(colonnes)

def compute(seqs_ligne, seqs_col):
	"""
	Resout le probleme avec PLNE
	"""
	N, M = len(seqs_ligne), len(seqs_col)

	lignes = list(range(N))
	colonnes = list(range(M))
	
	m = Model("mogplex")

	# Declaration et initialisation des variables de decision
	x = [[m.addVar(vtype=GRB.BINARY, name="x_%d%d" % (i, j)) for j in colonnes] for i in lignes]
	y = []
	z = []

	for i in lignes:
		tmp_y, tmp_z = [], []
		for j in colonnes:
			tmp_2_y, tmp_2_z = [], []
			for t in range(len(seqs_ligne[i])):
				if (j < sum(seqs_ligne[i][:t]) + len(seqs_ligne[i][:t])) or (j > M - (sum(seqs_ligne[i][t:]) + len(seqs_ligne[i][t+1:]))):
					tmp_2_y.append(0)
				else:
					tmp_2_y.append(m.addVar(vtype=GRB.BINARY, name="y_%d_%d%d" % (t+1, i, j)))
			for t in range(len(seqs_col[j])):
				if (i < sum(seqs_col[j][:t]) + len(seqs_col[j][:t])) or (i > N - (sum(seqs_col[j][t:]) + len(seqs_col[j][t+1:]))):
					tmp_2_z.append(0)
				else:
					tmp_2_z.append(m.addVar(vtype=GRB.BINARY, name="z_%d_%d%d" % (t+1, i, j)))
			tmp_y.append(tmp_2_y)
			tmp_z.append(tmp_2_z)
		y.append(tmp_y)
		z.append(tmp_z)

	# Maj du modele pour integrer les nouvelles variables
	m.update()
	
	# Formulation de  la fonction objectif
	obj = LinExpr();
	obj = 0
	for x_i in x:
		for x_i_j in x_i:
			obj += x_i_j
	
	# Definition de l'objectif
	m.setObjective(obj,GRB.MINIMIZE) 
	
	# Definition des contraintes
	for i in lignes:
		# Dans chaque ligne, il y a autant de cases noires que la somme de la taille des blocs
		m.addConstr(quicksum(seqs_ligne[i][t] for t in range(len(seqs_ligne[i]))) - quicksum(x[i][j] for j in colonnes) == 0)
		# Un bloc ne commence qu'a une seule case dans une ligne : somme des y_ijt = 1
		for t in range(len(seqs_ligne[i])):
			m.addConstr(quicksum(y[i][j][t] for j in colonnes) == 1)
			for j in colonnes:
				# x_ij est toujours superieur a y_ijt
				m.addConstr(x[i][j] - y[i][j][t] >= 0)
				# Si y_ijt = 1, alors la somme des x_ij des cases de j a j + st est egal a st * y_ijt
				m.addConstr(seqs_ligne[i][t]*y[i][j][t] - quicksum(x[i][k] for k in colonnes[j:(j+seqs_ligne[i][t])]) <= 0)
				# Si y_ijt vaut 1 alors la variable x_i(j+st) vaut 0
				b = j + seqs_ligne[i][t]
				if b < M:
					m.addConstr(y[i][j][t] + x[i][b] <= 1)
		# Si y_ijt = 1 alors les y_ijt des cases (i, 0) jusqu'a (i, j + st + 1) valent 0
		for t in range(len(seqs_ligne[i])-1):
			for j in colonnes:
				m.addConstr(y[i][j][t] + quicksum(y[i][k][t+1] for k in colonnes[:j + seqs_ligne[i][t] + 1]) <= 1)

	for j in colonnes:
		# Dans chaque colonne, il y a autant de cases noires que la somme de la taille des blocs
		m.addConstr(quicksum(seqs_col[j][t] for t in range(len(seqs_col[j]))) - quicksum(x[i][j] for i in lignes) == 0)
		# Un bloc ne commence qu'a une seule case dans une colonne : somme des z_ijt = 1
		for t in range(len(seqs_col[j])):
			m.addConstr(quicksum(z[i][j][t] for i in lignes) == 1)
			for i in lignes:
				# x_ij est toujours superieur a z_ijt
				m.addConstr(x[i][j] - z[i][j][t] >= 0)
				# Si z_ijt = 1, alors la somme des x_ij des cases de j a j + st est egal a st * z_ijt
				m.addConstr(seqs_col[j][t]*z[i][j][t] - quicksum(x[k][j] for k in lignes[i:(i+seqs_col[j][t])]) <= 0)
				# Si y_ijt vaut 1 alors la variable x_i(j+st) vaut 0
				b = i + seqs_col[j][t]
				if b < N:
					m.addConstr(z[i][j][t] + x[b][j] <= 1)
		# Si z_ijt = 1 alors les z_ijt des cases (i, j) jusqu'a (i + st + 1, j) valent 0
		for t in range(len(seqs_col[j])-1):
			for i in lignes:
				m.addConstr(z[i][j][t] + quicksum(z[k][j][t+1] for k in lignes[:i + seqs_col[j][t] + 1]) <= 1)

	# Resolution
	m.optimize()

	print ""      
	print 'Solution optimale:'
	print ""
	print 'Valeur de la fonction objectif :', m.objVal

	# Retourne la matrice de x
	aff_x = []
	for i, xi in enumerate(x):
		tmpx = []
		for j, xij in enumerate(xi):
			try:
				tmpx.append(xij.x)
			except:
				tmpx.append(xij)
		aff_x.append(tmpx)
	aff_x = np.array(aff_x)

	return aff_x


def main():
	for i in range(0):
		#if i == 9 or i == 16:
		#	continue
		seqs_ligne, seqs_col = read_file("instances/"+str(i)+".txt")
		debut = time.time()
		opt = compute(seqs_ligne, seqs_col)
		fin = time.time()
		print 'Temps d\'execution pour %d.txt'%i, fin-debut
		print '----------------------------------------'
		for a in range(len(seqs_ligne)):
			for b in range(len(seqs_col)):
				if opt[a][b] == 0:
					opt[a][b] = 2
		#np.save('%d.npy'%i, opt)
		plt.figure()
		plt.imshow(opt,'gray', interpolation='none')
		plt.show()
	for i in sys.argv[1:]:
		try:
			seqs_ligne, seqs_col = read_file(i)
			debut = time.time()
			opt = compute(seqs_ligne, seqs_col)
			fin = time.time()
			print 'Temps d\'execution pour %s'%i, fin-debut
			print '----------------------------------------'
			for a in range(len(seqs_ligne)):
				for b in range(len(seqs_col)):
					if opt[a][b] == 0:
						opt[a][b] = 2
			#np.save('%s.npy'%i, opt)
			plt.figure()
			plt.imshow(opt,'gray', interpolation='none')
			plt.show()
		except FileNotFoundError:
			print "Fichier %s non trouve"%i
			print '----------------------------------------'
			pass

if __name__ == '__main__':
	main()