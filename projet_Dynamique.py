#Projet MOGPL
#Binh Thanh LUONG - 3504859
#Melissa YAYA - 3673113

import numpy as np
import time
import matplotlib.pyplot as plt
import sys
import os.path
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

FREE = 0
BLACK = 1
WHITE = 2

def read_file (filename):
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
	return lignes,colonnes

def test_si_aucun(V ,j1, j2, val):
	"""
	Renvoie vrai si la valeur val n'apparait dans aucune des cases
	d'un vecteur V entre les cases j1 et j2 (comprises).
	"""
	for i in range(j1,j2+1):
		if V[i] == val:
			return False
	return True

def T (V, seq, j, l):
	"""
	Renvoie vrai s'il est possible de mettre une sous-sequence de s dans V
	"""
	if not l:
		return test_si_aucun(V,0,j,BLACK)
	if l == 1 and j == seq[l-1] - 1:
		return test_si_aucun(V,0,j,WHITE)
	if j <= seq[l-1] - 1:
		return False
	if V[j] == WHITE:
		return T (V, seq, j-1, l)
	else:
		if not test_si_aucun(V,j-seq[l-1]+1,j,WHITE):
			i = j
			while i >= 0:
				if V[i] == WHITE:
					break
				i -= 1
			return T(V,seq,i-1,l) and test_si_aucun(V,i,j,BLACK)
		else:
			if V[j-seq[l-1]] == BLACK:
				if V[j] == BLACK:
					return False
				else:
					return T(V,seq,j-1,l)
			else:
				if (V[j] == BLACK):
					return T(V,seq,j-seq[l-1]-1,l-1)
				else:
					return T(V,seq,j-seq[l-1]-1,l-1) or T(V,seq,j-1,l)

def T_ligne (M, seqs, i):
	"""
	Applique T sur la ligne i de la matrice M
	"""
	ligne = M[i]
	s = seqs[i]
	return T(ligne,seqs[i], len(ligne) - 1, len(s))

def T_col (M, seqs, j):
	"""
	Applique T sur la colonne j de la matrice M
	"""
	col = M[:,j]
	s = seqs[j]
	return T(col,seqs[j], len(col) - 1, len(s))

def propag_ligne (M, seqs, i):
	"""
	Colorie une ligne
	"""
	for j in range(len(M[0])):
		if not M[i][j]:
			M[i][j] = WHITE
			c1 = T_ligne(M,seqs,i)
			M[i][j] = BLACK
			c2 = T_ligne(M,seqs,i)
			M[i][j] = FREE
			if (not c1) and (not c2):
				return False
			if c1 and (not c2):
				M[i][j] = WHITE
			if (not c1) and c2:
				M[i][j] = BLACK
	return True

def propag_col (M, seqs, j):
	"""
	Colorie une colonne
	"""
	for i in range(len(M[:,0])):
		if not M[i][j]:
			M[i][j] = WHITE
			c1 = T_col(M,seqs,j)
			M[i][j] = BLACK
			c2 = T_col(M,seqs,j)
			M[i][j] = FREE
			if (not c1) and (not c2):
				return False
			if c1 and (not c2):
				M[i][j] = WHITE
			if (not c1) and c2:
				M[i][j] = BLACK
	return True

def propagation (Mat, seqs_ligne, seqs_col):
	"""
	Colorie une matrice
	"""
	N,M = len(seqs_ligne),len(seqs_col)
	lignesAVoir = [i for i in range(N)]
	colonnesAVoir = [j for j in range(M)]
	while lignesAVoir or colonnesAVoir:
		for i in lignesAVoir:
			if not propag_ligne(Mat,seqs_ligne,i):
				return False
			for j in range(M):
				if Mat[i][j] == FREE:
					if j not in colonnesAVoir:
						colonnesAVoir.append(j)
		lignesAVoir = []
		for j in colonnesAVoir:
			if not propag_col(Mat,seqs_col,j):
				return False
			for i in range(N):
				if Mat[i][j] == FREE:
					if i not in lignesAVoir:
						lignesAVoir.append(i)
		colonnesAVoir = []
	return True

def main():
	for i in range(0):
		l,c = read_file('instances/%d.txt'%i)
		m = np.zeros((len(l),len(c)), dtype=np.int)
		debut = time.time()
		propagation(m,l,c)
		fin = time.time()
		print 'Temps d\'execution pour %d.txt'%i, fin-debut
		#np.save('%s.npy'%i, m)
		plt.figure()
		plt.imshow(m,'gray', interpolation='none')
		plt.show()
	for i in sys.argv[1:]:
		try:
			l,c = read_file(i)
			m = np.zeros((len(l),len(c)), dtype=np.int)
			debut = time.time()
			propagation(m,l,c)
			fin = time.time()
			print 'Temps d\'execution pour %s'%i, fin-debut
			#np.save('%s.npy'%i, m)
			plt.figure()
			plt.imshow(m,'gray', interpolation='none')
			plt.show()
		except FileNotFoundError:
			print "Fichier %s non trouve"%i
			print '----------------------------------------'
			pass

if __name__ == '__main__':
	main()

