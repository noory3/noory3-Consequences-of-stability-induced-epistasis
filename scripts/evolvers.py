#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:10:07 2018

@author: N. Youssef

Includes functions needed to evolve along a tree
"""
import numpy as np
import helpers


############# BUILDING TREE #############
def tree_counter(T):
	Tips = T.count(',') + 1 
	Nodes = Tips - 1 
	return Tips, Nodes
        
def find_and_position(T, pattern):
	return [pos for pos, char in enumerate(T) if char == pattern]

def build_tree(T, D):
	"""
        Returns an array specifiing the child, branch length, parent 
	"""
	leftID = find_and_position(T, "(")
	rightID = find_and_position(T, ")")
	remID =sorted( list(set(list(range(0,len(T)))) - set(leftID+rightID))) #find the position for , and sequence numbers
	leafnodes = []
	for i in remID:
		if T[i] == "," or T[i] == ";":
			pass
		else:
			leafnodes.append(T[i])
	n = T.count(',') + 1  #number of sequences
	dets = np.nan*np.ones((2*(n) - 2,3)) 
	doneyet = 0 
	currentT = T
	counter = n 
	while doneyet == 0:
		leftID = list(find_and_position(currentT, "("))
		rightID = list(find_and_position(currentT, ")"))
		idx2 = min(rightID)
		leftID = [item for item in leftID if item <= idx2]
		idx1 = max(leftID)
		daughters= currentT[idx1+1:idx2].split(",")
		dets[int(daughters[0])] = [daughters[0], D[int(daughters[0])], counter]
		dets[int(daughters[1])] = [daughters[1], D[int(daughters[1])], counter]
		segL = currentT[0:idx1]
		segR = currentT[idx2+1:]
		newT = segL + str(counter) + segR
		#print(newT)
		currentT = newT
		counter = counter +1 
		#print(daughters)
		if len(list(find_and_position(currentT, "("))) == 0:
			doneyet = 1 
	#print(dets)
	return dets


############# Evolvers #############
def EvolveSeq(ParentSeq, BL, num_cores, Neff, GTR, ContactMapNS):
    '''
        Given a start sequence, returns the new sequence after time interval (BL).  
    '''
    Q = helpers.TransitionVector(ParentSeq, num_cores, Neff, GTR, ContactMapNS)
    site, newCodon, rate = helpers.NextSubstitution(Q)
    s = np.random.exponential(1/rate)
    if s > BL:
        return(ParentSeq)
    else:
       ChildSeq = list(ParentSeq)
       ChildSeq[site] = newCodon
       print("site: " + str(site) + "  newCodon: " + str(newCodon), flush = True)
       ChildSeq = EvolveSeq(ChildSeq, (BL-s), num_cores, Neff, GTR, ContactMapNS)
       return(ChildSeq)


def EvolveAlongTree(Tips, Nodes, dets, StartSeq, num_cores, Neff, GTR, ContactMapNS):
	Sequences = np.nan*np.ones((Tips + Nodes, len(StartSeq))) #array that contains all the sequences for nodes and tips 
	Sequences[len(Sequences)-1]  = StartSeq   #set root seq
	if dets[0,0] == 0:
		dets = np.flipud(dets)
	for branch in range(0, len(dets)):
		print("branch: " + str(branch), flush = True)
		parent = Sequences[int(dets[branch,2])] 
		BL = dets[branch,1]
		child = EvolveSeq(parent, BL, num_cores, Neff, GTR, ContactMapNS) 
		Sequences[int(dets[branch,0])] = child 
	return Sequences


############# WRITE TO OPUTPUT FILE ############# 
def _write_output_file(Sequences, Tips, output): 
	outfile = open( output , "w" )
	outfile.write(str(Tips) + "  " + str(len(Sequences[0])*3) + " \n")
	for k in range(0, Tips):
		seq = ''
		for n in range(0, len(Sequences[k])):
			seq += helpers.Codon[Sequences[k,n]]
		print('sequence' + str(k+1) + '\n' + seq, file = outfile)
	outfile.close()
