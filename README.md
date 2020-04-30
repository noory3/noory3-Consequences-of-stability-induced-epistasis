
# Consequences of Stability Induced Epistasis 
Repository for "Consequences of stability-induced epistasis for substitution rates." Noor Youssef, Ed Susko, and Joe Bielawski. 

# /Alignments
- contains natural protein alignments ("PDB_numbered.txt") and alignments simulated under three generative models (C-SI, S-SI, and S-SD) for three protein systems (PDB: 1QHW, 2PPN, and 1PEK). 

# /data
  simulation_parameters.txt 
   - contains the protein-specific tree, branch lengths, mutation parameters used for simulations
   - includes the equilibrated protein-specific root sequences  
    
  PDB_Ne2_branch_scale_factor.txt (PDB = 1qhw, 2ppn, and 1pek) 
  - contains the waiting time per 1000 substitution to be used to rescale simulation branch lengths 
  
# /Figures 
  contains all figures in the manuscript and supplementary material 
  
# /Results 
  /Inference
   - contains inference model results [M0, M3(k=2), M3(k=3), CLM3, BUSTED] for all simulated alignments under C-SI, S-SI, S-SD and three protein structures (1QHW, 2PPN, 1PEK). 
    
  /real_pdb
   - contains inference model results [M0, M3(k=2), M3(k=3), CLM3, BUSTED] for the respective real protein alignment (PDB = 1QHW, 2PPN, 1PEK)
  
  /site_specific_fitness
   - contains the site-specific fitness landscapes for each protein and each generative model. 
   - For C-SI simulations, the PDB_Ne2_C-SI folder contains the fitness landscape used for each simulated alignment (1, …50). 
   - For S-SI simulations, the PDB_Ne2_S-SD_AvgssFit folder contains the fitness landscape used for each simulated alignment (1, …, 50). 
   - For S-SD simulations, the PDB_Ne2_S-SD folder contains the fitness landscape used for each simulated alignment (1, …50).

  /dNdS
   - contains the nonsynonymous substitution rates (Kn) and the nonsynonymous mutation rates (Ln) used to calculate expected site-specific dN/dS for all simulations. 

  /post_prob_C_PDB_Ne2_GEN_MODEL_M3_k2 (PDB = 1QHW, 2PPN, or 1PEK and GEN_MODEL = C-SI, S-SI, S-SD and C = 1 or 2 ) 
  - contains the posterior probability for a site belonging to rate class C = 1 or 2 for simulation under the respective GEN_MODEL and protein structure (PDB)

# /scripts

/evolvers.py 
- contains functions necessary for evolving sequences along a phylogeny 

/helpers.py 
- contains general functions used to calculate free energy, and more.

/sequence_optimizer.py
- Starting at a random sequence, evolve using algorithm described in supp material (table S7) to reach a sequence with high fitness 
- The equilibrated sequences used are provided in “../data/simulation_parameters.txt”

/get_branch_scale_factor.py
- Evolve the equilibrated high fit sequence for 1000 substitution while keeping track of the waiting times. This value will be used to rescale branch length so that they have the desired interpretation as the expected number of substitutions per site.
- Scaling values are provided in “../data/PDB_Ne2_branch_scale_factor.txt” for PDB = 1qhw, 2ppn, and 1pek

/generate_S-SD.py
- generates alignments under S-SD model

/calc_avg_ssfitness_dNdS_S-SD.py
- calculates expected dN^h/dS^h for S-SD simulations
- calculates the average site-specific fitness landscape per alignment to be used for S-SI simulations 

/calc_dNdS_S-SI.py
- calculates expected dN^h/dS^h for S-SI simulations

/calc_dNdS_C-SI.py
- calculates expected dN^h/dS^h for C-SI simulations

/freq_to_fitness_converter_C-SI.py
- converts the c20 frequency profiles to protein-specific fitness vectors
/plotting_figureX.py 
- scripts to plot respective figure

/TableSX.py 
- scripts to get values presented in the supp tables 

# Consequences-of-stability-induced-epistasis
# noory3-Consequences-of-stability-induced-epistasis
