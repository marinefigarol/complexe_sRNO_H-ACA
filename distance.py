# -*- coding: utf-8 -*-

"""
@author: Marine FIGAROL & Christelle VIGUIER
"""

from math import sqrt

# ------------------------------- Distance CA --------------------------------- 

def dist_inter_atom(coord1, coord2) :
    
    """ Obj : Calcul de la distance interatomique entre les 2 carbones alpha.
        Input : coord1 liste des coordonnées du 1er atome
                coord2 list des coordonnées du 2nd atome
        Output : distance
    """
    
    # Nb de coordonées (ici on est en 3D, donc 3)
    taille = len(coord1)
    # Initialisation de la distance
    d = 0
    # Calcul de la distance inter-atomique
    for i in range(taille) :
        d += (coord1[i]-coord2[i])**2
    return(round(sqrt(d),4))
    

# -------------------------------- Barycentre ---------------------------------    
 
def barycentre(coord):
    
    """ Obj : Calcul du barycentre d'un résidu
        Input : coord liste des coordonnées du résidu
        Output : barycentre
    """
    
    taille = len(coord)                # Nb de coordonnées (ici, nombre d'atomes)
    # initialisation
    x = 0
    y = 0
    z = 0
    # calcul de la somme des coordonnées
    for i in range(taille) :
        x += coord[i][0]
        y += coord[i][1]
        z += coord[i][2]
        
    # coordonnées du barycentre
    x = x/taille
    y = y/taille
    z = z/taille 
    
    return([x,y,z]) 