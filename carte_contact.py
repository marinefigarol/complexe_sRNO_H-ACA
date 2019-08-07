# -*- coding: utf-8 -*-

"""
@author: Marine FIGAROL & Christelle VIGUIER
"""

# --------------------------------- Imports ----------------------------------- 


import pylab
from distance import *
import numpy as np
    

# ---------------------------- Carte contact ----------------------------------

def carte_contact(dic, domain) :

    """ Obj : Ce programme renvoie la carte de contact domaine-ARN associée
              à un fichier PDB. Calcul des distances sur le barycentre.   
        Input : dic le dictionnaire des conformations
                domain le domaine sur lequel on calcule les distances
        Output : Carte de contact résidu-résidu et matrice de distance
    """

    # On récupère la taille du domaine
    taille_dom = len(dic[domain]['res_list'])
    # Idem la taille de l'ARN
    taille_ARN = len(dic['B']['res_list'])
    
    # Création de la matrice de distance
    mat_dist = np.zeros((taille_dom, taille_ARN), dtype='f')
#        
    # On récupère le résidu du domaine
    for i in range(taille_dom) :
        residu1 = dic[domain]['res_list'][i]
        # Pour stocker les coordonnées du premier résidu
        coord1 = []
        
        # On récupère les coordonnées des atomes du résidu du domaine
        for atom1 in dic[domain][residu1]['atom_list'] :
            
            coord1.append(dic[domain][residu1][atom1]['coord'])
        
        # Caclul du barycentre du résidu du domaine
        bar1 = barycentre(coord1)
        
        # On récupère le résidu de l'ARN
        for j in range(taille_ARN) :
            residu2 = dic['B']['res_list'][j]
            # Pour stocker les coordonnées du résidu de l'ARN
            coord2 = []
            
            # On récupère les coordonnées des atomes du résidu de l'ARN
            for atom2 in dic['B'][residu2]['atom_list'] :
                
                coord2.append(dic['B'][residu2][atom2]['coord'])
                
            # Caclul du barycentre du résidu de l'ARN
            bar2 = barycentre(coord2)
            
            # Caclul de la distance entre les 2 résidus et stockage dans la matrice
            mat_dist[i][j] = dist_inter_atom(bar1, bar2)
    
    # On trace le heatmap
    pylab.pcolor(mat_dist)
    pylab.colorbar()
    pylab.ylabel('Domaine ' + domain)
    pylab.xlabel('ARN')
    pylab.title("Matrice de distance du domaine " + domain + " de la protéine "
    "contre l'ARN représenté sous forme de heatmap")
    pylab.show()