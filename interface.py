# -*- coding: utf-8 -*-

"""
@author: Marine FIGAROL & Christelle VIGUIER
"""

# --------------------------------- Imports ----------------------------------- 

import pylab
from distance import *
import numpy as np
from math import inf
from pdb import *


# --------------------------------- Interface --------------------------------- 

def interface(dic, seuil) :
    
    
    """ Obj : Savoir si un résidu appartient à l'interface entre la protéine
              et l'ARN. Calcul à partir des barycentres.
        Input : dic le dictionnaires du fichier pdb
                seuil le seuil d'appartenance à l'interface.
        Output : zone d'interface
    """
    
#==============================================================================

# Calcul des distances

#==============================================================================
    
    dom_prot = ['A1', 'A2', 'A3', 'A4']                         # Domaines de la protéine
    
    # Création d'un dictionnaire qui contiendra les résidus appartenant à l'interface
    res_interface = {}
    
    
    #Parcours des domaines de la protéine
    for domain in dom_prot :    
    
        # Stocke les résidus du domaine appartenant à l'interface
        res_interface[domain] = {}
        res_interface[domain]['list'] = [] 
    
        # Parcours des conformations
        for model in dic['model_list'] :
            
            # On récupère le résidu de la protéine
            for residu1 in dic[model][domain]['res_list'] :
                
                min = inf                           # On fixe la distance min à l'infini
                # Pour stocker les coordonnées du résidu de la protéine
                coord1 = []
                
                # On récupère les coordonnées des atomes du premier résidu
                for atom1 in dic[model][domain][residu1]['atom_list'] :
                    
                    coord1.append(dic[model][domain][residu1][atom1]['coord'])
                
                # Caclul du barycentre du résidu de la protéine
                bar1 = barycentre(coord1)
                
                # On récupère le résidu de l'ARN
                for residu2 in dic[model]['B']['res_list'] :
                    # Pour stocker les coordonnées du résidu de l'ARN
                    coord2 = []
                    
                    # On récupère les coordonnées des atomes du second résidu
                    for atom2 in dic[model]['B'][residu2]['atom_list'] :
                        
                        coord2.append(dic[model]['B'][residu2][atom2]['coord'])
                        
                    # Caclul du barycentre du résidu de l'ARN
                    bar2 = barycentre(coord2)
                    
                    # Caclul de la distance entre les 2 résidus
                    dist = dist_inter_atom(bar1, bar2)
                    #Recherche de la distance min entre le résidu et l'ARN
                    min = dist if dist < min else min
                
                # On regarde si le résidu appartient à l'interface
                if min < seuil :
                    # S'il n'appartenait pas déja à l'interface
                    if residu1 not in res_interface[domain]['list'] :
                        res_interface[domain]['list'].append(residu1)
                        
                        # Initialisation du nombre de fois où il appartient à l'interface
                        res_interface[domain][residu1] = 0                  
                        
                    res_interface[domain][residu1] += 1
    
    return(res_interface)
    
    

def printInterface(dicInterface) :
    
    """ Obj :  Produit un fichier contenant les fréqauences d'appartenance à
            l'interface 
        Input : dicRMSD le dictionnaire contenant les RMSD.
        Output : fichier contenant les RMSD.
    """ 
    
    # Création du fichier qui contiendra les RMSD
    fichier = open("interface.txt", "w")
    
    dom_prot = ['A1', 'A2', 'A3', 'A4']                         # Domaines de la protéine
    
    # Création des "colonnes" du fichier
    fichier.write("Domaine\t Residu n°\tFréquence")
    fichier.close()
    
#==============================================================================

# Calcul des fréquences et écriture dans le fichier

#==============================================================================
    
    # On complète maintenant le fichier
    
    fichier = open("interface.txt", "a")                              # Mode ajout
    
    # Parcours des domaines de la protéine
    for domain in dom_prot :
        # On parcours la liste des résidus qui appartiennent à l'interface
        for res in dicInterface[domain]['list'] :
            # Calcul de la fréquence
            dicInterface[domain][res] = dicInterface[domain][res]/500
            # Ajout dans le fichier
            fichier.write("\n\t" + domain + "\t\t " + str(res) + " \t " + str(dicInterface[domain][res]))
    
    fichier.close()
