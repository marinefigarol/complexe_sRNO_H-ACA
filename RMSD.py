# -*- coding: utf-8 -*-

"""
@author: Marine FIGAROL & Christelle VIGUIER
"""

# --------------------------------- Imports ----------------------------------- 

from distance import *
from math import sqrt
import pylab
import numpy as np

# ----------------------------- RMSD global -----------------------------------

def RMSD_prot (dprot1, dprot2) :
    
    """ Obj :  Calcul du RMSD global entre 2 protéines données. Calcul effectué
                sur les carbones alpha.
        Input : dprot1 et dprot2 les 2 dictionnaires des fichiers pdb.
        Output : RMSD global
    """
    
    somme_dist = 0                                  # Somme des distances entre les résidus des 2 protéines
    N = 0                                           # Pour avoir la taille de la protéine
    
    # Parcours des domaine de la protéine
    for domain in dprot1['domain_list'] :
        # Parcours des résidus du domaine
        for residu in dprot1[domain]['res_list'] :
            # Parcours des atomes du résidu
            for atom in dprot1[domain][residu]['atom_list'] :
                # Si on est sur un carbone alpha
                if atom == 'CA' :
                    # Récupération des coordonnées du carbone alpha des 2 protéines
                    coord1 = dprot1[domain][residu][atom]['coord']
                    coord2 = dprot2[domain][residu][atom]['coord']
                    
                    # Calcul de la distance entre les 2 CA
                    somme_dist += (dist_inter_atom(coord1, coord2))**2
                    N += 1
    
    # Calcul du RMSD global
    RMSD = sqrt(somme_dist/N)
    return (round(RMSD,4))
    


# ----------------------------- RMSD domaine ----------------------------------

def RMSD_domain (dprot1, dprot2, domain) :
    
    """ Obj :  Calcul du RMSD sur un domaine entre 2 protéines données. Calcul
                effectué sur les carbones alpha.
        Input : dprot1 et dprot2 les 2 dictionnaires des fichiers pdb ;
                domain le domaine sur lequel calculer le RMSD
        Output : RMSD du domaine
    """    
    
    somme_dist = 0                                  # Somme des distances entre les résidus des 2 protéines
    N = 0                                           # Pour avoir la taille de la protéine
    
    # Parcours des résidus du domaine désiré
    for residu in (dprot1[domain]['res_list']) :
        # Parcours des atomes du résidu
        for atom in dprot1[domain][residu]["atom_list"] :
           # Si on est sur un carbone alpha
            if atom == 'CA' :
                # Récupération des coordonnées du carbone alpha des 2 protéines
                coord1 = dprot1[domain][residu][atom]['coord']
                coord2 = dprot2[domain][residu][atom]['coord']
                    
                # Calcul de la distance entre les 2 CA
                somme_dist += (dist_inter_atom(coord1, coord2))**2
                N += 1

    # Calcul du RMSD du domaine
    RMSD = sqrt(somme_dist/N)
    return (round(RMSD,4))
    
    
    
# ---------------------------- Fichier RMSD -----------------------------------

def printRMSD(dicRMSD) :
    
    """ Obj :  Produit un fichier contenant les RMSD globaux et ceux des différents
               domaines pour les 500 conformations.
        Input : dicRMSD le dictionnaire contenant les RMSD.
        Output : fichier contenant les RMSD.
    """ 
    
    # Création du fichier qui contiendra les RMSD
    fichier = open("RMSD.txt", "w")
    
    # Création des "colonnes" du fichier
    fichier.write("Conformation\tRMSD global\tDomaine A1\tDomaine A2\tDomaine A3\tDomaine A4\tDomaine B")
    fichier.close()
    
    # On complète maintenant le fichier
    
    fichier = open("RMSD.txt", "a")                              # Mode ajout
    
    for model in dicRMSD['model_list'] :
        
        fichier.write("\n     "+str(model)+" \t\t")               # Numéro de la conformation
        fichier.write(str(dicRMSD[model]['global'])+" \t\t")
        fichier.write(str(dicRMSD[model]['domainA1'])+" \t\t")
        fichier.write(str(dicRMSD[model]['domainA2'])+" \t\t")
        fichier.write(str(dicRMSD[model]['domainA3'])+" \t\t")
        fichier.write(str(dicRMSD[model]['domainA4'])+" \t\t ")
        
    fichier.close()
    
    

# ---------------------------- RMSD croisé ------------------------------------ 
    
def RMSDcroise(dic, domain) :
    
    """ Obj :  Calcul du RMSD entre un domaine et ce même domaine sur toutes les
               conformations.
        Input : dic le dictionnaire des conformations
                domaine sur lequel on veut faire les calculs
        Output : heatmap des RMSD
    """
 
    # On crée une matrice qui contiendra les RMSD, pour ensuite afficher un heatmap
    taille = len(dic['model_list'])
    mat_RMSD = np.zeros((taille, taille), dtype='f')
    
    # Parcours des conformations 2 fois pour pouvoir avoir les variations au sein
    # du domaine au cours du temps
    for i in range(taille) :
        model1 = dic['model_list'][i]
        for j in range(i,taille) :
            model2 = dic['model_list'][j]
            
            somme_dist = 0                                  # Somme des distances entre les résidus des 2 protéines
            N = 0                                           # Pour avoir la taille de la protéine
            
            # Parcours des résidus du domaine désiré
            for residu in (dic[model1][domain]['res_list']) :
                # Parcours des atomes du résidu
                for atom in dic[model1][domain][residu]["atom_list"] :
                    # Si on est sur un carbone alpha
                    if atom == 'CA' :
                        # Récupération des coordonnées du carbone alpha des 2 atomes
                        coord1 = dic[model1][domain][residu][atom]['coord']
                        coord2 = dic[model2][domain][residu][atom]['coord']
                            
                        # Calcul de la distance entre les 2 CA
                        somme_dist += (dist_inter_atom(coord1, coord2))**2
                        N += 1
        
            # Calcul du RMSD du domaine que l'on affecte à la matrice
            mat_RMSD[i][j] = mat_RMSD[j][i] = sqrt(somme_dist/N)
    
    # Tracé du heatmap
    pylab.pcolor(mat_RMSD)
    pylab.colorbar()
    pylab.xlabel('Conformations du domaine ' + domain)
    pylab.ylabel('Conformations du domaine ' + domain)
    pylab.title("RMSD du domaine " + domain + " de la protéine représenté sous forme de heatmap")
    pylab.show()