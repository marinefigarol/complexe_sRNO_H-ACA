# -*- coding: utf-8 -*-

"""
@author: Marine FIGAROL & Christelle VIGUIER
"""

# --------------------------------- Imports ----------------------------------- 

from distance import *


# ------------------------------ Résidus clés ---------------------------------

def freq_res_cle(paires, dic, seuil) :
    
    """ Obj : Calcule la fréquence d'interaction des résidus clés.
        Input : paires les paires de résidus clés choisis
                dic le dictionnaire des conformations
                seuil le seuil en Angstrom qui définit un contact entre les résidus
        Output : fichier texte contenant les fréquences de contacts entre les
                 résidus clés.
    """
    
    # Création du fichier qui contiendra les fréquence de contact
    fichier = open("residus_cles.txt", "w")
    
    # Création des "colonnes" du fichier
    fichier.write("Residu n°\t Domaine\t Residu n°\t Domaine\t Fréquence")
    fichier.close()    
    
    
    # On complète maintenant le fichier 
    fichier = open("residus_cles.txt", "a")                              # Mode ajout
    
    # Parcours de la liste des résidus clés sélectionnés
    for residu in paires.keys() :
        # Initialisation de la fréquence de contact
        paires[residu]['freq'] = 0
        
        # Parcours des conformations pour pouvoir calculer la fréquence
        for model in dic['model_list'] :
        
            # Pour stocker les coordonnées du premier résidu
            coord1 = []
        
            # On récupère les coordonnées des atomes du premier résidu
            for atom1 in dic[model][paires[residu]['dom1']][residu]['atom_list'] :
            
                coord1.append(dic[model][paires[residu]['dom1']][residu][atom1]['coord'])
        
            # Caclul du barycentre du premier résidu
            bar1 = barycentre(coord1)
        
            # Pour stocker les coordonnées du second résidu
            coord2 = []
            
            # On récupère les coordonnées des atomes du second résidu
            for atom2 in dic[model][paires[residu]['dom2']][paires[residu]['res2']]['atom_list'] :
                
                coord2.append(dic[model][paires[residu]['dom2']][paires[residu]['res2']][atom2]['coord'])
                
            # Caclul du barycentre du second résidu
            bar2 = barycentre(coord2)
            
            # Caclul de la distance entre les 2 résidus
            dist = dist_inter_atom(bar1, bar2)
        
            # Si on est inférieur au seuil on a bien contact
            if dist <= seuil :
                paires[residu]['freq'] += 1
        
        # Ecriture dans le fichier
        fichier.write("\n\t" + str(residu) + "\t\t\t" + paires[residu]['dom1'] + "\t\t\t"
        + str(paires[residu]['res2']) + "\t\t\t" + paires[residu]['dom2'] + "\t\t\t" +
        str(paires[residu]['freq']/500))
    
    