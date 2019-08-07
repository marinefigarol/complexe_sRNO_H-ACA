# -*- coding: utf-8 -*-

"""
@author: Marine FIGAROL & Christelle VIGUIER
"""



# --------------------------------- Imports --------------------------------- 

from pdb import *
from RMSD import *
from interface import *
import matplotlib.pyplot as plt
from graphique import *
from carte_contact import *
from residu_cle import *


""" 
     Obj :  Répondre au problème posé dans le projet : quel est le rôle de
            chacun des domaines de la protéine dans l'interaction avec l'ARN ?
 """

#==============================================================================
#
# On parse le fichier contenant la structure de référence et on stocke ses
# informations dans un dictionnaire pour pouvoir l'analyser par la suite.
#
#==============================================================================

reference = parserPDB("pab21_structure_de_ref.pdb")

#==============================================================================
#
# On parse le fichier contenant les 500 structures et on stocke également ses
# informations dans un dictionnaire.
#
#==============================================================================

conformations = parserPDBmulti("pab21_prod_solute_500frames.pdb")

#==============================================================================

# Calcul des RMSD pour chacune des conformations : permet de mesurer la
# flexibilité des différents domaines de la protéine lors de l'interaction avec
# l'ARN. On stocke chacun des RMSD dans un dictionnaire.

#==============================================================================

RMSD = {}                                                                         # Stockage des RMSD
RMSD['model_list'] = []                                                           # Stockage de conformations
# Stockage des RMSD pou tracer les graphes
y_glob = []
y_A1 = []
y_A2 = []
y_A3 = []
y_A4 = []

for model in conformations['model_list'] :                                        # Parcours des conformations
    RMSD[model] = {}
    RMSD['model_list'].append(model)
    RMSD[model]['global'] = RMSD_prot(reference, conformations[model])            # Calcul du RMSD global
    RMSD[model]['domainA1'] = RMSD_domain(reference, conformations[model], 'A1')  # Calcul du RMSD sur le domaine A1
    RMSD[model]['domainA2'] = RMSD_domain(reference, conformations[model], 'A2')  # Calcul du RMSD sur le domaine A2
    RMSD[model]['domainA3'] = RMSD_domain(reference, conformations[model], 'A3')  # Calcul du RMSD sur le domaine A3
    RMSD[model]['domainA4'] = RMSD_domain(reference, conformations[model], 'A4')  # Calcul du RMSD sur le domaine A4
    #RMSD[model]['domainB'] = RMSD_domain(reference, conformations[model], 'B')    # Calcul du RMSD sur le domaine B    
    
    
    # Pour pouvoir tracer les graphes, on stocke les RMSD dans des listes
    y_glob.append(RMSD[model]['global'])
    y_A1.append(RMSD[model]['domainA1'])
    y_A2.append(RMSD[model]['domainA2'])
    y_A3.append(RMSD[model]['domainA3'])
    y_A4.append(RMSD[model]['domainA4'])

#==============================================================================

# Pour tracer les RMSD (en commentaire : on le trace uniquement pour le rapport)

#==============================================================================

#graphe_RMSD(conformations['model_list'], y_glob, 'global')
#graphe_RMSD(conformations['model_list'], y_A1, 'A1')
#graphe_RMSD(conformations['model_list'], y_A2, 'A2')
#graphe_RMSD(conformations['model_list'], y_A3, 'A3')
#graphe_RMSD(conformations['model_list'], y_A4, 'A4')


#==============================================================================

# On trace également les heatmaps de RMSD de quelques domaines pour le rapport

#==============================================================================

#RMSDcroise(conformations, 'A1')
#RMSDcroise(conformations, 'A2')
#RMSDcroise(conformations, 'A3')
#RMSDcroise(conformations, 'A4')


#==============================================================================

# On trace des heatmaps de matrice de distance des domaines du fichier de
# référence contre l'ARN pour le rapport

#==============================================================================

#carte_contact(reference, 'A1')
#carte_contact(reference, 'A2')
#carte_contact(reference, 'A3')
#carte_contact(reference, 'A4')


#==============================================================================

# On crée un fichier nommé "RMSD.txt" contenant les résultats

#==============================================================================

printRMSD(RMSD)



#==============================================================================

# Calcul de la fréquence d'appartenance à l'interface de chacun des résidus de
# la protéine tout au long de la dynamique.

#==============================================================================

# Recherche des résidus appartenant à l'interface, ie ayant un une distance de
# moins d'un certain seuil en Angstrom avec les résidus de l'ARN
res_interface = interface(conformations, 9)

# On crée un fichier nommé "interface.txt" contenant les résultats
printInterface(res_interface)


#==============================================================================

# Calcul de la fréquence de contact entre des résidus clés choisis à partir
# de l'analyse du mutant
 
#==============================================================================

# Liste des résidus clés
paires = {41 : {'dom1':'A4', 'res2': 32, 'dom2':'B'},
         100 :{'dom1':'A4', 'res2': 31, 'dom2':'B'},
         46 : {'dom1':'A4', 'res2': 25, 'dom2':'B'},
         34 : {'dom1':'A3', 'res2': 33, 'dom2':'B'},
         6 : {'dom1':'A3', 'res2': 63, 'dom2':'A4'},
         66 : {'dom1':'A4', 'res2': 41, 'dom2':'A3'},
         98 : {'dom1':'A4', 'res2': 30, 'dom2':'B'},
         26 : {'dom1':'A3', 'res2': 65, 'dom2':'A4'}}

# Calcul des fréquences de contact pour un seuil de 6 Angstrom et création d'un
# fichier nommé 'residus_cles.txt" contenant les résultats
freq_res_cle(paires, conformations, 9)