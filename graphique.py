# -*- coding: utf-8 -*-

"""
@author: Marine FIGAROL & Christelle VIGUIER
"""


# --------------------------------- Imports ----------------------------------- 

import matplotlib.pyplot as plt


# ------------------------------Graphique RMSD -------------------------------- 

def graphe_RMSD(x,y,domaine) :
    
    """ Obj :  Permet de tracer les graphiques de RMSD en fonction du temps.
        Input : x la liste des temps/conformations
                y la liste des RMSD
                domaine le domaine étudié (si global, écrire global)
        Output : Graphe RMSD
    """
    
    plt.plot(x,y, 'dodgerblue')
    plt.xlabel('Temps (en ns)')
    plt.ylabel('RMSD (en Angstrom)')
    plt.title('RMSD global de la protéine') if domaine == 'global' else plt.title("RMSD du domaine " + domaine + " de la protéine")
    plt.show()