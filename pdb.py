# -*- coding: utf-8 -*-

"""
@author: Marine FIGAROL & Christelle VIGUIER
"""

# --------------------------------- PDB simple --------------------------------

def parserPDB(fichier) :
    
    """ Obj : Retourne le dictionnaire d'un fichier pdb donné.
        Input : fichier le pdb dont on veut créer le dictionnaire (mettre son
                nom entre guillemets)
        Output : Dictionnaire de la fiche pdb.
    """
    
    data = open(fichier, "r")               # Ouverture du fichier
    dico_proteine = {}                      # Création du dictionnaire
    dico_proteine['domain_list'] = []       # Création de la liste des domaines

    line = data.readline()                  # Lecture de la 1ère ligne
    
    
    # Passe toutes les lignes ne commençant pas par ATOM
    while line[0:6] != "ATOM  " :
        line = data.readline()
        
    # Tant qu'on n'est pas sorti de la protéine
    while line != "" :
    
        if line[0:6] == "ATOM  " :
            # Récupération du domaine
            domain = line[72:76].strip()
            
            # On ne prend que le domaine B ou ceux commençant par A
            # (évite de prendre les molécules d'eau et autres)
            if domain[0] == "A" or domain == "B" :
                
                # Si le domaine n'est pas déja répertorié
                if domain not in dico_proteine['domain_list'] :
                    dico_proteine[domain] = {}
                    dico_proteine['domain_list'].append(domain)
                    dico_proteine[domain]['res_list'] = []
            
                num_residu = int(line[22:26])               # position résidu           
                
                aa = line[17:20].strip()                    # nom du résidu
                    
                # Création dico avec comme clé la position du résidu
                if num_residu not in dico_proteine[domain].keys() :
                    dico_proteine[domain]['res_list'].append(num_residu)
                    dico_proteine[domain][num_residu] = {}
                    dico_proteine[domain][num_residu]["AA"] = aa
                        
                    dico_proteine[domain][num_residu]['atom_list'] = []
                    
                    
                nomAtom = line[12:16].strip(' ')                                         # symbole de l'atome
                dico_proteine[domain][num_residu]['atom_list'].append(nomAtom)
                coor = [float(line[30:38]), float(line[38:46]), float(line[46:54])]   # liste des coordonnees x y z
                temp_fact = float(line[60:66])                                        # facteur de température
                    
                dico_proteine[domain][num_residu][nomAtom] = {}
                dico_proteine[domain][num_residu][nomAtom]["coord"] = coor
                dico_proteine[domain][num_residu][nomAtom]["temp"] = temp_fact
        
        # Passage à la ligne suivante
        line = data.readline()
    
    data.close()
    
    return(dico_proteine)



# ------------------------- PDB multiple --------------------------------------

def parserPDBmulti(fichier) :
    
    """ Obj : Retourne les dictionnaires d'un fichier multi-pdb donné.
        Input : fichier le fichier multi pdb (mettre son nom entre guillemets).
        Output : Dictionnaire de la fiche pdb.
    """
    
    data = open(fichier, "r")               # Ouverture du fichier
    
    # Création du dictionnaire contenant les infos de chacune des conformations
    dico_fichier = {}
    
    line = data.readline()                  # Lecture de la 1ère ligne

    dico_fichier['model_list'] = []         # Création liste des conformations

    # Tant qu'on est pas sorti du fichier
    while line != "" :
        
        # Nouvelle conformation
        if line[0:5] == "MODEL" :
            
            model = int(line[10:14])        # Récupération du numéro de la conformation
            dico_fichier[model] = {}        # Création du dictionnaire pour cette conformation
            dico_fichier['model_list'].append(model)
            
            # Création de la liste des domaines et celle des domaines correspondants
            dico_fichier[model]['domain_list'] = []
        
        
        # Sinon on parcourt la conformation pour créer le dictionnaire
        elif line[0:6] == "ATOM  " :
            
            domain = line[72:76].strip() 
            
            # On ne prend que le domaine B ou ceux commençant par A
            # (évite de prendre les molécules d'eau)
            if domain[0] == "A" or domain == "B" :
                
                # Si le domaine n'est pas déja répertorié
                if domain not in dico_fichier[model]['domain_list'] :
                    dico_fichier[model][domain] = {}
                    dico_fichier[model]['domain_list'].append(domain)
                    dico_fichier[model][domain]['res_list'] = []
            
            
                num_residu = "%s" % (line[22:26]).strip()        # position résidu
                num_residu = int(num_residu)                    # Conversion du string en entier
                
                aa = line[17:20].strip()                        # nom du résidu
                    
                # Création dico avec comme clé la position du résidu
                if num_residu not in dico_fichier[model][domain].keys() :
                    dico_fichier[model][domain]['res_list'].append(num_residu)
                    dico_fichier[model][domain][num_residu] = {}
                    dico_fichier[model][domain][num_residu]["AA"] = aa
                        
                    dico_fichier[model][domain][num_residu]['atom_list'] = []
                    
                    
                nomAtom = line[12:16].strip(' ')                                         # symbole de l'atome
                dico_fichier[model][domain][num_residu]['atom_list'].append(nomAtom)
                coor = [float(line[30:38]), float(line[38:46]), float(line[46:54])]     # liste des coordonnees x y z
                temp_fact = float(line[60:66])                                          # facteur de température
                    
                dico_fichier[model][domain][num_residu][nomAtom] = {}
                dico_fichier[model][domain][num_residu][nomAtom]["coord"] = coor
                dico_fichier[model][domain][num_residu][nomAtom]["temp"] = temp_fact
                
        
        # Passage à la ligne suivante
        line = data.readline()
                
    data.close()
    return(dico_fichier)