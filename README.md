# Etude du complexe sRNO H/ACA

Ce projet consiste à étudier le complexe sRNO H/ACA de Pyrococcus abyssi. Il a été réalisé dans le cadre de ma première année de master de bioinformatique. Le compte-rendu de ce projet est disponible dans le fichier *CR_python.pdf*. Les graphes générés sont également disponibles.

## Lancer le programme

Pour lancer ce programme, vous avez besoin des packages Python <code>matplotlib</code>, <code>pylab</code>, <code>numpy</code> et <code>math</code>.<br>
Une fois les fichiers téléchargés, placez-vous dans le dossier et lancer le programme :<br>
<pre>
<code>cd complexe_sRNO_H-ACA
python3 main.py</code>
</pre>

## Analyse des changements conformationnels globaux

La première étape consiste à écrire un programme permettant d'étudier si la protéine subit des changements conformationnels majeurs en solution. Une étude au cours du temps du RMSD a été réalisée (RMSD entre chacune des conformations sélectionnées au cours de la dynamique et la structure de référence). Le RMSD pour chacun des cinq domaines au cours de la dynamique a également été calculé afin de voir si des domaines sont plus flexibles que d'autres, c'est-à-dire si certains subissent plus de changements conformationnels. Le programme crée en sortie un fichier contenant pour la structure d'origine et pour chaque conformation (chaque fichier pdb de la dynamique) le RMSD global et celui de chacun des cinq domaines.

## Analyse des changements conformationnels locaux

La seconde partie consiste à calculer la fréquence d'appartenance à l'interface de chacun des résidusde la protéine tout au long de la dynamique. Pour cela, il faut identifier pour chaque conformation les résidus appartenant à l'interface et calculer la fréquence pour l'ensemble des conformations. La fréquence de contact entre différentes paires de résidus clés choisis est également calculée. Le seuil d'appartenance à cette interface a été fixée à 9 Å.

## Données

Le fichier *pab21_prob_500frames.pdb* contient les structures au format PDB des 500 conformations extraites de la dynamique moléculaire sélectionnées. Le fichier *pab21_structure_de_ref.pdb* contient la structure de référence de la protéine.
