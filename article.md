# Projet MegaBot

Le ***MEGABOT*** est un projet lancé par Julien Allali, il s'agit d'un robot quadrupède d'une taille plutot impressionante de 3m x 3m ,réalisé intégralement en métal et controlé à partir d'une manette de jeux-vidéos. L'entiereté du hardware et du software est open-source et disponible sur le Git d'Eirlab.

Une première version de l'électronique embarquée et du contrôle avait déja été réalisée par Julien et visionnable sur sa chaine personnelle ici (HYPERLINK). Elle posait des problèmes de vitesse et de fiabilité niveau électronique, le travail récent à donc été concentré sur ces deux points en y ajoutant une touche de modularité et de réparabilité. 

Du point de vue de l'électronique embarquée il s'agit principalement de distribution de la puissance et du controle des actionneurs. Concernant le premier point, la source d'énergie utilisée sont deux batteries Pb, elles sont reliées à un boiter de distribution qui permet d'alimenter les quatres pattes avec interruption et affichage de la puissance consommée sur un écran.

(INSERER PHOTO)

Pour le controle des actionneurs, la gestion de la puissance se fait indépendemment pour chaque patte et le tout est dirigé par un PC situé sur la plateforme centrale. Au sein de chaque boitier se situe deux cartes, une permmettant l'interfacage avec l'ordinateur et l'autre comportant l'électronique de puissance nécéssaire au déplacement des vérins. La carte de contrôle est pilotée par un microcontroleur STM32 permettant l'interfacage avec le PC et dotée d'un afficheur LCD rendant lisible les informations caractéristiques de chaque patte (adresse, élongation des vérins..). La carte de puissance comporte trois drivers moteur pour chaque vérin et deux convertisseurs buck qui alimentent le ventilateur et la logique. L'ensemble de l'électronique à été concue, réalisée et testée à Eirlab, beaucoup de temps à été investi sur la fiabilité, modularité et réparabilité du système, permmettant une bonne gestion des imprévus en protégeant à la fois l'humain et la machine.

(INSERER PHOTO DES DEUX CARTES ET DU BOITIER)

###(PARTIE INFO)

Pour la partie **informatique** de ce projet, un code fonctionnel et complet avait déjà été réalisé par Julien, avec plusieurs pistes d'amélioration.

Premièrement pour la partie **contrôle**, la cinématique indirecte était réalisée par une convergence de la cinématique directe, via la librairie *minimize*. En utilisant les matrices jacobiennes du déplacement des différents points de la patte, et en travaillant donc sur une discrétisation des mouvements en micro-déplacement, on obtient un résultat plus rapide. En y ajoutant un solveur pour contraindre par exemple les élongations maximales et minimales des vérins, on obtient une cinématique inverse correcte.
Un passage par une base absolue à aussi été proposé pour simplifier certaines phase du contrôle et du planning.

Un travaille de refonte de la partie **planning** a aussi été effectué pour correspondre à cette nouvelle vision "discrétisée" du mouvement. De plus l'utilisation du solveur de la librairie QPsolveur permet d'imposer des contraintes telles que la prise en compte du déplacement du centre de masse dans le triangle de sustentation.
Un des objectifs a été de réaliser une marche *courbe*, i.e. suivant des trajectoires circulaire, en faisant varier simplement le rayon pour aller plus ou moins droit.
Dans le but de pouvoir contrôler le MegaBot avec une manette de Xbox, une méthode de calcul de trajectoire à partir d'un vecteur unitaire de direction ainsi qu'un facteur de rotation a été implémentée.

Il resterait donc à combiner toutes ces améliorations pour réaliser une marche pour le MegaBot


L'équipe ayant contribuée au projet :
Julien ALLALI, Guillaume FORNES, Marc DUCLUSAUD, Victor BREGEON

Le projet est ouvert à de nouvelles participations, n’hésitez pas à nous contacter !