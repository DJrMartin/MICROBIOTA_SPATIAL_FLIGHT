# MICROBIOTA_SPATIAL_FLIGHT

L'immersion sèche est une technique où l'être humain se retrouve en suspension dans l'eau pendant 5 jours, mimant les effets physiologiques connus en absence de gravité.

## INUPUT :
- **Données de comptage de micro-organismes intestinales** de 16 personnes (humaines) avant et après une immersion sèche. Le séquençage est issu de la technologie 16S. (équipe DMEM, Univ. de Montpellier)
Sur ce jeu de données, les prélèvements biologiques ont été réalisés à J0 et J5.

- **Données de métalome** : ensemble des métaux contenu dans le sérum et les érithrocytes (équipe M2S, Univ de Rennes II)
Sur ce jeu de données, les prélèvements biologiques ont été réalisés à J-1 et J5.

- Données physiologiques : masse musculaire et flexibilité métabolique (équipe Univ de Strasbourg)
- Données de flexibilité métabolique.
- Données de protéomique au niveu du tissu sanguin.

## TRAITEMENT DES DONNEES :
Step 1 : Est-ce que l'échantillonnage (et donc la profondeur de séquençage de chaque échantillon) est suffissante pour avoir accès à toutes les informations ? Courbe de raréfaction.

Step 2 : Quel est la meilleure normalisation pour la différentes approches ? Total sum scaling, logratio additive, logratio centré, cumulative sum scaling ou autres ?

Step 3 : Analyse de la diversité bactérienne (Alpha et Beta) : richesse, Chao1, Simpson, Shannon, Bray-Curtis, Unifrac.

Step 4 : Analyse log poisson normée (ensemble des données) / réduction de dimension sur les données non zéros (lasso/ridge).

Step 5 : Approche sans a priori