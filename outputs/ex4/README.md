# Exercici 4

L'objectiu en aquest apartat és avaluar la qualitat de l'integrador numèric construit i de la funció "flux" també construida mitjançant les funcions amb condicions inicials: $$y' = 2y/t \; , \; \; y(1) = 1$$ amb solució $y(t) = t^2$ i del sistema: $$x' = y \; , \; \; y' = -x \; , \; \; \big(x(0), y(0)\big) = (1,0)$$

El codi associat és "scripts/ex4.py". No es demanen fitxers d'entrada ni es guarden resultats sinó que dins del mateix codi ja s'escull `t_pred`, el temps a predir i en executar el programa, es retorna amb `print()` el l'avaluació obtinguda i la comparació amb les solucions analítiques també avaluades en `t_pred`.