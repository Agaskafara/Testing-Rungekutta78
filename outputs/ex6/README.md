# Exercici 6

L'objectiu en aquest apartat és calcular i plotar les òrbites halo en el problema RTBP.

El codi associat es "scripts/ex6.py" i requereix:

* El path del fitxer amb les múltiples condicions inicials, passat com argument `--halos_input`. Per defecte, tenim el fitxer "outputs/ex6/condicions\_inicials/halos\_inp.txt".

* El path de la carpeta on es guarda l'output del programa "orbites.txt", passat com argument `--output_folder`. Per defecte, tenim la carpeta "outputs/ex6/orbites\_output".

Un cop hem executat el programa, obtenim el fixter "orbites.txt" i ja podem representar les òrbites halo a l'espai, l'objectiu de l'exercici. Per tal de generar-lo amb **gnuplot**, hem d'obrir un interpret de **gnuplot** en "ouputs/ex6" i executar les linies en el fitxer "outputs/ex6/comandes\_gnuplot.txt", o en "outputs/ex6/comandes\_gnuplot\_NP.txt", aquest últim sense els dos cossos primàris.

Podem veure els resultats obtinguts en la carpeta "outputs/ex6/retrats\_de\_fase".

A més, hem desenvolupat l'equació diferencial donat per RTBP amb **sagemath**, i amb el notebook del python, "rtbp\_diff\_eq.ipynb".
