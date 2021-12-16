# Exericici 1

L'objectiu d'aquest exercici i apartat és solucionar una EDO i alhora, calcular les equacions variacionals del flux.

Els programes associats són "scripts/s2\_ex1.py" i "scripts/utils/include_diffs.py".

Desglossem l'exercici en dues tàsques:

* Construir la funció "camp" que de les equacions diferencials: $$\vec{x}' = f(t, \vec{x}) \; , \; \vec{x}(t_0) = \vec{x}_0$$ $$A' = D_x f(t, \vec{x}) \cdot A \; , \; A(t_0) = I_n$$ donat $$\vec{y} = (\vec{x}, a_{1,1} , \dots, a_{1,n}, \dots, a_{n,1} , \dots, a_{n,n}) \in {\mathbb{R}}^{n^2 + n}$$ retorni $$ \big( f(t, \vec{x}), D_x f(t, \vec{x}) \cdot A \big)$$ on $a_{i,j}$ són els elements de la matriu $A$. Aquest procediment el realitzem de manera general per a tota funció dins "scripts/utils/include_diffs.py".

* Comprovar de manera pràctica que la matriu de variacions obtinguda coincideix amb l'aproximació: $$D_{\vec{x}}\phi(t ; t_0, \vec{x}_0) \approx \Big[\displaystyle\frac{\phi(t;t_0, \vec{x}_0 + \delta e_1) - \phi(t;t_0, \vec{x}_0)}{ \delta}, \dots, \displaystyle\frac{\phi(t;t_0, \vec{x}_0 + \delta e_n) - \phi(t;t_0, \vec{x}_0)}{ \delta}\Big]$$ Aquest procediment es realitza amb l'ajuda de la funció "camp" en "scripts/s2\_ex1.py" sobre el sistema diferencial: $$(x', y') = \big(\alpha (1 - r^2) x - y, x + \alpha (1 - r^2) y\big)$$ on $r^2 = x^2 + y^2$ i prenent $\alpha = 0.5$.
