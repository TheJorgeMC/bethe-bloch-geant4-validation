# Metodología de NIST PSTAR — nota para el reporte

**Fuente:** NIST Standard Reference Database 124, PSTAR
(`https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html`)

PSTAR calcula el poder de frenado y alcance de protones combinando dos
regímenes de energía con un tratamiento distinto en cada uno:

1. **Baja energía (T ≲ 1 MeV aprox., justo la región donde la teoría de
   Bethe pierde validez):** el poder de frenado electrónico se obtiene por
   **interpolación mediante spline cúbico** sobre datos experimentales
   recopilados y evaluados por ICRU (Report 49) y compilaciones asociadas
   (Andersen & Ziegler, Ziegler). En este régimen la teoría de Born-Bethe
   no es confiable porque la velocidad del protón es comparable a las
   velocidades orbitales de los electrones atómicos, por lo que se recurre
   a ajustes semiempíricos a los datos medidos.

2. **Alta energía (por encima del rango cubierto de forma confiable por
   los ajustes empíricos):** se usa la **fórmula de Bethe** (la misma
   forma funcional derivada de la teoría de Born a primer orden, con el
   logaritmo relativista y el prefactor 1/β²), incluyendo las
   correcciones estándar del formalismo ICRU 37 / ICRU 90:
   - **Corrección de capa (shell correction)** — relevante en el régimen
     donde la aproximación de Born empieza a fallar por la estructura
     atómica del absorbente (importante para Pb en casi todo el rango
     clínico, según ya se estableció en la derivación previa).
   - **Corrección de Barkas-Andersen** (término en z³, asimetría de carga).
   - **Corrección de Bloch** (interpolación entre el límite de Born y el
     límite clásico de Bohr, gobernada por el parámetro de Sommerfeld
     η = zα/β).

   El poder de frenado nuclear (STOP(n)) se calcula por separado con un
   modelo de colisiones elásticas núcleo-núcleo (potencial apantallado
   universal), dominante solo a energías muy bajas (keV) y despreciable
   en el rango clínico de 3–300 MeV.

3. **Empalme (matching):** las dos regiones se conectan de forma continua
   en energía y derivada en el punto de cruce, de modo que la tabla final
   entregada por PSTAR es una función suave sin discontinuidades, aun
   cuando internamente combina un ajuste spline (bajo T) con una fórmula
   analítica corregida (alto T).

**Relevancia para este proyecto:** en el rango clínico de protonterapia
(3–300 MeV, βγ ≳ 0.08), el protón ya está firmemente en el régimen donde
domina la fórmula de Bethe con correcciones ICRU 90, por lo que PSTAR
constituye una referencia adecuada y consistente con la formulación
teórica (Bethe-Bloch relativista + shell + Barkas/Bloch) que se está
validando en este trabajo — a diferencia del régimen de spline
experimental de bajísima energía, que queda fuera del rango de interés.

**Cita sugerida:**
Berger, M.J., Coursey, J.S., Zucker, M.A., and Chang, J. (2005),
*Stopping-Power and Range Tables for Electrons, Protons, and Helium Ions*
(NIST Standard Reference Database 124), National Institute of Standards
and Technology, Gaithersburg, MD. https://dx.doi.org/10.18434/T4NC7P
