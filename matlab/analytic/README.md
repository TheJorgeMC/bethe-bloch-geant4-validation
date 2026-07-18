# `matlab/analytic/` — Implementación Analítica de Bethe-Bloch (MATLAB)

Réplica en MATLAB de la implementación de Python (`python/analytic/analytic_solution.ipynb`),
verificada numéricamente contra ella (ver sección "Validación cruzada" abajo). Implementa
la ecuación de Bethe-Bloch relativista **sin correcciones** de capa, Barkas, Bloch, ni
efecto de densidad — ver `Derivacion_Completa_Bethe_Bloch.tex`, sección "Alcance del
Proyecto y Trabajo Futuro", para la justificación física completa.

No requiere ningún toolbox de MATLAB (ni Symbolic Math Toolbox, ni Statistics Toolbox).
Todo el código usa únicamente funciones base.

## Requisitos

- MATLAB (cualquier versión razonablemente reciente; se probó la lógica numérica en
  GNU Octave 8.4 como sustituto durante el desarrollo, pero **el uso previsto es MATLAB real**,
  ya que los scripts de salida usan `table`/`writetable`, que Octave no soporta).
- Ningún toolbox adicional.

## Cómo correrlo

1. Agrega esta carpeta al path de MATLAB (clic derecho → *Add to Path* → *Selected Folder*,
   o navega a `matlab/analytic/` en el explorador de archivos de MATLAB).
2. Corre, en este orden, desde la Command Window o el Editor (F5):

   ```matlab
   testBetheBloch()             % 7 tests básicos — correr primero
   testBetheBlochUncertainty()  % 7 tests de propagación de incertidumbre
   sanityCheck                  % chequeo puntual (protón 100 MeV en agua)
   generateCurves                % genera 4 CSV + 1 PNG (sin incertidumbre)
   generateCurvesWithUncertainty % genera 4 CSV + 1 PNG de 2 paneles (con incertidumbre)
   ```

   Los dos primeros son **funciones** (llevan paréntesis); los tres últimos son **scripts**
   (se corren solo con el nombre, sin paréntesis, o dando *Run* con el archivo abierto).

3. Los CSV y PNG se generan en `matlab/analytic/output/` (se crea automáticamente si no existe).

## Estructura de archivos

### Datos (constantes y materiales)

| Archivo | Contenido |
|---|---|
| `betheBlochConstants.m` | Constantes físicas: mₑc², Mₚc², K (PDG), factor eV→MeV, y las incertidumbres ΔK, Δmₑ (CODATA, despreciables) |
| `betheBlochWaterComposition.m` | Constituyentes del agua para la regla de Bragg: A(H), ΔA(H), A(O), ΔA(O) (IUPAC) y fracciones de masa (NIST STAR) |
| `betheBlochMaterials.m` | Struct con los 4 materiales (`water`, `aluminium`, `copper`, `lead`), cada uno con `Z_over_A`, `I_eV`, `Delta_I_eV`, `rho_g_cm3`, `Delta_rho_g_cm3`, `A_g_mol`, `Delta_A_g_mol`, `matno` — todo explícito, nada tratado como "despreciable" implícitamente |

### Física core

| Archivo | Contenido |
|---|---|
| `betheBlochBetaGamma.m` | β, γ relativistas a partir de la energía cinética |
| `betheBlochMassStoppingPower.m` | S(T), poder de frenado másico [MeV·cm²/g] — la ecuación de Bethe-Bloch sin correcciones |
| `betheBlochLinearStoppingPower.m` | S(T)·ρ, poder de frenado lineal [MeV/cm] |
| `betheBlochValidityCheck.m` | β, γ, βγ y bandera de validez (βγ ≥ 0.1) |

### Incertidumbre

| Archivo | Contenido |
|---|---|
| `betheBlochPartialDerivatives.m` | Las 4 derivadas parciales de S respecto a K, Z/A, mₑ, I — codificadas directamente (verificadas antes con SymPy en Python; no se usa Symbolic Math Toolbox aquí) |
| `betheBlochDeltaZA.m` | Propaga Δ(A) → Δ(Z/A), con la regla de Bragg completa para el agua |
| `betheBlochPropagateMassUncertainty.m` | ΔS por suma en cuadratura de los 4 términos; puede devolver el desglose término por término |
| `betheBlochPropagateLinearUncertainty.m` | ΔS_linear, agregando el término de Δρ vía regla del producto |

### Validación

| Archivo | Contenido |
|---|---|
| `testBetheBloch.m` | 7 tests básicos (límites de β/γ, positividad, monotonía, comparación NIST, escalado con ρ y z², bandera de validez) |
| `testBetheBlochUncertainty.m` | 7 tests de incertidumbre (K/mₑ despreciables, I domina, rango razonable, crecimiento a baja energía, desigualdad triangular, linealidad, y validación cruzada contra Kumazaki et al. 2007 — ver más abajo) |
| `sanityCheck.m` | Chequeo de un punto (100 MeV, agua) contra NIST PSTAR y contra el resultado exacto de Python |

### Generación de salida

| Archivo | Contenido |
|---|---|
| `betheBlochEnergyGrid.m` | Rejilla de 54 puntos (44 en 3–300 MeV clínico + 10 de extensión hasta 1000 MeV), idéntica a la de Python |
| `generateCurves.m` | CSV + PNG sin incertidumbre |
| `generateCurvesWithUncertainty.m` | CSV + PNG con incertidumbre (banda ±1σ y % relativo), más resumen del presupuesto de varianza en consola |

## Sobre el test de Kumazaki et al. (2007)

`test_water_5eV_shift_matches_literature_order_of_magnitude` en `testBetheBlochUncertainty.m`
valida el código contra un resultado publicado independiente, no solo contra su propia
consistencia interna:

> Kumazaki, Y., Akagi, T., Yanou, T., Suga, D., Hishikawa, Y., & Teshima, T. (2007).
> *Determination of the mean excitation energy of water from proton beam ranges.*
> Radiation Measurements, 42(10), 1683–1691. DOI: 10.1016/j.radmeas.2007.10.005

Los autores reportan que un cambio de I=75→80 eV en agua produce una diferencia de
0.8–1.2% en el stopping power entre 10–250 MeV. El test reproduce exactamente ese
experimento (construye agua con I=75 y con I=80, calcula S a 100 MeV con ambas, y verifica
que la diferencia caiga en ese mismo rango) — así que además de probar que el código no
tiene bugs internos, confirma que la física coincide con una fuente externa citable.

## Validación cruzada con Python

Todos los valores numéricos (materiales, S(T), ΔS(T), y el desglose de varianza K/ZA/mₑ/I)
se verificaron para coincidir con `python/analytic/analytic_solution.ipynb` dentro de error
de punto flotante (~10⁻⁶ relativo), en los 4 materiales y en múltiples puntos de energía
(incluyendo el resumen completo del presupuesto de incertidumbre a 100 MeV). Ver el
historial de chat del proyecto para el detalle de esa comparación si hace falta reproducirla.

## Alcance (recordatorio)

Esta implementación **no incluye** corrección de capa, Barkas, Bloch, ni efecto de densidad.
Se espera que subestime NIST PSTAR y Geant4, de forma creciente a baja energía y en
materiales de Z alto (plomo) — esto es un resultado documentado del proyecto, no un error.
Ver `Derivacion_Completa_Bethe_Bloch.tex` para la justificación física completa y el
razonamiento de por qué esas correcciones se dejan para una tesis dirigida futura.
