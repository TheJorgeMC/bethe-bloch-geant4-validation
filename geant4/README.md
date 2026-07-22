# slab — Validación de Bethe-Bloch con protones en Geant4

Aplicación Geant4 (11.x) que simula protones monoenergéticos atravesando un
slab absorbedor homogéneo, para comparar el poder de frenado simulado contra
la ecuación de Bethe-Bloch analítica (Fase 5 del proyecto). Esqueleto basado
en el estilo de `TestEm0`.

## Compilación

```bash
mkdir build && cd build
cmake -DGeant4_DIR=/ruta/a/geant4/lib/cmake/Geant4 ..
make -j$(nproc)
```

## Ejecución

```bash
./slab                          # interactivo (ejecuta macros/vis.mac)
./slab macros/run.mac           # batch, EM option4 (default)
./slab macros/run.mac 3         # batch, EM option3 (chequeo de sensibilidad)
./slab macros/energy_scan.mac   # barrido completo de energías
```

Alternativa al barrido con `foreach` (una macro archivable por punto):

```bash
python3 macros/generate_energy_scan.py
for m in run_scan/run_*MeV.mac; do ./slab "$m"; done
```

## Parámetros por macro

| Comando | Default | Significado |
|---|---|---|
| `/absorber/material <NIST>` | `G4_WATER` | material del slab |
| `/absorber/thickness <v> <u>` | 5 mm | espesor en z |
| `/absorber/sizeXY <v> <u>` | 10 cm | lado de la sección transversal |
| `/absorber/numberOfLayers <N>` | 1 | N>1 segmenta el slab (curva de Bragg) |
| `/physics/absorberCut <v> <u>` | 1 mm | corte de producción en `AbsorberRegion` |
| `/gun/energy <v> <u>` | 150 MeV | energía cinética del protón |
| `/analysis/setFileName <n>` | `slab_out` | archivo de salida (tras `/run/initialize`) |

## Salida

**Consola** (fin de run): media y rms (straggling) del depósito de energía,
balance por canales, dE/dx simulado "local" (~restringido) y "total"
(E_inc − E_salida), y los valores de referencia de `G4EmCalculator`:
dE/dx restringido (con el corte activo), dE/dx **no restringido** (el
comparable con Bethe-Bloch de libro de texto) y rango CSDA.

**Ntuple `slab`** (una fila por evento, MeV / mm):

| Columna | Significado |
|---|---|
| `E_incidente` | energía cinética inicial del primario |
| `E_deposit_primario` | depósito del primario en el slab (≈ dE/dx restringido × camino) |
| `E_deposit_secundarios` | depósito de rayos delta y su progenie dentro del slab |
| `E_transportada_fuera` | energía cinética que los secundarios sacan del slab |
| `E_salida_primario` | energía del primario al salir (−1 si no salió) |
| `longitud_traza` | camino real del primario dentro del slab |

Chequeo clave del diseño (sección 5.3): por evento,
`E_incidente ≈ E_deposit_primario + E_deposit_secundarios +
E_transportada_fuera + E_salida_primario`.

**Histograma `DepthEdep`**: energía depositada vs profundidad (curva de
Bragg si `numberOfLayers > 1`), un bin por capa.

### Formato de archivo

Por defecto se escribe **CSV**. En multi-hilo el ntuple CSV sale en archivos
por hilo (`*_nt_slab_t0.csv`, `_t1`, ...): concatenarlos en el post-análisis,
p. ej. `pandas.concat([pd.read_csv(f, comment='#', header=None) for f in files])`.
Para un único archivo fusionado, usar extensión `.root` en
`/analysis/setFileName` (p. ej. `dedx_150MeV.root`) y leerlo con `uproot`.

## Notas físicas

- El dE/dx continuo del primario es el **restringido** con `T_cut` fijado por
  `/physics/absorberCut` (traducido a energía material a material); depende de
  esa elección y no es una constante física. El **no restringido** se recupera
  sumando la energía de los deltas explícitos (columnas del ntuple) o bajando
  el corte hasta converger (p. ej. 0.01 mm).
- Sin física hadrónica (decisión de diseño): no hay reacciones nucleares
  inelásticas; los primarios "perdidos" solo pueden frenarse dentro del slab.
- Documentar siempre versión de Geant4, physics list (option3/4) y corte de
  producción junto a cada resultado.
