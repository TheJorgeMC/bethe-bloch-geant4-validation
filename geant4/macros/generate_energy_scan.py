#!/usr/bin/env python3
# ============================================================================
# generate_energy_scan.py — genera una familia de macros run_<E>MeV.mac a
# partir de una plantilla, como alternativa a energy_scan.mac (/control/foreach)
# cuando se quiere un macro archivable por punto de la curva (trazabilidad
# dato <-> condiciones, diseño seccion 7.2), o un espaciado de energias
# distinto (p.ej. logaritmico) sin editar macros a mano.
#
# Uso:
#   python3 generate_energy_scan.py            # energias por defecto
#   python3 generate_energy_scan.py --log 20   # 20 puntos log entre E_MIN y E_MAX
#
# Luego:  for m in run_scan/run_*MeV.mac; do ./slab "$m"; done
# ============================================================================
import argparse
import math
import pathlib

# --- Condiciones comunes del barrido (editar aqui, no en las macros) --------
MATERIAL = "G4_WATER"
THICKNESS = "5 mm"
SIZE_XY = "10 cm"
N_LAYERS = 50
CUT = "1 mm"
CUT_TAG = "cut1mm"   # aparece en el nombre del archivo de datos
N_EVENTS = 100000
ENERGIES_MEV = [50, 70, 100, 150, 200, 250, 500, 1000]  # default del diseño
E_MIN, E_MAX = 30.0, 1000.0  # rango para --log
OUT_DIR = pathlib.Path("run_scan")

TEMPLATE = """# Macro generada por generate_energy_scan.py — NO editar a mano.
# Punto del barrido: {energy:g} MeV, corte {cut} ({cut_tag})
/control/verbose 0
/run/verbose 1
/event/verbose 0
/tracking/verbose 0

/absorber/material {material}
/absorber/thickness {thickness}
/absorber/sizeXY {size_xy}
/absorber/numberOfLayers {n_layers}

/physics/absorberCut {cut}

/run/initialize

/gun/particle proton
/gun/energy {energy:g} MeV

/analysis/setFileName dedx_{energy:g}MeV_{cut_tag}

/run/beamOn {n_events}
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--log", type=int, metavar="N", default=0,
        help="usar N puntos con espaciado logaritmico entre "
             f"{E_MIN:g} y {E_MAX:g} MeV en vez de la lista por defecto")
    args = parser.parse_args()

    if args.log > 1:
        step = (math.log10(E_MAX) - math.log10(E_MIN)) / (args.log - 1)
        energies = [round(10 ** (math.log10(E_MIN) + i * step), 1)
                    for i in range(args.log)]
    else:
        energies = ENERGIES_MEV

    OUT_DIR.mkdir(exist_ok=True)
    for e in energies:
        macro = TEMPLATE.format(
            energy=e, material=MATERIAL, thickness=THICKNESS,
            size_xy=SIZE_XY, n_layers=N_LAYERS, cut=CUT,
            cut_tag=CUT_TAG, n_events=N_EVENTS)
        path = OUT_DIR / f"run_{e:g}MeV.mac"
        path.write_text(macro)
        print(f"escrito {path}")

    print(f"\n{len(energies)} macros en {OUT_DIR}/ — ejecutar con:")
    print(f'  for m in {OUT_DIR}/run_*MeV.mac; do ./slab "$m"; done')


if __name__ == "__main__":
    main()
