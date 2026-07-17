# Validación de la Ecuación de Bethe para Protones

🇬🇧 [English](README.md) | 🇪🇸 Español

Derivación teórica, implementación analítica y validación mediante simulación Monte Carlo de la fórmula de pérdida de energía (stopping power) de Bethe-Bloch para protones, contrastada con datos de referencia del NIST.

🚧 **Estado:** en desarrollo activo (etapa temprana).

## Resumen

Este proyecto valida la ecuación de Bethe-Bloch para protones a través de tres enfoques independientes:

1. **Teoría** — derivación completa de la ecuación de Bethe-Bloch (relativista, con correcciones de capa, Barkas, Bloch y densidad), implementada analíticamente en **Python y MATLAB** (con validación cruzada entre ambos).
2. **Simulación Monte Carlo** — simulación del transporte de partículas en **Geant4**.
3. **Datos de referencia** — tablas de poder de frenado del **NIST PSTAR**.

Los tres resultados se comparan entre sí en un rango de energías de **3–300 MeV**, correspondiente al rango clínicamente relevante en protonterapia (Paganetti, 2012), para cuatro materiales absorbentes que abarcan un amplio rango de número atómico: **agua, aluminio, cobre y plomo**. Se presta especial atención a las regiones de energía donde se espera que la aproximación de Bethe-Bloch comience a fallar (bajas energías, por debajo de βγ≈0.1, es decir T≲4.7 MeV para protones), lo cual es directamente relevante para aplicaciones de protonterapia cerca del pico de Bragg.

## Estructura del repositorio

```
theory/             Derivación teórica completa (LaTeX/Markdown)
python/
  analytic/         Implementación en Python de la ecuación de Bethe-Bloch
  analysis/         Comparación, gráficos y análisis de error (Python)
matlab/
  analytic/         Implementación en MATLAB de la ecuación de Bethe-Bloch
  analysis/         Comparación, gráficos y análisis de error (MATLAB)
geant4/             Código fuente de la simulación en Geant4
nist_data/          Tablas de referencia NIST PSTAR (convertidas a CSV)
figures/            Gráficos generados
report/             Reporte técnico final (PDF)
```

## Metodología

| Fuente | Método |
|---|---|
| Teoría | Fórmula analítica de Bethe-Bloch (Python + MATLAB, con validación cruzada) |
| Simulación | Monte Carlo en Geant4, haces de protones monoenergéticos a través de láminas delgadas de material absorbente |
| Referencia | Tablas de poder de frenado NIST PSTAR |

## Requisitos

- **Geant4** (se recomienda la versión LTS)
- **Python 3.x** — numpy, scipy, pandas, matplotlib
- **MATLAB** — instalación base, sin toolboxes especializados

Las instrucciones detalladas de configuración y ejecución se agregarán a medida que se complete cada módulo.

## Motivación

El poder de frenado de protones es central en la planeación de tratamientos de protonterapia, donde el perfil de dosis en profundidad (y la ubicación del pico de Bragg) depende directamente de qué tan bien se pueda predecir la pérdida de energía. El rango de energías de 3–300 MeV estudiado aquí corresponde al rango clínico utilizado en protonterapia (Paganetti, 2012), y fue elegido deliberadamente para abarcar tanto la región donde la teoría estándar de Bethe-Bloch es válida (βγ≳0.1) como la región cercana a su límite inferior de validez, donde las correcciones de capa y de Barkas se vuelven significativas.

## Referencias

- Paganetti, H. (2012). *Proton Therapy Physics*. CRC Press.
- Particle Data Group (Navas, S. et al.) (2024). "Passage of Particles Through Matter," *Prog. Theor. Exp. Phys.*
- Salvat, F. (2022). "Bethe stopping-power formula and its corrections," *Phys. Rev. A* 106, 032809.

## Licencia

Este proyecto está bajo la Licencia MIT — ver [LICENSE](LICENSE) para más detalles.

## Autor

[TheJorgeMC](https://github.com/TheJorgeMC)
