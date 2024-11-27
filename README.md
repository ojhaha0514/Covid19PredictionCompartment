
# SEIQRDVP Model for COVID-19 Impact Analysis

![MATLAB](https://img.shields.io/badge/MATLAB-R2021a-orange)
![COVID-19](https://img.shields.io/badge/COVID--19-SEIQRDVP%20Model-green)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)

This project implements the SEIQRDVP and SEIQRDV3P mathematical models to study the impact of COVID-19 variants, particularly the Omicron variant, in South Korea. The models are based on the paper titled *Mathematical Modeling of the Impact of Omicron Variant on the COVID-19 Situation in South Korea*.

## Overview

The project uses a modified SEIR model that incorporates vaccination levels and variant effects to predict the trajectory of COVID-19 cases. Key features of the model:
- **Incorporation of Variants**: Considers the weighted transmission rates of Delta and Omicron variants.
- **Vaccination Levels**: Simulates the effects of one-dose, two-dose, and booster vaccinations.
- **Comparative Performance**: Evaluates models like SEIR, SEIQR, SEIQRDVP, and SEIQRDV3P using RMSE metrics.

## Files in the Repository

- **`paper_seir.m`**: Main MATLAB script implementing the SEIQRDVP model.
- **`owid-covid-data.csv`**: Input data file containing COVID-19 case and vaccination statistics.

## Prerequisites

- MATLAB R2021a or later
- Required MATLAB toolboxes: Optimization Toolbox

## Usage

1. **Setup**:
   Place `paper_seir.m` and `owid-covid-data.csv` in the same directory.

2. **Run the Model**:
   Execute the main script in MATLAB:
   ```matlab
   paper_seir
   ```

3. **Output**:
   - Predicted daily cases over the test period.
   - RMSE metrics for different configurations (e.g., different vaccination levels or variant transmission rates).

## Results

Key findings from the implementation:
- The SEIQRDV3P model demonstrated a **27.4% improvement in RMSE** compared to the SEIQRDVP model.
- The Omicron variant significantly increased case predictions when compared to Delta variant-only scenarios.

## References

- Jooha Oh, Catherine Apio, and Taesung Park. *Mathematical Modeling of the Impact of Omicron Variant on the COVID-19 Situation in South Korea*. Genomics & Informatics, 2022. [doi.org/10.5808/gi.22025]{https://doi.org/10.5808/gi.22025}

## Authors

- **Jooha Oh**<sup>1*</sup>, **Catherine Apio**<sup>2</sup>
- Corresponding Author: [**Taesung Park**](mailto:tspark@stats.snu.ac.kr)<sup>1</sup>
- <sup>1</sup>Department of Statistics, Seoul National University, Seoul, Republic of Korea
- <sup>2</sup>Interdisciplinary Program in Bioinformatics, Seoul National University, Seoul, Republic of Korea



---

## Acknowledgments

- Publicly available data collected by Our World In Data were used in this study and can be found at following repository: [https://covid.ourworldindata.org]{https://covid.ourworldindata.org}.

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.


---
