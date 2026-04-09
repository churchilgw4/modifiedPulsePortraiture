# modifiedPulsePortraiture
We need two more packages than Tim Penucci's PulsePortraiture

1. pip install nautilus-sampler (https://nautilus-sampler.readthedocs.io/en/latest/guides/installation.html)
2. pip install pygeostat (https://www.ccgalberta.com/pygeostat/installation.html)

We need nautilus for bayesian mode of estimating toas and DM.
We need pygeostat for estimating skewness and kurtosis for phi and DM posteriors (will be helpful in outlier analysis).


# modifiedPulsePortraiture

A modified and extended version of Timothy T. Pennucci's [PulsePortraiture](https://github.com/pennucci/PulsePortraiture) package, adapted for wideband and multiband pulsar timing within the **Indian Pulsar Timing Array (InPTA)** pipeline. This package introduces Bayesian TOA/DM estimation via the **Maximum Likelihood with Analytic Noise (MLAN)** technique and adds support for simultaneous Band 3 + Band 5 (B35) observations at the uGMRT.

---

## Overview

`modifiedPulsePortraiture` simultaneously fits for:

- **Phase shifts / Times of Arrival (TOAs)**
- **Dispersion Measures (DMs)**
- **ν⁻⁴ delay parameters (GMs)**
- **Scattering timescales (τ) and scattering indices (α)**
- **Mean flux densities**

It extends the original PulsePortraiture with two fitting techniques:

| Technique | Description |
|-----------|-------------|
| **MLA** (Maximum Likelihood Analysis) | Pennucci et al. method — maximized over amplitude |
| **MLAN** (MLA + Nautilus sampler) | Abhimanyu et al. method — marginalized over amplitude and noise via Bayesian nested sampling |

Each technique supports two estimation modes:

| Mode | Flag | Description |
|------|------|-------------|
| **freq** | *(default)* | Frequentist optimization |
| **bayes** | `--bayes` | Bayesian posterior estimation |

The `_b35` variants of scripts are designed for simultaneous **Band 3 + Band 5** multiband fitting at the uGMRT.

---

## Repository Structure

```
modifiedPulsePortraiture/
├── mpplib.py                  # Core library (modified pplib)
├── mpplib_b35.py              # Core library for Band 3+5 simultaneous fitting
├── mpptoas.py                 # Main TOA/DM fitting script (MLA / MLAN)
├── mpptoas_b35.py             # Band 3+5 simultaneous TOA/DM fitting script
├── mpptoaslib.py              # TOA fitting support library
├── mpptoaslib_b35.py          # Band 3+5 TOA fitting support library
├── mpptoaslib_MLAN.py         # MLAN fitting support library
├── mpptoaslib_b35_MLAN.py     # Band 3+5 MLAN fitting support library
├── metafile_maker_cc_b35.py   # Utility: match Band 3 & Band 5 metafiles by MJD
├── mpp_gaussian_test.py       # Gaussian model testing utilities
├── mpp_telescope_codes.py     # Telescope/backend code definitions
└── setup.py                   # Package setup
```

---

## Dependencies

### From the original PulsePortraiture

Install [PulsePortraiture](https://github.com/pennucci/PulsePortraiture) and its dependencies first (PSRCHIVE, numpy, scipy, etc.).

### Additional packages required by modifiedPulsePortraiture

```bash
pip install nautilus-sampler
pip install pygeostat
```

- **[nautilus-sampler](https://nautilus-sampler.readthedocs.io/en/latest/guides/installation.html)** — Required for Bayesian (MLAN) mode via nested sampling.
- **[pygeostat](https://www.ccgalberta.com/pygeostat/installation.html)** — Required for estimating skewness and kurtosis of φ and DM posteriors (used in outlier analysis).

---

## Installation

```bash
git clone https://github.com/AvinashKumarPaladi/modifiedPulsePortraiture.git
cd modifiedPulsePortraiture
pip install .
```

Or install in editable/development mode:

```bash
pip install -e .
```

---

## Usage

### Step 1 — Prepare Metafiles

A **metafile** is a plain text file listing paths to your PSRCHIVE archive files, one per line. For Band 3+5 simultaneous fitting, use the provided utility to create time-matched metafiles:

```bash
python metafile_maker_cc_b35.py <band3.meta> <band5.meta>
```

This matches observations between Band 3 and Band 5 by MJD (within a ±0.1-day window) and outputs:
- `<band3.meta>.cc.b35`
- `<band5.meta>.cc.b35`

### Step 2 — Run TOA/DM Fitting

The four run modes below use pulsar **J2124-3358** as the example target. Replace paths with your own data accordingly.

---

#### B3 MLA — Band 3 only, frequentist

```bash
python3 pptoas.py \
  -d J2124-3358.DR2.run1/meta.files/J2124-3358.band3.200.meta.cc.b35 \
  -m J2124-3358.DR2.run1/template/J2124-3358.band3.200.spl \
  -o J2124-3358.DR2.run1/output/J2124-3358.B3.MLA_freq.tim
```

#### B3 MLA — Band 3 only, Bayesian

```bash
python3 pptoas.py \
  -d J2124-3358.DR2.run1/meta.files/J2124-3358.band3.200.meta.cc.b35 \
  -m J2124-3358.DR2.run1/template/J2124-3358.band3.200.spl \
  -o J2124-3358.DR2.run1/output/J2124-3358.B3.MLA_bayes.tim \
  --bayes
```

---

#### B3 MLAN — Band 3 only, frequentist

```bash
python3 pptoas.py \
  -d J2124-3358.DR2.run1/meta.files/J2124-3358.band3.200.meta.cc.b35 \
  -m J2124-3358.DR2.run1/template/J2124-3358.band3.200.spl \
  -o J2124-3358.DR2.run1/output/J2124-3358.B3.MLAN_freq.tim \
  --mlan
```

#### B3 MLAN — Band 3 only, Bayesian

```bash
python3 pptoas.py \
  -d J2124-3358.DR2.run1/meta.files/J2124-3358.band3.200.meta.cc.b35 \
  -m J2124-3358.DR2.run1/template/J2124-3358.band3.200.spl \
  -o J2124-3358.DR2.run1/output/J2124-3358.B3.MLAN_bayes.tim \
  --bayes --mlan
```

---

#### B35 MLA — Band 3 + Band 5 simultaneous, frequentist

```bash
python3 pptoas_b35.py \
  -a J2124-3358.DR2.run1/meta.files/J2124-3358.band3.200.meta.cc.b35 \
  -b J2124-3358.DR2.run1/meta.files/J2124-3358.band5.200.meta.cc.b35 \
  -c J2124-3358.DR2.run1/template/J2124-3358.band3.200.spl \
  -d J2124-3358.DR2.run1/template/J2124-3358.band5.200.spl \
  -o J2124-3358.DR2.run1/output/J2124-3358.B35.MLA_freq.tim
```

#### B35 MLA — Band 3 + Band 5 simultaneous, Bayesian

```bash
python3 pptoas_b35.py \
  -a J2124-3358.DR2.run1/meta.files/J2124-3358.band3.200.meta.cc.b35 \
  -b J2124-3358.DR2.run1/meta.files/J2124-3358.band5.200.meta.cc.b35 \
  -c J2124-3358.DR2.run1/template/J2124-3358.band3.200.spl \
  -d J2124-3358.DR2.run1/template/J2124-3358.band5.200.spl \
  -o J2124-3358.DR2.run1/output/J2124-3358.B35.MLA_bayes.tim \
  --bayes
```

---

#### B35 MLAN — Band 3 + Band 5 simultaneous, frequentist

```bash
python3 pptoas_b35.py \
  -a J2124-3358.DR2.run1/meta.files/J2124-3358.band3.200.meta.cc.b35 \
  -b J2124-3358.DR2.run1/meta.files/J2124-3358.band5.200.meta.cc.b35 \
  -c J2124-3358.DR2.run1/template/J2124-3358.band3.200.spl \
  -d J2124-3358.DR2.run1/template/J2124-3358.band5.200.spl \
  -o J2124-3358.DR2.run1/output/J2124-3358.B35.MLAN_freq.tim \
  --mlan
```

#### B35 MLAN — Band 3 + Band 5 simultaneous, Bayesian

```bash
python3 pptoas_b35.py \
  -a J2124-3358.DR2.run1/meta.files/J2124-3358.band3.200.meta.cc.b35 \
  -b J2124-3358.DR2.run1/meta.files/J2124-3358.band5.200.meta.cc.b35 \
  -c J2124-3358.DR2.run1/template/J2124-3358.band3.200.spl \
  -d J2124-3358.DR2.run1/template/J2124-3358.band5.200.spl \
  -o J2124-3358.DR2.run1/output/J2124-3358.B35.MLAN_bayes.tim \
  --bayes --mlan
```

---

### Command-Line Arguments Summary

| Script | Flag | Description |
|--------|------|-------------|
| `pptoas.py` | `-d` | Metafile listing archive paths |
| `pptoas.py` | `-m` | Template/model file (`.spl` spline or `.gmodel` Gaussian) |
| `pptoas.py` | `-o` | Output TOA file (`.tim`) |
| `pptoas.py` | `--bayes` | Enable Bayesian posterior estimation |
| `pptoas.py` | `--mlan` | Use MLAN (marginalized over amplitude & noise) |
| `pptoas_b35.py` | `-a` | Band 3 metafile |
| `pptoas_b35.py` | `-b` | Band 5 metafile |
| `pptoas_b35.py` | `-c` | Band 3 template |
| `pptoas_b35.py` | `-d` | Band 5 template |
| `pptoas_b35.py` | `-o` | Output TOA file (`.tim`) |
| `pptoas_b35.py` | `--bayes` | Enable Bayesian mode |
| `pptoas_b35.py` | `--mlan` | Use MLAN technique |

---

## Python API Usage

```python
from mpptoas import GetTOAs

# Initialize with a metafile and a model
gt = GetTOAs("mydata.meta", "mymodel.spl")

# Run wideband TOA/DM fitting
gt.get_TOAs(
    fit_DM=True,
    print_phase=True,
    print_flux=False,
    bayes=False,    # set True for Bayesian mode
    quiet=False
)

# Write TOAs to file
for toa in gt.TOA_list:
    toa.write_TOA(outfile="output.tim")
```

---

## Output

The scripts produce **IPTA-formatted TOA files** (`.tim`) with extended InPTA-specific flags:

| Flag | Description |
|------|-------------|
| `phi`, `phi_err` | Fitted phase and uncertainty |
| `DM`, `DM_err` | Fitted DM and uncertainty |
| `bandno` | uGMRT band identifier (3, 4, or 5) |
| `group`, `sys` | InPTA observation group and system identifier |
| `cycle` | `pre34`/`post34` or `pre36`/`post36` based on MJD |
| `snr` | Signal-to-noise ratio |
| `gof` | Goodness of fit (reduced χ²) |
| `technique` | Method used: `MLA_freq`, `MLA_bayes`, `MLAN_freq`, or `MLAN_bayes` |
| `log_z` | Log evidence (Bayesian mode only) |

---

## Results and Key Findings

Results are validated on pulsar **J2124-3358** from InPTA DR2 across all four technique combinations, comparing frequentist vs. Bayesian modes and B3 vs. B35 coverage.

**Note on plots:** PHI (phase) is plotted rather than raw TOA values, since the TOA range spans several years while the variations of interest are at the sub-microsecond level. Colour coding reflects SNR — yellow = high SNR, blue = low SNR.

**Key takeaways:**

1. DM and PHI from the frequentist and Bayesian modes are **mutually consistent**.
2. DM is **consistent** across MLA vs. MLAN and B3 vs. B3+5.
3. Uncertainties on DM and PHI (errDM, errPHI) are **consistent** between frequentist and Bayesian modes.
4. errDM and errPHI for **B3+5 are smaller than B3 alone**, as expected from the wider frequency leverage arm.
5. errDM and errPHI for **MLAN are larger than MLA** — consistent with results from Abhimanyu et al.
6. PHI values from MLA and MLAN differ because the **reference frequency is different** between the two techniques; points therefore appear scattered rather than along the diagonal.

---

## InPTA-Specific Notes

- Observatory/backend flags are automatically set for the **uGMRT** (upgraded Giant Metrewave Radio Telescope).
- Band assignments: Band 3 (300–500 MHz), Band 4 (525–1000 MHz), Band 5 (1260–1460 MHz).
- Observations are tagged `pre36`/`post36` based on MJD 58600, and `pre34`/`post34` based on MJD 58230.
- `metafile_maker_cc_b35.py` ensures time-matched Band 3 + Band 5 archive pairs for B35 fitting.

---

## Credits & Acknowledgements

- Original **PulsePortraiture** by [Timothy T. Pennucci](https://github.com/pennucci/PulsePortraiture) (with contributions by Scott M. Ransom and Paul B. Demorest).
- MLA technique: Pennucci et al.
- MLAN technique: Abhimanyu et al.
- Modified and extended for InPTA wideband timing by **Avinash Kumar Paladi** (avinashkumarpaladi@gmail.com).
- Uses [nautilus-sampler](https://nautilus-sampler.readthedocs.io) for nested sampling in Bayesian mode.
- Uses [pygeostat](https://www.ccgalberta.com/pygeostat/) for posterior shape statistics.

---

## Contact

For questions or issues, please open a GitHub issue or contact:  
**Avinash Kumar Paladi** — avinashkumarpaladi@gmail.com
