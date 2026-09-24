# Tested analysis environment for release preparation

Environment recorded during the 24 September 2026 release preflight on Linux:

| Component | Version |
| --- | --- |
| Python | 3.12.12 |
| NumPy | 2.4.2 |
| pandas | 3.0.1 |
| Matplotlib | 3.10.8 |
| Cooler | 0.10.4 |
| h5py | 3.16.0 |
| Biopython | 1.86 |
| PyArrow | 23.0.0 |

The four final-figure code bases depend directly on NumPy, pandas and Matplotlib; Figure 3 additionally uses Cooler/h5py. Other repository scripts may use Biopython. PyArrow was present in the tested environment but is not required directly by the final-figure scripts.

This environment record documents the tested working installation; it is not intended to claim that these exact versions are the only compatible versions.
