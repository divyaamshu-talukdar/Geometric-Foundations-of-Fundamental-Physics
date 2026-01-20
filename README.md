# Geometric Foundations of Fundamental Physics

This repository contains the computational implementation of the geometric framework presented in:

**"Geometric Foundations of Fundamental Physics: A Topological Approach to Standard Model Parameters"**  
*Divyaamshu Talukdar*  
European Physical Journal C (2026)

## Files

### 1. KelvinCell.py
**Purpose**: Quantum Field Theory on Kelvin Cell Lattice  
**Key Results**: 
- Demonstrates convergence of QFT on Kelvin lattice to continuum limit
- Convergence exponent: -3.47 ± 0.09 (superior to standard lattices)
- Verifies emergence of continuum physics from discrete geometry

### 2. GeometricFramework.py  
**Purpose**: Validation of Geometric Predictions  
**Key Results**:
- α⁻¹ = 137.035999 (4.78 ppm error)
- m_p/m_e = 1836.152673 (exact with correction)
- sin²θ_W = 0.230940 (0.14% error)
- α_s(M_Z) = 0.118174 (0.15% error)
- M_DM = 203.8 MeV (prediction)

### 3. VertexCorrectionDerivation.py
**Purpose**: Derivation of Geometric Correction Factor  
**Key Results**:
- f_corr = 0.9998835 from vertex curvature
- Essential for exact proton mass ratio
- Multi-step geometric derivation

## Requirements
- Python 3.8+
- NumPy 1.21+
- SciPy 1.7+
- Matplotlib 3.4+

## Installation
```bash
pip install numpy scipy matplotlib
