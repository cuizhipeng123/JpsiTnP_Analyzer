# JpsiTnP_Analyzer

This is a CMSSW-based analyzer module for measuring single muon efficiencies using the Tag-and-Probe method in J/ψ → μ⁺μ⁻ decays.  
Originally based on `MuonAnalyzer` for HEP data analysis.

## Structure

- `plugins/`: Main analyzer plugins for AOD and miniAOD.
- `python/`: Configuration fragments for reconstruction sequences.
- `test/`: CRAB configurations for TnP production.
- `data/`: Sample database JSONs for various data eras.
- `scripts/`: Utility tools for submission and analysis.

## Usage

This repository is intended to be used within a CMSSW environment (e.g., `CMSSW_13_0_13`).  
Typical steps:

```bash
cmsrel CMSSW_13_0_13
cd CMSSW_13_0_13/src
cmsenv
git clone git@github.com:cuizhipeng123/JpsiTnP_Analyzer.git
scram b -j 8
