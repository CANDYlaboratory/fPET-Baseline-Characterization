# fPET-Baseline-Characterization

This github repository is a companion to the paper "**On the analysis of functional PET (fPET)-FDG: baseline mischaracterization can introduce artifactual metabolic (de)activations**" (currently on [biorxiv](https://www.biorxiv.org/content/10.1101/2024.10.17.618550v1)).

The repository includes several functions used to generate the results presented in the paper, including `fPETregress`, `fPETsimulateTAC`, and `fPETgetBaseline`, along with sample data to test these functions with. While we are sharing this code primarily to allow others to replicate our results, we hope these functions may prove more generally useful to researchers interested in the analysis of fPET-FDG data.

To get started, try running the script `exampleScript.m` to see the functions used in context.

If you have any questions or comments regarding the paper or the code, you can email the corresponding authors: Sean Coursey (coursey.s@northeastern.edu) and Jingyuan Chen (JECHEN@mgh.harvard.edu).

-----

The data (`resting_CI_sample_subjects.mat`) is a randomly-chosen subset of the resting-state constant-infusion subjects from "Simultaneous BOLD-fMRI and constant infusion FDG-PET data of the resting human brain" [(Jamadar et al., _Sci Data_, 2020)](https://doi.org/10.1038/s41597-020-00699-5) whose data have been processed into 100-ROI parcellations [(Schaefer et al., _Cerebral Cortex_, 2018)](https://doi.org/10.1093/cercor/bhx179). The data was reconstructed with 16-second bins.

The provided constant infusion AIF (`cont_infusion_group_AIF.mat`) was derived from the above-mentioned resting-state dataset from Jamadar et al., and the bolus plus constant infusion AIF was derived from a dataset acquired at MGH. See our paper for details on the data and AIFs. Both are scale-less and normalized.
