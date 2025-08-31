# Pupillometric Production Effect Data

This repository contains the behavioral and pupillometric data associated with:

**Fawcett, Jonathan M., Roberts, Brady R. T., Willoughby, Hannah, Tiller, Jenny, Hourihan, Kathleen L., & MacLeod, Colin M.**  
*The Pupillometric Production Effect: Evidence for Enhanced Processing Preceding, During, and Following Production*  
_Preprint available at [SSRN](https://ssrn.com/abstract=5210001)_  
(Currently under revision at *Cognition*)

---

## Study Overview

The **production effect** refers to the well-established finding that words read aloud are remembered better than words read silently. While this effect has often been attributed to distinctive encoding via additional sensorimotor features, less work has examined attentional and motivational processes that differ between aloud and silent study.

Across **four experiments**, we used **pupillometry** to track attention and processing load during study. Participants read words aloud, silently, or (in some experiments) produced a control response (“check”). Using high-resolution eye-tracking, we measured changes in pupil dilation throughout study trials and linked these to subsequent recognition performance.

---

## Key Findings

- Replicated the **behavioral production effect**: better recognition for aloud words than silent words.  
- Identified a robust **pupillometric production effect**: larger pupil dilation during aloud trials than silent trials.  
- The magnitude of the pupillometric production effect correlated with the size of the behavioral production effect.  
- Control trials (saying “check”) produced pupil dilations similar to aloud trials but impaired memory for studied words, suggesting attentional engagement must be directed toward the **target word** to benefit memory.  
- Effects were observed **preceding, during, and following production**, consistent with a preparatory attentional component in addition to distinctive encoding.

---

## Data Description

The dataset includes:

- **Behavioral data**  
  - Item-level recognition judgments (hits, false alarms).  
  - Confidence ratings (in Experiments 2–4).  
  - Derived sensitivity measures (d′, familiarity, recollection estimates).  

- **Pupillometric data**  
  - Raw pupil size traces sampled via SR Research EyeLink 1000 Plus.  
  - Preprocessed signals (blink removal, interpolation, filtering, baseline correction, z-score normalization).  
  - Trial-wise averaged pupil dilation values (by condition).  

Full data, preprocessing scripts, and analysis code (a zip of the present directory) are also hosted on the **Open Science Framework (OSF): [https://osf.io/9mqhd/](https://osf.io/9mqhd/)**.

---

## Citation

If you use these data, please cite:

> Fawcett, J. M., Roberts, B. R. T., Willoughby, H., Tiller, J., Hourihan, K. L., & MacLeod, C. M.  
> The Pupillometric Production Effect: Evidence for Enhanced Processing Preceding, During, and Following Production.  
> *Cognition* (in revision). Preprint: [https://doi.org/10.2139/ssrn.5210001](https://doi.org/10.2139/ssrn.5210001)

A machine-readable citation file (`CITATION.cff`) is also included for GitHub integration.

---

## License

The dataset is licensed under **CC BY 4.0**. You are free to use, share, and adapt the materials provided that you give appropriate credit by citing the paper above.

---

## Contact

For questions about the dataset or analysis pipeline, please contact:  
**Jonathan M. Fawcett** – [Memorial University of Newfoundland]

If you use any of this analysis code in your own models, in whole or in part, please cite the above paper. If you find errors, please let me know.

---
