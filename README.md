# Acoustic-Beamforming-Advanced

MATLAB code for the baseline of "Learning an Interpretable End-to-End Network for Real-Time Acoustic Beamforming".

If you use the code, please cite our paper:
> [Liang, Hao and Zhou, Guanxing and Tu, Xiaotong and Jakobsson, Andreas and Ding, Xinghao and Huang, Yue, "Learning an Interpretable End-to-End Network for Real-Time Acoustic Beamforming", *Journal of Sound and Vibration*, 2024.](https://doi.org/10.1016/j.jsv.2024.118620 "https://doi.org/10.1016/j.jsv.2024.118620")

Also, the main code for "Learning an Interpretable End-to-End Network for Real-Time Acoustic Beamforming" is available at https://github.com/JoaquinChou/DAMAS_FISTA_Net or https://github.com/HauLiang/DAMAS-FISTA-Net (latest version).

## Preparation
- MATLAB >= R2022a

## Getting Started
This scan-frequency demo for acoustic imaging mainly contains the model-based acoustic beamforming methods.
> The previous version can be found at https://github.com/HauLiang/Acoustic-Beamforming-Methods.

As for the model-based methods:
* DAS: <br>
> [Van Veen, Barry D and Buckley, Kevin M, "Beamforming: A versatile approach to spatial filtering", *IEEE assp magazine*, 1988.](https://ieeexplore.ieee.org/abstract/document/665/ "https://ieeexplore.ieee.org/abstract/document/665/")
* MUSIC: <br>
> [Schmidt, Ralph, "Multiple emitter location and signal parameter estimation", *IEEE transactions on antennas and propagation*, 1986.](https://ieeexplore.ieee.org/abstract/document/1143830/ "https://ieeexplore.ieee.org/abstract/document/1143830/")
* DAMAS: <br>
> [Brooks, Thomas F and Humphreys, William M, "A deconvolution approach for the mapping of acoustic sources (DAMAS) determined from phased microphone arrays", *Journal of sound and vibration*, 2006.](https://www.sciencedirect.com/science/article/pii/S0022460X06000289 "https://www.sciencedirect.com/science/article/pii/S0022460X06000289")
* DAMAS2: <br>
> [Dougherty, Robert, "Extensions of DAMAS and benefits and limitations of deconvolution in beamforming", *11th AIAA/CEAS aeroacoustics conference*, 2005.](https://doi.org/10.2514/6.2005-2961 "https://doi.org/10.2514/6.2005-2961")
* DAMAS-FISTA: <br>
> [Liang, Hao and Zhou, Guanxing and Tu, Xiaotong and Jakobsson, Andreas and Ding, Xinghao and Huang, Yue, "Learning an Interpretable End-to-End Network for Real-Time Acoustic Beamforming", *The Journal of Sound and Vibration*, 2024.](https://doi.org/10.1016/j.jsv.2024.118620 "https://doi.org/10.1016/j.jsv.2024.118620")
* CLEAN-PSF: <br>
> [Högbom, JA, "Aperture synthesis with a non-regular distribution of interferometer baselines", *Astronomy and Astrophysics Supplement Series*, 1974.](http://adsabs.harvard.edu/pdf/1974A&AS...15..417H "http://adsabs.harvard.edu/pdf/1974A&AS...15..417H") 
* CLEAN-SC: <br>
> [Sijtsma, Pieter, "CLEAN based on spatial source coherence", *International journal of aeroacoustics*, 2007.](https://journals.sagepub.com/doi/abs/10.1260/147547207783359459 "https://journals.sagepub.com/doi/abs/10.1260/147547207783359459")
* FFT-NNLS: <br>
> [Ehrenfried, Klaus and Koop, Lars, "Comparison of iterative deconvolution algorithms for the mapping of acoustic sources", *AIAA journal*, 2007.](https://arc.aiaa.org/doi/abs/10.2514/1.26320?journalCode=aiaaj "https://arc.aiaa.org/doi/abs/10.2514/1.26320?journalCode=aiaaj")
* FFT-FISTA are available at https://github.com/dtu-act/EMBED: <br>
> [Lylloff, Oliver and Fernández-Grande, Efrén and Agerkvist, Finn and Hald, Jørgen and Tiana Roig, Elisabet and Andersen, Martin S. "Improving the efficiency of deconvolution algorithms for sound source localization". *The journal of the acoustical society of America*, 2015.](http://dx.doi.org/10.1121/1.4922516 "http://dx.doi.org/10.1121/1.4922516")
* FFT-DFISTA: <br>
> [Ding, Xinghao and Liang, Hao and Jakobsson, Andreas and Tu, Xiaotong and Huang, Yue. "High-Resolution Source Localization Exploiting the Sparsity of the Beamforming Map", *Signal Processing*, 2022.](https://www.sciencedirect.com/science/article/pii/S016516842100414X "https://www.sciencedirect.com/science/article/pii/S016516842100414X")

As for the deep network-based method:
* Acoustic-Net are available at https://github.com/JoaquinChou/Acousitc-Net: <br>
> [Zhou, Guanxing and Liang, Hao and Ding, Xinghao and Huang, Yue and Tu, Xiaotong and Abbas, Saqlain. "Acoustic-Net: A Novel Neural Network for Sound Localization and Quantification", in *19th Asia-Pacific Vibration Conference (APVC2021)*, 2022.](https://arxiv.org/abs/2203.16988 "https://arxiv.org/abs/2203.16988")

If you want to know more about acoustic imaging, please refer to the following papers:

* Fundamentals of acoustic beamforming:
> [de Santana, Leandro, "Fundamentals of Acoustic Beamforming", *Design and Operation of Aeroacoustic Wind Tunnel Tests for Group and Air Transport*, 2017.](https://www.sto.nato.int/publications/STO%20Educational%20Notes/STO-EN-AVT-287/EN-AVT-287-04.pdf "https://www.sto.nato.int/publications/STO%20Educational%20Notes/STO-EN-AVT-287/EN-AVT-287-04.pdf")
* A review of acoustic imaging methods:
> [Merino-Martínez, Roberto and Sijtsma, Pieter and Snellen, Mirjam and Ahlefeldt, Thomas and Antoni, Jerome and Bahr, Christopher J and Blacodon, Daniel and Ernst, Daniel and Finez, Arthur and Funke, Stefan and others, "A review of acoustic imaging methods using phased microphone arrays", *CEAS Aeronautical Journal*, 2019.](https://link.springer.com/article/10.1007/s13272-019-00383-4 "https://link.springer.com/article/10.1007/s13272-019-00383-4")

@ All rights are reserved by the authors.
