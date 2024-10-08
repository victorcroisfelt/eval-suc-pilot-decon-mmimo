# An Evaluation of Successive Pilot Decontamination in Massive MIMO

This is a research-oriented code package that is primarily intended to allow readers to replicate the results of the article mentioned below and also encourage and accelerate further research on this topic:

Victor Croisfelt Rodrigues and Taufik Abrão, "[An Evaluation of Successive Pilot Decontamination in Massive MIMO](http://www.uel.br/revistas/uel/index.php/semexatas/article/view/34450/24955)", Semina: Ciências Exatas e Tecnológicas, Londrina, Brazil, v. 39, n. 2, ago./dez. 2018. 

The package is based on the Matlab language and can, in fact, reproduce all the numerical results and figures discussed in the article. To contextualize, in the sequel, we present the abstract of the article and other important information.

I hope this content helps in your research and contributes to building the precepts behind open science. Remarkably, in order to boost the idea of open science and further drive the evolution of science, we also motivate you to share your published results to the public.

If you have any questions or if you have encountered any inconsistency, please do not hesitate to contact me via victorcroisfelt@gmail.com.

## Abstract
The demand for higher data rates can be satisfied by the spectral efficiency (SE) improvement offered by Massive multiple-input multiple-output (M-MIMO) systems. However, the pilot contamination remains as a fundamental issue to obtain the paramount SE in such systems. This propitiated the research of several methods to mitigate pilot contamination. One of these procedures is based on the coordination of the cells, culminating in proposals with multiple pilot training phases. This paper aims to expand the results of the original paper, whereby the concepts of large pilot training phases were offered. The evaluation of such method was conducted through more comprehensible numerical results, in which a large number of antennas were assumed and more rigorous SE expressions were used. The channel estimation approaches relying on multiple pilot training phases were considered cumbersome for implementation and an uninteresting solution to overcome pilot contamination; contradicting the results presented in the genuine paper.

## Content
From this code package, all figures reported in the above article can be simulated. To do this, you can run the script "simulationAllFigures.m", which, in turn, calls the files labeled with "function" in their names. Further details about each file can be found inside them.

## Acknowledgments
This work was supported in part by the Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq) of Brazil under Grant 102050/2018-0, and in part by the Universidade Estadual de Londrina (UEL), Paraná State Government under the PIBITI grant of the public notice 02/2018.

## Citing this Repository and License
This code is subject to the GPLv3 license. If you use any part of this repository for research, please consider citing our aforementioned work.

```bibtex
@article{Rodrigues_Abrão_2018,
  title={An evaluation of successive pilot decontamination in massive MIMO},
  volume={39},
  url={https://ojs.uel.br/revistas/uel/index.php/semexatas/article/view/34450},
  DOI={10.5433/1679-0375.2018v39n2p107},
  number={2},
  journal={Semina: Ciências Exatas e Tecnológicas},
  author={Rodrigues, Victor Croisfelt and Abrão, Taufik},
  year={2018},
  month={Dec.},
  pages={107–114}
}
