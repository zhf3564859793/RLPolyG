# RLPolyG
<!-- GETTING STARTED -->
## Getting Started

This is the code for [**De Novo Design of Polymers with Specified Properties Using Reinforcement Learning**](https://pubs.acs.org/doi/10.1021/acs.macromol.5c00427).

![image](https://github.com/zhf3564859793/RLPolyG/blob/main/Figure/Figure%201.jpg)


### Prerequisites

The version of python:

  ```sh
  python==3.9.19
  ```

### Installation

You can use pip to install the required packagess.

  ```sh
  pip install requirements.txt -r
  ```

<!-- USAGE EXAMPLES -->
## Usage

* Forward Model
  
  You can find the forward model in the Forward_model_revised_final.ipynb. The forward model predicts the target property (yield strength) of polymers based on their structural representation. It is a critical component of RLPolyG, enabling property evaluation during the reinforcement learning process.

* Inverse Model
  
  The inverse model aims to generate polymers with specified properties. You can find the inverse model in the RLPolyG_final.ipynb.
  1. Unbiased Generation Model: The inverse model is initially trained using the PI1M dataset to generate polymers without property bias. This model learns the underlying polymer structure distribution.
  2. Biased Generation Model: The ReLease algorithm is then applied to refine the model for property-biased generation. Reinforcement learning is used to guide the model toward generating polymers with the desired properties, such as high yield strength.

* Screening for Synthetic Accessibility and Degradability
  
  1. Synthetic Accessibility: You can find the related model in the SA score folder.
  2. Degradability: You can find the related model in the degradability prediction folder. We first collected degradability data from Yuan et al.’s work. Using this data, we constructed a random forest model to predict the degradability score of the generated polymers, which also took MFF as input features.


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [ReLeaSE](https://github.com/isayev/ReLeaSE)
* [Polymer-degradability-ranking](https://github.com/tsudalab/Polymer-degradability-ranking)
* [SA score](https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py)

## How to cite

```bibtex
@article{zhou2025novo,
  title = {De Novo Design of Polymers with Specified Properties Using Reinforcement Learning},
  author = {Zhou, Haifan and Fang, Yue and Li, Linyan and Liu, Pingwei and Gao, Hanyu},
  journal = {Macromolecules},
  year = {2025},
  publisher = {ACS Publications},
  doi = {10.1021/acs.macromol.5c00427},
  url = {https://doi.org/10.1021/acs.macromol.5c00427}
}
