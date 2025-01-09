# RLPolyG
<!-- GETTING STARTED -->
## Getting Started

This is the code for **De Novo Design of Polymers with Specified Properties Using Reinforcement Learning**.


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

* Forward model
  
You can find the forward model in the forward_model_0919.ipynb. The forward model predicts the target property (yield strength) of polymers based on their structural representation. It is a critical component of RLPolyG, enabling property evaluation during the reinforcement learning process.

* Inverse model
  
  The inverse model aims to generate polymers with specified properties.
  1. Unbiased Generation Model: The inverse model is initially trained using the PI1M dataset to generate polymers without property bias. This model learns the underlying polymer structure distribution.
  2. The ReLease algorithm is then applied to refine the model for property-biased generation. Reinforcement learning is used to guide the model toward generating polymers with the desired properties, such as   
     high or low yield strength.
