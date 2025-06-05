# ProKcat
Multimodal Regression for Enzyme Turnover Rates Prediction
This paper has been published in [[IJCAI 2025]](https://2025.ijcai.org/). This is the code.

## Data
BRENDA Release 2025.1 is now online.
This new release includes:
168 new EC Classes and 1620 updated EC Classes
6,857 new primary literature references
An updated metabolic pathway map featuring five new pathways: Glutathione-mediated detoxification, Curcuminoid biosynthesis, Monoterpenoid biosynthesis, Tropane alkaloid biosynthesis, Secologanin biosynthesis.
You can download the updated data in JSON and TXT formats [[here]](https://www.brenda-enzymes.org/download.php).

## Dependency
Pytorch (1.8.1+cu101)

Scikit-learn

esm

RDKit

BRENDApyrser

...

refer to the requirements.txt


## Thanks
Thanks for the work [[DLTKcat]](https://github.com/SizheQiu/DLTKcat). The data and baseline models are mainly obtained from this repository.

## Citation
If you find this project useful for your research, please use the following BibTeX entry.

```
@inproceedings{hu2025Multimodal,
  title={Multimodal Regression for Enzyme Turnover Rates Prediction.},
  author={Hu, Bozhen and Tan, Cheng and Li, Siyuan and Zheng, Jiangbin and Xia, Jun and Li, Stan Z.},
  booktitle={Thirty-fourth International Joint Conference on Artificial Intelligence (IJCAI 2025)},
  year={2025},
  organization={International Joint Conferences on Artificial Intelligence Organization}
}
@article{qiu2024dltkcat,
  title={DLTKcat: deep learning-based prediction of temperature-dependent enzyme turnover rates},
  author={Qiu, Sizhe and Zhao, Simiao and Yang, Aidong},
  journal={Briefings in Bioinformatics},
  volume={25},
  number={1},
  pages={bbad506},
  year={2024},
  publisher={Oxford University Press}
}
```

