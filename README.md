# Analysis-of-pan-age-PBMC-scRNAseq

The immune system undergoes progressive functional remodeling from neonatal stages to old age. Therefore, understanding how aging shapes immune cell function is vital for precise treatment of patients at different life stages. Here, we constructed the first transcriptomic atlas of immune cells encompassing human lifespan, ranging from newborns to supercentenarians, and comprehensively examined gene expression signatures involving cell signaling, metabolism, differentiation, and functions in all cell types to investigate immune aging changes.

This repository store all the necessary data and codes to reproduce the analysis and figures presented in this paper. Other big datasets are available at https://doi.org/10.5281/zenodo.13838466. 


![image](https://github.com/user-attachments/assets/a62c48bf-c140-4540-8822-e1b4c8812a9f)



# PHARE model

We developed a novel physiological age prediction model (PHARE) using PBMC single-cell transcriptome data from healthy individuals, leveraging a machine learning approach.

All cells received cell annotation based on Panage_data, and then the proportion of all immune cell types was calculated to train the PHARE model (https://github.com/cliffren/PHARE).

![image](https://github.com/user-attachments/assets/d70b990a-f53b-40f3-916c-ad8a58f2ecd4)

The web version of the PHARE tool is available at: https://xiazlab.org/phare/.

# REFERENCE 
Please cite this paper if you use our data or code.

Cangang Zhang, Tao Ren, Xiaofan Zhao, Yanhong Su, Qianhao Wang, Tianzhe Zhang, Boxiao He, Ling-Yun Wu, Lina Sun, Baojun Zhang, Zheng Xia. Biologically informed machine learning modeling of immune cells to reveal physiological and pathological aging process, bioRxiv 2024.04.01.587649; doi: https://doi.org/10.1101/2024.04.01.587649.

# Contact us
Please feel free to contact us at the following email addresses: cg_zhang2021@163.com; rentao@amss.ac.cn.
