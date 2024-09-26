# Analysis-of-pan-age-PBMC-scRNAseq

The immune system undergoes progressive functional remodeling from neonatal stages to old age. Therefore, understanding how aging shapes immune cell function is vital for precise treatment of patients at different life stages. Here, we constructed the first transcriptomic atlas of immune cells encompassing human lifespan, ranging from newborns to supercentenarians, and comprehensively examined gene expression signatures involving cell signaling, metabolism, differentiation, and functions in all cell types to investigate immune aging changes.


![image](https://github.com/user-attachments/assets/69899221-90ed-4856-9fda-3df130739a1e)


All the necessary data and code to reproduce the analysis and figures presented in this paper are available at https://doi.org/10.5281/zenodo.13838466. 

# PHARE model
All cells received cell annotation based on Panage_data, and then the proportion of all immune cell types was calculated to train the PHARE model (https://github.com/cliffren/PHARE).

We developed a novel physiological age prediction model (PHARE) using PBMC single-cell transcriptome data from healthy individuals, leveraging a machine learning approach.


![image](https://github.com/user-attachments/assets/d70b990a-f53b-40f3-916c-ad8a58f2ecd4)

The web version of the PHARE tool is available at: https://xiazlab.org/phare/

# REFERENCE 
Please cite this paper if you use our data or code.

Cangang Zhang, Tao Ren, Xiaofan Zhao, Yanhong Su, Qianhao Wang, Tianzhe Zhang, Boxiao He, Ling-Yun Wu, Lina Sun, Baojun Zhang, Zheng Xia. Biologically informed machine learning modeling of immune cells to reveal physiological and pathological aging process, bioRxiv 2024.04.01.587649; doi: https://doi.org/10.1101/2024.04.01.587649
