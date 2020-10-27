# LPI-SKF
This model is proposed to **predict lncRNA-protein interactions**, which can be summarized in four steps:  
**Firstly**, we downloaded experimentally verified lncRNA-protein interactions from the NPInter v2.0 database. 4158 lncRNA-protein interactions were selected according to the previous work, including 990 lncRNAs and 27 proteins. The interactions are shown in the Data/Interaction matrix.csv.  
**Secondly**, three different similarities of lncRNAs and three different similarities of proteins are calculated in different ways, respectively. Among these similarities, the lncRNA expression profile similarity, the protein pairwise sequence alignment similarity, and the sequence statistical feature similarity of both lncRNAs and proteins are shown in the Code/LPI-SKF/data.mat.  
**Thirdly**, the SKF approach is utilized to integrate the lncRNA similarities and protein similarities.   
**Finally**, the Laplacian Regularized Least Squares frame is applied to obtain our final prediction matrix. The code of SKF and LapRLS is shown in Code/LPI-SKF/LPI_SKF.m.  
**All data and code will be introduced in different folders**.  
