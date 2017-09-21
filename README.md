# EPLMI
Author: Yu-An Huang, Keith C.C. Chan1,Zhu-Hong You
Paper: Constructing Prediction Models from Expression Profiles for Large Scale lncRNA-miRNA Interaction Profiling

---------------------
File description
---------------------
Data:
  ----------------------------
  Expression_profile_data: 
    invalid_association_expression.m: The ids of known lncRNA-miRNA interactions without either of lncRNA and miRNA expression profiles.
    invalid_lnc_expression.m: The ids of lncRNAs without lncRNA expression profiles.
    invalid_mi_expression.m: The ids of lncRNAs without lncRNA expression profiles.
    Known_lncRNA_miRNA_association.m: The id pairs of known lncRNA-miRNA interactions.
    lnc_expression_similarity_matrix.m: The lncRNA-lncRNA similarity matrix of expression profiles.
    mi_expression_similarity_matrix.m: The miRNA-miRNA similarity matrix of expression profiles.
  ----------------------------
  Function_data:
    invalid_association_function.m: The ids of known lncRNA-miRNA interactions without either of lncRNA and miRNA function feature.
    invalid_lnc_function.m: The ids of lncRNAs without lncRNA function feature.
    invalid_mi_function.m: The ids of lncRNAs without lncRNA function feature.
    Known_lncRNA_miRNA_association.m: The id pairs of known lncRNA-miRNA interactions.
    lnc_function_similarity_matrix.m: The lncRNA-lncRNA similarity matrix of function.
    mi_function_similarity_matrix.m: The miRNA-miRNA similarity matrix of function.
  ----------------------------
  Sequence_data:
    invalid_lnc_seq.m: The ids of lncRNAs without lncRNA sequence information.
    Known_lncRNA_miRNA_association.m: The id pairs of known lncRNA-miRNA interactions.
    lnc_seq_similarity_matrix.m: The lncRNA-lncRNA similarity matrix of sequence.
    mi_seq_similarity_matrix.m: The miRNA-miRNA similarity matrix of sequence.
  ----------------------------
  lncRNA-miRNA_interaction:
    lncRNA_id.m The id list of lncRNA name.
    miRNA_id.m The id list of miRNA name.
    lncRNA-miRNA_idlist.m: The lncRNA-miRNA interaction data downloaded from lncRNASNP.
--------------------- 
Cross_validation:
  leave-one-out_CV:
    EPLMI_LOOCV.m: The codes for the implemention of leave-one-out cross validation.
    Plot_roc_curve: The codes for plotting the ROC curves based on the LOOCV results.
  ----------------------------
  5-fold_CV:
    main.m: The main program for the implemention of 5-fold cross validation.
    Generate_random_roder_5fold.m: The codes for generating the random sample lists.
    EPLMI5cv: The subprogram for the implemention of 5-fold cross validation.
    Position2AUC.m: The codes for plotting the ROC curves based on the ranks of testing samples.
  ----------------------------
  
