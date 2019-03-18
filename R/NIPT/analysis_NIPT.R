# 20190318

library(NIPTeR)

# setwd("~/workspace/project/NIPT")

# 读入参考组，正常情况下参考组应由自己实验室检测数据组成
NIPT_87_control_group <- readRDS(file="NIPTeR_cleaned_87_controls.rds")
# 校正CG值
bingc_control_group <- gc_correct(nipt_object=NIPT_87_control_group, method="bin")

# 读入受检数据，原始程序输入为bam文件
NIPT_raw_sample_bin <- readRDS(file="Trisomy21.rds")
# 校正CG值
bingc_sample <- gc_correct(nipt_object=NIPT_raw_sample_bin, method="bin")

# 卡方检验校正
NIPT_bin_chi_corrected_data <- chi_correct(nipt_sample=bingc_sample, nipt_control_group=bingc_control_group)
NIPT_bin_chi_corrected_sample <- NIPT_bin_chi_corrected_data$sample
NIPT_bin_chi_corrected_controls <- NIPT_bin_chi_corrected_data$control_group

# 计算样本均值
mean_match_sample <- mean(as.numeric(
  match_control_group(
    nipt_sample=NIPT_bin_chi_corrected_sample,
    nipt_control_group=NIPT_bin_chi_corrected_controls,
    n_of_samples=87,
    mode = "report"
    )
  ))


####### 计算标准Z-score ########
# 正常范围为[-3, 3]。当超出正常范围时，认为三体风险较高


# 计算13号染色体Z-score
z_score_13 <- calculate_z_score(nipt_sample=NIPT_bin_chi_corrected_sample,
                                nipt_control_group=NIPT_bin_chi_corrected_controls,
                                chromo_focus=13)
z_score_13$sample_Zscore

# 计算18号染色体Z-score
z_score_18 <- calculate_z_score(nipt_sample = NIPT_bin_chi_corrected_sample,
                                nipt_control_group=NIPT_bin_chi_corrected_controls,
                                chromo_focus=18)
z_score_18$sample_Zscore

# 计算21号染色体Z-score
z_score_21 <- calculate_z_score(nipt_sample = NIPT_bin_chi_corrected_sample,
                                nipt_control_group=NIPT_bin_chi_corrected_controls,
                                chromo_focus=21)
z_score_21$sample_Zscore




########### 计算归一化Z-score #############
# 正常范围为[-3, 3]。当超出正常范围时，认为三体风险较高


# 计算13号染色体Z-score
ncv_template_13 <- prepare_ncv(nipt_control_group=NIPT_bin_chi_corrected_controls,
                               chr_focus=13, max_elements=9,
                               use_test_train_set=F)
ncv_score_13 <- calculate_ncv_score(nipt_sample=NIPT_bin_chi_corrected_sample,
                                    ncv_template=ncv_template_13)
ncv_score_13$sample_score

# 计算18号染色体Z-score
ncv_template_18 <- prepare_ncv(nipt_control_group=NIPT_bin_chi_corrected_controls,
                               chr_focus=18, max_elements=9,
                               use_test_train_set=F)
ncv_score_18 <- calculate_ncv_score(nipt_sample=NIPT_bin_chi_corrected_sample,
                                    ncv_template=ncv_template_18)
ncv_score_18$sample_score

# 计算21号染色体Z-score
ncv_template_21 <- prepare_ncv(nipt_control_group=NIPT_bin_chi_corrected_controls,
                               chr_focus=21, max_elements=9,
                               use_test_train_set=F)
ncv_score_21 <- calculate_ncv_score(nipt_sample=NIPT_bin_chi_corrected_sample,
                                    ncv_template=ncv_template_21)
ncv_score_21$sample_score


########### 计算线性回归 Z-score #############
# 正常范围为[-3, 3]。当超出正常范围时，认为三体风险较高


# 计算13号染色体Z-score
RBZ_13 <- perform_regression(nipt_sample=NIPT_bin_chi_corrected_sample,
                             nipt_control_group=NIPT_bin_chi_corrected_controls, use_test_train_set=F,
                             chromo_focus=13)
RBZ_13$prediction_statistics

# 计算18号染色体Z-score
RBZ_18 <- perform_regression(nipt_sample=NIPT_bin_chi_corrected_sample,
                             nipt_control_group=NIPT_bin_chi_corrected_controls, use_test_train_set=F,
                             chromo_focus=18)
RBZ_18$prediction_statistics


# 计算21号染色体Z-score
RBZ_21 <- perform_regression(nipt_sample=NIPT_bin_chi_corrected_sample,
                             nipt_control_group=NIPT_bin_chi_corrected_controls, use_test_train_set=F,
                             chromo_focus=21)
RBZ_21$prediction_statistics
