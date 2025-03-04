library(dplyr)
library(rlang)
library(ratesci)

#结局指标是0和1赋值
f_result$pingjia <- as.numeric(f_result$pingjia)

#参数调试
# data=f_result
# treatment="ARM"
# strata="SITEID"
# response="pingjia"
# level=0.05
# AA="A"
# BB="B"

stratCI <- function(data, treatment, strata, response, level, AA, BB) {
  # 检查输入参数
  if (is.null(data)) {
    stop("ERROR: Input dataset not specified")
  }
  if (is.null(treatment)) {
    stop("ERROR: Treatment group must be specified")
  }
  if (is.null(strata)) {
    stop("ERROR: Strata must be specified")
  }
  if (is.null(response)) {
    stop("ERROR: Response must be specified")
  }
  
  # 提取所需的列
  temp <- data[, c(treatment, strata, response)]
  colnames(temp) <- c("treatment", "strata", "response")
  
  # 计算alpha
  alpha <- 1 - level / 2
  
  # 计算每个strata的riskdiff
  ready <- temp %>%
    group_by(strata) %>%
    summarise(
      n_trt1 = sum(treatment==AA, na.rm = TRUE),
      n_trt2 = sum(treatment==BB, na.rm = TRUE),
      prop1 = mean(response[treatment == AA], na.rm = TRUE),
      prop2 = mean(response[treatment == BB], na.rm = TRUE),
      r_diff = prop1 - prop2,
      ASE = sqrt((prop1 * (1 - prop1) / n_trt1) + (prop2 * (1 - prop2) / n_trt2))
    )
  
  # 计算权重
  ready$CMH_part <- (ready$n_trt1 * ready$n_trt2) / (ready$n_trt1 + ready$n_trt2)
  ready$invar_part <- 1 / ready$ASE^2
  #对于无穷的逆方差，赋值为0
  ready$invar_part[is.infinite(ready$invar_part)] <- 0
  ready$min_part1 <- ready$n_trt1 + ready$n_trt2
  ready$min_part2 <- ready$min_part1 * ready$r_diff
  ready$min_part3 <- ready$invar_part * ready$r_diff
  
  
  ready$CMH_sum <- sum(ready$CMH_part)
  ready$invar_sum <- sum(ready$invar_part)
  ready$min_sum1 <- sum(ready$min_part1)
  ready$min_sum2 <- sum(ready$min_part2)
  ready$min_sum3 <- sum(ready$min_part3)
  
  ready$min_alpha <- ready$r_diff*ready$invar_sum - ready$min_sum3;
  ready$min_beta <- ready$invar_part * (1 + ready$min_alpha*ready$min_sum2/ready$min_sum1);
  ready$min_part4 <- ready$r_diff*ready$min_beta;
  ready$min_part5 <- ready$r_diff*ready$min_alpha*ready$invar_part;
  
  ready$min_sum4 <- sum(ready$min_part4)
  ready$min_sum5 <- sum(ready$min_part5)
  
  # 计算权重
  ready$CMH_weight <- ready$CMH_part / ready$CMH_sum
  ready$INVAR_weight <- ready$invar_part / ready$invar_sum
  ready$MIN_weight <- ready$min_beta / ready$invar_sum - 
    (ready$min_alpha * ready$invar_part / (ready$invar_sum + ready$min_sum5)) * 
    (ready$min_sum4 / ready$invar_sum)
  
  # 读取权重数据
  w <- ready[, c("CMH_weight", "INVAR_weight", "MIN_weight")]
  # 横向合并 n_trt1 和 n_trt2
  n <- cbind(ready$n_trt1, ready$n_trt2)
  # 总样本量
  totaln <- sum(n)
  # 层数
  n_strata <- nrow(n)
  # 横向合并 prop1 和 prop2
  p <- cbind(ready$prop1, ready$prop2)
  # 计算两组的有效例数
  p1_1 <- ready$n_trt1 * ready$prop1
  p2_1 <- ready$n_trt2 * ready$prop2
  # 计算两组的总体例数
  n1 <- sum(ready$n_trt1)
  n2 <- sum(ready$n_trt2)
  # 计算两组的总体有效率
  p1 <- sum(p1_1) / n1
  p2 <- sum(p2_1) / n2
  
  #不进行校正的wlad结果
  Crude_Lower = (p1-p2) - qnorm(alpha) * sqrt(p1 * (1-p1)/n1 + p2 * (1-p2)/n2)
  Crude_Upper = (p1-p2) + qnorm(alpha) * sqrt(p1 * (1-p1)/n1 + p2 * (1-p2)/n2)
  #start to calculate Wald CI
  adj_rdiff= w * ready$r_diff
  wald_part1= apply(adj_rdiff,2,sum)
  wald_var = 1/ready$invar_part
  wald_var[is.infinite(wald_var)] <- 0  #去除无穷大的影响
  sq_w = w ^ 2
  adj_var = sq_w * wald_var
  wald_part2 = apply(adj_var,2,sum)
  #wald置信区间
  Wald_L <- wald_part1 - qnorm(alpha)*sqrt(wald_part2)
  Wald_U <-  wald_part1 + qnorm(alpha)*sqrt(wald_part2)
  Wald_Lower <- as.matrix(Wald_L)
  Wald_Upper <- as.matrix(Wald_U)
  
  ####### 计算Newcombe置信区间 #######
  #计算Z
  partz_a1 <- matrix(1, nrow = n_strata, ncol = 3)
  partz_a2 <- matrix(1, nrow = n_strata, ncol = 3)
  partz_b1 <- matrix(1, nrow = n_strata, ncol = 3)
  partz_b2 <- matrix(1, nrow = n_strata, ncol = 3)
  # 按层和列循环计算
  for (i in 1:n_strata) {
    for (j in 1:3) {
      partz_a1[i, j] <- as.numeric((sq_w[i, j] / ready$n_trt1[i]) * (p[i, 1] * (1 - p[i, 1])))
      partz_b1[i, j] <- as.numeric(w[i, j] * sqrt(p[i, 1] * (1 - p[i, 1]) / ready$n_trt1[i]))
      partz_a2[i, j] <- as.numeric((sq_w[i, j] / ready$n_trt2[i]) * (p[i, 2] * (1 - p[i, 2])))
      partz_b2[i, j] <- as.numeric(w[i, j] * sqrt(p[i, 2] * (1 - p[i, 2]) / ready$n_trt2[i]))
    }
  }
  z_value1 <- as.matrix(qnorm(alpha) * sqrt(apply(partz_a1,2,sum)) / apply(partz_b1,2,sum))
  z_value2 <- as.matrix(qnorm(alpha) * sqrt(apply(partz_a2,2,sum)) / apply(partz_b2,2,sum))
  z_value1[is.na(z_value1)] <- 0  #去除NA
  z_value2[is.na(z_value2)] <- 0
  
  #calculate L1 L2 U1 U2 based on stratified Wilson CI
  wilson_part1a <- matrix(0, nrow = n_strata, ncol = 3)
  wilson_part1b <- matrix(0, nrow = n_strata, ncol = 3)
  wilson_part2a <- matrix(0, nrow = n_strata, ncol = 3)
  wilson_part2b <- matrix(0, nrow = n_strata, ncol = 3)
  for (i in 1:n_strata) {
    for (j in 1:3) {
      wilson_part1a[i, j] <- (n[i, 1] * p[i, 1] + 0.5 * (z_value1[j, 1]^2)) / (n[i, 1] + z_value1[j, 1]^2)
      wilson_part1b[i, j] <- ((n[i, 1] * z_value1[j, 1]) / (n[i, 1] + z_value1[j, 1]^2)) *
        (sqrt(4 * n[i, 1] * p[i, 1] * (1 - p[i, 1]) + z_value1[j, 1]^2) / (2 * n[i, 1]))
      wilson_part2a[i, j] <- (n[i, 2] * p[i, 2] + 0.5 * (z_value2[j, 1]^2)) / (n[i, 2] + z_value2[j, 1]^2)
      wilson_part2b[i, j] <- ((n[i, 2] * z_value2[j, 1]) / (n[i, 2] + z_value2[j, 1]^2)) *
        (sqrt(4 * n[i, 2] * p[i, 2] * (1 - p[i, 2]) + z_value2[j, 1]^2) / (2 * n[i, 2]))
    }
  }
  wilson_L1 <- w * (wilson_part1a - wilson_part1b)
  wilson_U1 <- w * (wilson_part1a + wilson_part1b)
  wilson_L2 <- w * (wilson_part2a - wilson_part2b)
  wilson_U2 <- w * (wilson_part2a + wilson_part2b)
  
  newcombe_L1 <- apply(wilson_L1,2,sum)
  newcombe_U1 <- apply(wilson_U1,2,sum)
  newcombe_L2 <- apply(wilson_L2,2,sum)
  newcombe_U2 <- apply(wilson_U2,2,sum)
  
  lam1 <- sq_w / ready$n_trt1
  lam2 <- sq_w / ready$n_trt2
  lambda1 <- apply(lam1,2,sum)
  lambda2 <- apply(lam2,2,sum)
  
  NEWCOMBE_L <- wald_part1 - qnorm(alpha) * sqrt(lambda1 * newcombe_L1 * (1 - newcombe_L1) + lambda2 * newcombe_U2 * (1 - newcombe_U2))
  NEWCOMBE_U <- wald_part1 + qnorm(alpha) * sqrt(lambda2 * newcombe_L2 * (1 - newcombe_L2) + lambda1 * newcombe_U1 * (1 - newcombe_U1))
  
  Newcombe_Lower <- as.matrix(NEWCOMBE_L)
  Newcombe_Upper <- as.matrix(NEWCOMBE_U)
  
  
  #输出结果
  # 分层汇总表
  Strata_sum <- data.frame(
    Strata = 1:n_strata,
    n_trt1 = ready$n_trt1,
    n_trt2 = ready$n_trt2,
    prop1 = ready$prop1,
    prop2 = ready$prop2,
    r_diff = ready$r_diff,
    CMH_weight = ready$CMH_weight,
    invar_weight = ready$INVAR_weight,
    min_weight = ready$MIN_weight
  )
  
  # 未调整的置信区间
  Crude_NW <- moverci(x1 = sum(p1_1), n1 = n1, x2 = sum(p2_1), n2 = n2, type = "wilson")
  Crude_result <- data.frame(
    Crude_dif = sprintf("%0.4f",p1-p2),
    Crude_Lower = sprintf("%0.4f",Crude_Lower),
    Crude_Upper = sprintf("%0.4f",Crude_Upper),
    Crude_NewcombeCIL = sprintf("%0.4f",Crude_NW[1]),
    Crude_NewcombeCIU = sprintf("%0.4f",Crude_NW[3])
  )
  
  # 调整后的结果
  Adj_result <- data.frame(
    Adj_dif = sprintf("%0.4f",wald_part1),
    WaldCIL = sprintf("%0.4f",Wald_Lower),
    WaldCIU = sprintf("%0.4f",Wald_Upper),
    NewcombeCIL = sprintf("%0.4f",Newcombe_Lower),
    NewcombeCIU = sprintf("%0.4f",Newcombe_Upper)
  )
  
  return(list(Strata_sum=Strata_sum,Crude_result=Crude_result,Adj_result=Adj_result))
}

# stratCI(data=f_result, treatment=f_result$ARM, strata=f_result$SITEID, response=f_result$pingjia, level=0.05)
stratCI(data=f_result, treatment="ARM", strata="SITEID", response="pingjia", level=0.05,AA="A",BB="B")


fff <- read.csv("D:/心腔内超声/CMH率差/f_result1.csv")
fff$pingjia <- as.numeric(fff$pingjia)
stratCI(data=fff, treatment="ARM", strata="SITEID", response="pingjia", level=0.05,AA="A",BB="B")
