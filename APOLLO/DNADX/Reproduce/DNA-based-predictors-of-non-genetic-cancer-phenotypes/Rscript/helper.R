# calculate gene signatures
calc_signatures <- function(x, gmtFile, method="mean", scale=F){
  
  geneset.obj<- GSA.read.gmt(gmtFile)
  genenames<-row.names(x)
  np=length(geneset.obj$genesets)
  
  if(scale){xs=t(scale(t(x),center=T,scale=T))}else{xs<-x}
  
  if(method!="gsa"){
    val=matrix(NA,nrow=np,ncol=ncol(x))
    for(i in 1:np){
      gene.set=which(genenames %in% geneset.obj$genesets[[i]])
      gene.set=gene.set[!is.na(gene.set)]
      if(length(gene.set)>1){
        if(method=="mean"){
          val[i,]=colSums(xs[gene.set,,drop=F])/length(gene.set)
        }
        if(method=="median"){
          val[i,]=t(apply(xs[gene.set,,drop=F],2,median))
        }
        if(method=="pca"){
          y<-prcomp(as.matrix(xs[gene.set,,drop=F]))
          val[i,]=y$rotation[,1]
        }
      } else if (length(gene.set) == 1) {
        val[i,] <- unlist(xs[gene.set,])
      }
    }
  }else{
    val<-GSA.make.features(gsaObj,xs,geneset.obj$genesets,genenames)
  }
  dimnames(val)<-list(geneset.obj$geneset.names,dimnames(x)[[2]])
  return(val)
}

overlapSets<-function(x,y){
  # Ensure x and y are matrices to begin with
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(y)) y <- as.matrix(y)

  rn_x_orig <- rownames(x)
  rn_y_orig <- rownames(y)

  if (is.null(rn_x_orig) || is.null(rn_y_orig)) {
    warning("Input matrices to overlapSets must have rownames.")
    # Return empty matrices with the same number of columns
    return(list(x=x[0, , drop=FALSE], y=y[0, , drop=FALSE]))
  }

  # Ensure original rownames are character vectors for intersection
  rn_x_char <- as.character(rn_x_orig)
  rn_y_char <- as.character(rn_y_orig)
  
  common_rownames <- intersect(rn_x_char, rn_y_char)
  
  if (length(common_rownames) == 0) {
    warning("No common rownames in overlapSets.")
    return(list(x=x[0, , drop=FALSE], y=y[0, , drop=FALSE]))
  }
  
  # Filter x and y to include only common rownames, preserving matrix structure
  x_filtered <- x[common_rownames, , drop=FALSE]
  y_filtered <- y[common_rownames, , drop=FALSE]
  
  # Now sort the filtered matrices by their rownames.
  # Rownames of x_filtered and y_filtered should be common_rownames.
  
  # Get current rownames of filtered matrices and ensure they are characters
  current_rn_x <- as.character(rownames(x_filtered))
  current_rn_y <- as.character(rownames(y_filtered)) # Key for the error line

  if (length(current_rn_x) > 0) {
    x_sorted <- x_filtered[sort.list(current_rn_x), , drop=FALSE]
  } else {
    x_sorted <- x_filtered # Handles empty case (e.g. common_rownames was empty, though checked above)
  }

  if (length(current_rn_y) > 0) {
    y_sorted <- y_filtered[sort.list(current_rn_y), , drop=FALSE] # This was the error location
  } else {
    y_sorted <- y_filtered # Handles empty case
  }
  
  return(list(x=x_sorted, y=y_sorted))
}
assignDiffScore.dwd<-function(x,y){
  both<-overlapSets(x,y) 	# get the overlap of genes
  
  # Check if overlap is empty
  if (nrow(both$x) == 0 || ncol(both$x) == 0 || nrow(both$y) == 0 || ncol(both$y) == 0) {
    warning("Empty overlap in assignDiffScore.dwd, returning NAs")
    return(rep(NA, ncol(y))) # Return NAs for each sample in the original y (edata_matrix)
  }
  
  both$x<- apply(both$x,2,function(x){sign(x)*sqrt(x^2/sum(x^2))}) 	# set the distance to the origin to 1
  both$y<- apply(both$y,2,function(x){sign(x)*sqrt(x^2/sum(x^2))})	# set the distance to the origin to 1
  msproj<- apply(both$y,2,function(x,y){x%*%y},both$x[,1])	# project the samples on the MS-pL axis
  mlproj<- apply(both$y,2,function(x,y){x%*%y},both$x[,2])	# project the samples on the pL-mL axis
  diffScore<- mlproj - msproj 
  return( diffScore )	# return the point on the differentiation axis
}

GHI_RS <- function(edata) {
  getGene <- function(edata,gene_id){
    # Ensure gene_id is character for matching rownames
    gene_id_char <- as.character(gene_id)
    if (gene_id_char %in% rownames(edata)) {
      val <- edata[gene_id_char, ]
      return(unlist(val))
    } else {
      warning(paste("Gene", gene_id_char, "not found in edata for GHI_RS. Returning NA vector."))
      return(rep(NA, ncol(edata))) # Return NA for all samples if gene not found
    }
  }
  gene_ids <- c(2597,2990,60,7037,6175,2886,2064,596,5241,57758,2099,6790,4605,891,332,4288,4320,1515,968,2944,573)
  
  gene_expr_list <- list()
  for(gid in gene_ids){
    gene_expr_list[[as.character(gid)]] <- getGene(edata,gid)
  }

  safe_sum_or_na <- function(terms_list) {
    # Check if any term in the list is completely NA or has NA for all samples
    # This is for terms that are already vectors of per-sample values
    if(any(sapply(terms_list, function(term) all(is.na(term))))) {
      return(rep(NA, ncol(edata)))
    }
    # Ensure all terms are numeric and have the same length (ncol(edata))
    # If a term became shorter due to an issue, this might error or give bad results
    # For now, assume terms are correctly dimensioned or all NAs
    
    # Summing up, NA + number = NA. So this naturally propagates NAs per sample.
    res <- Reduce(`+`, terms_list)
    return(res)
  }

  safe_product_or_na <- function(coeff, term_vector) {
    # If the term_vector (per-sample values for a gene) is all NA, product is all NA.
    if(all(is.na(term_vector))) {
        return(rep(NA, ncol(edata)))
    }
    # Standard multiplication propagates NAs: NA * number = NA
    return(coeff * term_vector)
  }

  reference_genes_expr_for_avg <- list(
    gene_expr_list[["2597"]],
    gene_expr_list[["2990"]],
    gene_expr_list[["60"]],
    gene_expr_list[["7037"]],
    gene_expr_list[["6175"]]
  )
  
  if(any(sapply(reference_genes_expr_for_avg, function(g) all(is.na(g))))) {
    warning("At least one entire reference gene vector for GHI_RS average is NA or missing. GHI_RS will be NA for all samples.")
    return(rep(NA, ncol(edata)))
  }
  # Calculate average only over samples where all 5 ref genes are non-NA
  reference_Avg <- rowMeans(do.call(cbind, reference_genes_expr_for_avg), na.rm = FALSE) # if any ref gene is NA for a sample, avg is NA for that sample

  refNorm <- function(x_vector, ref_avg_vector) {
    if(all(is.na(x_vector)) || all(is.na(ref_avg_vector))) return(rep(NA, length(x_vector)))
    
    x_norm <- x_vector - ref_avg_vector # NA propagation here is per-sample
    
    # Handle per-sample normalization carefully
    # If a sample's x_norm is NA, it remains NA
    # For non-NA x_norm values for a sample:
    min_val_overall <- min(x_norm, na.rm = TRUE) # overall min of non-NA values
    max_val_overall <- max(x_norm, na.rm = TRUE) # overall max of non-NA values

    if (is.infinite(min_val_overall) || is.infinite(max_val_overall) || (max_val_overall - min_val_overall) == 0) {
        # This case means all non-NA x_norm values are identical, or no non-NA values exist.
        # If they are identical, their scaled value is 0 (or could be set to a fixed point like 7.5 if 0-15 is target)
        # To avoid NaN/Inf, and assuming 0 is a safe baseline for non-varying normalized data:
        x_norm[!is.na(x_norm)] <- 0 
        return(x_norm)
    }
    
    x_norm_scaled <- (x_norm - min_val_overall) * 15 / (max_val_overall - min_val_overall)
    return(x_norm_scaled)
  }
  
  norm_gene_expr_list <- list()
  for(gid_char in names(gene_expr_list)){
    norm_gene_expr_list[[gid_char]] <- refNorm(gene_expr_list[[gid_char]], reference_Avg)
  }

  GRB7_Group <- safe_sum_or_na(list(safe_product_or_na(0.9, norm_gene_expr_list[["2886"]]), 
                                   safe_product_or_na(0.1, norm_gene_expr_list[["2064"]])))
  if(!all(is.na(GRB7_Group))) { # Check if the GRB7_Group vector itself is not all NAs
    GRB7_Group_non_na_indices <- !is.na(GRB7_Group)
    GRB7_Group[GRB7_Group_non_na_indices & GRB7_Group[GRB7_Group_non_na_indices] < 8] <- 8
  }
  
  ER_Group <- safe_sum_or_na(list(
    norm_gene_expr_list[["596"]],
    safe_product_or_na(1.2, norm_gene_expr_list[["5241"]]),
    norm_gene_expr_list[["57758"]],
    safe_product_or_na(0.8, norm_gene_expr_list[["2099"]])
  )) / 4
  
  Prolif_Group <- safe_sum_or_na(list(
    norm_gene_expr_list[["6790"]],
    norm_gene_expr_list[["4605"]],
    norm_gene_expr_list[["891"]],
    norm_gene_expr_list[["332"]],
    norm_gene_expr_list[["4288"]]
  )) / 5
  if(!all(is.na(Prolif_Group))) {
      Prolif_Group_non_na_indices <- !is.na(Prolif_Group)
      Prolif_Group[Prolif_Group_non_na_indices & Prolif_Group[Prolif_Group_non_na_indices] < 6.5] <- 6.5
  }
  
  Invasion_Group <- safe_sum_or_na(list(norm_gene_expr_list[["4320"]], norm_gene_expr_list[["1515"]])) / 2
  
  CD68 <- norm_gene_expr_list[["968"]]
  GSTM1 <- norm_gene_expr_list[["2944"]]
  BAG1 <- norm_gene_expr_list[["573"]]
  
  RSU_components <- list(
    safe_product_or_na(0.47, GRB7_Group),
    safe_product_or_na(-0.34, ER_Group),
    safe_product_or_na(1.04, Prolif_Group),
    safe_product_or_na(0.10, Invasion_Group),
    safe_product_or_na(0.05, CD68),
    safe_product_or_na(-0.08, GSTM1),
    safe_product_or_na(-0.07, BAG1)
  )

  RSU <- safe_sum_or_na(RSU_components)
  RS <- 20*(RSU - 6.7) # This will also propagate NAs correctly
  return(RS)
}

# test association between a signature and a gene using Fisher's exact test
fisherTest <- function(score,CN) {
  
  topQ <- quantile(score,0.75)
  table <- matrix(rep(0,4),nrow = 2,ncol = 2)
  module <- ifelse(score>=topQ,"high","low")
  temp <- data.frame(module=module,CN=CN)
  
  table[1,1] <- length(which(temp$module=="high" & temp$CN=="mut"))
  table[1,2] <- length(which(temp$module=="high" & temp$CN=="wt"))
  table[2,1] <- length(which(temp$module=="low" & temp$CN=="mut"))
  table[2,2] <- length(which(temp$module=="low" & temp$CN=="wt"))
  
  fish <- fisher.test(table,alternative = 'greater')
  return(c(fish$p.value,fish$estimate))
}

# test associations between a signature and all genes, both fisher and spearman correlation/lm
sigCNTest <- function(score,CN_score,CN_gain,CN_loss) {
  # Ensure score is a numeric vector
  if(!is.numeric(score)) score <- as.numeric(score)

  num_samples <- length(score)
  min_finite_pairs <- 4 # Minimum for cor.test

  spearman_pos <- rep(NA, nrow(CN_score))
  spearman_neg <- rep(NA, nrow(CN_score))
  spearman_cor <- rep(NA, nrow(CN_score))
  
  lm_pos <- rep(NA, nrow(CN_score))
  lm_neg <- rep(NA, nrow(CN_score))
  beta_coeff <- rep(NA, nrow(CN_score))
  r.squared <- rep(NA, nrow(CN_score))
  
  gain_fisher_results <- matrix(NA, nrow = nrow(CN_score), ncol = 2, dimnames = list(rownames(CN_score), c("p.value", "estimate")))
  loss_fisher_results <- matrix(NA, nrow = nrow(CN_score), ncol = 2, dimnames = list(rownames(CN_score), c("p.value", "estimate")))

  for(j in 1:nrow(CN_score)){
    current_gene_name <- rownames(CN_score)[j]
    CN <- unname(unlist(CN_score[j,]))
    if(!is.numeric(CN)) CN <- as.numeric(CN)

    finite_pairs <- sum(is.finite(score) & is.finite(CN))
    
    if (finite_pairs < min_finite_pairs) {
      # Results for this gene remain NA for cor/lm
    } else {
      # Spearman rank correlation (suppress exact p-value calculation warnings for large N)
      spearman_pos[j] <- tryCatch(cor.test(score, CN, method = "spearman", alternative = "greater", exact = FALSE)$p.value, error = function(e) NA)
      spearman_neg[j] <- tryCatch(cor.test(score, CN, method = "spearman", alternative = "less", exact = FALSE)$p.value, error = function(e) NA)
      if (finite_pairs >= 2) { # spearman cor itself needs at least 2 pairs
          spearman_cor[j] <- tryCatch(cor(score, CN, method = "spearman", use = "pairwise.complete.obs"), error = function(e) NA)
      }
      
      # Linear model: Check for subtype_variable (THIS IS A CRITICAL DEPENDENCY)
      if (exists("subtype_variable", envir = .GlobalEnv) && 
          is.data.frame(subtype_variable) && 
          all(c("basal", "her2", "lumA", "lumB") %in% names(subtype_variable)) &&
          nrow(subtype_variable) == num_samples) {
        
        lm_data <- data.frame(score = score, CN = CN, 
                            basal = subtype_variable$basal, 
                            her2 = subtype_variable$her2, 
                            lumA = subtype_variable$lumA, 
                            lumB = subtype_variable$lumB)
        
        # Remove rows with NAs for lm explicitly to ensure model runs if possible
        lm_data_complete <- na.omit(lm_data)

        if(nrow(lm_data_complete) >= min_finite_pairs + 4) { # Need enough data after NA removal for lm with 5 predictors
            fit_res <- tryCatch(lm(score ~ CN + basal + her2 + lumA + lumB, data = lm_data_complete), error = function(e) NULL)
            if(!is.null(fit_res)){
                sum_fit <- summary(fit_res)
                if ("CN" %in% rownames(sum_fit$coefficients)) {
                    beta <- sum_fit$coefficients["CN", "Estimate"]
                    p_val <- sum_fit$coefficients["CN", "Pr(>|t|)"]
                    beta_coeff[j] <- beta
                    r.squared[j] <- sum_fit$r.squared
                    if(!is.na(beta) && !is.na(p_val)){
                        if(beta > 0) {
                            lm_pos[j] <- p_val / 2
                            lm_neg[j] <- 1 - (p_val / 2)
                        } else if(beta < 0) {
                            lm_pos[j] <- 1 - (p_val / 2)
                            lm_neg[j] <- p_val / 2
                        }
                    }
                }
            }
        } # else not enough data for lm after NA removal
      } else {
         # warning(paste0("Subtype_variable not available or misconfigured for gene: ", current_gene_name, ". Skipping lm."))
      }
    }
    
    # Fisher's exact test
    if(sum(is.finite(score)) > 0 && sum(is.finite(CN_gain[j,])) > 0 ){
        cn_gain_factor <- factor(ifelse(CN_gain[j,], "mut", "wt"), levels = c("mut", "wt"))
        fisher_res <- tryCatch(fisherTest(score, cn_gain_factor), error = function(e) { c(NA, NA) })
        gain_fisher_results[j, ] <- fisher_res
    } else {
        gain_fisher_results[j, ] <- c(NA,NA)
    }

    if(sum(is.finite(score)) > 0 && sum(is.finite(CN_loss[j,])) > 0){
        cn_loss_factor <- factor(ifelse(CN_loss[j,], "mut", "wt"), levels = c("mut", "wt"))
        fisher_res <- tryCatch(fisherTest(score, cn_loss_factor), error = function(e) { c(NA, NA) })
        loss_fisher_results[j, ] <- fisher_res
    } else {
        loss_fisher_results[j, ] <- c(NA,NA)
    }
  }
  
  p_values_df <- data.frame(
    spearman_pos = spearman_pos,
    spearman_neg = spearman_neg, 
    CN_gain = gain_fisher_results[,1], 
    CN_loss = loss_fisher_results[,1],
    lm_pos = lm_pos,
    lm_neg = lm_neg,
    spearman_correlation = spearman_cor,
    lm_beta = beta_coeff,
    lm_r_squared = r.squared,
    OR_gain = gain_fisher_results[,2],
    OR_loss = loss_fisher_results[,2]
  )
  rownames(p_values_df) <- rownames(CN_score)
  return(p_values_df)
}

calc_segments<-function(x, gmtFile, method="mean", scale=F, gsaObj=NA){
  
  geneset.obj<- GSA.read.gmt(gmtFile)
  genenames<-row.names(x)
  np=length(geneset.obj$genesets)
  
  if(scale){xs=t(scale(t(x),center=T,scale=T))}else{xs<-x}
  
  if(method!="gsa"){
    val=matrix(0,nrow=np,ncol=ncol(x))
    for(i in 1:np){
      gene.set=which(genenames %in% geneset.obj$genesets[[i]])
      gene.set=gene.set[!is.na(gene.set)]
      if(length(gene.set)>1){
        if(method=="mean"){
          val[i,]=colSums(xs[gene.set,,drop=F])/length(gene.set)
        }
        if(method=="median"){
          val[i,]=t(apply(xs[gene.set,,drop=F],2,median))
        }
        if(method=="pca"){
          y<-prcomp(as.matrix(xs[gene.set,,drop=F]))
          val[i,]=y$rotation[,1]
        }
      } else if (length(gene.set) == 1) {
        val[i,] <- unlist(xs[gene.set,])
      }
    }
  }else{
    val<-GSA.make.features(gsaObj,xs,geneset.obj$genesets,genenames)
  }
  
  dimnames(val)<-list(geneset.obj$geneset.names,dimnames(x)[[2]])
  return(val)
}

medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}

standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}

exp_wrap <- function(edata) {
  keep <- c('29126','1493','5133') # for PD1, PDL1 and CTLA4
  edata70 <- edata[rowSums(edata<=2)<(0.3*ncol(edata)),]
  for(i in 1:3) {
    if(!(keep[i] %in% rownames(edata70))) {
      edata70 <- rbind(edata70,edata[keep[i],])
    }
  }
  edata70[edata70<=2] <- 0  # missing data marked with 0
  edata70log2 <- log(edata70,base = 2)
  edata70log2[edata70log2=="-Inf"] <- 0
  
  exp <- medianCtr(edata70log2)
  exp <- standardize(exp)
  return(exp)
}

caret_wrap <- function(trainX,trainY,testX,testY,bi=T) {
  if(!bi) {
    # set cross validation resampling method
    train_control <- trainControl(method = 'LGOCV',number = 200,classProbs = F)
    
    # set alpha and lambda grid
    alpha <- seq(0.1,0.9,by=0.1)
    lambda <- list()
    for(i in 1:9) {
      init <- glmnet(trainX,trainY,alpha = alpha[i])
      lambda[[i]] <- init$lambda
    }
    lambda_min <- min(unlist(lapply(lambda,min)))
    lambda_max <- max(unlist(lapply(lambda,max)))
    
    tune_grid = expand.grid(alpha = seq(0.1,0.9,by=0.1), lambda = seq(lambda_min,lambda_max,length.out = 100))
    
    # train model
    glmnet_obj <- train(trainX, trainY, method = "glmnet", metric = "RMSE",
                        trControl = train_control,tuneGrid = tune_grid)
    }  else {
    trainY2 <- ifelse(trainY==1,'pos','neg')
    
    # set cross validation resampling method
    train_control <- trainControl(method = 'LGOCV',number = 200,classProbs = T)
    
    alpha <- seq(0.1,0.9,by=0.1)
    lambda <- list()
    for(i in 1:9) {
      init <- glmnet(trainX,trainY,alpha = alpha[i],family = 'binomial')
      lambda[[i]] <- init$lambda
    }
    lambda_min <- min(unlist(lapply(lambda,min)))
    lambda_max <- max(unlist(lapply(lambda,max)))
    
    tune_grid = expand.grid(alpha = seq(0.1,0.9,by=0.1), lambda = seq(lambda_min,lambda_max,length.out = 100))
    
    # train model
    glmnet_obj <- train(trainX, trainY2, method = "glmnet", metric = "Accuracy",
                        trControl = train_control,tuneGrid = tune_grid)
    }
  return(glmnet_obj)
}

plot_ROC <- function(perf1,perf2,a1,a2,main) {
  tiff(paste0(main,'.tiff'),width = 1.5,height = 1.5,units = 'in',res = 300)
  par(mai = c(0.2,0.2,0.05,0.05),cex.axis = 0.3)
  plot(perf1,col = 'red',lwd = 1.2,xlab = "",ylab = "",box.lwd=0.8,
       xaxis.xaxt = 'n',yaxis.yaxt = 'n')
  plot(perf2,col = 'darkblue',lwd = 1.2,add = T)
  abline(a=0,b=1,lwd = 0.8)
  axis(1,tck=(-0.02),lwd = 0.8)
  axis(2,tck=(-0.02),lwd = 0.8)
  mtext(side = 1,text = seq(0,1,by=0.2),at = seq(0,1,by=0.2),cex = 0.5,line = (-0.2))
  mtext(side = 2,text = seq(0,1,by=0.2),at = seq(0,1,by=0.2),cex = 0.5,line = 0.1)
  
  legend('bottomright',legend = c(paste0('AUC = ',a1),paste0('AUC = ',a2)), lty = c(1,1),lwd = c(1,1) ,
         col = c('red','darkblue'),cex = 0.6,bty = 'n')
  dev.off()
}

x_y <- function(beta,segment_anno) {
  x <- c()
  y <- c()
  for(i in 1:nrow(segment_anno)) {
    this_seg <- rownames(segment_anno)[i]
    this_start <- segment_anno[i,2]
    this_end <- segment_anno[i,3]
    this_x <- seq(this_start,this_end,by = 1)
    this_y <- rep(beta[this_seg],length(this_x))
    x <- c(x,this_x)
    y <- c(y,this_y)
  }
  return(data.frame(x=x,y=y))
}

# plot model feature for single signature
plot_seg_ss <- function(beta,main) {
  total_gene <- 24776
  vertical <- vertical_lines[-c(1,24)]
  text_pos <- vertical[1]/2
  for(i in 2:length(vertical)){
    thispos <- vertical[i] - (vertical[i] - vertical[i-1])/2
    text_pos <- rbind(text_pos,thispos)
  }
  text_pos <- rbind(text_pos,total_gene-(total_gene-vertical[22])/2)
  
  min_y <- min(beta)
  max_y <- max(beta)
  # beta whole arm
  index <- grep('wholearm',names(beta))
  beta_wholearm <- beta[index]
  beta <- beta[-index]
  
  pos_beta <- beta[beta>0]
  neg_beta <- beta[beta<0]
  
  pos_seg <- names(pos_beta)
  neg_seg <- names(neg_beta)
  
  pos_anno <- segment_anno[pos_seg,]
  neg_anno <- segment_anno[neg_seg,]
  
  # pos regions excluding whole arms
  pos_coor <- x_y(beta,pos_anno)
  neg_coor <- x_y(beta,neg_anno)
  
  x <- c(pos_coor$x,neg_coor$x)
  y <- c(pos_coor$y,neg_coor$y)
  color <- c(rep('red',nrow(pos_coor)),rep('darkblue',nrow(neg_coor)))
  
  # whole arm beta
  beta_pos_wholearm <- beta_wholearm[beta_wholearm>0]
  beta_neg_wholearm <- beta_wholearm[beta_wholearm<0]
  
  if( length(beta_pos_wholearm) == 0) {
    pos_wholearm_coor <- data.frame(x=0,y=0)
  } else {
    pos_wholearm_anno <- segment_anno[names(beta_pos_wholearm),]
    pos_wholearm_coor <- x_y(beta_wholearm,pos_wholearm_anno)
  }
  
  if(length(beta_neg_wholearm)==0) {
    neg_wholearm_coor <- data.frame(x=0,y=0)
  } else {
    neg_wholearm_anno <- segment_anno[names(beta_neg_wholearm),]
    neg_wholearm_coor <- x_y(beta_wholearm,neg_wholearm_anno)
  }
  x_wholearm <- c(pos_wholearm_coor$x,neg_wholearm_coor$x)
  y_wholearm <- c(pos_wholearm_coor$y,neg_wholearm_coor$y)
  color_wholearm <- c(rep('pink',nrow(pos_wholearm_coor)),rep('lightblue',nrow(neg_wholearm_coor)))
  par(cex.axis = 2)
  png(filename = paste(main,'.png',sep = ''),width = 28,height = 5,res = 72,units = 'in')
  par(cex.axis = 2,cex.lab = 2.5,mai=c(0.6,1.5,0.6,0.5))
  plot(x_wholearm,y_wholearm,type = 'h',col = color_wholearm,xlim=c(1,24776),ylim=c(min_y,max_y),lwd = 1,xaxs="i",xaxt = 'n',ylab = '',xlab = "")
  abline(v=vertical_lines,lwd = 1)
  abline(h = 0,lwd = 1)
  lines(x,y,xaxt = "n",yaxt = "n",col = color,type = 'h')
  mtext(c(1:20,"",22,'x'),side = 1,at = text_pos,line = 1.5,cex = 2.5)
  mtext('weight',side = 2,at = 0,line = 4,cex = 2.5)
  dev.off()
}

P_Plot <- function(p,main,vertical_lines = v,y1 = NULL,y2 = NULL,label_up = NULL,label_down = NULL) {
  pos <- p[,1]
  neg <- p[,2]
  gain <- p[,3]
  loss <- p[,4]
  total_gene <- nrow(p)
  
  # remove non-significant lines
  pos_index <- which(!(pos>2 & gain>2))
  neg_index <- which(!(neg>2 & loss>2))
  pos[pos_index] <- 0
  gain[pos_index] <- 0
  neg[neg_index] <- 0
  loss[neg_index] <- 0
  
  m1 <- (max(pos,gain))
  m2 <- (max(neg,loss))
  
  # reformat log_p values for better plot
  index <- which(gain>pos)
  pos2 <- rep(0,total_gene)
  pos2[index] <- pos[index]
  
  index <- which(loss > neg)
  neg2 <- rep(0,total_gene)
  neg2[index] <- neg[index]
  
  if (is.null(y1) & is.null(y2)) {
    # init plot to get suitable axis labels
    print(paste('m1 =',m1,'; m2 =',m2))
  } else {
    vertical <- vertical_lines[-c(1,24)]
    text_pos <- vertical[1]/2
    for(i in 2:length(vertical)){
      thispos <- vertical[i] - (vertical[i] - vertical[i-1])/2
      text_pos <- rbind(text_pos,thispos)
    }
    text_pos <- rbind(text_pos,total_gene-(total_gene-vertical[22])/2)
    
    png(filename = paste(main,'.png',sep = ''),width = 28,height = 5,res = 72,units = 'in')
    par(cex.axis = 2,font.lab = 2)
    #layout(matrix(c(rep(1,2),2),ncol = 1))
    
    layout(matrix(c(1,2),ncol = 1),heights = c(y1+y2/6,y2+y1/6))
    
    par(mai=c(0,1.5,0.6,0.5),lwd = 2)
    plot(c(1:total_gene),pos,type = "h", xaxt = "n",yaxt = "n",lwd=2,xlim = c(1,24776),xaxs="i",yaxs='i',bty = '7',
         ylim = c(0,y1),col = "red",ylab = NA,xlab = NA)
    lines(c(1:total_gene),gain,xaxt = "n",yaxt = "n",col = "orange",lwd=2,type = 'h')
    lines(c(1:total_gene),pos2,xaxt = "n",yaxt = "n",col = "red",lwd=2,type = 'h')
    axis(2,labels = label_up,at = label_up,las = 1)
    abline(v = vertical_lines,lwd=1)
    abline(h = 2,lwd = 1,lty = 2)
    
    par(mai = c(0.6,1.5,0,0.5),lwd = 2)
    plot(c(1:total_gene),neg,type = "h", xaxt = "n",yaxt = "n",lwd=2,xlim = c(1,24776),xaxs="i",yaxs='i',bty = 'u',
         ylim = c(y2,0),col = "darkblue",ylab = NA,xlab = NA)
    lines(c(1:total_gene),loss,xaxt = "n",yaxt = "n",col = "cornflowerblue",lwd=2,type = 'h')
    lines(c(1:total_gene),neg2,xaxt = "n",yaxt = "n",col = "darkblue",lwd=2,type = 'h')
    abline(v = vertical_lines,lwd=1)
    abline(h = 2,lwd = 1,lty = 2)
    axis(2,labels = label_down,at = label_down,las = 1)
    
    mtext(c(1:20,"",22,'x'),side = 1,at = text_pos,line = 1.5,cex = 2.5)
    #mtext(0,side = 2,at = 0,line = 1,cex = 2)
    mtext(expression("-Log"[10]*'q'),side = 2,at = 0,line = 4,cex = 2.5)
    
    legend("bottom",fill = c("red","blue","orange","cornflowerblue"),
           legend = c(expression("-log"[10]*'q'[REGRESSION]*'(positive correlation)'),
                      expression("-log"[10]*'q'[REGRESSION]*'(negative correlation)'),
                      expression("-log"[10]*'q'[FISHER]*'(gain)'),
                      expression("-log"[10]*'q'[FISHER]*'(loss)')),bty = 'n',
           cex = 1, ncol = 2 )
    dev.off()
    
  }
}
