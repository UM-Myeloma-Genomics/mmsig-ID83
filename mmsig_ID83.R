#MMsig modification for ID83 signatures

#load in test data, require ID83 matrix format
###for easy matrix creation you can use SigProfiler matrix generator output .all file

setwd("~/mmsig_ID83/")
mtx.def<- read.delim("test_matrix.txt")

#load reference signature file 
cosmID<- read.table("~/COSMIC_v3.2_ID_GRCh37.txt", sep = "\t", stringsAsFactors = F, header = T)
colnames(cosmID)[1]<- "Type"

#make user-supplied signature catalogue for fitting
indel_list<- list()
for(i in (2:ncol(mtx.def))){
  
  signature_pos<- c("ID1", "ID2", "ID3", "ID8", "ID9") 
  indel_list[[paste(colnames(mtx.def)[i])]]<- c(paste(colnames(mtx.def)[i]) , signature_pos)
}    

#load in mmsig for modification #here 

library(mmsig) 

source("~/mmsig-master/R/vcfcleaner.R")
source("~/mmsig-master/R/refcheck.R")
source("~/mmsig-master/R/fit_signatures.R")
source("~/mmsig-master/R/spit.R")
source("~/mmsig-master/R/spat.R")
source("~/mmsig-master/R/em_signatures.R")
source("~/mmsig-master/R/bootstrap.R")

###generate variables
samples.muts<- mtx.def[,-1]
tcons.defn <- t(cosmID[, 2:ncol(cosmID)])
mutlist = cosmID[,1] 
consigts.defn <- sapply(cosmID[, 2:ncol(cosmID)],  
                        as.numeric) 
samples <- colnames(samples.muts)

rownames(consigts.defn) <- mutlist
ref_signatures <- colnames(consigts.defn)
sample.sigt.profs<- indel_list
sigt.profs <- sample.sigt.profs #need this arg for user-defined sig catalogue

#####constants
cos_sim_threshold = 0.01
force_include<- c("ID1", "ID2") #must include clock-like sigs
dbg = FALSE

####################
###Alter em_signatures to take 83 class matrix
####################
em_signatures.ind<-
  function(sigts.defn, mut.freqs, max.iter, dbg) {
    nof.sigts.defn <- ncol(sigts.defn)
    alpha = stats::runif(nof.sigts.defn); alpha=alpha/sum(alpha) 
    for (i in 1:max.iter) {
      contr = t(array(alpha, dim=c(nof.sigts.defn,83))) * sigts.defn
      probs = contr/array(rowSums(contr), dim=dim(contr))
      probs[is.na(probs)] = 0
      probs = probs * mut.freqs
      old_alpha = alpha
      alpha = colSums(probs)/sum(probs)
      if (sum(abs(alpha-old_alpha))<1e-5) {
        break
      }
    }
    spit(dbg, "em: exit iteration: %d", i)
    
    return( alpha )
  }
####################
#rewrite bootstrap for indels
#####################

samples.muts<- mtx.def 
rownames(samples.muts)<- samples.muts[,1]
samples.muts<- samples.muts[,-1] 

bootstrap_fit_signatures_ID <- function(samples.muts,
                                        consigts.defn,
                                        sigt.profs,
                                        cos_sim_threshold,
                                        force_include,
                                        iterations = 1000){
  
  # Setup
  samples <- names(samples.muts)
  classes <- rownames(samples.muts)
  
  # List to populate with signatures
  mutSigs <- list()
  
  # Setup progress bar if desired
  #pbar <- create_progress_bar('text')
  #pbar$init(length(samples))
  
  for(i in 1:length(samples)){
    # Loop through samples, generating a data frame of signature contributions for each
    sub <- as.integer(samples.muts[classes,i])
    total <- sum(sub)
    # sample new 83-classes profiles from the multinomial distribution
    bootMat <- data.frame(rmultinom(n = iterations, size = total, prob = sub/total))
    row.names(bootMat) <- classes
    
    # prepare the signatures to fit for each sample
    sig.prof <- list()
    
    for(s in 1:ncol(bootMat)){
      sig.prof[[names(bootMat)[s]]] <- sigt.profs[[samples[i]]]
    }
    
    ### Run mmSig
    sig_out <- fit_signatures.ind(samples.muts=bootMat,
                                  consigts.defn=consigts.defn,
                                  sigt.profs=sig.prof,
                                  cos_sim_threshold=cos_sim_threshold,
                                  force_include=force_include,
                                  dbg=FALSE)
    
    mutSigs[[i]] <- sig_out
    #pbar$step()
  }
  
  names(mutSigs) <- samples
  
  # Generate final summary data frame
  
  ## Summary statistics   ####here
  my_summary <- function(x){
    c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
  }
  
  mutSigsSummary <- list()
  for(i in 1:length(mutSigs)){
    s <- names(mutSigs)[i]
    temp <- mutSigs[[i]]
    out <- data.frame(t(sapply(temp[names(temp) != "mutations"], my_summary)))
    names(out) <- c('mean', 'CI025', 'CI975')
    out$signature <- row.names(out)
    out$sample <- s
    out <- out[c('sample', 'signature', 'mean', 'CI025', 'CI975')]
    mutSigsSummary[[i]] <- out
  }
  
  mutSigsSummary <- bind_rows(mutSigsSummary)
  
  return(mutSigsSummary)
}

####################
#Alter fit_signatures to take 83 class matrix
####################
fit_signatures.ind<-
  function(samples.muts,
           consigts.defn,
           sigt.profs,
           cos_sim_threshold,
           force_include,
           dbg=dbg) {
    
    #start running line by line here
    max.em.iter=2000
    consigt.names <- colnames(consigts.defn)
    samples <- colnames(samples.muts)
    
    # Mutational signature fitting procedure for each individual sample  
    sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
    rownames(sigt.fraction) <- consigt.names
    colnames(sigt.fraction) <- samples
    
    for (j in 1:length(samples)) {
      sample.mut.freqs = as.numeric(samples.muts[,j])
      sample.mut.freqs[is.na(sample.mut.freqs)] = 0
      sample.sigts <- unique(sigt.profs[[ samples[j] ]])
      sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
      sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
      spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
      alpha <- em_signatures.ind(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      spat(dbg, "alpha", alpha)
      sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
      sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]
      
      if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
      
      spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
      reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs)
      sample.cos.sim.meas <- MutationalPatterns::cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
      spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
      
      rem.alpha <- sampleAlpha                     # holds the final result
      rem.sample.consigts.defn <- sample.consigts.defn
      spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
      reducing = TRUE
      
            # Signature profile shrinkage by cosine similarity (removing signatures that are not necessary to explain profile)
      while (reducing) {
        spat(dbg, "in the while, rem.alpha: ", rem.alpha)
        cosReduction <- NULL
        rem.names <- setdiff(names(rem.alpha), force_include)
        if(length(rem.names) == 0){ ## Avoiding script crash when only the forced signatures are present.
          spit(dbg, "removed all signatures except forced inclusion: exiting while...")
          break
        }
        
        for(c in rem.names){
          spit(dbg, "doing c: %s", c)
          red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
          red.alpha <- em_signatures.ind(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
          red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
          red.cos.sim.meas <- MutationalPatterns::cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
          cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
        }
        names(cosReduction) <- rem.names
        if (min(cosReduction) < cos_sim_threshold) {
          spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
          rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
          rem.alpha <-  em_signatures.ind(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
          reducing = TRUE
        }
        else {
          spit(dbg, "exiting while...")
          reducing = FALSE
        }
      }
      
   
      
      spit(dbg,"... while exited")
      rem.alpha.names <- names(rem.alpha)
      
      for (n in 1:length(consigt.names)) {
        if (consigt.names[n] %in% rem.alpha.names) {
          sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
        }
        else {
          sigt.fraction[n,j] <- 0
        }
      }
    }
    spat(dbg, "sigt.fraction", sigt.fraction)
    tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
    colsums.samples.muts <- colSums(samples.muts)
    sig <- cbind(tdf.sigt.fraction, "indels"=colsums.samples.muts)
    
    
    return(sig)
  }

###end of new function definition


############
###this separates signature and bootstrapping, but can skip in favor of the complete function following
#############
output <- list()

samples.muts.nozero<-samples.muts[,colSums(samples.muts)!=0] #need to do this to remove 0 samples which throw an error
samples.muts<- samples.muts.nozero

sigfit <- fit_signatures.ind(samples.muts = samples.muts, consigts.defn = consigts.defn, 
                             sigt.profs = sigt.profs, cos_sim_threshold = 0.01, 
                             force_include = c("ID1", "ID2"), dbg = FALSE)


bootID<- bootstrap_fit_signatures_ID(samples.muts = samples.muts.nozero, consigts.defn = consigts.defn, 
                                     sigt.profs = sigt.profs, cos_sim_threshold = 0.1, 
                                     force_include = c("ID1", "ID2"), iterations = 1000) #removed progress bar from function, it still works, 0 indel samples gone

#####
######
######

################
#####define complete fitting function
###############
mm_fit_indels<-  function (muts.input, sig.input, input.format = "vcf", sample.sigt.profs = NULL, 
                           bootstrap = FALSE, iterations = 1000, strandbias = FALSE, 
                           refcheck = TRUE, cos_sim_threshold = 0.001, force_include = c("ID1", 
                                                                                         "ID2"), dbg = FALSE) 
{
  "\n\n    "
  options(scipen = 999)
  consigts.defn <- sig.input
  samples.muts <- samples.muts #using sigprofiler matrix for conversion
  tcons.defn <- tcons.defn
  mutlist = mutlist
  #consigts.defn <- sapply(consigts.defn[, 1:ncol(consigts.defn)], 
  #                       as.numeric)
  # rownames(consigts.defn) <- mutlist
  ref_signatures <- colnames(consigts.defn)
  samples <- colnames(samples.muts)
  if (is.list(sample.sigt.profs)) {
    spit(dbg, "using mm signature profiles from input argument")
    sigt.profs <- sample.sigt.profs
  }
  else {
    spit(dbg, "defaulting to use all signatures in the provided reference")
    mm.sigts <- ref_signatures
    sigt.profs <- list()
    for (i in 1:length(samples)) {
      sigt.profs[[samples[i]]] <- mm.sigts
    }
  }
  output <- list()
  sigfit <- fit_signatures.ind(samples.muts = samples.muts, consigts.defn = consigts.defn, 
                               sigt.profs = sigt.profs, cos_sim_threshold = cos_sim_threshold, 
                               force_include = force_include, dbg = dbg)
  output$estimate <- sigfit
  if (bootstrap) {
    sig_est <- sigfit
    sig_est$sample <- row.names(sig_est)
    sig_est <- melt(sig_est[names(sig_est) != "indels"], 
                    id.vars = "sample", variable.name = "signature", 
                    value.name = "estimate")
    sig_est$signature <- as.character(sig_est$signature)
    sigboot <- bootstrap_fit_signatures_ID(samples.muts = samples.muts, 
                                           consigts.defn = consigts.defn, sigt.profs = sigt.profs, 
                                           iterations = iterations, cos_sim_threshold = cos_sim_threshold, 
                                           force_include = force_include)
    sigboot <- left_join(sigboot, sig_est, by = c("sample", 
                                                  "signature"))
    output$bootstrap <- sigboot
  }
  if (strandbias & input.format == "vcf") {
    strand_bias_out <- getStrandBias(muts.input)
    output$strand_bias_all_3nt <- strand_bias_out$all_3nt
    output$strand_bias_mm1 <- strand_bias_out$mm1
    output$strand_bias_SBS35 <- strand_bias_out$SBS35
  }
  else if (strandbias & input.format != "vcf") {
    warning("Transcriptional strand bias cannot be estimated from 96 classes input. \n Please provide vcf-like input format.")
  }
  output$mutmatrix <- samples.muts
  return(output)
}

#####run on test data

indelfit<- mm_fit_indels(muts.input = samples.muts, sig.input = consigts.defn, sample.sigt.profs = sigt.profs, 
                         bootstrap = TRUE, iterations = 1000, strandbias = FALSE, 
                         cos_sim_threshold = 0.01, force_include = c("ID1", "ID2"), dbg = FALSE) 
