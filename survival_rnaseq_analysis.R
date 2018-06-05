library(survival)
library(limma)

# input counts, and filter genes whose expression is zero in more than half the samples
rna <- as.matrix(read.table("KICH.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header=T, row.names=1, sep="\t"))
rna <- rna[-1,]
x <- t(apply(rna,1,as.numeric))
r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
remove <- which(r > dim(rna)[2]*0.5)
rna <- rna[-remove,]

# get the index of the normal/control samples
n_index <- which(substr(colnames(rna),14,14) == "1")
t_index <- which(substr(colnames(rna),14,14) == "0")

# voom normalization
cond <- factor(ifelse(seq(1, dim(rna)[2],1) %in% t_index, 1,  0))
d <- model.matrix(~1 + cond)
x <- t(apply(rna, 1, as.numeric))
rna_vm <- voom(x, d, plot=FALSE)$E
colnames(rna_vm) <- gsub("\\.","-",substr(colnames(rna),1,12))

# zscore scaling against normals
mean_n <- rowMeans(rna_vm[, n_index])
sd_n <- apply(rna_vm[, n_index], 1, sd)
tumor <- rna_vm[,t_index]
z_rna <- matrix(nrow=nrow(tumor), ncol=ncol(tumor))
colnames(z_rna) <- colnames(rna_vm[,t_index])
rownames(z_rna) <- rownames(rna_vm[,t_index])
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,"\\|"))[[1]])
for(i in 1:nrow(tumor)){
 for(j in 1:ncol(tumor)){
  z_rna[i,j] <- (tumor[i, j] - mean_n[i]) / sd_n[i]
 }
}


clinical <- t(read.table('KIPAN.merged_only_clinical_clin_format.txt',header=T, row.names=1, sep='\t'))
clinical <- as.data.frame(clinical)
clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
sum(clinical$IDs %in% colnames(z_rna))
ind_keep <- grep('days_to_new_tumor_event_after_initial_treatment',colnames(clinical))

new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- min(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,'NA')
  }
}

# do the same to death
ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c('new_tumor_days', 'death_days', 'followUp_days')

# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                                   as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}

# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}


# create vector for death censoring
table(clinical$patient.vital_status)
# alive dead
# 372   161

all_clin$death_event <- ifelse(clinical$patient.vital_status == 'alive', 0,1)

#finally add row.names to clinical
rownames(all_clin) <- clinical$IDs

write.table(z_rna,"KICH.txt",sep="\t")
write.table(all_clin,"all_clin.txt",sep="\t")

