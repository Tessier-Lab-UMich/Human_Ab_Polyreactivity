getwd()
setwd("/Users/htdennis/Desktop/rstudio")
install.packages("doBy")
install.packages("readxl")
library("readxl")
library(doBy)
library(dplyr)
library(randomForest)

#
# Data cleanup
#
drop_empty_columns <- function(data)
{
  L_names <- grep("^L", names(data), value=TRUE)
  H_names <- grep("^H", names(data), value=TRUE)
  
  col_filt <- c()
  for(L in L_names)
  {
    t_noins <- nrow(data) - length(which(data[[L]] == '-'))
    
    if(t_noins == 0)
    {
      col_filt <- c(col_filt, L)
    }
  }
  
  for(H in H_names)
  {
    t_noins <- nrow(data) - length(which(data[[H]] == '-'))
    
    if(t_noins == 0)
    {
      col_filt <- c(col_filt, H)
    }
  }
  
  data <- data[, !names(data) %in% col_filt]
  
  return(data)
}

#
# Normalizing SASA (A^2) values
#
get_norm_sasa <- function()
{
  SASA <- c(64.7809
            , 210.02
            , 95.2439
            , 143.924
            , 186.7
            , 23.1338
            , 149.451
            , 151.242
            , 177.366
            , 139.524
            , 164.674
            , 111.533
            , 147.855
            , 229.619
            , 200.306
            , 113.187
            , 110.209
            , 81.2159
            , 111.597
            , 124.237
            , 0);
  
  AA <- c('A','R','C','E','F','G','H','I','K','L','M','P','Q','W','Y','N','D','S','T','V','-')
  
  return(data.frame(AA=AA, Norm_SASA=SASA))
}

#
# Definition of CDR residues
#
cdr_residues <- function(hc='H', lc='L')
{
  Residue <- c("L24","L25","L26","L27","L28","L29","L30","L30A","L30B","L30C","L30D","L30E","L30F","L31","L32","L33","L34",
               "L50","L51","L52","L52A","L52B","L52C","L52D","L52E","L53","L54","L55","L56",
               "L89","L90","L91","L92","L93","L94","L95","L95A","L95B","L95C","L95D","L95E","L96","L97",
               "H26","H27","H28","H29","H30","H31","H31A","H31B","H32","H33","H34","H35",
               "H50","H51","H52","H52A","H52B","H52C","H53","H54","H55","H56","H57","H58","H59","H60","H61","H62","H63","H64","H65",
               "H93","H94","H95","H96","H97","H98","H99","H99A","H99B","H99C","H99D","H99E","H99F","H99G","H99H","H99I","H99J","H99K","H99L","H100N","H100O","H100P","H100Q","H100R","H100S","H100T","H100U","H100V","H100W","H100X","H100Y","H100Z","H100","H101","H102")
  
  Residue[1:43] <- sapply(Residue[1:43], function(x) paste(lc,substr(x,2,100),sep=''))
  Residue[44:109] <- sapply(Residue[44:109], function(x) paste(hc,substr(x,2,100),sep=''))
  
  d <- data.frame(Residue=Residue)
  d$Residue <- as.character(d$Residue)
  
  d$CDR <- 'H3'
  d$CDR[1:17] <- 'L1'
  d$CDR[18:29] <- 'L2'
  d$CDR[30:43] <- 'L3'
  d$CDR[44:55] <- 'H1'
  d$CDR[56:74] <- 'H2'
  
  return(d)
}

#
# Read toolkit sequence output and prepare in a form needed to
# predict SASA
#
combine_vh_vl_data <- function(vh, vl)
{
  VH <- read.table(vh, header=T, colClasses='character')
  VL <- read.table(vl, header=T, colClasses='character')
  
  names(VH)[3:5] <- c('h1_len','h2_len','h3_len')
  names(VL)[3:5] <- c('l1_len','l2_len','l3_len')
  
  sequence <- merge(VH, VL, by="Name")
  sequence$h1_len <- as.integer(sequence$h1_len)
  sequence$h2_len <- as.integer(sequence$h2_len)
  sequence$h3_len <- as.integer(sequence$h3_len)
  
  sequence$l1_len <- as.integer(sequence$l1_len)
  sequence$l2_len <- as.integer(sequence$l2_len)
  sequence$l3_len <- as.integer(sequence$l3_len)
  
  sequence$VHF <- substr(sequence$VH, 1, 3)
  sequence$VLF <- substr(sequence$VL, 1, 3)
  
  sequence$Kappa <- FALSE
  sequence$Kappa[grep("VK", sequence$VL)] <- TRUE
  sequence$germline.species <- 'human'
  
  if('PID.VH' %in% names(sequence))
  {
    sequence$PID.VH <- as.numeric(as.character(sequence$PID.VH))
    sequence$PID.VL <- as.numeric(as.character(sequence$PID.VL))
  }
  
  names(sequence)[1] <- "PreCADName"
  
  sequence <- drop_empty_columns(sequence)
  
  return(sequence)
}

#
# Get the new data column types in sync with the types used to build the models
# Also identify those PreCADs that cannot be used since they have sequence positions not
# in the data used to fit the models
#
sync_data_classes <- function(sequence, fit_sequence)
{
  LH_names_in_new <- names(sequence)[c(grep("^L", names(sequence)),grep("^H", names(sequence)))]
  LH_names_in_fit <- names(fit_sequence)[c(grep("^L", names(fit_sequence)),grep("^H", names(fit_sequence)))]
  
  addn_names_in_fit <- setdiff(LH_names_in_fit, LH_names_in_new)
  
  common_LH_names <- intersect(LH_names_in_new, LH_names_in_fit)
  for(N in common_LH_names)
  {
    sequence[[N]] <- factor(as.character(sequence[[N]]), levels=levels(fit_sequence$AA))
    sequence[[N]] <- as.numeric(sequence[[N]])
  }
  
  for(N in addn_names_in_fit)
  {
    sequence[[N]] <- 1	# Corresponding to insertions
  }
  
  new_names <- setdiff(LH_names_in_new, LH_names_in_fit)
  omitted_precads <- c()
  for(N in new_names)
  {
    omitted_precads <- unique(c(omitted_precads, sequence$PreCADName[sequence[[N]]!='-']))
  }
  
  sequence <- subset(sequence, !PreCADName %in% omitted_precads)
  sequence <- sequence[, !names(sequence) %in% new_names]
  
  if(length(new_names)>0)
  {
    print(paste('Following PreCADS have positions not in fit data ', paste(omitted_precads, collapse=',')))
    print(paste('Following positions seen in PreCADS without models', paste(new_names, collapse=',')))
    print(paste('No SASA predictions will be made for these.'))
  }
  
  for(N in intersect(names(sequence), names(fit_sequence)))
  {
    if(is.factor(fit_sequence[[N]]))
    {
      sequence[[N]] <- factor(as.character(sequence[[N]]), levels=levels(fit_sequence[[N]]))
    }
    else if(is.character(fit_sequence[[N]]))
    {
      sequence[[N]] <- as.character(sequence[[N]])
    }
    else if(is.numeric(fit_sequence[[N]]))
    {
      sequence[[N]] <- as.numeric(as.character(sequence[[N]]))
    }
    else if(is.logical(fit_sequence[[N]]))
    {
      sequence[[N]] <- as.logical(as.numeric(sequence[[N]]))
    }
  }
  
  return(list(sequence=sequence, missing=new_names, omitted_precads=omitted_precads))
}

#
# Add the individual position amino-acid and residue information to the sequence data
# Optionally, only do it for the CDRs
#
reshape_sequence_data <- function(sequence, ignore_G=TRUE, cdr_only)
{
  require(dplyr)
  LH <- names(sequence)[c(grep("^L", names(sequence)),grep("^H", names(sequence)))]
  
  if(cdr_only == TRUE)
  {
    cdrs <- cdr_residues()
    cdr_resi <- cdrs$Residue
    
    #
    # Double check if any LH positions are missing. This is likely to
    # happen in H3s or L3s of lengths not captured in the function above.
    # Shouldn't really happen since H3s upto length 35 are supported
    #
    H3 <- grep("H99", LH, value=TRUE)
    H3 <- c(H3, grep("H100", LH, value=TRUE))
    L3 <- grep("L95", LH, value=TRUE)
    
    cdr_resi <- unique(c(cdr_resi, H3, L3))
    LH <- LH[LH %in% cdr_resi]
  }
  
  # Replicate the data frame as many times as there are residues
  r <- bind_rows(replicate(length(LH), sequence, simplify=FALSE))
  
  # Add the AA_enc infornation by unwrapping the columns by row
  r$AA_enc <- sapply(apply(sequence[, LH], 2, unlist), unlist)
  
  # Add the residue names
  r$Residue <- unlist(sapply(LH, function(x){ rep(x, nrow(sequence))}, simplify=FALSE))
  
  # Remove positions that are insertions or glycines
  # 1 => '-'
  # 2 => 'G'
  if(ignore_G==TRUE)
  {
    r <- subset(r, AA_enc > 1)
  }
  else
  {
    r <- subset(r, AA_enc > 0)
  }
  
  return(r)
}

#
# Apply the SASA models to an input sequence to predict the SASA
# - input_sequence : New sequences for which SASA predictions are needed. Data needs to be in the same format
#    as the input data provided for training. Exmaple of such data is present in pdb_training_sequence.txt file
# - models : Randomforest models built using the build_models function
# - set_of_positions : Positions for which a prediction is needed
# - name_col : Name of the column containing the id of the sequence
#
apply_models_input_synced <- function(input_sequence, outcome='Fractional.SASA.predict', ignore_G=TRUE, models, set_of_positions=NULL, name_col='PreCADName', sasa_template=NULL)
{
  require(randomForest)
  
  # Get intersection between model positions and input sequence positions
  # Predictions will be made for these positions, since they are the ones for
  # which a model has been passed
  if(is.null(set_of_positions))
  {
    set_of_positions <- unique(input_sequence$Residue)
  }
  
  missing_positions <- setdiff(set_of_positions, names(models))
  
  set_of_positions <- intersect(set_of_positions, names(models))
  
  if(length(set_of_positions) == 0)
  {
    print('No valid positions found for predictions')
    return(data.frame())
  }
  
  input_sequence <- subset(input_sequence, Residue %in% set_of_positions)
  
  input_sequence_gly <- NULL
  if(ignore_G == TRUE)
  {
    input_sequence_gly <- subset(input_sequence, AA_enc == 2)
    input_sequence <- subset(input_sequence, AA_enc > 2)
  }
  
  predictions <- data.frame()
  
  # Predict for one position at a time over all input sequences
  for(pos in set_of_positions)
  {
    print(paste('Predicting for position', pos))
    
    if(!is.null(models[[pos]]))
    {
      data <- subset(input_sequence, Residue == pos)
      
      if(nrow(data) > 0)
      {
        data[[outcome]] <- predict(models[[pos]], newdata=data)
        
        # Set to zero for GLY]
        if(ignore_G == TRUE)
        {
          data[[outcome]][data$AA_enc <= 2] <- 0
        }
        
        # Update predictions data frame
        predictions <- rbind(predictions, data[,c('Residue', name_col, outcome, 'AA_enc')])
      }
    }
    else
    {
      missing_positions <- unique(c(missing_positions, pos))
    }
  }
  
  if(!is.null(input_sequence_gly) & nrow(input_sequence_gly))
  {
    gly <- input_sequence_gly[,c('Residue',name_col,'AA_enc')]
    gly[[outcome]] <- 0
    
    predictions <- rbind(predictions, gly)
  }
  
  # Set to NA for positions that are present but for which no models were found
  if(length(missing_positions))
  {
    missing <- NULL
    if(0 & !is.null(sasa_template))
    {
      fit_sasa_nonG <- subset(sasa_template, AA!='G' & Residue %in% missing_positions)
      for_missing <- data.frame(tapply(fit_sasa_nonG$Fractional.SASA, fit_sasa_nonG$Residue, mean))
      names(for_missing) <- c(outcome)
      for_missing$Residue <- rownames(for_missing)
      
      missing <- subset(input_sequence, Residue %in% missing_positions)
      missing <- merge(missing, for_missing, all.x=TRUE)
      
    }
    else
    {
      missing <- subset(input_sequence, Residue %in% missing_positions)
      
      if(nrow(missing) > 0)
        missing[[outcome]] <- NA
    }
    
    if(!is.null(missing) & nrow(missing) > 0)
    {
      missing[[outcome]][missing$AA_enc <= 2] <- 0
      predictions <- rbind(predictions, missing[,c('Residue', name_col, outcome, 'AA_enc')])
    }
  }
  
  # Put back AA column for user friendliness
  if(nrow(predictions))
  {
    trans <- data.frame(AA=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'))
    trans$AA_enc <- factor(trans$AA, levels=c("-","G","A","S","C","D","P","T","N","V","L","E","Q","H","I","M","K","F","Y","R","W"))
    trans$AA_enc <- as.numeric(trans$AA_enc)
    predictions <- merge(predictions, trans, by='AA_enc', suffixes=c('.P',''))
    
    # Add normalizing SASA
    norm <- get_norm_sasa()
    predictions <- merge(predictions, norm, sort=FALSE)
  }
  
  return(list(Pred=predictions, missing=missing_positions))
}

#
# SASA predictions
#
master_predict_sasa <- function(VH, VL, models, sasa_template, cdr_only=TRUE, test=FALSE, ...)
{
  # Read numbered sequences and germline info
  sequence <- combine_vh_vl_data(VH, VL)
  
  if(test == TRUE)
  {
    sequence <- sequence[1:min(nrow(sequence),10),]
  }
  
  # Change and update data-types needed by the prediction function
  sync <- sync_data_classes(sequence, sasa_template)
  
  # Reshape the data to get one line per position for which prediction is desired
  reshaped_data <- reshape_sequence_data(sync$sequence, cdr_only=cdr_only)
  
  # Predict the SASA
  sasa_pred <- apply_models_input_synced(reshaped_data, models=models, sasa_template=sasa_template, ...)
  
  # SASA wide format
  sasa_wide <- reshape(sasa_pred$Pred[, c('PreCADName','Residue', 'Fractional.SASA.predict')], direction="wide", idvar='PreCADName', timevar='Residue')
  sasa_wide[is.na(sasa_wide)] <- 0
  names(sasa_wide)[2:length(sasa_wide)] <- sapply(names(sasa_wide)[2:length(sasa_wide)], function(x) strsplit(x,'.',fixed=TRUE)[[1]][4])
  
  # Original sequence
  sequence <- subset(sequence, PreCADName %in% sync$sequence$PreCADName)
  
  return(list(SASA=sasa_pred, omitted_precads=sync$omitted_precads, sequence=sync$sequence, sasa_wide=sasa_wide, orig_sequence=sequence))
}


###################################################################
v2 <- load('aug_2017_sasa_used_for_pnas.bin') ## Model A
sasa_prediction <- master_predict_sasa(VH='cncn_VH.txt', VL='cncn_VL.txt', models=models, cdr_only=FALSE, sasa_template=sasa_template)
title<-read.csv(file='position.csv',header=FALSE)
#
# Calculate chemical flags and similarity to clinical mAbs using chemical space
#
ResultsMatrix<-function(sasa_prediction,title){
  orderSASA<-data.frame(matrix(ncol=dim(sasa_prediction$sasa_wide)[2],nrow=dim(sasa_prediction$sasa_wide)[1]))
  col_name2<-rep(NA,dim(title)[2])
  col_name<-colnames(sasa_prediction$sasa_wide)
  ## order SASA and remove NA in column names
  for (i in 1:dim(title)[2]){
    if (any(col_name==title[i][1,1])){
      orderSASA[,i]<-sasa_prediction$sasa_wide[,as.character(title[i][1,1])]
      col_name2[i]<-as.character(title[i][1,1])
    }else{
      orderSASA[,i]<-list(rep(0,dim(sasa_prediction$sasa_wide)[1]))
      col_name2[i]<-c("NA")
    }
  }
  colnames(orderSASA)<-col_name2
  
  ## Generate the df includes all VH and VL with the same position as predicted SASA values
  comb_seq<-function(VH,VL)
  {
    sequence<-combine_vh_vl_data(VH,VL)
  }
  orderseq<-comb_seq(VH='cncn_VL.txt', VL='cncn_VH.txt')
  orderseq<-orderseq[ , -which(names(orderseq) %in% c("VH","h1_len","h2_len","h3_len","PID_VH","canonical","VL","l1_len","l2_len","l3_len","PID_VL","VHF","VLF","Kappa","germline.species"))]
  orderSASA<-orderSASA[, grep("^(NA)", names(orderSASA), value = TRUE, invert = TRUE)]
  orderSASA<-orderBy(~PreCADName,orderSASA)
  rownames(orderSASA)<-1:nrow(orderSASA)
  orderseq<-orderBy(~PreCADName,orderseq)
  rownames(orderseq)<-1:nrow(orderseq)
  
  ## Define CDR regions 
  p1<-c("H26","H50","H93","L24","L50","L89","H1","L1","L1")
  p2<-c("H35","H65","H102","L34","L56","L97","H113","L107","H113")
  #list(SASA=sasa_pred, omitted_precads=sync$omitted_precads, sequence=sync$sequence, sasa_wide=sasa_wide, orig_sequence=sequence)
  num1<-rep(NA,length(p1))
  for (i in 1:length(p1)){
    num1[i]=which( colnames(orderseq)==p1[i] )
  }
  num2<-rep(NA,length(p2))
  for (i in 1:length(p2)){
    num2[i]=which( colnames(orderseq)==p2[i] )
  }
  ## separate the SASA and sequence in to sub-region and save in a list
  d1<-orderseq[,num1[1]:num2[1]]
  d2<-orderseq[,num1[2]:num2[2]]
  d3<-orderseq[,num1[3]:num2[3]]
  d4<-orderseq[,num1[4]:num2[4]]
  d5<-orderseq[,num1[5]:num2[5]]
  d6<-orderseq[,num1[6]:num2[6]]
  d7<-orderseq[,num1[7]:num2[7]]
  d8<-orderseq[,num1[8]:num2[8]]
  d9<-orderseq[,num1[9]:num2[9]]
  d1_sasa<-orderSASA[,num1[1]:num2[1]]
  d2_sasa<-orderSASA[,num1[2]:num2[2]]
  d3_sasa<-orderSASA[,num1[3]:num2[3]]
  d4_sasa<-orderSASA[,num1[4]:num2[4]]
  d5_sasa<-orderSASA[,num1[5]:num2[5]]
  d6_sasa<-orderSASA[,num1[6]:num2[6]]
  d7_sasa<-orderSASA[,num1[7]:num2[7]]
  d8_sasa<-orderSASA[,num1[8]:num2[8]]
  d9_sasa<-orderSASA[,num1[9]:num2[9]]
  Seq.list<-list(H1=d1,H2=d2,H3=d3,L1=d4,L2=d5,L3=d6,VH=d7,VL=d8,Fv=d9)
  SASA.list<-list(H1=d1_sasa,H2=d2_sasa,H3=d3_sasa,L1=d4_sasa,L2=d5_sasa,L3=d6_sasa,VH=d7_sasa,VL=d8_sasa,Fv=d9_sasa)
  
  ## Extract features from Fv and CDR region
  ## Calculate chemcial flags values and flag numbers for each mAbs
  Region<-c("Fv","VH","VH","VH","H1H2H3L1L2L3","H1H2H3L2L3","H1H3L1L2L3","H1H3L2","H1H3L3","H2H3L1L2L3","H2L1","H2L1L3")
  AA_rule<-c("DFLMNPT","ADEHILMNPQ","GKMQTWY","DEHLMNQRSW","DEFITVW","HLMPQRWY","AFHPRWY","HIKMPQRVWY","FHKLPQRWY","ADENTV","DEHMNPTY","EPQRVWY")
  flagvalue<-c(19.2,12.2,23.2,19.8,4.7,5,4.6,2.2,4,4.9,3.9,2.8)
  sign<-c("<","<",">","<","<",">",">",">",">","<","<",">")
  min_num<-c(1,2,4,5,10,11)
  f_sub<-data.frame(matrix(NA,nrow=dim(orderseq)[1],ncol=12))
  for (num in 1:dim(orderSASA)[1]){
    for (i in 1:12){
      r=nchar(Region[i])/2
      f_subsub<-data.frame(matrix(NA,nrow=1,ncol=r))
      for (r1 in 1:r){
        total_sub=data.frame()
        Region.name<-as.character(substr(Region[i],2*r1-1,2*r1))
        Re_sub<-Seq.list[Region.name][[1]][num,]
        for (j in 1:nchar(AA_rule[i])){
          test_v3<-data.frame(which(Re_sub==substr(AA_rule[i],j,j),arr.ind=TRUE))
          test_v4<-test_v3[test_v3[,"row"]==1,]
          total_sub<-rbind(total_sub,test_v4["col"])
        }
        test_x2<-data.frame()
        test_x2<-data.frame(SASA.list[Region.name][[1]][num,as.list(total_sub)$col])
        if (dim(test_x2)[1]==0){
          f_subsub[,r1]=0}else{
            f_subsub[,r1]=sum(test_x2[test_x2>0.1])
          }
      }
      f_sub[num,i]=sum(f_subsub)
    }
    f_sub[num,3]<-f_sub[num,3]+sum(nchar(Seq.list$VH[num,][Seq.list$VH[num,]=="G"]))
  }
  
  flag<-data.frame(matrix(NA,nrow=dim(orderSASA)[1],ncol=12))
  min_num<-c(1,2,4,5,10,11)
  for (i in 1:12){
    if (i %in% min_num){
      flag[,i]<-f_sub[,i]<flagvalue[i] 
    }else{
      flag[,i]<-f_sub[,i]>flagvalue[i]
    }
  }
  
  flagnum<-data.frame(matrix(NA,nrow=dim(orderSASA)[1],ncol=1))
  for (i in 1:dim(orderSASA)[1]){
    flagnum[i,]=length(flag[i,][flag[i,]==TRUE])
  }
 
  ## Calculate chemical space for each mAb
  Region_simi<-c("H1","H2","H3","H1H2H3","L1","L2","L3","L1L2L3","H1H2H3L1L2L3","VH","VL","Fv","Framework")
  Range_AA<-c("ACDEFGHIKLMNPQRSTVWY")
  name<-orderSASA[,1]
  check.limit<-setNames(replicate(dim(orderSASA)[1],as.data.frame(matrix(,13,22)),simplify=FALSE),name)
  
  countchar<-function(char,s){
    s2<-gsub(char,"",s)
    return(sum(nchar(s))-sum(nchar(s2)))
  }
  
  for (num in 1:dim(orderSASA)[1]){
    mAb.list<-as.character(name[num])
    colnames(check.limit[mAb.list][[1]])<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","length","netcharge")
    rownames(check.limit[mAb.list][[1]])<-Region_simi
    for (j in 1:nchar(Range_AA)){
      for (i in 1:12){
        r=nchar(Region_simi[i])/2
        SASA.AA<-data.frame(matrix(NA,nrow=1,ncol=r))
        length.residues<-data.frame(matrix(NA,r,1))
        netcharge<-data.frame(matrix(NA,r,1))
        G_AA<-data.frame(matrix(NA,r,1))
        for (r1 in 1:r){
          list.name<-as.character(substr(Region_simi[i],2*r1-1,2*r1))
          Re_simi_sub<-Seq.list[list.name][[1]][num,]
          test_v3<-data.frame(which(Re_simi_sub==substr(Range_AA,j,j),arr.ind=TRUE))
          test_v4<-test_v3[test_v3[,"row"]==1,]
          test_x2<-data.frame(SASA.list[list.name][[1]][num,as.list(test_v4["col"])$col])
          if (dim(test_x2)[2]==0){
            SASA.AA[,r1]=0}else{
              SASA.AA[,r1]=sum(test_x2)
            }
          length.residues[r1,]<-sum(nchar(gsub("-","",Seq.list[list.name][[1]][num,])))
          char<-Seq.list[list.name][[1]][num,]
          netcharge[r1,]<-countchar("R",char)+countchar("K",char)+0.1*countchar("H",char)-countchar("D",char)-countchar("E",char)
          G_AA[r1,]<-countchar("G",char)
        }
        check.limit[mAb.list][[1]][i,j]<-sum(SASA.AA)
        check.limit[mAb.list][[1]][i,21]<-sum(length.residues)
        check.limit[mAb.list][[1]][i,22]<-sum(netcharge)
        check.limit[mAb.list][[1]][i,6]<-sum(G_AA)
      }
      check.limit[mAb.list][[1]][13,j]<-check.limit[mAb.list][[1]][12,j]-check.limit[mAb.list][[1]][9,j]
      check.limit[mAb.list][[1]][13,21]<-check.limit[mAb.list][[1]][12,21]-check.limit[mAb.list][[1]][9,21]
      check.limit[mAb.list][[1]][13,22]<-check.limit[mAb.list][[1]][12,22]-check.limit[mAb.list][[1]][9,22]
    }
  }
  
  ## load maximum and minimum limits for clinical mAbs
  total.num<-dim(orderSASA)[1][1]
  CheSpace_max<-read_excel("ChemicalSpace.xlsx",sheet = "maximum")
  CheSpace_min<-read_excel("ChemicalSpace.xlsx",sheet = "minimum") 
  CheSpace_max<-CheSpace_max[,2:23]
  CheSpace_min<-CheSpace_min[,2:23]
  ## Results matrix for chemical flags
  Results<-data.frame(matrix(NA,nrow=as.numeric(dim(orderSASA)[1]),ncol=3))
  Results[,1]<-name
  Results[,2]<-flagnum
  for (i in 1:dim(orderSASA)[1]){
    print(i)
    Results[i,3]<-(1-(sum(check.limit[i][[1]]>CheSpace_max)+sum(check.limit[i][[1]]<CheSpace_min))/(13*22))*100
    if (Results[i,3]<99){
      Results[i,4]=c("warning")} else {
        Results[i,4]=c("normal")
      }
    }
  
  colnames(Results)<-c("mAb","chemicalflags","Similarity(%)","comments")
  
  ## Results details for each mAb (failed chemical rules)
  Results.details<-setNames(replicate(dim(orderSASA)[1],as.data.frame(matrix(,12,4)),simplify=FALSE),name)
  Region<-c("Fv","VH","VH","VH","H1H2H3L1L2L3","H1H2H3L2L3","H1H3L1L2L3","H1H3L2","H1H3L3","H2H3L1L2L3","H2L1","H2L1L3")
  AA_rule<-c("DFLMNPT","ADEHILMNPQ","GKMQTWY","DEHLMNQRSW","DEFITVW","HLMPQRWY","AFHPRWY","HIKMPQRVWY","FHKLPQRWY","ADENTV","DEHMNPTY","EPQRVWY")
  flagvalue<-c(19.2,12.2,23.2,19.8,4.7,5,4.6,2.2,4,4.9,3.9,2.8)
  sign<-c("<","<",">","<","<",">",">",">",">","<","<",">")
  #as.numeric(dim(orderSASA)[1]
  for (i in 1:as.numeric(dim(orderSASA)[1])){
    print(i)
    mAb.list<-as.character(name[i])
    row.num<-length(Region[flag[i,]==TRUE])
    Results.details[mAb.list][[1]]<-data.frame(matrix(NA,row.num,4))
    Results.details[mAb.list][[1]][,1]<-data.frame(Region[flag[i,]==TRUE])
    Results.details[mAb.list][[1]][,2]<-data.frame(AA_rule[flag[i,]==TRUE])
    Results.details[mAb.list][[1]][,3]<-data.frame(sign[flag[i,]==TRUE])
    Results.details[mAb.list][[1]][,4]<-data.frame(flagvalue[flag[i,]==TRUE])
    colnames(Results.details[mAb.list][[1]])<-c("Region","AAcomp","Sign","flag.value")
  }
  return(list(Results.details,Results))
  
}

#
## Get results 
## Results: chemical flag number, similarity
## Results.details: failed chemical rules for each mAb
#
Results.summary<-ResultsMatrix(sasa_prediction,title)
Results<-Results.summary[[1]]
Results.details<-Results.summary[[2]]


normal_result_flag2=normal_result[normal_result$chemicalflags==2,]
normal_result_flag3=normal_result[normal_result$chemicalflags==3,]
normal_result_flag4=normal_result[normal_result$chemicalflags==4,]
normal_result_flag5=normal_result[normal_result$chemicalflags==5,]
normal_result_flag6=normal_result[normal_result$chemicalflags==6,]
normal_result_flag7=normal_result[normal_result$chemicalflags==7,]
normal_result_flag8=normal_result[normal_result$chemicalflags==8,]
normal_result_flag9=normal_result[normal_result$chemicalflags==9,]
normal_result_flag10=normal_result[normal_result$chemicalflags==10,]
normal_result_flag11=normal_result[normal_result$chemicalflags==11,]
normal_result_flag12=normal_result[normal_result$chemicalflags==12,]

all_result_flag0=Results.details[Results.details$chemicalflags==0,]
all_result_flag1=Results.details[Results.details$chemicalflags==1,]
all_result_flag2=Results.details[Results.details$chemicalflags==2,]
all_result_flag3=Results.details[Results.details$chemicalflags==3,]
all_result_flag4=Results.details[Results.details$chemicalflags==4,]
all_result_flag5=Results.details[Results.details$chemicalflags==5,]
all_result_flag6=Results.details[Results.details$chemicalflags==6,]
all_result_flag7=Results.details[Results.details$chemicalflags==7,]
all_result_flag8=Results.details[Results.details$chemicalflags==8,]
all_result_flag9=Results.details[Results.details$chemicalflags==9,]
all_result_flag10=Results.details[Results.details$chemicalflags==10,]
all_result_flag11=Results.details[Results.details$chemicalflags==11,]
all_result_flag12=Results.details[Results.details$chemicalflags==12,]




