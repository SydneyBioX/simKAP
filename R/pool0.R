#' national_selection_pool
#'
#' this is the Australia current national recipient selection pool
#'
#' @param recip_matrix a data.frame, this is the input recipient matrix
#' @param donor_matrix a data.frame, this is the input donor matrix
#' @param whichdonor a categorical donor id, this is the incoming kidney from one donor
#' @param bloodGroup_strict_match a boolean value, TRUE/FALSE, this is to define whether we want the strict match for bloodgroup or not
#' @param AB_priority a boolean value, TRUE/FALSE, this is to define whether AB group can get from other blood group or not
#' @return a data.frame contains eligible recipients
#'
#' @examples
#'
#' data("rawdata", package = "simKAP")
#' rawdata$donor_rank = rank(rawdata[,'donor_kdri']) / nrow(rawdata)
#'
#'
#' selectionpool1.2 <- selection_default(rawdata,rawdata,1,
#' bloodGroup_strict_match = TRUE,AB_priority = TRUE);
#'
#' @export
#'

selection_default <- function(recip_matrix, donor_matrix,
                             whichdonor, bloodGroup_strict_match = TRUE,
                             AB_priority = TRUE){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (is.na(whichdonor)) stop("no kidney coming in")
    if (any(! c("recip_blgroup","recip_liststartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("donor_blgroup","tx_date") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")


    # Step 1) Eligible sub-recipients are those who started waiting before the donor became available
    tx_date = donor_matrix[whichdonor, "tx_date"]
    #print(tx_date)
    # Check blood group
    bg_match <- blood_group_match(donor_matrix$donor_blgroup[whichdonor],
                                recip_matrix$recip_blgroup,
                                strict_match = bloodGroup_strict_match,
                                AB_priority = AB_priority)
    eligible_recipM = recip_matrix[recip_matrix$recip_liststartdate < tx_date &
                                       bg_match,]
    # here, every recipt is eligible
    sub_recipM = eligible_recipM

    # Step 2) rank donor by donor kdri
    #linMap <- function(x, from, to)(x - min(na.omit(x))) / max(na.omit(x) - min(na.omit(x))) * (to - from) + from
    #donor_kdri = recip_donor_matrix$donor_kdri
    #donor_kdri=linMap(donor_kdri, 0,  1)
    #recip_donor_matrix$donor_rank=rank(recip_donor_matrix[,'donor_kdri']) / nrow(recip_donor_matrix) #from small to large, eg,rank(c(1,2,3))/3 returns 0.3333333 0.6666667 1.0000000
    # Step 3) get the rank of this donor
    #donor_percentile = recip_donor_matrix[whichdonor,"donor_rank"]

    return(sub_recipM)
}




#' us_selection_pool (updated recip epts for the new tx_date)
#'
#' this is the us inspired recipient selection pool
#' @param recip_matrix a data.frame, this is the input recipient matrix
#' @param donor_matrix a data.frame, this is the input donor matrix
#' @param whichdonor a categorical donor id, this is the incoming kidney from one donor
#' @param bloodGroup_strict_match a boolean value, TRUE/FALSE, this is to define whether we want the strict match for bloodgroup or not
#' @param AB_priority a boolean value, TRUE/FALSE, this is to define whether AB group can get from other blood group or not
#' @param threshold_number a numerical value between 0 and 1, this is the threshold number we would like to use for matching
#' @return a data.frame contains eligible recipients
#'
#' @examples
#' data("rawdata", package = "simKAP")
#' rawdata$donor_rank = rank(rawdata[,'donor_kdri']) / nrow(rawdata)
#'
#' rawdata$recip_tx_date=rawdata$tx_date
#' selectionpool2.2 <- selection_corisk(rawdata,rawdata,threshold_number = 0.1,
#' whichdonor=1,bloodGroup_strict_match = TRUE,AB_priority = TRUE);
#' @export
#'
#'
#'
#'



selection_corisk <- function(recip_matrix, donor_matrix,
                             threshold_number = 0.1, whichdonor,
                             bloodGroup_strict_match = TRUE,
                             AB_priority = TRUE){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (is.na(whichdonor)) stop("no kidney coming in")
    if (any(! c("recip_blgroup","recip_liststartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("donor_blgroup","tx_date") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")

    # Step 1) Eligible sub-recipients are those who started waiting before the donor became available
    tx_date = donor_matrix[whichdonor,"tx_date"]
    # Check blood group
    bg_match <- blood_group_match(donor_matrix$donor_blgroup[whichdonor],
                                recip_matrix$recip_blgroup,
                                strict_match = bloodGroup_strict_match,
                                AB_priority = AB_priority)
    eligible_recipM = recip_matrix[recip_matrix$recip_liststartdate < tx_date &
                                       bg_match,]

    # # Step 2) rank donor by donor kdri
    #
    # donor_kdri = donor_matrix$donor_kdri
    # donor_kdri = linMap(donor_kdri, 0,  1)
    # donor_matrix$donor_rank = rank(donor_matrix[,'donor_kdri']) / nrow(donor_matrix) #from small to large, eg,rank(c(1,2,3))/3 returns 0.3333333 0.6666667 1.0000000
    # # Step 3) rank recip by recip epts

    # the au epts is updated accorading to tx_date
    eligible_recipM$recip_age=
        as.numeric(tx_date-eligible_recipM$recip_tx_date+eligible_recipM$recip_age*365)/365#changed to the correct age 20210504
    eligible_recipM$recip_waittime=as.numeric(tx_date-eligible_recipM$recip_rrtstartdate)
    eligible_recipM$recip_epts <-(0.049*as.numeric(eligible_recipM$recip_age>25)*(eligible_recipM$recip_age-25)
                                  +0.493*(ifelse(eligible_recipM$recip_graftno == 1, 0, 1))
                                  +0.287*log(eligible_recipM$recip_waittime/365.25+1)
                                  +0.598*as.numeric(eligible_recipM$recip_waittime<365.25))

    recip_percentile = rank(eligible_recipM[,'recip_epts']) / nrow(eligible_recipM)
    eligible_recipM$recip_epts=recip_percentile #changed to the current rank 20210504
    # Step 4) going to search for the top a percent and the bottom 1-a percent, i.e. top matches top, bottom matches bottom
    #donor_percentile = donor_matrix[whichdonor, "donor_rank"]
    donor_percentile = donor_matrix[whichdonor, "donor_rank"]
    donor_percentile=donor_percentile[1] #we can have two because of two kidneys, same number
    #print(donor_percentile)
    #print(threshold_number)
    if (is.na(donor_percentile)){ sub_recipM=eligible_recipM } else if(donor_percentile < threshold_number){
        index = recip_percentile < threshold_number
        sub_recipM = eligible_recipM[index,]
    } else {
        index = recip_percentile >= threshold_number
        sub_recipM = eligible_recipM[index,]
    }

    return(sub_recipM)
}


#' au_slicingwindow_selection_pool (updated recip epts for the new tx_date)
#'
#' this is the au with threshold matching selection pool
#'
#' @param recip_matrix a data.frame, this is the input recipient matrix
#' @param donor_matrix a data.frame, this is the input donor matrix
#' @param whichdonor a categorical donor id, this is the incoming kidney from one donor
#' @param bloodGroup_strict_match a boolean value, TRUE/FALSE, this is to define whether we want the strict match for bloodgroup or not
#' @param AB_priority a boolean value, TRUE/FALSE, this is to define whether AB group can get from other blood group or not
#' @param threshold_number a numerical value between 0 and 1, this is the threshold number we would like to use for matching
#' @return a data.frame contains eligible recipients
#'
#' @examples
#' data("rawdata", package = "simKAP")
#' rawdata$donor_rank = rank(rawdata[,'donor_kdri']) / nrow(rawdata)
#'
#' rawdata$recip_tx_date=rawdata$tx_date
#' selectionpool1.1 <- selection_irisk(rawdata,rawdata,threshold_number = 0.1,
#' whichdonor=1,bloodGroup_strict_match = TRUE,AB_priority = TRUE);
#' @export
#'


selection_irisk = function(recip_matrix,
                            donor_matrix,
                            threshold_number = 0.1,
                            whichdonor,
                            bloodGroup_strict_match = TRUE,
                            AB_priority = TRUE){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (is.na(whichdonor)) stop("no kidney coming in")
    if (any(! c("recip_blgroup","recip_liststartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("donor_blgroup","tx_date") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")

    # Step 1) Eligible sub-recipients are those who started waiting before the donor became available
    tx_date = donor_matrix[whichdonor, "tx_date"]
    # Check blood group
    bg_match <- blood_group_match(donor_matrix$donor_blgroup[whichdonor],
                                recip_matrix$recip_blgroup,
                                strict_match = bloodGroup_strict_match,
                                AB_priority = AB_priority)
    eligible_recipM = recip_matrix[recip_matrix$recip_liststartdate < tx_date &
                                       bg_match,]
    # # Step 2) rank donor by donor kdri
    #
    # donor_kdri = donor_matrix$donor_kdri
    # donor_kdri = linMap(donor_kdri, 0,  1)
    # donor_matrix$donor_rank = rank(donor_matrix[,'donor_kdri']) / nrow(donor_matrix) #from small to large, eg,rank(c(1,2,3))/3 returns 0.3333333 0.6666667 1.0000000

    # Step 3) going to search between this donor_rank - threshold_number and this donor_rank + threshold_number
    donor_percentile = donor_matrix[whichdonor, "donor_rank"]
    upperlimit = min(donor_percentile +  threshold_number, 1)
    lowerlimit = max(donor_percentile -  threshold_number, 0)

    # Step 4) Looking at all eligible recipt (define as everyone that is on the registory), rank recipient by recip epts first

    # the au epts is updated accorading to tx_date
    eligible_recipM$recip_age=
        as.numeric(tx_date-eligible_recipM$recip_tx_date+eligible_recipM$recip_age*365)/365#changed to the correct age 20210504
    eligible_recipM$recip_waittime=as.numeric(tx_date-eligible_recipM$recip_rrtstartdate)
    eligible_recipM$recip_epts <-(0.049*as.numeric(eligible_recipM$recip_age>25)*(eligible_recipM$recip_age-25)
                                  +0.493*(ifelse(eligible_recipM$recip_graftno == 1, 0, 1))
                                  +0.287*log(eligible_recipM$recip_waittime/365.25+1)
                                  +0.598*as.numeric(eligible_recipM$recip_waittime<365.25))

    recip_percentile = rank(eligible_recipM[, 'recip_epts']) / length(eligible_recipM[, 'recip_epts'])
    eligible_recipM$recip_epts=recip_percentile #changed to the current rank 20210504
    index = (recip_percentile < upperlimit) & (recip_percentile > lowerlimit)
    sub_recipM = eligible_recipM[index,]

    return(sub_recipM)
}




#' age based only match with <35 matched to >35, 35-65 and then >65 (updated recip epts for the new tx_date)
#'
#'
#' @param recip_matrix a data.frame, this is the input recipient matrix
#' @param donor_matrix a data.frame, this is the input donor matrix
#' @param whichdonor a categorical donor id, this is the incoming kidney from one donor
#' @param bloodGroup_strict_match a boolean value, TRUE/FALSE, this is to define whether we want the strict match for bloodgroup or not
#' @param AB_priority a boolean value, TRUE/FALSE, this is to define whether AB group can get from other blood group or not
#' @return a data.frame contains eligible recipients
#'
#' @examples
#' data("rawdata", package = "simKAP")
#' rawdata$donor_rank = rank(rawdata[,'donor_kdri']) / nrow(rawdata)
#'
#' rawdata$recip_tx_date=rawdata$tx_date
#' selectionpool2.2 <- selection_corisk_age(rawdata,rawdata,
#' whichdonor=1,bloodGroup_strict_match = TRUE,AB_priority = TRUE);
#' @export
#'
#'
#'
#'



selection_corisk_age <- function(recip_matrix, donor_matrix,
                                 whichdonor,
                                 bloodGroup_strict_match = TRUE,
                                 AB_priority = TRUE){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (is.na(whichdonor)) stop("no kidney coming in")
    if (any(! c("recip_blgroup","recip_liststartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("donor_blgroup","tx_date") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")

    # Step 1) Eligible sub-recipients are those who started waiting before the donor became available
    tx_date = donor_matrix[whichdonor,"tx_date"]
    # Check blood group
    bg_match <- blood_group_match(donor_matrix$donor_blgroup[whichdonor],
                                  recip_matrix$recip_blgroup,
                                  strict_match = bloodGroup_strict_match,
                                  AB_priority = AB_priority)
    eligible_recipM = recip_matrix[recip_matrix$recip_liststartdate < tx_date &
                                       bg_match,]

    # # Step 2) calculate the updated age
    eligible_recipM$recip_age=
        as.numeric(tx_date-eligible_recipM$recip_tx_date+eligible_recipM$recip_age*365)/365 #changed to the correct age 20210504
    eligible_recipM$recip_waittime=as.numeric(tx_date-eligible_recipM$recip_rrtstartdate)
    eligible_recipM$recip_epts <-(0.049*as.numeric(eligible_recipM$recip_age>25)*(eligible_recipM$recip_age-25)
                                  +0.493*(ifelse(eligible_recipM$recip_graftno == 1, 0, 1))
                                  +0.287*log(eligible_recipM$recip_waittime/365.25+1)
                                  +0.598*as.numeric(eligible_recipM$recip_waittime<365.25))

    recip_percentile = rank(eligible_recipM[,'recip_epts']) / nrow(eligible_recipM)
    eligible_recipM$recip_epts=recip_percentile #changed to the current rank 20210504
    # Step 3) match within the age group
    #print(whichdonor)
    donor_percentile = donor_matrix[whichdonor, "donor_age"]
    donor_percentile=donor_percentile[1] #we can have two because of two kidneys, same number
    #print("original donor age")
    #print(donor_percentile)
    #print(threshold_number)
    recip_age=eligible_recipM$recip_age
    if (is.na(donor_percentile)){ sub_recipM=eligible_recipM }
    else if(donor_percentile < 35){
        #print("1")
        index = recip_age < 35
        sub_recipM = eligible_recipM[index,]
    } else if(donor_percentile >65){
        #print("2")
        index = recip_age > 65
        sub_recipM = eligible_recipM[index,]}
    else{
        #print("3")
        index = recip_age <=65 & recip_age >=35
        sub_recipM = eligible_recipM[index,]
    }

    return(sub_recipM)
}


#' us_selection_pool (updated recip epts for the new tx_date) with flexible calculation of epts, kdri
#'
#' this is the us inspired recipient selection pool with flexible calculation of epts, kdri
#' @param recip_matrix a data.frame, this is the input recipient matrix
#' @param donor_matrix a data.frame, this is the input donor matrix
#' @param whichdonor a categorical donor id, this is the incoming kidney from one donor
#' @param bloodGroup_strict_match a boolean value, TRUE/FALSE, this is to define whether we want the strict match for bloodgroup or not
#' @param AB_priority a boolean value, TRUE/FALSE, this is to define whether AB group can get from other blood group or not
#' @param threshold_number a numerical value between 0 and 1, this is the threshold number we would like to use for matching
#' @param recip_variable_names a categorical vector contains variables used in the calculation
#' @param donor_variable_names a categorical vector contains variables used in the calculation
#' @param historical_data a data.frame, this is data used to train the model, notice that this must contatin the variable "binary_outcome" which indicates the binary response variable used in the training data
#' @return a data.frame contains eligible recipients
#'
#' @examples
#' data("rawdata", package = "simKAP")
#' rawdata$donor_rank = rank(rawdata[,'donor_kdri']) / nrow(rawdata)
#'
#' rawdata$recip_tx_date=rawdata$tx_date
#' recip_variable_names=c("recip_age","recip_diabetes")
#' donor_variable_names=c("donor_age","donor_diabetes")
#' historical_data=rawdata
#' historical_data$binary_outcome=as.vector(ifelse(historical_data$tx_status==0&historical_data$tx_period>2920,0,ifelse(historical_data$tx_status==1& historical_data$tx_period<5*365,1,2))) #8yr 5yr
#' historical_data=historical_data[-which(historical_data$binary_outcome==2),]
#' historical_data$binary_outcome=as.factor(as.character(historical_data$binary_outcome))
#'
#' selectionpool2.2 <- selection_corisk_flexible(rawdata,rawdata,threshold_number = 0.1,
#' whichdonor=1,bloodGroup_strict_match = TRUE,AB_priority = TRUE,recip_variable_names,donor_variable_names,historical_data);
#' @export
#'
#'
#'
#'



selection_corisk_flexible <- function(recip_matrix, donor_matrix,
                             threshold_number = 0.1, whichdonor,
                             bloodGroup_strict_match = TRUE,
                             AB_priority = TRUE,recip_variable_names,donor_variable_names,historical_data){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (is.na(whichdonor)) stop("no kidney coming in")
    if (any(! c("recip_blgroup","recip_liststartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("donor_blgroup","tx_date") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")


    #donor model
    xnam <- paste(donor_variable_names, sep="")
    form=as.formula(paste("binary_outcome~ ", paste(xnam, collapse= "+")))
    model <- glm(form,family=binomial(link='logit'),data=historical_data)
    fitted.results <- predict(model,newdata=donor_matrix,type='response')
    donor_kdri_new=fitted.results
    donor_matrix$donor_kdri=donor_kdri_new
    # Step 1) Eligible sub-recipients are those who started waiting before the donor became available
    tx_date = donor_matrix[whichdonor,"tx_date"]
    # Check blood group
    bg_match <- blood_group_match(donor_matrix$donor_blgroup[whichdonor],
                                  recip_matrix$recip_blgroup,
                                  strict_match = bloodGroup_strict_match,
                                  AB_priority = AB_priority)
    eligible_recipM = recip_matrix[recip_matrix$recip_liststartdate < tx_date &
                                       bg_match,]

    # # Step 2) rank donor by donor kdri
    #
    # donor_kdri = donor_matrix$donor_kdri
    # donor_kdri = linMap(donor_kdri, 0,  1)
    # donor_matrix$donor_rank = rank(donor_matrix[,'donor_kdri']) / nrow(donor_matrix) #from small to large, eg,rank(c(1,2,3))/3 returns 0.3333333 0.6666667 1.0000000
    # # Step 3) rank recip by recip epts

    # the au epts is updated accorading to tx_date
    eligible_recipM$recip_age=
        as.numeric(tx_date-eligible_recipM$recip_tx_date+eligible_recipM$recip_age*365)/365#changed to the correct age 20210504
    eligible_recipM$recip_waittime=as.numeric(tx_date-eligible_recipM$recip_rrtstartdate)
    # eligible_recipM$recip_epts <-(0.049*as.numeric(eligible_recipM$recip_age>25)*(eligible_recipM$recip_age-25)
    #                               +0.493*(ifelse(eligible_recipM$recip_graftno == 1, 0, 1))
    #                               +0.287*log(eligible_recipM$recip_waittime/365.25+1)
    #                               +0.598*as.numeric(eligible_recipM$recip_waittime<365.25))
    #recip model
    xnam <- paste(recip_variable_names, sep="")
    form=as.formula(paste("binary_outcome~ ", paste(xnam, collapse= "+")))
    model <- glm(form,family=binomial(link='logit'),data=historical_data)
    #print(dim(eligible_recipM))
    #print(summary(eligible_recipM$recip_age))
    if(dim(eligible_recipM)[1]!=0){
    fitted.results <- predict(model,newdata=eligible_recipM,type='response')
    recip_epts_new=fitted.results
    eligible_recipM$recip_epts=recip_epts_new}
    recip_percentile = rank(eligible_recipM[,'recip_epts']) / nrow(eligible_recipM)
    eligible_recipM$recip_epts=recip_percentile
    #print(recip_percentile)
    # Step 4) going to search for the top a percent and the bottom 1-a percent, i.e. top matches top, bottom matches bottom
    donor_matrix$donor_rank = rank(donor_matrix[,'donor_kdri']) / nrow(donor_matrix)
    donor_percentile = donor_matrix[whichdonor,"donor_rank"]
    donor_percentile=donor_percentile[1] #we can have two because of two kidneys, same number
    #print(donor_percentile)
    #print(threshold_number)
    if (is.na(donor_percentile)){ sub_recipM=eligible_recipM } else if(donor_percentile < threshold_number){
        index = recip_percentile < threshold_number
        #print(index)
        sub_recipM = eligible_recipM[index,]
    } else {
        index = recip_percentile >= threshold_number
        sub_recipM = eligible_recipM[index,]
    }

    return(sub_recipM)
}


#' au_slicingwindow_selection_pool (updated recip epts for the new tx_date) with flexible calculation of epts, kdri
#'
#' this is the au with threshold matching selection pool with flexible calculation of epts, kdri
#'
#' @param recip_matrix a data.frame, this is the input recipient matrix
#' @param donor_matrix a data.frame, this is the input donor matrix
#' @param whichdonor a categorical donor id, this is the incoming kidney from one donor
#' @param bloodGroup_strict_match a boolean value, TRUE/FALSE, this is to define whether we want the strict match for bloodgroup or not
#' @param AB_priority a boolean value, TRUE/FALSE, this is to define whether AB group can get from other blood group or not
#' @param threshold_number a numerical value between 0 and 1, this is the threshold number we would like to use for matching
#' @param recip_variable_names a categorical vector contains variables used in the calculation
#' @param donor_variable_names a categorical vector contains variables used in the calculation
#' @param historical_data a data.frame, this is data used to train the model, notice that this must contatin the variable "binary_outcome" which indicates the binary response variable used in the training data
#' @return a data.frame contains eligible recipients
#'
#' @examples
#' data("rawdata", package = "simKAP")
#' rawdata$donor_rank = rank(rawdata[,'donor_kdri']) / nrow(rawdata)
#'
#' rawdata$recip_tx_date=rawdata$tx_date
#' recip_variable_names=c("recip_age","recip_diabetes")
#' donor_variable_names=c("donor_age","donor_diabetes")
#' historical_data=rawdata
#' historical_data$binary_outcome=as.vector(ifelse(historical_data$tx_status==0&historical_data$tx_period>2920,0,ifelse(historical_data$tx_status==1& historical_data$tx_period<5*365,1,2))) #8yr 5yr
#' historical_data=historical_data[-which(historical_data$binary_outcome==2),]
#' historical_data$binary_outcome=as.factor(as.character(historical_data$binary_outcome))
#'
#' selectionpool2.2 <- selection_irisk_flexible(rawdata,rawdata,threshold_number = 0.1,
#' whichdonor=1,bloodGroup_strict_match = TRUE,AB_priority = TRUE,recip_variable_names,donor_variable_names,historical_data);
#' @export
#'
#'


selection_irisk_flexible = function(recip_matrix,
                           donor_matrix,
                           threshold_number = 0.1,
                           whichdonor,
                           bloodGroup_strict_match = TRUE,
                           AB_priority = TRUE,recip_variable_names,donor_variable_names,historical_data){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (is.na(whichdonor)) stop("no kidney coming in")
    if (any(! c("recip_blgroup","recip_liststartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("donor_blgroup","tx_date") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")

    #donor model
    xnam <- paste(donor_variable_names, sep="")
    form=as.formula(paste("binary_outcome~ ", paste(xnam, collapse= "+")))
    model <- glm(form,family=binomial(link='logit'),data=historical_data)
    fitted.results <- predict(model,newdata=donor_matrix,type='response')
    donor_kdri_new=fitted.results
    donor_matrix$donor_kdri=donor_kdri_new
    # Step 1) Eligible sub-recipients are those who started waiting before the donor became available
    tx_date = donor_matrix[whichdonor, "tx_date"]
    # Check blood group
    bg_match <- blood_group_match(donor_matrix$donor_blgroup[whichdonor],
                                  recip_matrix$recip_blgroup,
                                  strict_match = bloodGroup_strict_match,
                                  AB_priority = AB_priority)
    eligible_recipM = recip_matrix[recip_matrix$recip_liststartdate < tx_date &
                                       bg_match,]
    # # Step 2) rank donor by donor kdri
    #
    # donor_kdri = donor_matrix$donor_kdri
    # donor_kdri = linMap(donor_kdri, 0,  1)
    # donor_matrix$donor_rank = rank(donor_matrix[,'donor_kdri']) / nrow(donor_matrix) #from small to large, eg,rank(c(1,2,3))/3 returns 0.3333333 0.6666667 1.0000000

    # Step 3) going to search between this donor_rank - threshold_number and this donor_rank + threshold_number
    donor_matrix$donor_rank = rank(donor_matrix[,'donor_kdri']) / nrow(donor_matrix)
    donor_percentile = donor_matrix[whichdonor, "donor_rank"]
    donor_percentile=donor_percentile[1]
    #print(donor_percentile)
    upperlimit = min(donor_percentile +  threshold_number, 1)
    lowerlimit = max(donor_percentile -  threshold_number, 0)
    #print(c(lowerlimit,upperlimit))
    # Step 4) Looking at all eligible recipt (define as everyone that is on the registory), rank recipient by recip epts first
    # the au epts is updated accorading to tx_date
    eligible_recipM$recip_age=
        as.numeric(tx_date-eligible_recipM$recip_tx_date+eligible_recipM$recip_age*365)/365#changed to the correct age 20210504
    eligible_recipM$recip_waittime=as.numeric(tx_date-eligible_recipM$recip_rrtstartdate)
    # eligible_recipM$recip_epts <-(0.049*as.numeric(eligible_recipM$recip_age>25)*(eligible_recipM$recip_age-25)
    #                               +0.493*(ifelse(eligible_recipM$recip_graftno == 1, 0, 1))
    #                               +0.287*log(eligible_recipM$recip_waittime/365.25+1)
    #                               +0.598*as.numeric(eligible_recipM$recip_waittime<365.25))
    #recip model
    xnam <- paste(recip_variable_names, sep="")
    form=as.formula(paste("binary_outcome~ ", paste(xnam, collapse= "+")))
    model <- glm(form,family=binomial(link='logit'),data=historical_data)
    if(dim(eligible_recipM)[1]!=0){
    fitted.results <- predict(model,newdata=eligible_recipM,type='response')
    recip_epts_new=fitted.results
    eligible_recipM$recip_epts=recip_epts_new}
    recip_percentile = rank(eligible_recipM[,'recip_epts']) / nrow(eligible_recipM)
    eligible_recipM$recip_epts=recip_percentile
    #print(recip_percentile)
    #eligible_recipM$recip_rank=recip_percentile
    index = (recip_percentile < upperlimit) & (recip_percentile > lowerlimit)
    sub_recipM = eligible_recipM[index,]

    return(sub_recipM)
}




#' get user defined epts and kdri
#'
#' this is the function where users can feed the variables in to get their own "epts" and "kdri" as risk vectors
#'
#' @param recip_matrix a data.frame, this is the input recipient matrix
#' @param donor_matrix a data.frame, this is the input donor matrix
#' @param historical_data a data.frame, this is data used to train the model, notice that this must contatin the variable "binary_outcome" which indicates the binary response variable used in the training data
#' @param donor_variable_names a categorical vector contains variables used in the calculation
#' @param recip_variable_names a categorical vector contains variables used in the calculation
#' @return a list of length two contains recip_epts_new and donor_kdri_new
#'
#' @examples
#' data("rawdata", package = "simKAP")
#' rawdata$donor_rank = rank(rawdata[,'donor_kdri']) / nrow(rawdata)
#'
#' rawdata$recip_tx_date=rawdata$tx_date
#' recip_variable_names=c("recip_age","recip_diabetes")
#' donor_variable_names=c("donor_age","donor_diabetes")
#' historical_data=rawdata
#' historical_data$binary_outcome=as.vector(ifelse(historical_data$tx_status==0&historical_data$tx_period>2920,0,ifelse(historical_data$tx_status==1& historical_data$tx_period<5*365,1,2))) #8yr 5yr
#' historical_data=historical_data[-which(historical_data$binary_outcome==2),]
#' historical_data$binary_outcome=as.factor(as.character(historical_data$binary_outcome))
#' epts_kdri_new <- epts_kdri_new(rawdata,rawdata,historical_data,recip_variable_names,donor_variable_names);
#' @export
#'


epts_kdri_new = function(recip_matrix,
                                    donor_matrix,historical_data,recip_variable_names,donor_variable_names){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    #recip model
    xnam <- paste(recip_variable_names, sep="")
    form=as.formula(paste("binary_outcome~ ", paste(xnam, collapse= "+")))
    model <- glm(form,family=binomial(link='logit'),data=historical_data)
    fitted.results <- predict(model,newdata=recip_matrix,type='response')
    recip_epts_new=fitted.results
    #donor model
    xnam <- paste(donor_variable_names, sep="")
    form=as.formula(paste("binary_outcome~ ", paste(xnam, collapse= "+")))
    model <- glm(form,family=binomial(link='logit'),data=historical_data)
    fitted.results <- predict(model,newdata=donor_matrix,type='response')
    donor_kdri_new=fitted.results
    return(list(recip_epts_new,donor_kdri_new))
}










