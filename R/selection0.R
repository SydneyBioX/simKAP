#' selection
#'
#' this is the max selection function
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features, this is the input donor matrix
#' @param allocation_score a vector, this is the score for each pair
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @param graft_number a numerical value, this is the graft number
#' @return a matrix contains allo_score, which is the alocation score for each new pair
#'
#' @examples
#' data("rawdata", package = "simKAP")
#' selected <- select_max(allocation_score=allocation_national(rawdata,rawdata,rawdata),
#' graft_number = 1);
#' @export





select_max <- function(recip_matrix = NULL,
                       donor_matrix = NULL,
                       allocation_score,
                       HLA = NULL,
                       graft_number = 1
) {
    #if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    #if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    #if (is.na(allocation_score)) stop("no allocation score")
    #if (any(! c("recip_id") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")

    threshold <- sort(allocation_score, decreasing = TRUE)[graft_number]
    select_id <- which(allocation_score >= threshold)
    if (length(select_id) >= graft_number) {
        select_id <- select_id[seq_len(graft_number)]
        select_recip_id <- recip_matrix$recip_id[select_id]
        return(select_recip_id)
    } else {
        return(NULL)
    }
}



#' selection
#'
#' this is the sdm3
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features, this is the input donor matrix
#' @param allocation_score a vector, this is the score for each pair
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @param graft_number a numerical value, this is the graft number
#' @param parameter_a an interger hyperparameter to tune
#' @param parameter_b an interger hyperparameter to tune
#' @param parameter_c an interger hyperparameter to tune
#' @param quality a boolean value, TRUE/FALSE, quality measure indicator
#' @param age a boolean value, TRUE/FALSE, age indicator
#' @param state a boolean value, TRUE/FALSE, state indicator
#' @param bloodgroup a boolean value, TRUE/FALSE, blood group indicator
#' @param accept_prob a numerical value between 0 and 1, which is a user defined minimum acceptence probability, then for different recipient, we randomly sample from this value to 1, the default value is NULL, which means we give an accept rate based on HLA matchability
#' @param yearly_scale a vector representing yearly scale, this allows the user to define different accpetence rate for different years
#'
#'
#' @return a data.frame with matched recipient
#'
#' @examples
#' data("rawdata", package = "simKAP")
#' data("newdata", package = "simKAP")
#' HLA_matrix <- hla_match(raw_recip_matrix_subset, raw_donor_matrix[1,])
#'
#' selected <- dm_national_choice(raw_recip_matrix_subset,raw_donor_matrix[1,],
#' allocation_score=allocation_national(raw_recip_matrix_subset,raw_donor_matrix[1,],HLA_matrix),
#' HLA_matrix,1,1,1,1,TRUE,TRUE,TRUE,TRUE,accept_prob = NULL, yearly_scale=rep(1, 100));
#' @export



## Respect the ordering by national allocation score
## include the EPTS and KDRI
## include age difference
## ignore blood matching

dm_national_choice <- function(recip_matrix = NULL,
                              donor_matrix = NULL,
                              allocation_score = NULL,
                              HLA = NULL,
                              graft_number = 1,
                              parameter_a = 1, parameter_b = 1,
                              parameter_c = 1,
                              quality = TRUE, age = TRUE, state = TRUE, bloodgroup = TRUE,
                              accept_prob = NULL, yearly_scale=rep(1, 100)
){


    HLAall = as.numeric(HLA$tx_misa) + as.numeric(HLA$tx_misb) + as.numeric(HLA$tx_misdr)

    recip_matrix = cbind(recip_matrix, allocation_score)
    recip_matrix = recip_matrix[with(recip_matrix, order(-allocation_score)), ]


    if(nrow(recip_matrix) == 0 ) {
        return(NULL)
    }


    donor_state <- donor_matrix[,"donor_state"]
    donor_state <- switch(as.character(donor_state),
                          NSW = "NSW", ACT = "NSW",
                          VIC = "VIC", TAS = "VIC",
                          SA = "SA", NT = "SA",
                          WA = "WA", QLD = "QLD")

    accept = 0
    i = 1 #start at patient index 1
    recip = c()

    #print(c("graft",graft_number))

    while( accept != graft_number ){



        ########### adding the yearly factor ##############

        days <- recip_matrix[i,"recip_sim_listdate"] - min(recip_matrix$recip_sim_listdate)
        #print(c("days", days))
        year = round(days/365+1)
        #print(c("year", year))

        ####################


        recip_age = round((donor_matrix[,"tx_date"]-recip_matrix[i,"recip_birthdate"])/365)
        age_difference = donor_matrix[,"donor_age"] - recip_age


        if (!is.null(accept_prob)){
            #print(c("yearly_scale", yearly_scale[year]))
            p_accept = sample(accept_prob,1)*yearly_scale[year]
            #print(c(p_accept, "p_Accept state"))
        }

        if (is.null (accept_prob)){
            pra_val = recip_matrix[i,"recip_pra"]/100
            hla_val = 2*HLA[i,"tx_misdr"] + HLA[i,"tx_misa"] + HLA[i,"tx_misb"]


            #print(c("yearly_scale", yearly_scale[year]))
            p_accept = ((1-pra_val)*parameter_a)*(1-hla_val/50*parameter_b)*yearly_scale[year]
            #print(c(p_accept, "p_Accept"))
        }


        p_accept = min(1, p_accept)


        if(i == nrow(recip_matrix) ){

            return(NULL)

        }else{


            if(quality & recip_matrix[i, "recip_epts"] < donor_matrix[,"donor_kdri"] ){

                samp = rbinom(1,1,p_accept)

                if (samp == 1) {
                    accept = accept + 1
                    recip = cbind(recip, i)}

            } else if (state & recip_matrix[i, "recip_state"]
                       %in% donor_state ) {

                accept = accept + 1
                recip = cbind(recip, i)

            } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB") {

                accept = accept + 1
                recip = cbind(recip, i)


            } else if (age & recip_age < 18
                       & age_difference < 30 ) {

                accept = accept + 1
                recip = cbind(recip, i)

            }else{

                samp = rbinom(1,1,p_accept*parameter_c)


                if (samp == 1) {
                    accept = accept + 1
                    recip = cbind(recip, i)

                }

            }
        }

        i = i + 1

    }


    #print(recip)
    match_recip = as.vector(recip_matrix[recip,"recip_id"])
    #print(c("matched", match_recip))


    return(match_recip)

}





#' selection
#'
#' this is the sdm1
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features, this is the input donor matrix
#' @param allocation_score a vector, this is the score for each pair
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @param graft_number a numerical value, this is the graft number
#' @param parameter_a an interger hyperparameter to tune
#' @param parameter_b an interger hyperparameter to tune
#' @param parameter_c an interger hyperparameter to tune
#' @param quality a boolean value, TRUE/FALSE, quality measure indicator
#' @param age a boolean value, TRUE/FALSE, age indicator
#' @param state a boolean value, TRUE/FALSE, state indicator
#' @param bloodgroup a boolean value, TRUE/FALSE, blood group indicator
#' @param accept_prob a numerical value between 0 and 1, which is a user defined minimum acceptence probability, then for different recipient, we randomly sample from this value to 1, the default value is NULL, which means we give an accept rate based on HLA matchability
#' @param yearly_scale a vector representing yearly scale, this allows the user to define different accpetence rate for different years
#'
#'
#' @return a data.frame with matched recipient
#' @examples
#' data("newdata", package = "simKAP")
#' HLA_matrix <- hla_match(raw_recip_matrix_subset, raw_donor_matrix[1,])
#'
#' selected <- dm_national_formula(raw_recip_matrix_subset,raw_donor_matrix[1,],
#' allocation_score=allocation_national(raw_recip_matrix_subset,raw_donor_matrix[1,],HLA_matrix),
#' HLA_matrix,1,1,1,1,TRUE,TRUE,TRUE,TRUE,accept_prob = NULL, yearly_scale=rep(1, 100));
#' @export



## Respect the ordering by national allocation score
## include the state balance
## include the EPTS and KDRI
## include age difference
## ignore blood matching

dm_national_formula <- function(recip_matrix = NULL,
                               donor_matrix = NULL,
                               allocation_score = NULL,
                               HLA = NULL,
                               graft_number = 1,
                               # graft_number = num_kidney_remain,
                               parameter_a = 1, parameter_b = 1,parameter_c = 1,
                               quality = TRUE, age = TRUE, state = TRUE,
                               bloodgroup = TRUE, accept_prob = NULL, yearly_scale=rep(1, 100)
){

    HLAall = as.numeric(HLA$tx_misa) + as.numeric(HLA$tx_misb) + as.numeric(HLA$tx_misdr)

    recip_matrix = cbind(recip_matrix, allocation_score)
    recip_matrix = recip_matrix[with(recip_matrix, order(-allocation_score)), ]


    if(nrow(recip_matrix) == 0 ) {
        return(NULL)
    }


    donor_state <- donor_matrix[,"donor_state"]
    donor_state <- switch(as.character(donor_state),
                          NSW = "NSW", ACT = "NSW",
                          VIC = "VIC", TAS = "VIC",
                          SA = "SA", NT = "SA",
                          WA = "WA", QLD = "QLD")

    accept = 0
    i = 1 #start at patient index 1
    recip = c()

    #print(c("graft",graft_number))

    while( accept != graft_number ){



        ########### adding the yearly factor ##############

        days <- recip_matrix[i,"recip_sim_listdate"] - min(recip_matrix$recip_sim_listdate)
        #print(c("days", days))
        year = round(days/365+1)
        #print(c("year", year))

        ####################


        recip_age = round((donor_matrix[,"tx_date"]-recip_matrix[i,"recip_birthdate"])/365)
        age_difference = donor_matrix[,"donor_age"] - recip_age


        if (!is.null(accept_prob)){
            # print(c("yearly_scale", yearly_scale[year]))
            p_accept = sample(accept_prob,1)*yearly_scale[year]
            #print(c(p_accept, "p_Accept state"))
        }

        if (is.null (accept_prob)){
            pra_val = recip_matrix[i,"recip_pra"]/100
            hla_val = 2*HLA[i,"tx_misdr"] + HLA[i,"tx_misa"] + HLA[i,"tx_misb"]


            #print(c("yearly_scale", yearly_scale[year]))
            p_accept = ((1-pra_val)*parameter_a)*(1-hla_val/50*parameter_b)*yearly_scale[year]
            # print(c(p_accept, "p_Accept"))
        }


        p_accept = min(1, p_accept)


        if(i == nrow(recip_matrix)){

            return(NULL)

        }else{

            if(recip_matrix[i,"recip_pra"] > 50
               & HLAall[i] == 0){

                accept = accept + 1
                recip = cbind(recip, i)

            }else if (recip_matrix[i,"recip_pra"] > 80
                      & HLAall[i] == 1){

                if (quality & recip_matrix[i,"recip_epts"]
                    < donor_matrix[,"donor_kdri"]) {

                    samp = rbinom(1,1,0.85)

                    if( samp == 1 ){
                        accept = accept + 1
                        recip = cbind(recip, i)

                        # this is a second chance at acceptance if it was rejected due to quality reasons
                    } else if(state & recip_matrix[i, "recip_state"]
                              %in% donor_state){

                        samp = rbinom(1,1,0.85)

                        if ( samp == 1 ) {
                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (age & recip_age < 18
                                   & age_difference < 30 ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        }
                    }

                } else {

                    samp = rbinom(1,1,0.75)

                    if( samp == 1 ){
                        accept = accept + 1
                        recip = cbind(recip, i)

                    } else if(state & recip_matrix[i, "recip_state"] %in% donor_state){

                        samp = rbinom(1,1,0.80)

                        if ( samp == 1 ) {
                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (age & recip_age < 18
                                   & age_difference < 30 ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        }
                    }
                }

            } else if (recip_matrix[i,"recip_pra"] <= 50 &
                       HLAall[i] == 0 ){

                if (quality & recip_matrix[i,"recip_epts"] < donor_matrix[,"donor_kdri"]) {

                    accept = accept + 1
                    recip = cbind(recip, i)

                    # this is for when epts > kdri
                } else if(state & recip_matrix[i, "recip_state"] %in% donor_state){

                    samp = rbinom(1,1,0.80)

                    if ( samp == 1 ) {
                        accept = accept + 1
                        recip = cbind(recip, i)

                    } else if ( age & recip_age < 18
                                & age_difference < 30 ) {
                        accept = accept + 1
                        recip = cbind(recip, i)

                    } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                        accept = accept + 1
                        recip = cbind(recip, i)

                    }
                }


            } else if (recip_matrix[i,"recip_pra"] <= 50 &
                       HLAall[i] == 1 ){


                if (quality & recip_matrix[i,"recip_epts"] < donor_matrix[,"donor_kdri"]) {

                    samp = rbinom(1,1,0.85)

                    if (samp == 1){

                        accept = accept + 1
                        recip = cbind(recip, i)

                        # this is for when epts > kdri
                    } else if(state & recip_matrix[i, "recip_state"] %in% donor_state){

                        samp = rbinom(1,1,0.80)

                        if ( samp == 1 ) {
                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if ( age & recip_age < 18
                                    & age_difference < 30 ) {
                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        }
                    }

                } else {

                    samp = rbinom(1,1,0.75)

                    if( samp == 1 ){
                        accept = accept + 1
                        recip = cbind(recip, i)

                    } else if(state & recip_matrix[i, "recip_state"] %in% donor_state){

                        samp = rbinom(1,1,0.80)

                        if ( samp == 1 ) {
                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (age & recip_age < 18
                                   & age_difference < 30 ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        }
                    }
                }

                # all other HLA combinations
            } else {

                if (quality & recip_matrix[i,"recip_epts"] < donor_matrix[,"donor_kdri"]) {

                    samp = rbinom(1,1,p_accept)

                    if( samp == 1 ){
                        accept = accept + 1
                        recip = cbind(recip, i)

                    } else if(state & recip_matrix[i, "recip_state"]
                              %in% donor_state){

                        samp = rbinom(1,1,p_accept)

                        if ( samp == 1 ) {
                            accept = accept + 1
                            recip = cbind(recip, i)
                        } else if (age & recip_age < 18
                                   & age_difference < 30 ) {
                            accept = accept + 1
                            recip = cbind(recip, i)

                        }else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        }
                    }

                    # for epts > kdri
                } else {

                    samp = rbinom(1,1, p_accept)

                    if( samp == 1 ){
                        accept = accept + 1
                        recip = cbind(recip, i)

                    } else if(state & recip_matrix[i, "recip_state"]
                              %in% donor_state){

                        samp = rbinom(1,1,0.85)

                        if ( samp == 1 ) {
                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (age & recip_age < 18
                                   & age_difference < 30 ) {

                            accept = accept + 1
                            recip = cbind(recip, i)

                        } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                            accept = accept + 1
                            recip = cbind(recip, i)
                        }
                    }
                }
            }
        }
        i = i+1
    }

    #print(c("recip",recip))
    match_recip = as.vector(recip_matrix[recip,"recip_id"])

    #print out the position on the list

    return(match_recip)

}





#' selection
#'
#' this is the sdm2
#'
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features, this is the input donor matrix
#' @param allocation_score a vector, this is the score for each pair
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @param graft_number a numerical value, this is the graft number
#' @param parameter_a an integer hyperparameter to tune
#' @param parameter_b an integer hyperparameter to tune
#' @param top_matches an integer, top match parameter
#' @param age a boolean value, TRUE/FALSE, age indicator
#' @param bloodgroup a boolean value, TRUE/FALSE, blood group indicator
#' @param parameter_state an integer, state parameter
#'
#'
#' @return a data.frame with matched recipient
#' @examples
#' data("newdata", package = "simKAP")
#' HLA_matrix <- hla_match(raw_recip_matrix_subset, raw_donor_matrix[1,])
#'
#' selected <- dm_national_top(raw_recip_matrix_subset,raw_donor_matrix[1,],
#' allocation_score=allocation_national(raw_recip_matrix_subset,raw_donor_matrix[1,],HLA_matrix),
#' HLA_matrix,1,1,1,1.5,1.5,20,15);
#' @export


##Take the top 15 - look to see which has desirable qualities
##include the state balance
##include the EPTS and KDRI

dm_national_top <- function(recip_matrix = NULL,
                           donor_matrix = NULL,
                           allocation_score = NULL,
                           HLA = NULL,
                           graft_number = 1,
                           # graft_number = num_kidney_remain,
                           parameter_a = 1, parameter_b = 1,
                           bloodgroup = 1.5,
                           age = 1.5,
                           parameter_state = 20,
                           top_matches = 15
                           # parameters can be set to 0
){
    #if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    #if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    #if (dim(HLA)[1]==0) stop("empty HLA matrix")
    #if (is.na(allocation_score)) stop("no allocation score")
    if (any(! c("recip_birthdate","recip_pra","recip_epts","recip_state","recip_id") %in% colnames(recip_matrix))) {
        stop("recipient matrix wrong column names")
    }
    if (any(! c("donor_state","tx_date","donor_age","donor_kdri","donor_blgroup") %in% colnames(donor_matrix))) {
        stop("donor matrix wrong column names")
    }
    if (any(! c("tx_misa","tx_misb","tx_misdr") %in% colnames(HLA))) {
        stop("HLA matrix wrong column names")
    }



    # print("national")

    HLAall = as.numeric(HLA$tx_misa) + as.numeric(HLA$tx_misb) + as.numeric(HLA$tx_misdr)

    recip_matrix = cbind(recip_matrix, allocation_score)
    recip_matrix = recip_matrix[with(recip_matrix, order(-allocation_score)), ]

    donor_state <- donor_matrix[,"donor_state"]
    donor_state <- switch(as.character(donor_state),
                          NSW = "NSW", ACT = "NSW",
                          VIC = "VIC", TAS = "VIC",
                          SA = "SA", NT = "SA",
                          WA = "WA", QLD = "QLD")


    if(nrow(recip_matrix) == 0){

        return(NULL)

    }else{

        min_nrow = min(top_matches, nrow(recip_matrix))
        recip_matrix = recip_matrix[c(1:min_nrow),]

        blood <- ifelse(donor_matrix[,"donor_blgroup"] %in% "AB", bloodgroup, 1 )

        recip_matrix$rank_score = rank(recip_matrix$allocation_score)*blood
        recip_matrix$rank_epts = rank(-recip_matrix$recip_epts)
        recip_matrix$rank_state = ifelse(recip_matrix$recip_state== donor_state,parameter_state,0)

        index <- which(as.numeric(as.character(recip_matrix$recip_age)) < 18)
        recip_matrix[index,]$rank_score <- recip_matrix[index,]$rank_score*age

        recip_matrix$rank_sum = recip_matrix$rank_score*parameter_a +
            recip_matrix$rank_epts*parameter_b +
            recip_matrix$rank_state

        ranked_recip = recip_matrix[with(recip_matrix, order(-rank_sum)), ]
        num <- min(nrow(ranked_recip),graft_number)
        match_recip = ranked_recip[c(1:num), "recip_id"]


        #i_rank = paste(i,as.vector(recip_matrix[i, "rank_sum"]), sep="-")
        # print(match_recip)
        # print(graft_number)
        # print( ranked_recip[c(1:graft_number), "rank_sum"])

        return(c(match_recip))

    }
}



#' selection
#'
#' this is the sdm1
#'
#' @return matched recipient
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features, this is the input donor matrix
#' @param allocation_score a vector, this is the score for each pair
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @param graft_number a numerical value, this is the graft number
#' @param parameter_state_a an interger hyperparameter to tune
#' @param parameter_state_b an interger hyperparameter to tune
#' @param parameter_state_c an interger hyperparameter to tune
#' @param quality a boolean value, TRUE/FALSE, quality measure indicator
#' @param age a boolean value, TRUE/FALSE, age indicator
#' @param bloodgroup a boolean value, TRUE/FALSE, blood group indicator
#' @param accept_prob a numerical value between 0 and 1, acceptence probability
#' @param yearly_scale a vector representing yearly scale
#'
#'
#' @return a data.frame with matched recipient
#' @examples
#' data("newdata", package = "simKAP")
#' HLA_matrix <- hla_match(raw_recip_matrix_subset, raw_donor_matrix[1,])
#'
#' selected <- dm_state_choice(raw_recip_matrix_subset,raw_donor_matrix[1,],
#' allocation_score=allocation_national(raw_recip_matrix_subset,raw_donor_matrix[1,],HLA_matrix),
#' HLA_matrix,1,1,1,1,TRUE,TRUE,TRUE,accept_prob = NULL, yearly_scale=rep(1, 100));
#' @export
#'



## Respect the ordering by national allocation score
## include the EPTS and KDRI
## include age difference
## ignore blood matching

dm_state_choice <- function(recip_matrix = NULL,
                           donor_matrix = NULL,
                           allocation_score = NULL,
                           HLA = NULL,
                           graft_number = 1,
                           parameter_state_a = 1,
                           parameter_state_b = 1,
                           parameter_state_c =1,
                           quality = TRUE, age = TRUE,
                           bloodgroup = TRUE,
                           accept_prob = NULL, yearly_scale=rep(1, 100)
){

    recip_matrix = cbind(recip_matrix, allocation_score)
    recip_matrix = recip_matrix[with(recip_matrix, order(-allocation_score)), ]
    print("State")

    if( nrow( recip_matrix) == 0 ) {
        return(NULL)
    }

    accept = 0
    i = 1
    recip = c()

    while(accept != graft_number ){


        ########### adding the yearly factor ##############

        days <- recip_matrix[i,"recip_sim_listdate"] - min(recip_matrix$recip_sim_listdate)
        #print(c("days", days))
        year = round(days/365+1)
        #print(c("year", year))

        ####################

        recip_age = round((donor_matrix[,"tx_date"]-recip_matrix[i,"recip_birthdate"])/365)
        age_difference = donor_matrix[,"donor_age"] - recip_age

        if (!is.null(accept_prob) ){
            #print(c("yearly_scale", yearly_scale[year]))
            p_accept = sample(accept_prob,1)*yearly_scale[year]
            #print(c(p_accept, "p_Accept state"))
        }

        if (is.null (accept_prob)){
            pra_val = recip_matrix[i,"recip_pra"]/100
            hla_val = 2*HLA[i,"tx_misdr"] + HLA[i,"tx_misa"] + HLA[i,"tx_misb"]
            #print(c("yearly_scale", yearly_scale[year]))
            p_accept = ((1-pra_val)*parameter_state_a)*(1-hla_val/50*parameter_state_b)*yearly_scale[year]
            #print(c(p_accept, "p_Accept"))
        }

        p_accept = min(1, p_accept)

        if(i == nrow(recip_matrix) ){

            return(NULL)

        }else{


            if(quality & recip_matrix[i, "recip_epts"] < donor_matrix[,"donor_kdri"] ){

                samp = rbinom(1,1,p_accept)

                if (samp == 1) {
                    accept = accept + 1
                    recip = cbind(recip, i)}

            } else if (age & recip_age < 18 & age_difference < 30 ) {

                accept = accept + 1
                recip = cbind(recip, i)

            } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB") {

                accept = accept + 1
                recip = cbind(recip, i)

            } else{

                samp = rbinom(1,1,p_accept*parameter_state_c)

                if (samp == 1) {
                    accept = accept + 1
                    recip = cbind(recip, i)

                }
            }
        }

        i = i + 1

    }


    match_recip = as.vector(recip_matrix[recip,"recip_id"])

    #output the position on the list some other way
    #output 0 or 1 recipient

    return(match_recip)

}

#' selection
#'
#' this is the sdm1
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features, this is the input donor matrix
#' @param allocation_score a vector, this is the score for each pair
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @param graft_number a numerical value, this is the graft number
#' @param parameter_state_a an integer hyperparameter to tune
#' @param parameter_state_b an integer hyperparameter to tune
#' @param parameter_state_c an integer hyperparameter to tune
#' @param age a boolean value, TRUE/FALSE, age indicator
#' @param bloodgroup a boolean value, TRUE/FALSE, blood group indicator
#' @param accept_prob a numerical value between 0 and 1, acceptence probability
#' @param yearly_scale a vector, yearly scale
#' @param quality a boolean value, TRUE/FALSE, quality indicator
#'
#' @return a data.frame with matched recipient
#' @examples
#' data("newdata", package = "simKAP")
#' HLA_matrix <- hla_match(raw_recip_matrix_subset, raw_donor_matrix[1,])
#'
#' selected <- dm_state_formula(raw_recip_matrix_subset,raw_donor_matrix[1,],
#' allocation_score=allocation_national(raw_recip_matrix_subset,raw_donor_matrix[1,],HLA_matrix),
#' HLA_matrix,1,1,1,TRUE,TRUE,TRUE,accept_prob = NULL, yearly_scale=rep(1, 100));
#' @export
#'


## Respect the ordering by national allocation score
## include the EPTS and KDRI
## include age difference
## ignore blood matching

dm_state_formula <- function(recip_matrix = NULL,
                            donor_matrix = NULL,
                            allocation_score = NULL,
                            HLA = NULL,
                            graft_number = 1,
                            parameter_state_a = 1, parameter_state_b = 1,
                            parameter_state_c = 1,
                            quality = TRUE, age = TRUE, bloodgroup = TRUE,
                            accept_prob = NULL, yearly_scale=rep(1, 100)
){

    recip_matrix = cbind(recip_matrix, allocation_score)
    recip_matrix = recip_matrix[with(recip_matrix, order(-allocation_score)), ]
    #print("State")

    if( nrow( recip_matrix) == 0 ) {
        return(NULL)
    }

    accept = 0
    i = 1
    recip = c()

    while(accept != graft_number ){


        ########### adding the yearly factor ##############

        days <- recip_matrix[i,"recip_sim_listdate"] - min(recip_matrix$recip_sim_listdate)
        #print(c("days", days))
        year = round(days/365+1)
        #print(c("year", year))

        ####################

        recip_age = round((donor_matrix[,"tx_date"]-recip_matrix[i,"recip_birthdate"])/365)
        age_difference = donor_matrix[,"donor_age"] - recip_age

        if (!is.null(accept_prob) ){
            #print(c("yearly_scale", yearly_scale[year]))
            p_accept = sample(accept_prob,1)*yearly_scale[year]
            #print(c(p_accept, "p_Accept state"))
        }

        if (is.null (accept_prob)){
            pra_val = recip_matrix[i,"recip_pra"]/100
            hla_val = 2*HLA[i,"tx_misdr"] + HLA[i,"tx_misa"] + HLA[i,"tx_misb"]
            #print(c("yearly_scale", yearly_scale[year]))
            p_accept = ((1-pra_val)*parameter_state_a)*(1-hla_val/50*parameter_state_b)*yearly_scale[year]
            #print(c(p_accept, "p_Accept"))
        }

        p_accept = min(1, p_accept)


        if(i == nrow(recip_matrix) ){

            return(NULL)

        }else{


            if( quality
                & recip_matrix[i, "recip_epts"] < donor_matrix[,"donor_kdri"] ){

                samp = rbinom(1,1,p_accept)


                if ( samp == 1 ) {
                    accept = accept + 1
                    recip = cbind(recip, i)

                } else if ( age
                            & recip_age < 18
                            & age_difference < 30 ) {

                    accept = accept + 1
                    recip = cbind(recip, i)

                } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                    accept = accept + 1
                    recip = cbind(recip, i)

                }

            }else{

                samp = rbinom(1,1,p_accept)


                if ( samp == 1 ) {
                    accept = accept + 1
                    recip = cbind(recip, i)

                } else if ( age
                            & recip_age < 18
                            & age_difference < 30 ) {

                    accept = accept + 1
                    recip = cbind( recip, i )

                } else if (bloodgroup & donor_matrix[,"donor_blgroup"] %in% "AB" ) {

                    accept = accept + 1
                    recip = cbind(recip, i)

                }
            }
        }

        i = i + 1

    }


    match_recip = as.vector(recip_matrix[recip,"recip_id"])

    #output the position on the list some other way
    #output 0 or 1 recipient

    return(match_recip)

}




#' selection
#'
#' this is the sdm2
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features, this is the input donor matrix
#' @param allocation_score a vector, this is the score for each pair
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @param graft_number a numerical value, this is the graft number
#' @param parameter_state_a an integer hyperparameter to tune
#' @param parameter_state_b an integer hyperparameter to tune
#' @param top_state_matches an integer, top match parameter
#' @param age a boolean value, TRUE/FALSE, age indicator
#' @param bloodgroup a boolean value, TRUE/FALSE, blood group indicator
#' @param parameter_state an integer, state parameter
#'
#'
#' @return a data.frame with matched recipient
#'
#' @examples
#' data("newdata", package = "simKAP")
#' HLA_matrix <- hla_match(raw_recip_matrix_subset, raw_donor_matrix[1,])
#'
#' selected <- dm_state_top(raw_recip_matrix_subset,raw_donor_matrix[1,],
#' allocation_score=allocation_national(raw_recip_matrix_subset,raw_donor_matrix[1,],HLA_matrix),
#' HLA_matrix,1,1,1.5,1.5,15,10);
#'
#' @export
#'


dm_state_top <- function(recip_matrix = NULL, #this matrix will only include people from same state
                        donor_matrix = NULL,
                        allocation_score = NULL,
                        HLA = NULL,
                        graft_number = 1,
                        parameter_state_a = 1, parameter_state_b = 1,
                        bloodgroup = 1.5,
                        age = 1.5,
                        top_state_matches = 15,
                        parameter_state = 10
                        #parameters can be set to 0

){
    #if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    #if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    #if (dim(HLA)[1]==0) stop("empty HLA matrix")
    #if (is.na(allocation_score)) stop("no allocation score")
    if (any(! c("recip_birthdate","recip_pra","recip_epts","recip_state","recip_id") %in% colnames(recip_matrix))) {
        stop("recipient matrix wrong column names")
    }
    if (any(! c("donor_state","tx_date","donor_age","donor_kdri","donor_blgroup") %in% colnames(donor_matrix))) {
        stop("donor matrix wrong column names")
    }
    if (any(! c("tx_misa","tx_misb","tx_misdr") %in% colnames(HLA))) {
        stop("HLA matrix wrong column names")
    }



    #print("state")

    HLAall = as.numeric(HLA$tx_misa) + as.numeric(HLA$tx_misb) + as.numeric(HLA$tx_misdr)

    recip_matrix = cbind(recip_matrix, allocation_score)
    recip_matrix = recip_matrix[with(recip_matrix, order(-allocation_score)), ]

    donor_state <- donor_matrix[,"donor_state"]
    donor_state <- switch(as.character(donor_state),
                          NSW = "NSW", ACT = "NSW",
                          VIC = "VIC", TAS = "VIC",
                          SA = "SA", NT = "SA",
                          WA = "WA", QLD = "QLD")


    if(nrow(recip_matrix) == 0){

        return(NULL)

    }else{

        min_nrow = min(top_state_matches, nrow(recip_matrix))
        recip_matrix = recip_matrix[c(1:min_nrow),]

        blood <- ifelse(donor_matrix[,"donor_blgroup"] %in% "AB", bloodgroup, 1 )

        recip_matrix$rank_score = rank(recip_matrix$allocation_score)*blood
        recip_matrix$rank_epts = rank(-recip_matrix$recip_epts)
        recip_matrix$rank_state = ifelse(recip_matrix$recip_state == donor_state,parameter_state,0)

        index <- which(as.numeric(as.character(recip_matrix$recip_age)) < 18)
        recip_matrix[index,]$rank_score <- recip_matrix[index,]$rank_score*age

        recip_matrix$rank_sum = recip_matrix$rank_score*parameter_state_a +
            recip_matrix$rank_epts*parameter_state_b +
            recip_matrix$rank_state

        ranked_recip = recip_matrix[with(recip_matrix, order(-rank_sum)), ]
        num <- min(nrow(ranked_recip),graft_number)
        match_recip = ranked_recip[c(1:num), "recip_id"]


        #i_rank = paste(i,as.vector(recip_matrix[i, "rank_sum"]), sep="-")
        # print(match_recip)
        # print(graft_number)
        # print( ranked_recip[c(1:graft_number), "rank_sum"])

        return(c(match_recip))
    }

}


