
#' state algorithm overall
#'
#' this is the Australia current state allocation algorithm
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features,this is the input donor matrix
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @param state a categorical variable representing state names, specify which state
#' @return a matrix contains allo_score, which is the alocation score for each new pair
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#'
#'
#' state_score <- australia_state_algorithm(rawdata,rawdata,rawdata,"NSW");
#' @export
#'
#'





australia_state_algorithm <- function(recip_matrix,
                                      donor_matrix,
                                      HLA,
                                      state = NULL) {
    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (dim(HLA)[1]==0) stop("empty HLA matrix")
    if (any(! c("recip_state","recip_age","recip_liststartdate","recip_blgroup") %in% colnames(recip_matrix))) {
        stop("recipient matrix wrong column names")
    }
    if (any(! c("tx_date","donor_state","donor_blgroup") %in% colnames(donor_matrix))) {
        stop("donor matrix wrong column names")
    }


    if (is.null(state)) {
        state <- donor_matrix$donor_state
    }

    state <- match.arg(toupper(state),
                       c("NSW", "ACT", "NT", "QLD", "SA", "TAS", "VIC", "WA"))

    state_algorithm_FUN <- switch(state,
                                  NSW = NSWACT_allocation,
                                  ACT = NSWACT_allocation,
                                  NT = SA_allocation,
                                  QLD = QLD_allocation,
                                  SA = SA_allocation,
                                  TAS = VICTAS_allocation,
                                  VIC = VICTAS_allocation,
                                  WA = WA_allocation)


    allocation_scores <- state_algorithm_FUN(recip_matrix,
                                             donor_matrix,
                                             HLA)
    return(allocation_scores)


}


#' Australia_state_selectionpool
#'
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with recipient features, this is the input donor matrix
#' @param whichdonor a categorical donor id, this is the incoming kidney from one donor
#' @param bloodGroup_strict_match a boolean value, TRUE/FALSE, this is to define whether we want the strict match for bloodgroup or not
#' @return a data.frame representing recipient pool under the state allocation tule
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#' state_pool <- australia_state_selection(rawdata,rawdata,"D948826",
#' bloodGroup_strict_match = TRUE);
#' @export


australia_state_selection <- function(recip_matrix,
                                          donor_matrix,
                                          whichdonor,
                                          bloodGroup_strict_match = TRUE){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (any(! c("recip_state","recip_age","recip_liststartdate","recip_blgroup") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("tx_date","donor_state","donor_blgroup") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")


    # Step 1) Eligible sub-recipients are those who started waiting before the donor became available
    transplantdate <- donor_matrix[whichdonor, "tx_date"]
    donor_state <- as.character(donor_matrix[whichdonor, "donor_state"])
    recipient_state <- as.character(recip_matrix$recip_state)

    donor_state <- switch(donor_state,
                          NSW = "NSW", ACT = "NSW",
                          VIC = "VIC", TAS = "VIC",
                          SA = "SA", NT = "SA",
                          WA = "WA", QLD = "QLD")

    recipient_state <- sapply(recipient_state, function(x) switch(as.character(x),
                                                                  NSW = "NSW", ACT = "NSW",
                                                                  VIC = "VIC", TAS = "VIC",
                                                                  SA = "SA", NT = "SA",
                                                                  WA = "WA", QLD = "QLD"))

    bg_match <- blood_group_match(donor_matrix$donor_blgroup[whichdonor],
                                recip_matrix$recip_blgroup,
                                strict_match = bloodGroup_strict_match)

    eligible_recipM <- recip_matrix[recip_matrix$recip_liststartdate < transplantdate &
                                        recipient_state == donor_state &
                                        bg_match, ]


    return(eligible_recipM)
}



#' state algorithm NSWACT
#'
#' this is the Australia current state allocation algorithm
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features,this is the input donor matrix
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @return a matrix contains allo_score, which is the alocation score for each new pair
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#'
#'
#' state_score <- NSWACT_allocation(rawdata,rawdata,rawdata);
#' @export

















NSWACT_allocation <- function(recip_matrix,
                              donor_matrix,
                              HLA) {
    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (dim(HLA)[1]==0) stop("empty HLA matrix")
    if (any(! c("recip_state","recip_age","recip_liststartdate","recip_blgroup","recip_rrtstartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("tx_date","donor_state","donor_blgroup") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")


    HLA_A <- HLA$tx_misa
    HLA_B <- HLA$tx_misb
    HLA_DR <- HLA$tx_misdr

    HLAall = HLA_A + HLA_B + HLA_DR

    allo_score = c(rep(0, length(HLAall)))

    index <- HLA_DR == 0
    allo_score[index] <- 50000000
    allo_score[index] <- allo_score[index] - 1000000*(HLA_A[index] + HLA_B[index])

    index <- recip_matrix[, "recip_age"] < 18
    allo_score[index] <- allo_score[index] + 100000

    waiting_month <- as.numeric(donor_matrix$tx_date - recip_matrix$recip_rrtstartdate)/365*12
    index <- waiting_month > 0
    allo_score[index] <- allo_score[index] + waiting_month[index]*100

    index_waiting <- allo_score < 48000000

    if (sum(index_waiting) > 0) {
        allo_score[index_waiting] <- 40000000



        index <- recip_matrix[index_waiting, "recip_age"] < 18

        if (sum(index) > 0) {
            allo_score[index_waiting][index] <- allo_score[index_waiting][index] + 100000
        }


        index <- waiting_month[index_waiting] > 0

        if (sum(index) > 0) {
            allo_score[index_waiting][index] <- allo_score[index_waiting][index] +
                waiting_month[index_waiting][index]*100
        }

    }


    return(allo_score)

}





#' state algorithm QLD
#'
#' this is the Australia current state allocation algorithm
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features,this is the input donor matrix
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @return a matrix contains allo_score, which is the alocation score for each new pair
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#'
#'
#' state_score <- QLD_allocation(rawdata,rawdata,rawdata);
#' @export




QLD_allocation <- function(recip_matrix,
                           donor_matrix,
                           HLA) {
    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (dim(HLA)[1]==0) stop("empty HLA matrix")
    if (any(! c("recip_state","recip_age","recip_liststartdate","recip_blgroup","recip_rrtstartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("tx_date","donor_state","donor_blgroup") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")


    HLA_A <- HLA$tx_misa
    HLA_B <- HLA$tx_misb
    HLA_DR <- HLA$tx_misdr

    HLAall = HLA_A + HLA_B + HLA_DR

    allo_score <- c(rep(50000000, length(HLAall)))



    allo_score <- allo_score - 1000000*HLAall


    waiting_month <- as.numeric(donor_matrix$tx_date - recip_matrix$recip_rrtstartdate)/365*12
    index <- waiting_month > 0
    allo_score[index] <- allo_score[index] + waiting_month[index]*100

    index_waiting <- allo_score < 48000000

    if (sum(index_waiting) > 0) {
        allo_score[index_waiting] <- 40000000



        index <- waiting_month[index_waiting] > 0
        allo_score[index_waiting][index] <- allo_score[index_waiting][index] +
            waiting_month[index_waiting][index]*100
    }


    return(allo_score)

}






#' state algorithm SA
#'
#' this is the Australia current state allocation algorithm
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features,this is the input donor matrix
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @return a matrix contains allo_score, which is the alocation score for each new pair
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#'
#'
#' state_score <- SA_allocation(rawdata,rawdata,rawdata);
#' @export
#'
#'



SA_allocation <- function(recip_matrix,
                          donor_matrix,
                          HLA) {
    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (dim(HLA)[1]==0) stop("empty HLA matrix")
    if (any(! c("recip_state","recip_age","recip_liststartdate","recip_blgroup","recip_rrtstartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("tx_date","donor_state","donor_blgroup") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")


    HLA_A <- HLA$tx_misa
    HLA_B <- HLA$tx_misb
    HLA_DR <- HLA$tx_misdr
    HLA = HLA_A + HLA_B + HLA_DR

    allo_score <- c(rep(30000000, length(HLA))) #allocation score


    #STATE HLA:
    allo_score <- allo_score - HLA * 10000000

    # RESET
    index <- HLA > 3
    allo_score[index] <- 0


    waiting_month <- as.numeric(donor_matrix$tx_date - recip_matrix$recip_rrtstartdate)/365*12
    index <- waiting_month > 0
    allo_score[index] <- allo_score[index] + waiting_month[index]*1



    return(allo_score)

}






#' state algorithm VICTAS
#'
#' this is the Australia current state allocation algorithm
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features,this is the input donor matrix
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @return a matrix contains allo_score, which is the alocation score for each new pair
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#'
#'
#' state_score <- VICTAS_allocation(rawdata,rawdata,rawdata);
#' @export
#'
#'
#'



VICTAS_allocation <- function(recip_matrix,
                              donor_matrix,
                              HLA) {
    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (dim(HLA)[1]==0) stop("empty HLA matrix")
    if (any(! c("recip_state","recip_age","recip_liststartdate","recip_blgroup","recip_rrtstartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("tx_date","donor_state","donor_blgroup") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")


    HLA_A <- HLA$tx_misa
    HLA_B <- HLA$tx_misb
    HLA_DR <- HLA$tx_misdr
    HLA = HLA_A + HLA_B + HLA_DR

    allo_score <- c(rep(40000000, length(HLA_A))) #allocation score
    allo_score <- allo_score - 20000000 * (HLA_B + HLA_DR)

    index <- recip_matrix[, "recip_age"] < 18
    allo_score[index] <- allo_score[index] + 100000

    index <- HLA_B + HLA_DR > 2
    allo_score[index] <- 0

    waiting_month <- as.numeric(donor_matrix$tx_date - recip_matrix$recip_rrtstartdate)/365*12
    index <- waiting_month > 0
    allo_score[index] <- allo_score[index] + waiting_month[index]*1

    return(allo_score)

}






#' state algorithm WA
#'
#' this is the Australia current state allocation algorithm
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features,this is the input donor matrix
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @return a matrix contains allo_score, which is the alocation score for each new pair
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#'
#'
#' state_score <- WA_allocation(rawdata,rawdata,rawdata);
#' @export
#'
#'
#'



WA_allocation <- function(recip_matrix,
                          donor_matrix,
                          HLA) {
    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (dim(HLA)[1]==0) stop("empty HLA matrix")
    if (any(! c("recip_state","recip_age","recip_liststartdate","recip_blgroup","recip_rrtstartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("tx_date","donor_state","donor_blgroup") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")


    HLA_A <- HLA$tx_misa
    HLA_B <- HLA$tx_misb
    HLA_DR <- HLA$tx_misdr
    HLA <- HLA_A + HLA_B + HLA_DR

    allo_score <- c(rep(40000000, length(HLA))) #allocation score

    allo_score <- allo_score - HLA_A * 3000000 - HLA_B * 3000000 - HLA_DR * 5000000


    waiting_month <- as.numeric(donor_matrix$tx_date - recip_matrix$recip_rrtstartdate)/365*12
    index <- waiting_month > 0
    allo_score[index] <- allo_score[index] + waiting_month[index]*100000

    #Homozygous at HLA-DR
    index <- (recip_matrix$recip_dr1 == recip_matrix$recip_dr2) &
        waiting_month > 60
    allo_score[index] <- allo_score[index] + 5000000

    return(allo_score)

}




#' state balance function
#'
#' this is the function to obtain the number of kidneys allocated to another state
#' @param SimResults a data.frame, this is the current simulation results
#' @param recip_matrix a data.frame, this is the recipient matrix
#' @param donor_matrix a data.frame, this is the donor vector
#'
#' @return a vector, which contains all the kidneys allocated to another state
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#' data("newdata", package = "KidneyAllocation")
#'
#' state_balance <- state_balance(rawdata,raw_recip_matrix_subset,raw_donor_matrix[1,]);
#' @export
#'
#'


###################################
# STATE-BALANCE ALLOCATED KIDNEYS #
###################################

state_balance = function(SimResults,
                        recip_matrix = NULL,
                        donor_matrix = NULL
) {


    #count the number of kidney allocated in each state

    # print(lubridate::month(SimResults$tx_date))
    index <- lubridate::year(SimResults$tx_date) == lubridate::year(donor_matrix$tx_date)

    recipient_state = SimResults[index,]$recip_state
    recipient_state = factor(recipient_state, levels = c("NSW", "VIC", "QLD", "SA", "WA" ))
    donor_state = SimResults[index,]$donor_state
    donor_state = factor(donor_state, levels = c("NSW", "VIC", "QLD", "SA", "WA" ))

    match_table = table(recipient_state, donor_state)
    recip_number = rowSums(match_table)
    donor_number = colSums(match_table)
    excess = as.matrix((donor_number - recip_number)[which(donor_number - recip_number > 0)])

    if(length(excess) == 1){
        excess_vector = apply(excess , 1, function(x) rep(rownames(excess), excess[,1]))
    }else{
        excess_vector = apply(excess , 1, function(x) rep(rownames(excess), excess[,1]))[,1]}

    return(excess_vector)

}


