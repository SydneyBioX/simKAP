
#' allocation_national
#'
#' this is the Australia current national allocation algorithm
#' @param recip_matrix a data.frame with recipient features, this is the input recipient matrix
#' @param donor_matrix a data.frame with donor features,this is the input donor matrix
#' @param HLA a matrix, this is the calculated HLA mismatching counts matrix
#' @return a matrix contains allo_score, which is the alocation score for each new pair
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#'
#'
#' allocation_score <- allocation_national(rawdata,rawdata,rawdata);
#' @export



# score function national

allocation_national <- function(recip_matrix, donor_matrix, HLA){

    if (dim(recip_matrix)[1]==0) stop("empty recipient matrix")
    if (dim(donor_matrix)[1]==0) stop("empty donor matrix")
    if (dim(HLA)[1]==0) stop("empty HLA matrix")
    if (any(! c("recip_pra","recip_age","recip_rrtstartdate") %in% colnames(recip_matrix))) stop("recipient matrix wrong column names")
    if (any(! c("tx_date") %in% colnames(donor_matrix))) stop("donor matrix wrong column names")
    if (any(! c("tx_misa","tx_misb","tx_misdr") %in% colnames(HLA))) stop("HLA matrix wrong column names")




    HLAall = as.numeric(HLA$tx_misa) + as.numeric(HLA$tx_misb) + as.numeric(HLA$tx_misdr)

    changed = c(rep(0 , length(HLAall)))
    allo_score = c(rep(0, length(HLAall)))

    index <- HLAall == 0 & recip_matrix$recip_pra >= 50 & changed == 0
    allo_score[index] <- 60000000
    changed[index] <- 1

    index <-  HLAall == 1 & recip_matrix$recip_pra > 80 & changed  == 0
    allo_score[index] <- 59000000
    changed[index] <- 1



    index <- HLAall == 2 & recip_matrix$recip_pra > 80 & changed == 0
    allo_score[index] <- 58000000
    changed[index] <- 1


    index <- HLAall == 0 & recip_matrix$recip_pra < 50 & changed == 0
    allo_score[index] <- 57000000
    changed[index] <- 1

    index <- (HLA$tx_misdr == 0 & (HLA$tx_misa + HLA$tx_misb == 1) &
                  recip_matrix$recip_pra <= 80 & changed == 0)
    allo_score[index] <-  56000000
    changed[index] <- 1


    index <- (HLA$tx_misdr == 0 & (HLA$tx_misa + HLA$tx_misb) == 2  &
                  recip_matrix$recip_pra <= 80 & changed == 0)
    allo_score[index] <-  55000000
    changed[index] <- 1


    index <- changed == 0
    #allo_score[index] <- 54000000
    allo_score[index] <- 0
    changed[index] <- 1



    index <- recip_matrix$recip_age < 18
    allo_score[index]  = allo_score[index] + 30000

    waiting_time <- as.numeric(donor_matrix$tx_date - recip_matrix$recip_rrtstartdate)/365*12

    #index <- recip_matrix$recip_waittime/30 > 0
    #allo_score[index] = allo_score[index] + recip_matrix$recip_waittime[index]*1

    index <- waiting_time > 0
    allo_score[index] = allo_score[index] + waiting_time[index]*1

    return(allo_score)
}




