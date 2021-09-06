#' wl function
#'
#' this is the dynamic wl function
#' @param yearly_scale_in a numeirc vector between 0 and 1, this is the input yearly scale rate
#' @param yearly_scale_out a numeirc vector between 0 and 1, this is the input yearly scale rate
#' @param rate_in a numeirc value between 0 and 1, this is the rate
#' @param rate_out a numeirc value between 0 and 1, this is the rate
#' @param donor_matrix a data.frame, this is the input donor matrix
#' @param IN a boolean value, TRUE/FALSE, whether to use the user defined rate_in or not
#' @param OUT a boolean value, TRUE/FALSE, whether to use the user defined rate_out or not
#' @param waitlist_Risk a boolean value, TRUE/FALSE, whether to have this or not
#'
#' @return wl data which is a matrix for all recipients on the waiting list
#'
#' @examples
#' data("rawdata", package = "simKAP")
#'
#'
#' wl <- dynamic_waitlist(rawdata,yearly_scale_in = NULL,yearly_scale_out = NULL,rate_in = 0.7,
#' IN= TRUE, rate_out = 0.4,waitlist_Risk = TRUE);
#'
#' @import stats
#'
#'
#'
#' @export




dynamic_waitlist <- function(donor_matrix,
                             yearly_scale_in = NULL,
                             yearly_scale_out = NULL,
                             rate_in = 0.7,
                             rate_out = 0.4,
                             IN = FALSE,
                             OUT = FALSE,
                             waitlist_Risk = TRUE
){
    if (IN==FALSE & OUT == FALSE) stop("wrong input parameter")
    if (IN==TRUE & OUT == TRUE) stop("wrong input parameter")
    first_tx_date <- min(donor_matrix$tx_date) #"2006-07-02"
    final_tx_date <- max(donor_matrix$tx_date) #"2017-12-31"

    tx_date_years <- as.numeric(final_tx_date - first_tx_date)/365 #all together, 11 years
    tx_date_ranges <- split(sort(donor_matrix$tx_date),
                            cut(seq_along(sort(donor_matrix$tx_date)), tx_date_years, labels = FALSE))
# group all into 11 groups, around 540 transplant dates in each group
    if (is.null(yearly_scale_in) & is.null(yearly_scale_out)){
        yearly_scale_in = c(sort(1 + stats::rexp(length(tx_date_ranges), 1/0.1))) #a vector of length 11, generated randomly from the exp distribution
        yearly_scale_out = c(sort(stats::runif(length(tx_date_ranges), 0.3, 1), decreasing = TRUE)) #a vector of length 11,
        #generated randomly from the uniform distribution and also ordered from the largest to the smallest
        #between 0.3 and 1
    }


    if (IN ){

        rate_parameter = rate_in
        yearly_scale_rate = yearly_scale_in
        #print(c("in", rate_parameter, yearly_scale_rate))

    } else if (OUT){
        #print("out")
        rate_parameter = rate_out
        yearly_scale_rate = yearly_scale_out

    } else {
        print("Please specify rate parameter for waitlist")
    }

    wl_list = NULL

    for (i in seq(1:length(tx_date_ranges))){
        date_range <- tx_date_ranges[[i]] #for example, for those kidneys transplanted in the "first year"
        rate = rate_parameter*yearly_scale_rate[i] #get a random increase number of the wl,
        #using the rate_param with the randomly generated rate vector for each year (this could be the same rate for each year,
        #such as a vector of length 11 with 1.1 as element, which means, the final rate is 0.7*1.1 all the time)
        diff_tx_date = as.numeric(max(date_range) - min(date_range))
        #print(diff_tx_date), the number of total transplantation
        jumps_number <- stats::rpois(1, lambda = rate * diff_tx_date) #the added number of people on the wl
        jumps_time <- round(stats::runif(n = jumps_number, min = 0, max = diff_tx_date)) #randomly generate the added time of those added people
        #this has to be within the period ("first year")
        jump_sum <- table(jumps_time)
        wl_list_sub <-  rep(0, diff_tx_date)
        names(wl_list_sub) <- c(1:diff_tx_date)
        index <- which(names(wl_list_sub) %in% names(jump_sum))
        #extract those above selected patients on the wl to make up those added patients for the fake wl
        wl_list_sub[index] <- jump_sum[1:length(index)]
        names(wl_list_sub)<- seq(min(date_range) + 1, max(date_range), 1)
        wl_list = c(wl_list, wl_list_sub)
    }

    waitlist_Risk = waitlist_Risk

    return(list(wl_list, waitlist_Risk))
}
