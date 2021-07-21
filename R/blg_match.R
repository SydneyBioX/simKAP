#' blood group match
#'
#' this is the blood group match function
#' @param donor_bl a character value,this is the donor blood group
#' @param recipient_bl a character value,this is the irecipient blood group
#' @param strict_match a boolean value, TRUE/FALSE, this is the strict match indicator
#' @param AB_priority a boolean value, TRUE/FALSE, this is the AB priority indicator
#' @return a data.frame with matched recipient and donor
#'
#' @examples
#' data("rawdata", package = "KidneyAllocation")
#'
#'
#' blgmatch <- blood_group_match(rawdata$donor_blgroup,rawdata$recip_blgroup,TRUE,TRUE);
#' @export



# blg match




blood_group_match <- function(donor_bl,
                            recipient_bl,
                            strict_match = TRUE,
                            AB_priority = TRUE) {
    if (length(donor_bl)==0) stop("empty donor vector")
    if (length(recipient_bl)==0) stop("empty recipient vector")
    if (any(!unique(donor_bl) %in%  c("A","B","AB","O",NA) )) stop("wrong blood type")
    if (any(!unique(recipient_bl) %in%  c("A","B","AB","O",NA) )) stop("wrong blood type")


    recipient_bl <- as.character(recipient_bl)
    donor_bl <- as.character(donor_bl)



    if (strict_match) {

        if (AB_priority) {
            res <- recipient_bl %in%  c(donor_bl, "AB")

        } else {

            res <- recipient_bl %in%  donor_bl
        }

    } else {
        if (donor_bl == "A") {
            res <- recipient_bl %in% c("A", "AB")
        }

        if (donor_bl == "B") {
            res <- recipient_bl %in% c("B", "AB")
        }

        if (donor_bl == "O") {
            res <- recipient_bl %in% c("A", "B", "O", "AB")
        }

        if (donor_bl == "AB") {
            res <- recipient_bl %in% c("AB")
        }
    }

    return(res)
}

#
# bloodGroupMatch("A", c("A", "B", "B", "A", "AB", "O"))
# bloodGroupMatch("B", "A")
# bloodGroupMatch("B", "B")
# bloodGroupMatch("B", "AB")
#
#
# bloodGroupMatch("B", "AB", strict_match = FALSE)
# bloodGroupMatch("O", "AB", strict_match = FALSE)
# bloodGroupMatch("AB", "O", strict_match = FALSE)
# bloodGroupMatch("A", c("A", "B", "B", "A", "AB", "O"), strict_match = FALSE)
# bloodGroupMatch("O", c("A", "B", "B", "A", "AB", "O"), strict_match = FALSE)
#
