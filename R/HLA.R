
#' HLA functions
#'
#' this is the function for calculating HLA
#' @param Rmatrix a matrix, this is the recipient matrix
#' @param Dvector a matrix, this is the donor matrix
#' @return HLA mismatching matrix
#'
#' @examples
#' data("rawdata", package = "simKAP")
#'
#' HLA_matrix <- hla_match(rawdata, rawdata);
#'
#' @export




hla_match <- function(Rmatrix, Dvector){
    if (dim(Rmatrix)[1]==0) stop("empty recipient matrix")
    if (length(Dvector)==0) stop("empty donor vector")

    if (any(! c("recip_a1","recip_a2","recip_b1","recip_b2","recip_dr1","recip_dr2") %in% colnames(Rmatrix))) stop("recipient matrix wrong column names")
    if (any(! c("donor_a1","donor_a2","donor_b1","donor_b2","donor_dr1","donor_dr2") %in% names(Dvector))) stop("donor vector wrong  names")



    rraw = cbind(Rmatrix,Dvector)

    #misa
    col1 = c(2,9,24,28)
    col2 = c(203,23,2403,68)
    col3 = c(210,24,2403,69)
    col4 = c(25,26,34,66)
    col5 = c(29,30,31,32,33,74)
    vec1 = rraw$donor_a1
    vec2 = rraw$donor_a2
    vec3 = rraw$recip_a1
    vec4 = rraw$recip_a2


        A1 = A2 = misa = rep(0, nrow(rraw))
        for (i in 1:nrow(rraw)) {
            #print(i)
            if (vec3[i] %in% col1) {

                res <- ((col2[which(col1 == vec3[i])]) != vec1[i]) &
                    ((col3[which(col1 == vec3[i])]) != vec1[i]) &
                    ((col2[which(col1 == vec3[i])]) != vec2[i]) &
                    ((col3[which(col1 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) &
                    (vec2[i] != vec3[i])

                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                A1[i] <- unique(res)
            }else if ((vec3[i] %in% col2)) {

                res <- ((col1[which(col2 == vec3[i])]) != vec1[i]) &
                    ((col1[which(col2 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) & (vec2[i] != vec3[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                A1[i] <- unique(res)

            } else if ((vec3[i] %in% col3)) {

                res <- ((col1[which(col3 == vec3[i])]) != vec1[i]) &
                    ((col1[which(col3 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) &
                    (vec2[i] != vec3[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                A1[i] <- unique(res)

            }else if ((vec3[i] %in% col4)) {

                A1[i] = (10 != vec1[i]) & (10 != vec2[i]) &
                    (vec1[i] != vec3[i]) & (vec2[i] != vec3[i])

            }else if ((vec3[i] %in% col5)) {

                A1[i] = (19 != vec1[i]) &  (19 != vec2[i]) &
                    (vec1[i] != vec3[i]) & (vec2[i] != vec3[i])

            }else if (vec3[i] == 10) {

                A1[i] = !(vec1[i] %in% col4) &  !(vec2[i] %in% col4) &
                    (vec1[i] != 10) &  (vec2[i] != 10)

            }else if (vec3[i] == 19) {

                A1[i] = !(vec1[i] %in% col5) &  !(vec2[i] %in% col5) &
                    (vec1[i] != 19) &  (vec2[i] != 19)

            }else {

                A1[i] = (vec3[i] != vec1[i]) &  (vec3[i] != vec2[i])

            }
        }
        for (i in 1:nrow(rraw)) {
            #print(i)
            if (vec4[i] %in% col1) {

                res <- ((col2[which(col1 == vec4[i])]) != vec1[i]) &
                    ((col3[which(col1 == vec4[i])]) != vec1[i]) &
                    ((col2[which(col1 == vec4[i])]) != vec2[i]) &
                    ((col3[which(col1 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                A2[i] <- unique(res)

            } else if ((vec4[i] %in% col2)) {

                res <- ((col1[which(col2 == vec4[i])]) != vec1[i]) &
                    ((col1[which(col2 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                A2[i] <- unique(res)

            } else if ((vec4[i] %in% col3)) {

                res <- ((col1[which(col3 == vec4[i])]) != vec1[i]) &
                    ((col1[which(col3 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                A2[i] <- unique(res)

            }else if ((vec4[i] %in% col4)) {

                A2[i] = (10 != vec1[i]) & (10 != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])

            }else if ((vec4[i] %in% col5)) {

                A2[i] = (19 != vec1[i]) & (19 != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])

            }else if (vec4[i] == 10) {

                A2[i] = !(vec1[i] %in% col4) &  !(vec2[i] %in% col4) &
                    (vec1[i] != 10) &  (vec2[i] != 10)

            }else if (vec4[i] == 19) {

                A2[i] = !(vec1[i] %in% col5) & !(vec2[i] %in% col5) &
                    (vec1[i] != 19) &  (vec2[i] != 19)

            }else{
                A2[i] = (vec4[i] != vec1[i]) & (vec4[i] != vec2[i])
            }
        }

        misa = A1 + A2
        for (i in 1:nrow(rraw)) {
            if ((vec3[i] == vec4[i]) & (misa[i] > 1)) {
                misa[i] = 1
            }
        }


#misb
        col1 = c(5,7,12,14,16,17,21,4005,27,39,40,51,70)
        col2 = c(51,703,44,64,38,57,49,49,2708,3901,60,5102,71)
        col3 = c(52,703,45,65,39,58,50,50,2708,3902,61,5103,72)
        col4 = c(62,63,75,76,77)
        col5 = c(54,55,56)

        vec1 = rraw$donor_b1
        vec2 = rraw$donor_b2
        vec3 = rraw$recip_b1
        vec4 = rraw$recip_b2

        B1 = B2 = misb = rep(0, nrow(rraw))
        for (i in 1:nrow(rraw)) {
            #print(i)
            if (vec3[i] %in% col1) {

                res <- ((col2[which(col1 == vec3[i])]) != vec1[i]) &
                    ((col3[which(col1 == vec3[i])]) != vec1[i]) &
                    ((col2[which(col1 == vec3[i])]) != vec2[i]) &
                    ((col3[which(col1 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) &
                    (vec2[i] != vec3[i])

                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                B1[i] <- unique(res)
            }else if ((vec3[i] %in% col2)) {

                res <- ((col1[which(col2 == vec3[i])]) != vec1[i]) &
                    ((col1[which(col2 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) & (vec2[i] != vec3[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                B1[i] <- unique(res)

            } else if ((vec3[i] %in% col3)) {

                res <- ((col1[which(col3 == vec3[i])]) != vec1[i]) &
                    ((col1[which(col3 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) &
                    (vec2[i] != vec3[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                B1[i] <- unique(res)

            }else if ((vec3[i] %in% col4)) {

                B1[i] = (15 != vec1[i]) &  (15 != vec2[i]) &
                    (vec1[i] != vec3[i]) & (vec2[i] != vec3[i])

            }else if ((vec3[i] %in% col5)) {

                B1[i] = (22 != vec1[i]) &  (22 != vec2[i]) &
                    (vec1[i] != vec3[i]) & (vec2[i] != vec3[i])

            }else if (vec3[i] == 15) {

                B1[i] = !(vec1[i] %in% col4) &  !(vec2[i] %in% col4) &
                    (vec1[i] != 15) &  (vec2[i] != 15)

            }else if (vec3[i] == 22) {

                B1[i] = !(vec1[i] %in% col5) &  !(vec2[i] %in% col5) &
                    (vec1[i] != 22) &  (vec2[i] != 22)

            }else {

                B1[i] = (vec3[i] != vec1[i]) &  (vec3[i] != vec2[i])

            }
        }
        for (i in 1:nrow(rraw)) {
            #print(i)
            if (vec4[i] %in% col1) {

                res <- ((col2[which(col1 == vec4[i])]) != vec1[i]) &
                    ((col3[which(col1 == vec4[i])]) != vec1[i]) &
                    ((col2[which(col1 == vec4[i])]) != vec2[i]) &
                    ((col3[which(col1 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                B2[i] <- unique(res)

            } else if ((vec4[i] %in% col2)) {

                res <- ((col1[which(col2 == vec4[i])]) != vec1[i]) &
                    ((col1[which(col2 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                B2[i] <- unique(res)

            } else if ((vec4[i] %in% col3)) {

                res <- ((col1[which(col3 == vec4[i])]) != vec1[i]) &
                    ((col1[which(col3 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
                if (length(res) > 1 & sum(res) == 1) {stop("check")}
                B2[i] <- unique(res)

            }else if ((vec4[i] %in% col4)) {

                B2[i] = (15 != vec1[i]) & (15 != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])

            }else if ((vec4[i] %in% col5)) {

                B2[i] = (22 != vec1[i]) & (22 != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])

            }else if (vec4[i] == 15) {

                B2[i] = !(vec1[i] %in% col4) &  !(vec2[i] %in% col4) &
                    (vec1[i] != 15) &  (vec2[i] != 15)

            }else if (vec4[i] == 22) {

                B2[i] = !(vec1[i] %in% col5) & !(vec2[i] %in% col5) &
                    (vec1[i] != 22) &  (vec2[i] != 22)

            }else{
                B2[i] = (vec4[i] != vec1[i]) & (vec4[i] != vec2[i])
            }
        }

        misb = B1 + B2
        for (i in 1:nrow(rraw)) {
            if ((vec3[i] == vec4[i]) & (misb[i] > 1)) {
                misb[i] = 1
            }
        }





#misdr

        col1 = c(1,2,3,5,6,14)
        col2 = c(103,15,17,11,13,1403)
        col3 = c(103,16,18,12,14,1404)
        vec3 = rraw$donor_dr1
        vec4 = rraw$donor_dr2
        vec1 = rraw$recip_dr1
        vec2 = rraw$recip_dr2



        DR1 = DR2 = misdr = rep(0, nrow(rraw))
        for (i in 1:nrow(rraw)) {
            #print(i)
            if (vec3[i] %in% col1) {
                DR1[i] = ((col2[which(col1 == vec3[i])]) != vec1[i]) &
                    ((col3[which(col1 == vec3[i])]) != vec1[i]) &
                    ((col2[which(col1 == vec3[i])]) != vec2[i]) &
                    ((col3[which(col1 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) &
                    (vec2[i] != vec3[i])
                # print(paste("dr1=",DR1[i]))
            } else if ((vec3[i] %in% col2)) {
                DR1[i] = ((col1[which(col2 == vec3[i])]) != vec1[i]) &
                    ((col1[which(col2 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) &
                    (vec2[i] != vec3[i])
            } else if ((vec3[i] %in% col3)) {
                DR1[i] = ((col1[which(col3 == vec3[i])]) != vec1[i]) &
                    ((col1[which(col3 == vec3[i])]) != vec2[i]) &
                    (vec1[i] != vec3[i]) &
                    (vec2[i] != vec3[i])
            } else {
                DR1[i] = (vec3[i] != vec1[i]) &  (vec3[i] != vec2[i])
            }
        }


        for (i in 1:nrow(rraw)) {
            if (vec4[i] %in% col1) {
                DR2[i] = ((col2[which(col1 == vec4[i])]) != vec1[i]) &
                    ((col3[which(col1 == vec4[i])]) != vec1[i]) &
                    ((col2[which(col1 == vec4[i])]) != vec2[i]) &
                    ((col3[which(col1 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
            } else if ((vec4[i] %in% col2)) {
                DR2[i] = ((col1[which(col2 == vec4[i])]) != vec1[i]) &
                    ((col1[which(col2 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
            }else if ((vec4[i] %in% col3)) {
                DR2[i] = ((col1[which(col3 == vec4[i])]) != vec1[i]) &
                    ((col1[which(col3 == vec4[i])]) != vec2[i]) &
                    (vec1[i] != vec4[i]) &
                    (vec2[i] != vec4[i])
            }else{
                DR2[i] = (vec4[i] != vec1[i]) &  (vec4[i] != vec2[i])
            }
        }

        misdr = DR1 + DR2

        for (i in 1:nrow(rraw)) {
            if ((vec3[i] == vec4[i]) & (misdr[i] > 1)) {
                misdr[i] = 1
            }
        }




    HLAmat = data.frame(tx_misa = misa, tx_misb = misb, tx_misdr = misdr)
    colnames(HLAmat) <- c("tx_misa", "tx_misb", "tx_misdr")
    return(HLAmat)
}

