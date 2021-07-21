#' simulation function
#'
#' this is the simulation function
#' @param recip_matrix a matrix, indicating the recipient matrix
#' @param donor_matrix a matrix, indicating the donor matrix
#' @param algorithm_FUN a function, this is the allocation algorithm
#' @param eligible_FUN a function, this is the eligible recipient pool function
#' @param matching_FUN a function, this is the selection function
#' @param waitlist_FUN a function,this is the waiting list function
#' @param national_algorithm_threshold an integer value, gives the threshold number to run national-state algorithm
#' @param state_algorithm a boolean value, TRUE/FALSE, whether to have state algorithm or not
#' @param state_algorithm_FUN a function, state algorithm function
#' @param state_eligible_FUN a function, state eligible function
#' @param state_matching_FUN a function, state matching function
#' @param resampleN an integer value, number of simulation results
#' @param num_donor an integer value, number of donor
#' @param num_recip an integer value, number of recipients
#' @param is_parallel a boolean value, TRUE/FALSE, to run paralel or not
#' @param ncores an integer value, number of cores
#' @param verbose a boolean value, TRUE/FALSE, verbose indicator
#' @param waitlist_arg a list,arguments for the wl
#' @param eligible_arg a list,arguments for the eligible pool
#' @param state_eligible_arg a list,arguments for the state eligible pool
#' @param matching_arg a list,arguments for matching functions
#' @param state_matching_arg a list,arguments for state matching functions
#' @param state_balance a boolean value, TRUE/FALSE, state balance adjustment indicator
#' @param dynamic_waitlist a boolean value, TRUE/FALSE, whether to have dynamic wl or not
#' @param National_II a boolean value, TRUE/FALSE, national algorithm round two indicator
#'
#'
#'
#'
#' @return a list with 2 objects, a data.frame with matched recipient and donor features, another data.frame with the discarded kidneys if any
#'
#' @examples
#' data("newdata", package = "KidneyAllocation")
#'
#'
#'simulation_result <-run_simulation(raw_recip_matrix_subset, raw_donor_matrix,
#'algorithm_FUN = allocation_national,eligible_FUN = selection_default,
#'matching_FUN = dm_national_formula,waitlist_FUN = dynamic_waitlist,
#'state_algorithm = TRUE,national_algorithm_threshold = 54000000,
#'state_algorithm_FUN = australia_state_algorithm,
#'state_eligible_FUN = australia_state_selection,
#'state_matching_FUN = dm_state_formula,eligible_arg = list(AB_priority = TRUE),
#'resampleN = 1,
#'  waitlist_arg = list(waitlist_Risk = TRUE,rate_in = 0.70,rate_out = 0.40),,
#'verbose = FALSE,num_donor = 800,num_recip = 300,
#'state_balance = TRUE,dynamic_waitlist = TRUE, National_II = TRUE);
#'
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom lubridate month year
#'
#' @export

# This function is the poisson process

# This function will rerun national algorithm relaxing the strict blood group match
# if a kidney can not be matched with the same blood group
# Need to think of a way to generalise this

run_simulation <- function(recip_matrix, donor_matrix,
                             algorithm_FUN, eligible_FUN, matching_FUN,
                             waitlist_FUN,
                             national_algorithm_threshold = 54000000,
                             state_algorithm = TRUE,
                             state_algorithm_FUN = NULL,
                             state_eligible_FUN = NULL,
                             state_matching_FUN = NULL,
                             resampleN = 20,
                             num_donor = 800,
                             num_recip = 300,
                             is_parallel = FALSE,
                             ncores = 2,
                             verbose = TRUE,
                             waitlist_arg = list(waitlist_Risk = TRUE,
                                                 rate_in = 0.70,
                                                 rate_out = 0.40,
                                                 yearly_scale_in = c(sort(1+rexp(length(tx_date_ranges), 1/0.1))),
                                                 yearly_scale_out = c(sort(runif(length(tx_date_ranges), 0.3, 1), decreasing = TRUE))),
                             eligible_arg = list(AB_priority = TRUE),
                             state_eligible_arg = list(),
                             matching_arg = list(),
                             state_matching_arg = list(),
                             state_balance = TRUE,
                             dynamic_waitlist = TRUE,
                             National_II = TRUE) {

    SimResults = NULL
    discard_kidney=NULL
    current_nrow = 0
    statebalance_kidney = NULL
    excess_vector = NULL
    select_sb_idx <- NULL
    dynamics_matrix <- NULL

    # function to rbind each dataframe within a list , when using foreach
    # (because the simulation returns 2 dataframes )
    comb <- function(...) {
        mapply('rbind', ..., SIMPLIFY=FALSE)
    }

    if (any(!colnames(donor_matrix)%in% c("donor_source","donor_state","donor_blgroup",
                                            "donor_a1", "donor_a2","donor_b1","donor_b2"
                                            ,"donor_dr1" ,"donor_dr2"
                                             ,"donor_dq1" ,"donor_dq2","donor_height"
                                             ,"donor_weight","donor_eth_code" ,"donor_eth_country"
                                             ,"donor_eth"  ,"donor_diabetes", "donor_hypertension"
                                             ,"donor_smoker", "donor_id","donor_death"
                                             ,"donor_dcd"  ,"donor_creatinine" ,"donor_age"
                                             ,"donor_sex"  ,"donor_kdri", "tx_date"
                                             ,"num_kidney" ))) stop("wrong column names")
    if (any(!colnames(recip_matrix)%in% c("recip_id" , "recip_graftno", "recip_pra"
                                           ,"recip_blgroup" ,"recip_a1","recip_a2"
                                           ,"recip_b1" , "recip_b2","recip_dr1"
                                           ,"recip_dr2", "recip_dq1","recip_dq2"
                                           ,"recip_state_initial","recip_state_current", "recip_age_rrtstart"
                                           , "recip_sex","recip_eth_detailed" ,"recip_eth"
                                           , "recip_primrenal", "recip_biopsy" , "recip_creatinine"
                                           ,"recip_height","recip_weight" ,"recip_smoker"
                                           ,"recip_rrtstartdate", "recip_deathdate", "recip_state"
                                           , "recip_age", "recip_waittime","recip_lung"
                                           ,"recip_coronary", "recip_pvd","recip_cvd"
                                           ,"recip_diabetes", "recip_comorb_date"  ,"recip_cancer"
                                           ,"recip_liststartdate", "recip_waitstatus", "recip_listtime"
                                           ,"recip_epts", "riskscore", "recip_algorithm"
                                           ,"recip_id_original", "recip_original_listdate" ,"recip_original_rrtstartdate"
                                           ,"recip_tx_date", "recip_birthdate" ,"pra_group"
                                           ,"recip_age_group"))) stop("wrong column names")

    if (nrow(donor_matrix) < num_donor) {
        warning("Number of donor in donor_matrix is less than num_donor,
                 reset the number")
        num_donor <- nrow(donor_matrix)
    }


    if (is_parallel) {

        #print("here")

        cl <-  makeCluster(ncores)
        registerDoSNOW(cl)

        # set up progress bar
        opts <- list(progress = function(n) {
             setTxtProgressBar(
               txtProgressBar(min = 1, max = resampleN, style = 3),
                n
            )
        })


        # this will run the simulation in parallel
        finalList <-  foreach(i=1:resampleN,
                                      .combine='comb', .multicombine=TRUE ,
                                      .verbose = T,
                                      # .export = ls(globalenv()),

                                      .export = c(
                                                  'allocation_national',

                                                  'australia_state_algorithm',
                                                  'australia_state_selection',

                                                  'blood_group_atch', 'comb', 'dm_national_choice', 'dm_national_formula', 'dm_national_top', 'dm_state_choice', 'dm_state_formula',
                                                  'dm_state_top',    'hla_match',
                                                  'linMap',
                                                  'NSWACT_allocation',
                                                  'QLD_allocation',
                                                  'raw_donor_matrix', 'raw_recip_matrix',
                                                  'raw_recip_matrix_subset',  'recip_sample_list', 'run_simulation',
                                                  'SA_allocation', 'select_max', 'selection_default', 'selection_corisk',
                                                   'selection_irisk',   'state_balance',
                                                  'VICTAS_allocation', 'WA_allocation') ,
                                      .options.snow = opts,
                                      .errorhandling = 'remove'
                                #       .packages=c("doSNOW", "doParallel")
        ) %dopar% {


            sim_donor_matrix = donor_matrix[sample(1:nrow(donor_matrix), num_donor), ]

            # sort donor by time
            sim_donor_matrix = sim_donor_matrix[order(sim_donor_matrix$tx_date), ]

            # for the dynamic waiting list, we start with the second donor to change the wl,
            # we get the added people from outside the current eligible people on the list
            # we kick off people from the current list
            # therefore, we have two check function below with the final result
            # update_recip_matrix <- rbind(update_recip_matrix, new_recip_matrix)

            if (dynamic_waitlist){
                #this is to get the update_recip_matrix (eligible recipients on the list for this donor)
                #in_wl_list:number to add
                #out_wl_list:number to kick

                tx_year=  format(as.Date(sim_donor_matrix$tx_date), format = "%Y")
                #print(c("transplants per year",table(tx_year)))

                # Sample an initial recipient matrix
                first_tx_date <- min(sim_donor_matrix$tx_date)
                final_tx_date <- max(sim_donor_matrix$tx_date)


                update_recip_matrix <- recip_matrix[recip_matrix$recip_liststartdate <= first_tx_date, ]

                idx <- sample(nrow(update_recip_matrix), min(num_recip, nrow(update_recip_matrix)))
                update_recip_matrix <- update_recip_matrix[idx, ]
                update_recip_matrix$recip_sim_listdate <- update_recip_matrix$recip_liststartdate
                #print(c("2006 waitlist",nrow(update_recip_matrix),  nrow(sim_donor_matrix)))

                #print(nrow(update_recip_matrix))
                #print(summary(update_recip_matrix$recip_waittime/365))

                # Decided the number of patient go in and out for the whole process
                # using poisson process - new arrivals into list & non-transplant departures from list
                tx_date_years <- as.numeric(final_tx_date - first_tx_date)/365
                #print(c("tx_date_years", tx_date_years))
                tx_date_ranges <- split(sort(sim_donor_matrix$tx_date),
                                        cut(seq_along(sort(sim_donor_matrix$tx_date)), tx_date_years, labels = FALSE))


                in_wl_list <- do.call(waitlist_FUN,
                                      args = append(list(donor_matrix = sim_donor_matrix,
                                                         IN = TRUE),
                                                    waitlist_arg))

                in_wl_list <- in_wl_list[[1]]
                out_wl_list <- do.call(waitlist_FUN,
                                       args = append(list(donor_matrix = sim_donor_matrix,
                                                          OUT = TRUE),
                                                     waitlist_arg))
                waitlist_Risk <- out_wl_list[[2]] #this is logic, to kick those risked or not
                #print(c("RISK", waitlist_Risk))
                out_wl_list <- out_wl_list[[1]]


                #print(c("in", sum(in_wl_list[1:365]), sum(in_wl_list)/tx_date_years, sum(in_wl_list)))
                #print(c("out", sum(out_wl_list[1:365]), sum(out_wl_list)/tx_date_years, sum(out_wl_list)))

            } else {

                update_recip_matrix = recip_matrix
            }
            #output here is update_recip_matrix
            for (d in 1:num_donor) {

                if (verbose && d %% 100 == 0) {
                    cat(d, "...")
                }

                if (dynamic_waitlist) {

                     # this is the real step to get the recip on the list
                    if (d > 1) {
                        #start with the second donor

                        dynamics_vector <- NULL


                        current_tx_date <-sim_donor_matrix$tx_date[d]
                        #print(c("current tx_date", class(current_tx_date),current_tx_date))

                        dynamics_vector<- c(as.character(sim_donor_matrix$tx_date[d]), sim_donor_matrix$donor_id[d], sim_donor_matrix$num_kidney[d])
                        #print(c("dynamics_vector",dynamics_vector ))

                        last_tx_date <- sim_donor_matrix$tx_date[d - 1]#the first transplant date

                        original_score <- allocation_national(update_recip_matrix,sim_donor_matrix[d,],
                                                              hla_match(update_recip_matrix,sim_donor_matrix[d,]))
                        #print(table(original_score < 54000000))#notice here, update_recip_matrix is the updated one from the previous if

                        in_wl_list_tmp <- in_wl_list[names(in_wl_list) <= current_tx_date &
                                                         names(in_wl_list) > last_tx_date]
                        # print(in_wl_list_tmp)#notice here, in_wl_list is obtained from previous if
                        num_in <- sum(in_wl_list_tmp)

                        num_out <- sum(out_wl_list[names(out_wl_list) <= current_tx_date &
                                                       names(out_wl_list) > last_tx_date])

                        include_idx <- which(recip_matrix$recip_liststartdate <= current_tx_date &
                                                 !recip_matrix$recip_id %in% update_recip_matrix$recip_id)
                        #those included cant be those current on the list, but need to be before the current donor

                        #print("Before ")
                        #print(nrow(update_recip_matrix))

                        dynamics_vector <- c(dynamics_vector, nrow(update_recip_matrix))
                        #print(c("dynamics_vector",dynamics_vector ))

                        #print(c(num_out, "out"))

                        if (length(include_idx) > 0 ) {
                            new_recip_matrix <- recip_matrix[include_idx, ]
                            #print(c(num_in, "in", nrow(new_recip_matrix)))



                            new_recip_matrix <- new_recip_matrix[sample(nrow(new_recip_matrix),
                                                                        min(num_in,nrow(new_recip_matrix))), ]#sample to get the added recipients

                            if (num_in > 0){
                                #get a new dialysis start date,subtract the length of days from the current donor's transplant date
                                new_recip_matrix$recip_sim_listdate <- rep(names(in_wl_list_tmp), in_wl_list_tmp)[1:nrow(new_recip_matrix)]
                                time_rrt_to_list <- as.numeric(new_recip_matrix$recip_liststartdate -  new_recip_matrix$recip_rrtstartdate)
                                new_recip_matrix$recip_rrtstartdate <- as.Date(as.Date(new_recip_matrix$recip_sim_listdate) - time_rrt_to_list) #have to overwrite for other functions

                            } #output the new_recip_matrix

                            #### Here we need to insert the kicking off mechanism - survival
                              # we kick off people from the current list: update_recip_matrix
                            if (waitlist_Risk){
                                #print("RISK RUNNING")
                                excl_recip_matrix <- update_recip_matrix[order(update_recip_matrix$riskscore, decreasing = TRUE),]
                                excl_recip_matrix <- excl_recip_matrix[1:min(num_out,nrow(excl_recip_matrix)),]

                                idx_exclude <- which(update_recip_matrix$recip_id %in% excl_recip_matrix$recip_id)

                            } else {
                                idx_exclude <- sample(nrow(update_recip_matrix), min(num_out,nrow(update_recip_matrix) ))

                            }

                            # cat("Number exclude")
                            # print(length(idx_exclude))
                            # cat("Number include")
                            # print(nrow(new_recip_matrix))

                            dynamics_vector <- c(dynamics_vector,length(idx_exclude),
                                                 nrow(new_recip_matrix))

                            if (length(idx_exclude) > 0) {
                                update_recip_matrix <- update_recip_matrix[-idx_exclude, ]
                            }  #output the update_recip_matrix


                            update_recip_matrix <- rbind(update_recip_matrix, new_recip_matrix)
                            #print(c("After",nrow(update_recip_matrix))) #combine then to get the current wl

                            if(nrow(update_recip_matrix) == 0){
                                print("ERROR: There are not enough people on the waitlist.
                                          Please specify larger number.")
                            }

                            dynamics_vector <- c(dynamics_vector,nrow(update_recip_matrix))
                            #print(c("dynamics_vector",dynamics_vector ))

                            dynamics_matrix = as.data.frame(rbind(dynamics_matrix,dynamics_vector))
                            colnames(dynamics_matrix) <- c("date", "donor_id", "num_kidney", "before", "exclude", "include", "after")
                        }

                    }
                }

                ##################################################
                ############## State balance  ####################
                ##################################################

                # Factors in one donor two kidney
                if(state_balance & is.null(state_algorithm_FUN)){

                    print("WARNING: Can only run state balance mechanism if state algorithm is specified")

                }

                if (state_balance & state_algorithm){

                    # print(c("excess",length(excess_vector)))

                    if(lubridate::month(sim_donor_matrix[d, ]$tx_date) == 12 & ! is.null(SimResults)){

                        # triggers the counting mechanism
                        # Counting is run once per year: in december of year donor is donating
                        if(length(which(
                            lubridate::month(SimResults$tx_date) == 12 &
                            lubridate::year(SimResults$tx_date) %in%
                            lubridate::year(sim_donor_matrix[d, ]$tx_date))) == 0 ){

                            excess_vector = c(excess_vector, state_balance(SimResults,
                                                                          recip_matrix = update_recip_matrix,
                                                                          donor_matrix = sim_donor_matrix[d, ]))

                            # Note: does not balance state imbalance caused by december donors
                            # May need update_excess_vector if state balance not achieved before next dec

                        }
                    }# goal is to count and output excess_vector

                    SB_state = NULL

                    # this only runs if donor_state is state that received too many kidney
                    # + excess state exist
                    if (length(excess_vector) != 0 &
                        !(as.character(sim_donor_matrix[d,]$donor_state) %in% excess_vector)){

                        num_kidney <- sim_donor_matrix$num_kidney[d]
                        num_send <- min(length(excess_vector), num_kidney)
                        SB_state = excess_vector[1:num_send] #new state for kidney
                        newstate = paste(c(SB_state), collapse = " / ")

                        for (SBK in SB_state){
                            excess_vector = excess_vector[-match(SBK,excess_vector)]
                        }

                    }# goal is to output SB_state

                    # if is.na(SB_state) is TRUE, continue usual simulation
                    if(!is.null(SB_state) & sum(!is.na(SB_state)) > 0){

                        # CASE1: If 1 kidney from donor
                        # CASE2/3: If 2 kidney from donor, same/ different SB_state

                        select_statebalance_idx = matrix( ncol=5, nrow = 0)
                        for (SBK in SB_state){

                            index <- which(as.character(update_recip_matrix$recip_state) %in% SBK
                                           & as.Date(update_recip_matrix$recip_liststartdate) < as.Date(sim_donor_matrix[d,]$tx_date)
                                           & update_recip_matrix$recip_blgroup %in% sim_donor_matrix[d,]$donor_blgroup)

                            eligible_recip_matrix_statebalance = update_recip_matrix[index,]


                            if (nrow(eligible_recip_matrix_statebalance) > 0) {

                                HLA_matrix_statebalance <- hla_match(eligible_recip_matrix_statebalance, sim_donor_matrix[d, ])
                                rownames(HLA_matrix_statebalance) <- eligible_recip_matrix_statebalance$recip_id

                                statebalance_recip_score <- state_algorithm_FUN(recip_matrix = eligible_recip_matrix_statebalance,
                                                                                donor_matrix = sim_donor_matrix[d, ],
                                                                                HLA = HLA_matrix_statebalance,
                                                                                state = SBK)

                                names(statebalance_recip_score) <- eligible_recip_matrix_statebalance$recip_id

                                if(length(SB_state) == 2 & SB_state[1] %in% SB_state[2]){
                                    num_need = 2
                                }else{num_need =1}

                                select_sb_idx <- do.call(state_matching_FUN,
                                                         args = append(list(recip_matrix = eligible_recip_matrix_statebalance,
                                                                            donor_matrix = sim_donor_matrix[d, ],
                                                                            allocation_score = statebalance_recip_score,
                                                                            HLA = HLA_matrix_statebalance,
                                                                            graft_number = num_need),
                                                                       state_matching_arg))

                                if(is.null(select_sb_idx)) next

                                statebalance_idx = as.matrix(cbind(select_sb_idx,
                                                                   score = statebalance_recip_score[select_sb_idx],
                                                                   HLA_matrix_statebalance[select_sb_idx,]), ncol=5)

                                select_statebalance_idx = as.matrix(rbind(select_statebalance_idx,
                                                                          statebalance_idx), ncol=5)

                                if (dim(select_statebalance_idx)[1] == length(SB_state)) break

                            }
                        }

                        # if NULL then next donor, donor is skipped (still in the system)
                        if(is.null(select_sb_idx)) next
                        select_statebalance_idx = unique(select_statebalance_idx)

                        statebalance_recip_id = paste(c(select_statebalance_idx[,1]), collapse = " / ")

                        newstate_vector = cbind(sim_donor_matrix[d, ],
                                                num_send = num_send,
                                                newstate = newstate,
                                                recip_id_sb = statebalance_recip_id)

                        statebalance_kidney=rbind(statebalance_kidney,newstate_vector)

                        allo_alg <- "State"


                        if (length(select_statebalance_idx[,1]) >= 1) {

                            sb_recip = select_statebalance_idx[,1]

                            simVector <- cbind(update_recip_matrix[update_recip_matrix$recip_id %in% sb_recip, ],
                                               sim_donor_matrix[d, ],
                                               score = select_statebalance_idx[select_sb_idx %in% sb_recip,2],
                                               tx_misa = select_statebalance_idx[select_sb_idx %in% sb_recip,3],
                                               tx_misb = select_statebalance_idx[select_sb_idx %in% sb_recip,4],
                                               tx_misdr = select_statebalance_idx[select_sb_idx %in% sb_recip,5],
                                               algorithm = allo_alg,
                                               recip_waittime_new=as.Date(sim_donor_matrix[d,"tx_date"])-
                                                   as.Date(update_recip_matrix[update_recip_matrix$recip_id %in%sb_recip,"recip_rrtstartdate"]))

                            SimResults  <- rbind(SimResults, simVector)
                            update_recip_matrix <- update_recip_matrix[!update_recip_matrix$recip_id %in%
                                                                           simVector$recip_id,]

                        }

                        #print(excess_vector)
                        if(length(excess_vector) != 0) next

                    } # statebalance function

                } # state_balance = TRUE


                # Select eligible recipients
                eligible_recip_matrix <- do.call(eligible_FUN,
                                                 args = append(list(recip_matrix = update_recip_matrix,
                                                                    donor_matrix = sim_donor_matrix,
                                                                    whichdonor = d),
                                                               eligible_arg))

                num_kidney <- sim_donor_matrix$num_kidney[d]

                # cat("Number of eligible recipient")
                # print(nrow(eligible_recip_matrix))

                if (nrow(eligible_recip_matrix) > 0) {

                    # Calculate HLA
                    HLA_matrix <- hla_match(eligible_recip_matrix,
                                                 sim_donor_matrix[d, ])

                    rownames(HLA_matrix) <- eligible_recip_matrix$recip_id
                    # Calculate Socre
                    recip_score <- algorithm_FUN(recip_matrix = eligible_recip_matrix,
                                                 donor_matrix = sim_donor_matrix[d, ],
                                                 HLA = HLA_matrix)

                    names(recip_score) <- eligible_recip_matrix$recip_id

                    index_national <- na.omit(recip_score) >= national_algorithm_threshold
                    national_scores <- recip_score[index_national]
                    # print(c("length national score",length(national_scores)))
                    # print(c("length recip score",length(recip_score)))

                    # to check how many kidneys are going to the state algorithm

                    max_recip_score <- sort(recip_score,
                                            decreasing = TRUE)[seq_len(num_kidney)]

                    eligible_recip_matrix_national <-
                        eligible_recip_matrix[which(eligible_recip_matrix$recip_id %in% names(national_scores)),]

                    print(dim(eligible_recip_matrix_national))

                    if (state_algorithm &
                        !is.null(state_algorithm_FUN)) {
                        num_kidney_national <- sum(na.omit(max_recip_score) >=
                                                       national_algorithm_threshold)
                        num_kidney_state <- num_kidney - num_kidney_national
                    } else {
                        num_kidney_national <- num_kidney
                    }



                    select_state_idx <- NULL
                    select_national_idx <- NULL
                    if (num_kidney_national != 0) {
                        # Select matched recipient

                        select_national_idx <- do.call(matching_FUN,
                                                       args = append(list(recip_matrix = eligible_recip_matrix_national,
                                                                          donor_matrix = sim_donor_matrix[d, ],
                                                                          allocation_score = national_scores,
                                                                          HLA = HLA_matrix,
                                                                          graft_number = num_kidney_national),
                                                                     matching_arg))
                        # print(select_national_idx)
                        allo_alg <- "National"

                        if (!is.null(select_national_idx) & length(select_national_idx) >= 1) {

                            simVector <- cbind(eligible_recip_matrix[eligible_recip_matrix$recip_id %in%
                                                                         select_national_idx, ],
                                               sim_donor_matrix[d, ],
                                               score = recip_score[select_national_idx],
                                               HLA_matrix[select_national_idx, ],
                                               algorithm = allo_alg,
                                               recip_waittime_new=sim_donor_matrix[d,"tx_date"]-eligible_recip_matrix[eligible_recip_matrix$recip_id %in%
                                                                                                                          select_national_idx,"recip_rrtstartdate"])

                            SimResults  <- rbind(SimResults, simVector)
                            update_recip_matrix <- update_recip_matrix[!update_recip_matrix$recip_id %in%
                                                                           simVector$recip_id,]

                        }

                        # update the number of the state kidney to be allocated
                        # length(select_national_idx) is the number of kidney donated in the national
                        num_kidney_state <- num_kidney - length(select_national_idx)

                    }


                    if (num_kidney_state != 0 &
                        state_algorithm &
                        !is.null(state_algorithm_FUN)) {

                        # Run state algorithm
                        donor_state <- as.character(sim_donor_matrix$donor_state)[d]
                        #
                        if (verbose) {
                            cat("run state algorithm:", donor_state, "\n")
                            cat("Blood group:",
                                sim_donor_matrix$donor_blgroup[d], "\n")
                        }

                        # Reselect the eligible recipient
                        eligible_recip_matrix_state <- do.call(state_eligible_FUN,
                                                               args = append(list(recip_matrix = update_recip_matrix,
                                                                                  donor_matrix = sim_donor_matrix,
                                                                                  whichdonor = d),
                                                                             state_eligible_arg))


                        # print(nrow(eligible_recip_matrix))
                        if (nrow(eligible_recip_matrix_state) > 0) {
                            HLA_matrix_state <- hla_match(eligible_recip_matrix_state,
                                                               sim_donor_matrix[d, ])
                            rownames(HLA_matrix_state) <- eligible_recip_matrix_state$recip_id

                            state_recip_score <- state_algorithm_FUN(recip_matrix = eligible_recip_matrix_state,
                                                                     donor_matrix = sim_donor_matrix[d, ],
                                                                     HLA = HLA_matrix_state)
                            names(state_recip_score) <- eligible_recip_matrix_state$recip_id

                            # Select matched recipient

                            select_state_idx <- do.call(state_matching_FUN,
                                                        args = append(list(recip_matrix = eligible_recip_matrix_state,
                                                                           donor_matrix = sim_donor_matrix[d, ],
                                                                           allocation_score = state_recip_score,
                                                                           HLA = HLA_matrix_state,
                                                                           graft_number = num_kidney_state),
                                                                      state_matching_arg))

                            #print(select_state_idx)
                            allo_alg <- "State"
                            if (!is.null(select_state_idx) & length(select_state_idx) >= 1) {

                                simVector <- cbind(eligible_recip_matrix_state[eligible_recip_matrix_state$recip_id %in%
                                                                                   select_state_idx, ],
                                                   sim_donor_matrix[d, ],
                                                   score = state_recip_score[select_state_idx],
                                                   HLA_matrix_state[select_state_idx, ],
                                                   algorithm = allo_alg,
                                                   recip_waittime_new=sim_donor_matrix[d,"tx_date"]-
                                                       eligible_recip_matrix_state[eligible_recip_matrix_state$recip_id %in%select_state_idx,"recip_rrtstartdate"])
                                SimResults  <- rbind(SimResults, simVector)
                                update_recip_matrix <- update_recip_matrix[!update_recip_matrix$recip_id %in%
                                                                               simVector$recip_id,]
                            }

                        }

                    }

                    select_idx <- c(select_national_idx,
                                    select_state_idx)

                    # print(c("nat/state", select_idx))

                } else {
                    # If there is no eligible pool in national, disregard the donor
                    select_idx <- NULL
                }


                select_national_idx2 <- NULL

                if(National_II){


                    if (length(select_idx) < num_kidney) {
                        #print("rerun national algorithm...")

                        num_kidney_remain <- num_kidney - length(select_idx)

                        # reselect the pool with bloodGroup_strict_match = FALSE
                        eligible_recip_matrix <- do.call(eligible_FUN,
                                                         args = append(list(recip_matrix = update_recip_matrix,
                                                                            donor_matrix = sim_donor_matrix,
                                                                            whichdonor = d,
                                                                            bloodGroup_strict_match = FALSE),
                                                                       eligible_arg))


                        if (nrow(eligible_recip_matrix) > 0) {

                            # Calculate HLA
                            HLA_matrix <- hla_match(eligible_recip_matrix,
                                                         sim_donor_matrix[d, ])

                            rownames(HLA_matrix) <- eligible_recip_matrix$recip_id
                            # Calculate Socre
                            recip_score <- algorithm_FUN(recip_matrix = eligible_recip_matrix,
                                                         donor_matrix = sim_donor_matrix[d, ],
                                                         HLA = HLA_matrix)

                            names(recip_score) <- eligible_recip_matrix$recip_id

                            max_recip_score <- sort(recip_score,
                                                    decreasing = TRUE)[seq_len(num_kidney_remain)]

                            index_national <- na.omit(recip_score) >= national_algorithm_threshold
                            national_scores <- recip_score[index_national]

                            eligible_recip_matrix_national <-
                                eligible_recip_matrix[which(eligible_recip_matrix$recip_id %in% names(national_scores)),]

                            select_national_idx2 <- do.call(matching_FUN,
                                                            args = append(list(recip_matrix = eligible_recip_matrix_national,
                                                                               donor_matrix = sim_donor_matrix[d, ],
                                                                               allocation_score = national_scores,
                                                                               HLA = HLA_matrix,
                                                                               graft_number = num_kidney_remain),
                                                                          matching_arg))
                            #print(c("nat2", select_national_idx2))

                            allo_alg <- "National_round2"

                            if (!is.null(select_national_idx2) & length(select_national_idx2) >= 1) {

                                simVector <- cbind(eligible_recip_matrix[eligible_recip_matrix$recip_id %in%
                                                                             select_national_idx2, ],
                                                   sim_donor_matrix[d, ],
                                                   score = recip_score[select_national_idx2],
                                                   HLA_matrix[select_national_idx2, ],
                                                   algorithm = allo_alg,
                                                   recip_waittime_new=sim_donor_matrix[d,"tx_date"]-
                                                       eligible_recip_matrix[eligible_recip_matrix$recip_id %in%select_national_idx2,"recip_rrtstartdate"])

                                SimResults  <- rbind(SimResults, simVector)
                                update_recip_matrix <- update_recip_matrix[!update_recip_matrix$recip_id %in%
                                                                               simVector$recip_id,]

                            }

                        }
                    }

                }


                if (is.null(select_national_idx2) & is.null(select_idx) |
                    num_kidney > length(select_national_idx2) + length(select_idx)) {
                    #print("kidney got discarded :( ")
                    discard_vector=sim_donor_matrix[d, ]
                    discard_kidney=rbind(discard_kidney,discard_vector)
                }


                # cat("Number of kidney assigned:",
                #     length(select_idx),
                #     "/",
                #     num_kidney,
                #     "\n"
                # )


                current_nrow <- nrow(SimResults)

                # if(DEBUG) print(simVector)
            }## end d


            list(SimResults = SimResults,
                 discard_kidney = discard_kidney,
                 statebalance_kidney = statebalance_kidney,
                 dynamics_matrix = dynamics_matrix)



        }


        SimResults = finalList$SimResults
        discard_kidney = finalList$discard_kidney
        statebalance_kidney = finalList$statebalance_kidney
        dynamics_matrix = finalList$dynamics_matrix


    } else {#i think this is duplication without the run parallel

        for (r in 1:resampleN) {
            # set.seed(r)

            if (verbose) {
                cat("\n")
                print(paste("Resample = ", r))
            }

            sim_donor_matrix = donor_matrix[sample(1:nrow(donor_matrix), num_donor), ]

            # sort donor by time
            sim_donor_matrix = sim_donor_matrix[order(sim_donor_matrix$tx_date), ]

            if (dynamic_waitlist){

                #######################################################
                ############## Dynamic waiting list ###################
                #######################################################

                tx_year=  format(as.Date(sim_donor_matrix$tx_date), format = "%Y")
                #print(c("transplants per year",table(tx_year)))

                # Sample an initial recipient matrix
                first_tx_date <- min(sim_donor_matrix$tx_date)
                final_tx_date <- max(sim_donor_matrix$tx_date)


                update_recip_matrix <- recip_matrix[recip_matrix$recip_liststartdate <= first_tx_date, ]

                idx <- sample(nrow(update_recip_matrix), min(num_recip, nrow(update_recip_matrix)))
                update_recip_matrix <- update_recip_matrix[idx, ]
                update_recip_matrix$recip_sim_listdate <- update_recip_matrix$recip_liststartdate
                #print(c("2006 waitlist",nrow(update_recip_matrix),  nrow(sim_donor_matrix)))

                #print(nrow(update_recip_matrix))
                #print(summary(update_recip_matrix$recip_waittime/365))

                # Decided the number of patient go in and out for the whole process
                # using poisson process - new arrivals into list & non-transplant departures from list
                tx_date_years <- as.numeric(final_tx_date - first_tx_date)/365
                #print(c("tx_date_years", tx_date_years))
                tx_date_ranges <- split(sort(sim_donor_matrix$tx_date),
                                        cut(seq_along(sort(sim_donor_matrix$tx_date)), tx_date_years, labels = FALSE))


                in_wl_list <- do.call(waitlist_FUN,
                                      args = append(list(donor_matrix = sim_donor_matrix,
                                                         IN = TRUE),
                                                    waitlist_arg))

                in_wl_list <- in_wl_list[[1]]
                out_wl_list <- do.call(waitlist_FUN,
                                       args = append(list(donor_matrix = sim_donor_matrix,
                                                          OUT = TRUE),
                                                     waitlist_arg))
                waitlist_Risk <- out_wl_list[[2]]
                #print(c("RISK", waitlist_Risk))
                out_wl_list <- out_wl_list[[1]]


                # print(c("in", sum(in_wl_list[1:365]), sum(in_wl_list)/tx_date_years, sum(in_wl_list)))
                # print(c("out", sum(out_wl_list[1:365]), sum(out_wl_list)/tx_date_years, sum(out_wl_list)))

            } else {

                update_recip_matrix = recip_matrix
            }

            for (d in 1:num_donor) {

                if (verbose && d %% 100 == 0) {
                    cat(d, "...")
                }

                if (dynamic_waitlist) {

                    ########################################################
                    ############## Dynamic waiting list ####################
                    ########################################################

                    if (d > 1) {

                        dynamics_vector <- NULL


                        current_tx_date <-sim_donor_matrix$tx_date[d]
                        #print(c("current tx_date", class(current_tx_date),current_tx_date))

                        dynamics_vector<- c(as.character(sim_donor_matrix$tx_date[d]), sim_donor_matrix$donor_id[d], sim_donor_matrix$num_kidney[d])
                        #print(c("dynamics_vector",dynamics_vector ))

                        last_tx_date <- sim_donor_matrix$tx_date[d - 1]

                        original_score <- allocation_national(update_recip_matrix,sim_donor_matrix[d,],
                                                              hla_match(update_recip_matrix,sim_donor_matrix[d,]))
                        #print(table(original_score < 54000000))

                        in_wl_list_tmp <- in_wl_list[names(in_wl_list) <= current_tx_date &
                                                         names(in_wl_list) > last_tx_date]
                        # print(in_wl_list_tmp)
                        num_in <- sum(in_wl_list_tmp)

                        num_out <- sum(out_wl_list[names(out_wl_list) <= current_tx_date &
                                                       names(out_wl_list) > last_tx_date])

                        include_idx <- which(recip_matrix$recip_liststartdate <= current_tx_date &
                                                 !recip_matrix$recip_id %in% update_recip_matrix$recip_id)

                        #print("Before ")
                        # print(nrow(update_recip_matrix))

                        dynamics_vector <- c(dynamics_vector, nrow(update_recip_matrix))
                        #print(c("dynamics_vector",dynamics_vector ))

                        #print(c(num_out, "out"))

                        if (length(include_idx) > 0 ) {
                            new_recip_matrix <- recip_matrix[include_idx, ]
                            #print(c(num_in, "in", nrow(new_recip_matrix)))



                            new_recip_matrix <- new_recip_matrix[sample(nrow(new_recip_matrix),
                                                                        min(num_in,nrow(new_recip_matrix))), ]

                            if (num_in > 0){
                                new_recip_matrix$recip_sim_listdate <- rep(names(in_wl_list_tmp), in_wl_list_tmp)[1:nrow(new_recip_matrix)]
                                time_rrt_to_list <- as.numeric(new_recip_matrix$recip_liststartdate -  new_recip_matrix$recip_rrtstartdate)
                                new_recip_matrix$recip_rrtstartdate <- as.Date(as.Date(new_recip_matrix$recip_sim_listdate) - time_rrt_to_list) #have to overwrite for other functions

                            }

                            #### Here we need to insert the kicking off mechanism - survival

                            if (waitlist_Risk){
                                #print("RISK RUNNING")
                                excl_recip_matrix <- update_recip_matrix[order(update_recip_matrix$riskscore, decreasing = TRUE),]
                                excl_recip_matrix <- excl_recip_matrix[1:min(num_out,nrow(excl_recip_matrix)),]

                                idx_exclude <- which(update_recip_matrix$recip_id %in% excl_recip_matrix$recip_id)

                            } else {
                                idx_exclude <- sample(nrow(update_recip_matrix), min(num_out,nrow(update_recip_matrix) ))

                            }

                            # cat("Number exclude")
                            # print(length(idx_exclude))
                            # cat("Number include")
                            # print(nrow(new_recip_matrix))

                            dynamics_vector <- c(dynamics_vector,length(idx_exclude),
                                                 nrow(new_recip_matrix))

                            if (length(idx_exclude) > 0) {
                                update_recip_matrix <- update_recip_matrix[-idx_exclude, ]
                            }


                            update_recip_matrix <- rbind(update_recip_matrix, new_recip_matrix)
                            #print(c("After",nrow(update_recip_matrix)))

                            if(nrow(update_recip_matrix) == 0){
                                print("ERROR: There are not enough people on the waitlist.
                                          Please specify larger number.")
                            }

                            dynamics_vector <- c(dynamics_vector,nrow(update_recip_matrix))
                            #print(c("dynamics_vector",dynamics_vector ))

                            dynamics_matrix = as.data.frame(rbind(dynamics_matrix,dynamics_vector))
                            colnames(dynamics_matrix) <- c("date", "donor_id", "num_kidney", "before", "exclude", "include", "after")
                        }

                    }
                }

                ##################################################
                ############## State balance  ####################
                ##################################################

                # Factors in one donor two kidney
                if(state_balance & is.null(state_algorithm_FUN)){

                    print("WARNING: Can only run state balance mechanism if state algorithm is specified")

                }

                if (state_balance & state_algorithm){

                    # print(c("excess",length(excess_vector)))

                    if(lubridate::month(sim_donor_matrix[d, ]$tx_date) == 12){

                        # triggers the counting mechanism
                        # Counting is run once per year: in december of year donor is donating
                        if(length(which(
                            lubridate::month(SimResults$tx_date) == 12 &
                            lubridate::year(SimResults$tx_date) %in%
                            lubridate::year(sim_donor_matrix[d, ]$tx_date))) == 0 ){

                            excess_vector = c(excess_vector, state_balance(SimResults,
                                                                          recip_matrix = update_recip_matrix,
                                                                          donor_matrix = sim_donor_matrix[d, ]))

                            # Note: does not balance state imbalance caused by december donors
                            # May need update_excess_vector if state balance not achieved before next dec

                        }
                    }# goal is to count and output excess_vector

                    SB_state = NULL

                    # this only runs if donor_state is state that received too many kidney
                    # + excess state exist
                    if (length(excess_vector) != 0 &
                        !(as.character(sim_donor_matrix[d,]$donor_state) %in% excess_vector)){

                        num_kidney <- sim_donor_matrix$num_kidney[d]
                        num_send <- min(length(excess_vector), num_kidney)
                        SB_state = excess_vector[1:num_send] #new state for kidney
                        newstate = paste(c(SB_state), collapse = " / ")

                        for (SBK in SB_state){
                            excess_vector = excess_vector[-match(SBK,excess_vector)]
                        }

                    }# goal is to output SB_state

                    # if is.na(SB_state) is TRUE, continue usual simulation
                    if(!is.null(SB_state) & sum(!is.na(SB_state)) > 0){

                        # CASE1: If 1 kidney from donor
                        # CASE2/3: If 2 kidney from donor, same/ different SB_state

                        select_statebalance_idx = matrix(, ncol=5, nrow=0)
                        for (SBK in SB_state){

                            index <- which(as.character(update_recip_matrix$recip_state) %in% SBK
                                           & as.Date(update_recip_matrix$recip_liststartdate) < as.Date(sim_donor_matrix[d,]$tx_date)
                                           & update_recip_matrix$recip_blgroup %in% sim_donor_matrix[d,]$donor_blgroup)

                            eligible_recip_matrix_statebalance = update_recip_matrix[index,]


                            if (nrow(eligible_recip_matrix_statebalance) > 0) {

                                HLA_matrix_statebalance <- hla_match(eligible_recip_matrix_statebalance, sim_donor_matrix[d, ])
                                rownames(HLA_matrix_statebalance) <- eligible_recip_matrix_statebalance$recip_id

                                statebalance_recip_score <- state_algorithm_FUN(recip_matrix = eligible_recip_matrix_statebalance,
                                                                                donor_matrix = sim_donor_matrix[d, ],
                                                                                HLA = HLA_matrix_statebalance,
                                                                                state = SBK)

                                names(statebalance_recip_score) <- eligible_recip_matrix_statebalance$recip_id

                                if(length(SB_state) == 2 & SB_state[1] %in% SB_state[2]){
                                    num_need = 2
                                }else{num_need =1}

                                select_sb_idx <- do.call(state_matching_FUN,
                                                         args = append(list(recip_matrix = eligible_recip_matrix_statebalance,
                                                                            donor_matrix = sim_donor_matrix[d, ],
                                                                            allocation_score = statebalance_recip_score,
                                                                            HLA = HLA_matrix_statebalance,
                                                                            graft_number = num_need),
                                                                       state_matching_arg))

                                if(is.null(select_sb_idx)) next

                                statebalance_idx = as.matrix(cbind(select_sb_idx,
                                                                   score = statebalance_recip_score[select_sb_idx],
                                                                   HLA_matrix_statebalance[select_sb_idx,]), ncol=5)

                                select_statebalance_idx = as.matrix(rbind(select_statebalance_idx,
                                                                          statebalance_idx), ncol=5)

                                if (dim(select_statebalance_idx)[1] == length(SB_state)) break

                            }
                        }

                        # if NULL then next donor, donor is skipped (still in the system)
                        if(is.null(select_sb_idx)) next
                        select_statebalance_idx = unique(select_statebalance_idx)

                        statebalance_recip_id = paste(c(select_statebalance_idx[,1]), collapse = " / ")

                        newstate_vector = cbind(sim_donor_matrix[d, ],
                                                num_send = num_send,
                                                newstate = newstate,
                                                recip_id_sb = statebalance_recip_id)

                        statebalance_kidney=rbind(statebalance_kidney,newstate_vector)

                        allo_alg <- "State"


                        if (length(select_statebalance_idx[,1]) >= 1) {

                            sb_recip = select_statebalance_idx[,1]

                            simVector <- cbind(update_recip_matrix[update_recip_matrix$recip_id %in% sb_recip, ],
                                               sim_donor_matrix[d, ],
                                               score = select_statebalance_idx[select_sb_idx %in% sb_recip,2],
                                               tx_misa = select_statebalance_idx[select_sb_idx %in% sb_recip,3],
                                               tx_misb = select_statebalance_idx[select_sb_idx %in% sb_recip,4],
                                               tx_misdr = select_statebalance_idx[select_sb_idx %in% sb_recip,5],
                                               algorithm = allo_alg,
                                               recip_waittime_new=as.Date(sim_donor_matrix[d,"tx_date"])-
                                                   as.Date(update_recip_matrix[update_recip_matrix$recip_id %in%sb_recip,"recip_rrtstartdate"]))

                            SimResults  <- rbind(SimResults, simVector)
                            update_recip_matrix <- update_recip_matrix[!update_recip_matrix$recip_id %in%
                                                                           simVector$recip_id,]

                        }

                        #print(excess_vector)
                        if(length(excess_vector) != 0) next

                    } # statebalance function

                } # state_balance = TRUE


                # Select eligible recipients
                eligible_recip_matrix <- do.call(eligible_FUN,
                                                 args = append(list(recip_matrix = update_recip_matrix,
                                                                    donor_matrix = sim_donor_matrix,
                                                                    whichdonor = d),
                                                               eligible_arg))

                num_kidney <- sim_donor_matrix$num_kidney[d]

                # cat("Number of eligible recipient")
                # print(nrow(eligible_recip_matrix))

                if (nrow(eligible_recip_matrix) > 0) {

                    # Calculate HLA
                    HLA_matrix <- hla_match(eligible_recip_matrix,
                                                 sim_donor_matrix[d, ])

                    rownames(HLA_matrix) <- eligible_recip_matrix$recip_id
                    # Calculate Socre
                    recip_score <- algorithm_FUN(recip_matrix = eligible_recip_matrix,
                                                 donor_matrix = sim_donor_matrix[d, ],
                                                 HLA = HLA_matrix)

                    names(recip_score) <- eligible_recip_matrix$recip_id

                    index_national <- na.omit(recip_score) >= national_algorithm_threshold
                    national_scores <- recip_score[index_national]
                    #print(c("length national score",length(national_scores)))
                    #print(c("length recip score",length(recip_score)))

                    # to check how many kidneys are going to the state algorithm

                    max_recip_score <- sort(recip_score,
                                            decreasing = TRUE)[seq_len(num_kidney)]

                    eligible_recip_matrix_national <-
                        eligible_recip_matrix[which(eligible_recip_matrix$recip_id %in% names(national_scores)),]

                    #print(dim(eligible_recip_matrix_national))

                    if (state_algorithm &
                        !is.null(state_algorithm_FUN)) {
                        num_kidney_national <- sum(na.omit(max_recip_score) >=
                                                       national_algorithm_threshold)
                        num_kidney_state <- num_kidney - num_kidney_national
                    } else {
                        num_kidney_national <- num_kidney
                    }



                    select_state_idx <- NULL
                    select_national_idx <- NULL
                    if (num_kidney_national != 0) {
                        # Select matched recipient

                        select_national_idx <- do.call(matching_FUN,
                                                       args = append(list(recip_matrix = eligible_recip_matrix_national,
                                                                          donor_matrix = sim_donor_matrix[d, ],
                                                                          allocation_score = national_scores,
                                                                          HLA = HLA_matrix,
                                                                          graft_number = num_kidney_national),
                                                                     matching_arg))
                        # print(select_national_idx)
                        allo_alg <- "National"

                        if (!is.null(select_national_idx) & length(select_national_idx) >= 1) {

                            simVector <- cbind(eligible_recip_matrix[eligible_recip_matrix$recip_id %in%
                                                                         select_national_idx, ],
                                               sim_donor_matrix[d, ],
                                               score = recip_score[select_national_idx],
                                               HLA_matrix[select_national_idx, ],
                                               algorithm = allo_alg,
                                               recip_waittime_new=sim_donor_matrix[d,"tx_date"]-eligible_recip_matrix[eligible_recip_matrix$recip_id %in%
                                                                                                                          select_national_idx,"recip_rrtstartdate"])

                            SimResults  <- rbind(SimResults, simVector)
                            update_recip_matrix <- update_recip_matrix[!update_recip_matrix$recip_id %in%
                                                                           simVector$recip_id,]

                        }

                        # update the number of the state kidney to be allocated
                        # length(select_national_idx) is the number of kidney donated in the national
                        num_kidney_state <- num_kidney - length(select_national_idx)

                    }


                    if (num_kidney_state != 0 &
                        state_algorithm &
                        !is.null(state_algorithm_FUN)) {

                        # Run state algorithm
                        donor_state <- as.character(sim_donor_matrix$donor_state)[d]
                        #
                        if (verbose) {
                            cat("run state algorithm:", donor_state, "\n")
                            cat("Blood group:",
                                sim_donor_matrix$donor_blgroup[d], "\n")
                        }

                        # Reselect the eligible recipient
                        eligible_recip_matrix_state <- do.call(state_eligible_FUN,
                                                               args = append(list(recip_matrix = update_recip_matrix,
                                                                                  donor_matrix = sim_donor_matrix,
                                                                                  whichdonor = d),
                                                                             state_eligible_arg))


                        # print(nrow(eligible_recip_matrix))
                        if (nrow(eligible_recip_matrix_state) > 0) {
                            HLA_matrix_state <- hla_match(eligible_recip_matrix_state,
                                                               sim_donor_matrix[d, ])
                            rownames(HLA_matrix_state) <- eligible_recip_matrix_state$recip_id

                            state_recip_score <- state_algorithm_FUN(recip_matrix = eligible_recip_matrix_state,
                                                                     donor_matrix = sim_donor_matrix[d, ],
                                                                     HLA = HLA_matrix_state)
                            names(state_recip_score) <- eligible_recip_matrix_state$recip_id

                            # Select matched recipient

                            select_state_idx <- do.call(state_matching_FUN,
                                                        args = append(list(recip_matrix = eligible_recip_matrix_state,
                                                                           donor_matrix = sim_donor_matrix[d, ],
                                                                           allocation_score = state_recip_score,
                                                                           HLA = HLA_matrix_state,
                                                                           graft_number = num_kidney_state),
                                                                      state_matching_arg))

                            #print(select_state_idx)
                            allo_alg <- "State"
                            if (!is.null(select_state_idx) & length(select_state_idx) >= 1) {

                                simVector <- cbind(eligible_recip_matrix_state[eligible_recip_matrix_state$recip_id %in%
                                                                                   select_state_idx, ],
                                                   sim_donor_matrix[d, ],
                                                   score = state_recip_score[select_state_idx],
                                                   HLA_matrix_state[select_state_idx, ],
                                                   algorithm = allo_alg,
                                                   recip_waittime_new=sim_donor_matrix[d,"tx_date"]-
                                                       eligible_recip_matrix_state[eligible_recip_matrix_state$recip_id %in%select_state_idx,"recip_rrtstartdate"])
                                SimResults  <- rbind(SimResults, simVector)
                                update_recip_matrix <- update_recip_matrix[!update_recip_matrix$recip_id %in%
                                                                               simVector$recip_id,]
                            }

                        }

                    }

                    select_idx <- c(select_national_idx,
                                    select_state_idx)

                    # print(c("nat/state", select_idx))

                } else {
                    # If there is no eligible pool in national, disregard the donor
                    select_idx <- NULL
                }


                select_national_idx2 <- NULL

                if(National_II){


                    if (length(select_idx) < num_kidney) {
                        #print("rerun national algorithm...")

                        num_kidney_remain <- num_kidney - length(select_idx)

                        # reselect the pool with bloodGroup_strict_match = FALSE
                        eligible_recip_matrix <- do.call(eligible_FUN,
                                                         args = append(list(recip_matrix = update_recip_matrix,
                                                                            donor_matrix = sim_donor_matrix,
                                                                            whichdonor = d,
                                                                            bloodGroup_strict_match = FALSE),
                                                                       eligible_arg))


                        if (nrow(eligible_recip_matrix) > 0) {

                            # Calculate HLA
                            HLA_matrix <- hla_match(eligible_recip_matrix,
                                                         sim_donor_matrix[d, ])

                            rownames(HLA_matrix) <- eligible_recip_matrix$recip_id
                            # Calculate Socre
                            recip_score <- algorithm_FUN(recip_matrix = eligible_recip_matrix,
                                                         donor_matrix = sim_donor_matrix[d, ],
                                                         HLA = HLA_matrix)

                            names(recip_score) <- eligible_recip_matrix$recip_id

                            max_recip_score <- sort(recip_score,
                                                    decreasing = TRUE)[seq_len(num_kidney_remain)]

                            index_national <- na.omit(recip_score) >= national_algorithm_threshold
                            national_scores <- recip_score[index_national]

                            eligible_recip_matrix_national <-
                                eligible_recip_matrix[which(eligible_recip_matrix$recip_id %in% names(national_scores)),]

                            select_national_idx2 <- do.call(matching_FUN,
                                                            args = append(list(recip_matrix = eligible_recip_matrix_national,
                                                                               donor_matrix = sim_donor_matrix[d, ],
                                                                               allocation_score = national_scores,
                                                                               HLA = HLA_matrix,
                                                                               graft_number = num_kidney_remain),
                                                                          matching_arg))
                            #print(c("nat2", select_national_idx2))

                            allo_alg <- "National_round2"

                            if (!is.null(select_national_idx2) & length(select_national_idx2) >= 1) {

                                simVector <- cbind(eligible_recip_matrix[eligible_recip_matrix$recip_id %in%
                                                                             select_national_idx2, ],
                                                   sim_donor_matrix[d, ],
                                                   score = recip_score[select_national_idx2],
                                                   HLA_matrix[select_national_idx2, ],
                                                   algorithm = allo_alg,
                                                   recip_waittime_new=sim_donor_matrix[d,"tx_date"]-
                                                       eligible_recip_matrix[eligible_recip_matrix$recip_id %in%select_national_idx2,"recip_rrtstartdate"])

                                SimResults  <- rbind(SimResults, simVector)
                                update_recip_matrix <- update_recip_matrix[!update_recip_matrix$recip_id %in%
                                                                               simVector$recip_id,]

                            }

                        }
                    }

                }


                if (is.null(select_national_idx2) & is.null(select_idx) |
                    num_kidney > length(select_national_idx2) + length(select_idx)) {
                    #print("kidney got discarded :( ")
                    discard_vector=sim_donor_matrix[d, ]
                    discard_kidney=rbind(discard_kidney,discard_vector)
                }


                # cat("Number of kidney assigned:",
                #     length(select_idx),
                #     "/",
                #     num_kidney,
                #     "\n"
                # )


                current_nrow <- nrow(SimResults)

                # if(DEBUG) print(simVector)
            }## end d
        }## end r
    }

    #####################################

    if(state_balance & state_algorithm & dynamic_waitlist){
        return(list(SimResults,discard_kidney, statebalance_kidney, dynamics_matrix))}

    if(state_balance & state_algorithm & !(dynamic_waitlist)){
        return(list(SimResults,discard_kidney, statebalance_kidney))}

    if(!(state_balance & state_algorithm) & dynamic_waitlist){
        return(list(SimResults,discard_kidney, dynamics_matrix))}

    ######################################

    return(list(SimResults,discard_kidney))
}
