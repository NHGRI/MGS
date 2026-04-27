#!/usr/bin/env Rscript

SCRIPT_NAME <- "my_script.R"

# disable strings as factors, but re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))

# for parallel processing
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(iterators))

# for profiling
# suppressPackageStartupMessages(library(pryr))
# suppressPackageStartupMessages(library(profvis))


#' Combines dataframes in a list
#'
#' @importFrom data.table rbindlist
rbindlist_df <- function(...) {
    return(data.table::rbindlist(list(...), use.names = TRUE, fill = TRUE))
}


#' Command line interface wrapper
#'
#' @importFrom optparse make_option
#' @importFrom optparse OptionParser
#' @importFrom optparse parse_args
command_line_interface <- function() {
    optionList <- list(
        optparse::make_option(c("-f", "--file"),
            type = "character",
            help = paste0(
                "Input file.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--flag"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = paste0(
                "My flag.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--threads"),
            type = "numeric",
            default = 1,
            help = paste0(
                "Number of threads to use.",
                " [default %default]"
            )
        ),

        optparse::make_option(c("--verbose"),
            type = "logical",
            action = "store_true",
            default = FALSE,
            help = paste0(
                "Verbose mode (write extra info to std.err).",
                " [default %default]"
            )
        )
    )

    parser <- optparse::OptionParser(
        usage = "%prog",
        option_list = optionList,
        description = paste0(
            "Runs an analysis."
        )
    )

    # a hack to fix a bug in optparse that won"t let you use positional args
    # if you also have non-boolean optional args:
    getOptionStrings <- function(parserObj) {
        optionStrings <- character()
        for (item in parserObj@options) {
            optionStrings <- append(optionStrings,
                                    c(item@short_flag, item@long_flag))
        }
        optionStrings
    }
    optStrings <- getOptionStrings(parser)
    arguments <- optparse::parse_args(parser, positional_arguments = TRUE)

    # read in the parameters
    param <- list()
    for (i in names(arguments$options)) {
        param[[i]] <- arguments$options[[i]]
    }

    param[["df"]] <- read.csv(
        param[["file"]],
        header = F,
        stringsAsFactors = F
    )

    return(run_analysis(param))
}


#' Run analysis function
#'
#' @importFrom data.table fread
#' @importFrom parallel detectCores
#' @importFrom parallel makeForkCluster
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom doMC registerDoMC
#' @importFrom iterators isplit
#' @importFrom parallel stopCluster
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom plyr rename
run_analysis <- function(param) {

    # set the parallel backend: fork, snow, domc
    param[["parallel_type"]] <- "domc"

    # print out parameters
    if (param[["verbose"]]) {
        message(paste0("Parameters [", SCRIPT_NAME, "]:"))
        for (i in names(param)) {
            message("\t", i, " = ", param[[i]])
        }
    }

    # for parallel processing with foreach
    # more info here: https://ljdursi.github.io/beyond-single-core-R
    time_init <- Sys.time()
    if (param[["threads"]] > parallel::detectCores()) {
        warning(paste0(
            "User requested more threads than available. Setting number of",
            " threads to doParallel::detectCores() - 1, that is",
            " ", parallel::detectCores() - 1, "."
        ))
        param[["threads"]] <- parallel::detectCores() - 1
    }
    if (param[["parallel_type"]] == "fork") {
        sink("/dev/null") # use sink to grab annoying output when outfile = ""
        # (v1) Forking = copy the R session in its current state. Fast
        # because only copies objects if they are modified; however it takes
        # much more memory as files are copied
        cl <- parallel::makeForkCluster(param[["threads"]], outfile = "")
        sink()
        doParallel::registerDoParallel(cl)
    } else if (param[["parallel_type"]] == "snow") {
        # since snow could be on any machine, set outfile to /dev/null rather
        # than the parent process (unlink fork)
        # (v2) Spawn new subprocess: benefits can be done on any machine,
        # new process will not have unneeded data
        if (param[["verbose"]]) {
            message(paste0(
                "[", SCRIPT_NAME, "]:\t",
                " Running parallel with snow. Any std.out or std.err (apart",
                " from a failed exit) will be lost."
            ))
        }
        cl <- parallel::makeCluster(param[["threads"]], outfile = "/dev/null")
        parallel::clusterExport(cl, c(
            "SCRIPT_NAME", "param"
        ))
        doSNOW::registerDoSNOW(cl)
    } else if (param[["parallel_type"]] == "domc") {
        # (v3) Shared memory process
        # doesn't use persistent workers like the snow package or ths
        # now-derived functions in the parallel package
        doMC::registerDoMC(param[["threads"]])
    } else {
        stop(paste0("[", SCRIPT_NAME, "]:\tERROR invalid parallel_type."))
    }

    df <- param[["df"]]
    df$iter_var <- paste(
        df$V1,
        df$V2,
        sep = "::"
    )

    # iterate over each index variant
    #
    # note that if using foreach and combining the results into a dataframe,
    # errors will not stop execution, but will not necessarily be returned to
    # the user. Will show up as:
    # Error in { : task 1 failed - "replacement has 1 row, data has 0"
    #
    df_list_results <- foreach::foreach(
        iter_info = iterators::iter(df, by = "row"), # iterator row
        #iter_info = iterators::isplit(df, df$iter_var), # iterator split
        .combine = rbindlist_dflist,
        .inorder = FALSE,
        .multicombine = TRUE,
        .errorhandling = "stop"
    ) %dopar% {
        #message(paste0("(iter_e)", iter_info$my_column, "\n"))

        # return list
        df_list_tmp <- list(some_result1, some_result2)

        return(df_list_tmp)
    } # end outer foreach loop

    # stop the cluster
    if (param[["parallel_type"]] != "domc") {
        parallel::stopCluster(cl)
    }
    time_end <- Sys.time()
    if (param[["verbose"]]) {
        message(paste0(
            "\nForeach loop execution time", " [", SCRIPT_NAME, "]:\t",
            difftime(time_end, time_init, units = c("hours")),
            " hours."
        ))
    }

    # add param to return list of results
    df_list_results[["param"]] <- param

    return(df_list_results)
}


main <- function() {
    # run analysis
    run_time <- system.time(df_results <- command_line_interface())
    message(paste0(
        "Analysis execution time", " [", SCRIPT_NAME, "]:\t",
        run_time[["elapsed"]]/3600, # proc.time sec to hours
        " hours."
    ))
    # write results to std.out
    # write.table(
    #     df_results[["my_result_df"]],
    #     "",
    #     row.names = FALSE,
    #     col.names = TRUE,
    #     quote = FALSE,
    #     sep = "\t",
    #     na = ""
    # )
    # save the my_result_df results
    out_file_base <- df_results[["param"]][["output_file_tag"]]
    gzfh <- gzfile(
        paste(out_file_base, ".tsv.gz", sep = ""),
        "w",
        compression = 9
    )
    write.table(
        df_results[["my_result_df"]],
        gzfh,
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t",
        na = ""
    )
    close(gzfh)

    return(0)
}


dev <- function() {
    # for profiling
    suppressPackageStartupMessages(library(profvis))

    base <- paste0("profile-", gsub(" ", "_", Sys.time()))
    # run the analysis to understand memory usage
    prof <- profvis::profvis(command_line_interface())
    saveRDS(prof, paste0(base, ".Rds.gz"), compress = TRUE)
    #prof <- readRDS(paste0(base, ".Rds.gz"))
    #print(prof)

    # only run if snow
    prof <- snow::snow.time(run_analysis())
    pdf(file = paste0(base, ".pdf"), height = 5, width = 6)
    print(plot(prof))
    dev.off()
}


test <- function() {
    # a test function
}


main()
#dev()
#test()
