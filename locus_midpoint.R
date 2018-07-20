## Function to determine rounded midpoints of GRanges loci
locMidpoint <- function (x, ...)
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- which(strand(x) == "+" | strand(x) == "*")
    on_plus_start <- start(x)[on_plus]
    on_plus_end <- end(x)[on_plus]
    start(x)[on_plus] <- round( on_plus_start + ((on_plus_end - on_plus_start)/2) )
    end(x)[on_plus] <- round( on_plus_start + ((on_plus_end - on_plus_start)/2) )

    on_minus <- which(strand(x) == "-")
    on_minus_start <- end(x)[on_minus]
    on_minus_end <- start(x)[on_minus]
    start(x)[on_minus] <- round( on_minus_end + ((on_minus_start - on_minus_end)/2) ) 
    end(x)[on_minus] <- round( on_minus_end + ((on_minus_start - on_minus_end)/2) )

    x
}

## Function to determine rounded midpoints of GRanges loci
## and to include 2 flanking bases (midpoint-1 bp, midpoint, midpoint+1bp)
locMidpointFlank2bp <- function (x, ...)
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- which(strand(x) == "+" | strand(x) == "*")
    on_plus_start <- start(x)[on_plus]
    on_plus_end <- end(x)[on_plus]
    start(x)[on_plus] <- round( on_plus_start + ((on_plus_end - on_plus_start)/2) ) - 1
    end(x)[on_plus] <- round( on_plus_start + ((on_plus_end - on_plus_start)/2) ) + 1

    on_minus <- which(strand(x) == "-")
    on_minus_start <- end(x)[on_minus]
    on_minus_end <- start(x)[on_minus]
    start(x)[on_minus] <- round( on_minus_end + ((on_minus_start - on_minus_end)/2) ) - 1
    end(x)[on_minus] <- round( on_minus_end + ((on_minus_start - on_minus_end)/2) ) + 1

    x
}

## and to include 500 flanking bases (midpoint-1 bp, midpoint, midpoint+1bp)
locMidpointFlank <- function (x, leftFlank, rightFlank, ...)
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- which(strand(x) == "+" | strand(x) == "*")
    on_plus_start <- start(x)[on_plus]
    on_plus_end <- end(x)[on_plus]
    start(x)[on_plus] <- round( on_plus_start + ((on_plus_end - on_plus_start)/2) ) - leftFlank
    end(x)[on_plus] <- round( on_plus_start + ((on_plus_end - on_plus_start)/2) ) + rightFlank

    on_minus <- which(strand(x) == "-")
    on_minus_start <- end(x)[on_minus]
    on_minus_end <- start(x)[on_minus]
    start(x)[on_minus] <- round( on_minus_end + ((on_minus_start - on_minus_end)/2) ) - rightFlank
    end(x)[on_minus] <- round( on_minus_end + ((on_minus_start - on_minus_end)/2) ) + leftFlank

    x
}

