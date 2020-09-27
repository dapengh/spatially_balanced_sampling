packages <- c("spsurvey", "BalancedSampling", "parallel", "foreach", "doParallel", "doRNG", "SDraw")
lapply(packages, library, character.only = TRUE)

numLevels <- function(samplesize, shift.grid, startlev, maxlev, sfobject) {
  # As necessary, export sfobject to the parallel processing cluster
  
  if(!all(sf::st_geometry_type(sfobject) %in% c("POINT", "MULTIPOINT"))) {
    cl <- getDefaultCluster()
    clusterExport(cl, "sfobject", envir=environment())
  }
  
  # Determine the number of levels for hierarchical randomization
  
  if(is.null(startlev)) {
    nlev <- ceiling(logb(samplesize, 4))
    if(nlev == 0) {
      nlev <- 1
    }
  } else {
    nlev <- startlev
  }
  cel.wt <- 99999
  celmax.ind <- 0
  sint <- 1
  if(shift.grid) {
    roff.x <- runif(1, 0, 1)
    roff.y <- runif(1, 0, 1)
  }
  bbox <- sf::st_bbox(sfobject)
  grid.extent <- max(bbox["xmax"] - bbox["xmin"], bbox["ymax"] - bbox["ymin"])
  temp <- 0.04*grid.extent
  grid.xmin <- bbox["xmin"] - temp
  grid.ymin <- bbox["ymin"] - temp
  grid.extent <- 1.08*grid.extent
  grid.xmax <- grid.xmin + grid.extent
  grid.ymax <- grid.ymin + grid.extent
  while (any(cel.wt/sint > 1) && celmax.ind < 2 && nlev <= maxlev) {
    cat( "Current number of levels:", nlev, "\n");
    celmax <- max(cel.wt)
    nlv2 <- 2^nlev
    dx <- dy <- grid.extent/nlv2
    xc <- seq(grid.xmin, grid.xmax, length=nlv2+1)
    yc <- seq(grid.ymin, grid.ymax, length=nlv2+1)
    if(shift.grid) {
      xc <- rep(xc, nlv2+1) + (roff.x * dx)
      yc <- rep(yc, rep(nlv2+1, nlv2+1)) + (roff.y * dy)
    } else {
      xc <- rep(xc, nlv2+1)
      yc <- rep(yc, rep(nlv2+1, nlv2+1))
    }
    
    # Determine total inclusion probability for each grid cell and, as necessary,
    # adjust the indicator for whether maximum of the total inclusion probabilities
    # is changing
    
    cel.wt <- spsurvey::cellWeight(xc, yc, dx, dy, sfobject)
    if(max(cel.wt) == celmax) {
      celmax.ind <- celmax.ind + 1
      if(celmax.ind == 2) {
        warning("\nSince the maximum value of total inclusion probability for the grid cells was \nnot changing, the algorithm for determining the number of levels for \nhierarchical randomization was terminated.\n")
      }
    }
    
    # Adjust sampling interval and number of hierarchical levels
    
    sint <- sum(cel.wt)/samplesize
    if(nlev == maxlev) {
      nlev <- nlev + 1
    } else {
      nlev <- nlev + max(1, ceiling(logb(cel.wt[cel.wt > 0]/sint, 4)))
      if(nlev > maxlev) {
        nlev <- maxlev
      }
    }
  }
  
  #  Print the final number of levels
  
  cat( "Final number of levels:", nlev-1, "\n");
  
  # Create the output list
  
  rslt <- list(nlev=nlev, xc=xc, yc=yc, dx=dx, dy=dy, cel.wt=cel.wt, sint=sint)
  
  # Return the list
  
  return(rslt)
}

grtspts <- function(ptsframe, samplesize = 100, SiteBegin = 1,
                    shift.grid = TRUE, do.sample = TRUE, startlev = NULL, maxlev = 11) {
  
  # Determine the number of levels for hierarchical randomization
  
  temp <- numLevels(samplesize, shift.grid,
                    startlev, maxlev, ptsframe)
  nlev <- temp$nlev
  dx <- temp$dx
  dy <- temp$dy
  xc <- temp$xc
  yc <- temp$yc
  cel.wt <- temp$cel.wt
  sint <- temp$sint
  
  # Assign the final number of levels
  
  endlev <- nlev - 1
  
  # Remove cells with zero weight
  
  indx <- cel.wt > 0
  xc <- xc[indx]
  yc <- yc[indx]
  cel.wt <- cel.wt[indx]
  
  # Construct the hierarchical address for all cells
  
  hadr <- spsurvey::constructAddr(xc, yc, dx, dy, nlev)
  
  # Construct randomized hierarchical addresses
  
  ranhadr <- spsurvey::ranho(hadr)
  
  # Determine order of the randomized hierarchical addresses
  
  rord <- order(ranhadr)
  
  if(do.sample) {
    
    # Select grid cells that get a sample point
    
    rstrt <- runif(1, 0, sint)
    ttl.wt <- c(0, cumsum(cel.wt[rord]))
    idx <- ceiling((ttl.wt - rstrt)/sint)
    smpdx <- spsurvey::pickGridCells(samplesize, idx)
    rdx <- rord[smpdx]
    n.cells <- length(unique(rdx))
    if(length(rdx) > n.cells) {
      temp <- sum(sapply(split(rdx, rdx), length) > 1)
      warning(paste("\nOf the ", n.cells, " grid cells from which sample points were selected,\n", temp, " (", round(100*temp/n.cells, 1), "%) of the cells contained more than one sample point.\n", sep=""))
    }
    
    # Pick sample point(s) in selected cells
    
    id <- spsurvey::pickFiniteSamplePoints(rdx, xc, yc, dx, dy, ptsframe)
    rho <- ptsframe[match(id, ptsframe$id), ]
    
  } else {
    
    # Pick all points in the frame
    
    id <- selectframe(rord, xc, yc, dx, dy, ptsframe)
    rho <- ptsframe[match(id, ptsframe$id), ]
  }
  
  # Construct sample hierarchical address
  
  np <- nrow(rho)
  nlev <- max(1, trunc(logb(np,4)))
  ifelse(np == 4^nlev, nlev, nlev <- nlev + 1)
  ad <- matrix(0, 4^nlev, nlev)
  rv4 <- 0:3
  pwr4 <- 4.^(0.:(nlev - 1.))
  for(i in 1:nlev)
    ad[, i] <- rep(rep(rv4, rep(pwr4[i], 4.)),pwr4[nlev]/pwr4[i])
  rho4 <- as.vector(ad%*%matrix(rev(pwr4), nlev, 1))
  
  # Place sample in reverse hierarchical order
  
  rho <- rho[unique(floor(rho4 * np/4^nlev)) + 1.,]
  
  # Create desired attributes
  
  rho$siteID <- SiteBegin - 1 + 1:nrow(rho)
  temp <- sf::st_coordinates(rho)
  rho$xcoord <- temp[,1]
  rho$ycoord <- temp[,2]
  rho$wgt <- 1/rho$mdm
  
  # Create the output sf object
  
  rho <- subset(rho, select = c("siteID", "id", "xcoord", "ycoord", "mdcaty",
                                "wgt"))
  row.names(rho) <- 1:nrow(rho)
  
  # Assign the final number of levels as an attribute of the output object
  
  attr(rho, "nlev") <- endlev
  
  # Return the sample
  
  rho
}

grts <- function(design, DesignID = "Site", SiteBegin = 1, type.frame = NULL,
                 src.frame = "shapefile", in.shape = NULL, sf.object = NULL, sp.object = NULL,
                 att.frame = NULL, id = NULL, xcoord = NULL, ycoord = NULL, stratum = NULL,
                 mdcaty = NULL, startlev = NULL, maxlev = 11, maxtry = NULL, shift.grid =
                   TRUE, do.sample = rep(TRUE, length(design)), shapefile = TRUE, prjfilename =
                   NULL, out.shape = "sample.shp") {
  
  # Ensure that a design list is provided
  
  if(is.null(design))
    stop("\nA design list must be provided.")
  
  # Ensure that the design list is named and determine strata names from the
  # design list
  
  strata.names <- names(design)
  if(is.null(strata.names)) {
    if(length(design) > 1) {
      stop("\nThe design list must be named.")
    } else {
      warning("\nSince the single stratum specified in the design list was not named, \n\"None\" will be used for the stratum name.\n")
      strata.names <- "None"
      names(design) <- strata.names
    }
  }
  
  # Ensure that type.frame contains a valid value
  
  if(is.null(type.frame))
    stop("\nA value must be provided for argument type.frame.")
  temp <- match(type.frame, c("finite", "linear", "area"), nomatch=0)
  if(temp == 0)
    stop(paste("\nThe value provided for argument type.frame, \"", type.frame, "\", is not a valid value.", sep=""))
  
  # Ensure that src.frame contains a valid value
  
  temp <- match(src.frame, c("sf.object", "shapefile", "sp.object", "att.frame"),
                nomatch=0)
  if(temp == 0)
    stop(paste("\nThe value provided for argument src.frame, \"", src.frame, "\", is not a valid value.", sep=""))
  
  # If src.frame equals "shapefile", then create an sf object from the shapefile
  
  if(src.frame == "shapefile") {
    if(is.null(shapefile))
      stop("\nA shapefile name is required when the value provided for argument src.frame \nequals \"shapefile\".")
    nc <- nchar(in.shape)
    if(substr(in.shape, nc-3, nc) != ".shp") {
      if(substr(in.shape, nc-3, nc-3) == ".") {
        in.shape <- paste(substr(in.shape, 1, nc-4), ".shp", sep="")
      } else {
        in.shape <- paste(in.shape, ".shp", sep="")
      }
    }
    sf.object <- st_read(in.shape, quiet = TRUE)
  }
  
  # If src,frame equals "sf.object", ensure that an sf object was provided
  
  if(src.frame == "sf.object") {
    if(is.null(sf.object))
      stop("\nAn sf package object is required when the value provided for argument src.frame \nequals \"sf.object\".")
  }
  
  # If src.frame equals "sp.object", then create an sf object from the sp object
  
  if(src.frame == "sp.object") {
    if(is.null(sp.object))
      stop("\nAn sp package object is required when the value provided for argument src.frame \nequals \"sp.object\".")
    sf.object <- sf::st_as_sf(sp.object)
  }
  
  # If src.frame equals "att.frame", ensure that type.frame equals "finite",
  # ensure that a data frame object is assigned to argument att.frame, and create
  # an sf object from att.frame
  
  if(src.frame == "att.frame") {
    if(type.frame != "finite") {
      stop(paste("\nThe value provided for argument type.frame must equal \"finite\" when argument \nsrc.frame equals \"att.frame\"  The value provided for argument type.frame was \n\"", type.frame, "\".", sep=""))
    }
    if(is.null(att.frame)) {
      stop(paste("\nA data frame object must be assigned to argument att.frame when argument\nsrc.frame equals \"att.frame\"."))
    }
    if(is.null(xcoord) | is.null(ycoord)) {
      stop(paste("\nValues must be provided for arguments xcoord and ycoord when argument src.frame \nequals \"att.frame\"."))
    }
    if(!(all(c(xcoord, ycoord) %in% names(att.frame)))) {
      stop(paste("\nThe values provided for arguments xcoord and ycoord do not occur among the \nnames for att.frame."))
    }
    sf.object <- sf::st_as_sf(att.frame, coords = c(xcoord, ycoord))
  }
  
  # If src.frame does not equal "att.frame" and att.frame is not NULL, create an
  # sf object composed of att.frame and the geometry column from sf.object
  
  if(src.frame != "att.frame" & !is.null(att.frame)) {
    geom <- st_geometry(sf.object)
    sf.object <- sf::st_set_geometry(att.frame, geom)
  }
  
  # Ensure that the class attribute for sf.object contains only the values "sf"
  # and "data.frame"
  
  class(sf.object) <- c("sf", "data.frame")
  
  # Ensure that the geometry types for sf.object are consistent
  
  temp <- sf::st_geometry_type(sf.object)
  tst <- all(temp %in% c("POINT", "MULTIPOINT")) |
    all(temp %in% c("LINESTRING", "MULTILINESTRING")) |
    all(temp %in% c("POLYGON", "MULTIPOLYGON"))
  if(!tst) {
    stop(paste("\nThe geometry types for the survey frame object passed to function grts: \n\"", unique(sf::st_geometry_type(sf.object)), "\" are not consistent.", sep=""))
  }
  
  # Create ID values
  
  id <- "id"
  sf.object$id <- 1:nrow(sf.object)
  
  # If stratum equals NULL, ensure that the design list specifies a single stratum
  # and add an attribute named "stratum" to sf.object.  Otherwise, ensure that the
  # name provided for stratum identifies an attribute in sf.object.
  
  if(is.null(stratum)) {
    if(length(strata.names) > 1)
      stop("\nThe attribute in sf.object that identifies stratum membership was not provided \nand the design list specifies more than one stratum.")
    stratum <- "stratum"
    sf.object$stratum <- factor(rep(strata.names, nrow(sf.object)))
  } else {
    temp <- match(stratum, names(sf.object), nomatch=0)
    if(temp == 0)
      stop(paste("\nThe value provided for the attribute in sf.object that identifies stratum \nmembership for each feature, \"", stratum, "\", does not occur among the \nattributes in sf.object.", sep=""))
  }
  
  # Ensure that the stratum attribute in sf.object is a factor
  
  if(!is.factor(sf.object$stratum))
    sf.object[, stratum] <- as.factor(sf.object[, stratum, drop = TRUE])
  
  # Check whether strata names from the design list occur among the values for the
  # stratum attribute in sf.object
  
  temp <- match(strata.names, levels(sf.object[, stratum, drop = TRUE]),
                nomatch=0)
  if(any(temp == 0)) {
    temp.str <- vecprint(strata.names[temp == 0])
    stop(paste("\nThe following strata names in the design list do not occur among the strata \nnames in the stratum attribute in sf.object:\n", temp.str, sep=""))
  }
  
  # If seltype is not "Equal" for every stratum, then do the following: (1) ensure
  # that mdcaty is not NULL and (2) ensure that the name provided for mdcaty
  # identifies an attribute in sf.object
  
  seltype.ind <- FALSE
  for(s in strata.names) {
    if(design[[s]]$seltype != "Equal") {
      seltype.ind <- TRUE
    }
  }
  if(seltype.ind) {
    if(is.null(mdcaty))
      stop(paste("\nThe name of the attribute in sf.object that identifies the unequal probability \ncategory for each feature must be provided.", sep=""))
    temp <- match(mdcaty, names(sf.object), nomatch=0)
    if(temp == 0)
      stop(paste("\nThe value provided for the attribute in sf.object that identifies the unequal \nprobability category for each feature, \"", mdcaty, "\", does not occur among \nthe attributes in sf.object.", sep=""))
  }
  
  # Ensure that startlev and maxlev are valid and compatible values
  
  if(!is.null(startlev)) {
    if(startlev < 1)
      stop("\nThe value for startlev cannot be less than 1")
    if(startlev > 11)
      stop("\nThe value for startlev cannot be greater than 11")
    if(maxlev < 1)
      stop("\nThe value for maxlev cannot be less than 1")
    if(maxlev > 11)
      stop("\nThe value for maxlev cannot be greater than 11")
    if(startlev > maxlev)
      stop("\nThe value for startlev cannot be greater than the value for maxlev")
  } else {
    if(maxlev < 1)
      stop("\nThe value for maxlev cannot be less than 1")
    if(maxlev > 11)
      stop("\nThe value for maxlev cannot be greater than 11")
  }
  
  # As necessary, initialize parallel processing
  
  if(type.frame != "finite") {
    ncore <- detectCores()
    tst <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if(nzchar(tst)) {
      ncore <- ifelse(ncore == 1L, 1L, 2L)
    } else {
      temp <- as.integer(floor(0.1 * ncore))
      ncore <- ifelse(ncore == 1L, 1L, ifelse(temp == 0, ncore - 1L, ncore - temp))
    }
    cl <- makeCluster(ncore, methods=FALSE)
    invisible(clusterEvalQ(cl, library(sf)))
    setDefaultCluster(cl)
  }
  
  # Begin the section for a finite population (discrete points)
  
  if(type.frame == "finite") {
    
    first <- TRUE
    SiteBegin <- SiteBegin
    
    # Ensure that do.sample is the correct length and is named
    
    if(length(do.sample) > 1) {
      if(length(do.sample) != length(design))
        stop("\nArgument do.sample must be the same length as the design list.")
      if(is.null(names(do.sample))) {
        names(do.sample) <- strata.names
      } else {
        temp <- match(names(do.sample), strata.names, nomatch=0)
        if(any(temp) == 0)
          temp.str <- vecprint(names(do.sample)[temp == 0])
        stop(paste("\nThe following names in do.sample do not occur among the names in design:\n", temp.str, sep=""))
      }
    } else if(is.null(names(do.sample))) {
      names(do.sample) <- strata.names
    }
    
    # Begin the loop for strata
    
    for(s in strata.names) {
      
      cat(paste("\nStratum:", s, "\n"))
      
      # Create the sample frame
      
      temp <- sf.object[, stratum, drop = TRUE] == s
      grtspts.ind <- TRUE
      if(sum(temp) == 0) {
        warning(paste("\nThe stratum column in the attributes data frame contains no values that match \nthe stratum named \"", s, "\" in the design list.\n", sep=""))
        next
      } else if(sum(temp) == 1) {
        warning(paste("\nThe stratum column in the attributes data frame contains a single value that \nmatches the stratum named \"", s, "\" in the design list. \nThe sample for this stratum will be composed of a single point.\n", sep=""))
        grtspts.ind <- FALSE
      }
      
      sframe <- subset(sf.object, temp)
      if(design[[s]]$seltype == "Equal") {
        sframe$mdcaty <- "Equal"
      } else if(design[[s]]$seltype == "Unequal") {
        sframe$mdcaty <- factor(sframe[, mdcaty, drop = TRUE])
      } else if(design[[s]]$seltype == "Continuous") {
        sframe$mdcaty <- sframe[, mdcaty, drop = TRUE]
      } else {
        stop(paste("\nThe value provided for the type of random selection, \"", design[[s]]$seltype, "\", \nfor stratum \"", s, "\" is not valid.", sep=""))
      }
      
      # If seltype is not "Equal", ensure that mdcaty contains valid values
      
      if(design[[s]]$seltype == "Unequal") {
        if(any(is.na(sframe$mdcaty)))
          stop(paste("\nMissing values were detected among the unequal probability category values for \nstratum \"", s, "\".", sep=""))
      } else if(design[[s]]$seltype == "Continuous") {
        if(any(is.na(sframe$mdcaty)))
          stop(paste("\nMissing values were detected among the unequal probability category values for \nstratum \"", s, "\".", sep=""))
        if(!is.numeric(sframe$mdcaty))
          stop(paste("\nThe type of random selection for stratum \"", s, "\" is \"Continuous\", \nbut the unequal probability category values are not numeric.", sep=""))
        if(any(sframe$mdcaty < 0))
          stop(paste("\nNonpositive values were detected among the unequal probability category values \nfor stratum \"", s, "\".", sep=""))
      }
      
      # If seltype is "Unequal", ensure that caty.n is provided and that the names
      # in caty.n are included amont the levels of mdcaty
      
      if(design[[s]]$seltype == "Unequal") {
        if(is.null(design[[s]]$caty.n))
          stop(paste("The type of random selection was set to \"Unequal\", but caty.n was not \nprovided for stratum \"", s, "\".", sep=""))
        temp <- match(names(design[[s]]$caty.n),
                      levels(as.factor(sframe$mdcaty)), nomatch=0)
        if(any(temp == 0)) {
          temp.str <- vecprint(names(design[[s]]$caty.n)[temp == 0])
          stop(paste("\nThe following names in caty.n for stratum \"", s, "\" do not occur \namong the levels of the mdcaty variable in att.frame:\n", temp.str, sep=""))
        }
      }
      
      # Ensure that panel and caty.n contain valid values
      
      if(!is.numeric(design[[s]]$panel))
        stop(paste(" The design list must contain numeric values in the panel argument for \nstratum \"", s, "\".\n", sep=""))
      design[[s]]$panel <- round(design[[s]]$panel)
      design[[s]]$panel <- design[[s]]$panel[design[[s]]$panel > 0]
      if(length(design[[s]]$panel) == 0)
        stop(paste(" The design list does not not contain any valid values of the panel \nargument for stratum \"", s, "\".\n", sep=""))
      
      if(design[[s]]$seltype == "Unequal") {
        if(!is.numeric(design[[s]]$caty.n))
          stop(paste(" The design list must contain numeric values in the caty.n argument for \nstratum \"", s, "\".\n", sep=""))
        design[[s]]$caty.n <- round(design[[s]]$caty.n)
        design[[s]]$caty.n <- design[[s]]$caty.n[design[[s]]$caty.n > 0]
        if(length(design[[s]]$caty.n) == 0)
          stop(paste(" The design list does not not contain any valid values of the caty.n \nargument for stratum \"", s, "\".\n", sep=""))
      }
      
      # As necessary, remove rows from sframe that have values of mdcaty which are not
      # included among the names in caty.n
      
      if(design[[s]]$seltype == "Unequal") {
        temp <- sframe$mdcaty %in% names(design[[s]]$caty.n)
        if(any(!temp)) {
          sframe <- sframe[temp,]
        }
      }
      
      # Determine overall sample size for the stratum
      
      if(is.null(design[[s]]$over))
        design[[s]]$over <- 0
      if(design[[s]]$seltype != "Unequal") {
        samplesize <- sum(design[[s]]$panel)
        n.desired <- sum(samplesize, design[[s]]$over)
      } else {
        if(sum(design[[s]]$panel) != sum(design[[s]]$caty.n))
          stop("\nThe sum of panel sample sizes does not equal sum of caty.n sample sizes")
        samplesize <- sum(design[[s]]$caty.n)
        if(design[[s]]$over == 0) {
          n.desired <- design[[s]]$caty.n
        } else {
          over.n <- design[[s]]$over * design[[s]]$caty.n /
            sum(design[[s]]$caty.n)
          if(any(over.n != floor(over.n)))
            warning(paste("\nOversample size is not proportional to category sample sizes for stratum\n\"", s, "\".\n", sep=""))
          n.desired <- design[[s]]$caty.n + ceiling(over.n)
        }
      }
      
      # Calculate mdm - inclusion probabilities
      
      if(design[[s]]$seltype == "Equal")
        sframe$mdm <- spsurvey::mdmpts(sframe$mdcaty, c(Equal=n.desired))
      else if(design[[s]]$seltype == "Unequal")
        sframe$mdm <- spsurvey::mdmpts(sframe$mdcaty, n.desired)
      else
        sframe$mdm <- n.desired * sframe$mdcaty / sum(sframe$mdcaty)
      
      # Select the sample
      
      sf::st_agr(sframe) <- "constant"
      if(grtspts.ind) {
        stmp <- grtspts(sframe, sum(n.desired), SiteBegin,
                        shift.grid, do.sample[s], startlev, maxlev)
      } else {
        stmp <- sframe
        stmp$siteID <- SiteBegin
        stmp$wgt <- 1/sframe$mdm
        stmp <- subset(stmp, select = c("siteID", "id", "mdcaty", "wgt"))
        row.names(stmp) <- 1
        attr(stmp, "nlev") <- NA
      }
      
      # Determine whether the realized sample size is less than the desired size
      
      if(nrow(stmp) < sum(n.desired))
        warning(paste("\nThe size of the selected sample was less than the desired size for stratum\n\"", s, "\".\n", sep=""))
      
      # Add the stratum variable
      
      stmp$stratum <- as.factor(rep(s,nrow(stmp)))
      
      # Add panel and oversample structure
      
      stmp$panel <- as.character(rep("OverSamp",nrow(stmp)))
      n.panel <- length(design[[s]]$panel)
      if(nrow(stmp) < samplesize) {
        n.short <- samplesize - nrow(stmp)
        n.temp <- n.short / n.panel
        if(n.temp != floor(n.temp)) {
          n.temp <- c(ceiling(n.temp), rep(floor(n.temp), n.panel-1))
          i <- 1
          while(sum(n.temp) != n.short) {
            i <- i+1
            n.temp[i] <- n.temp[i] + 1
          }
        }
        np <- c(0, cumsum(design[[s]]$panel - n.temp))
      } else {
        np <- c(0, cumsum(design[[s]]$panel))
      }
      for(i in 1:n.panel)
        stmp$panel[(np[i]+1):np[i+1]] <- names(design[[s]]$panel[i])
      
      # If an oversample is present or the realized sample size is less than the
      # desired size, then adjust the weights
      
      if(design[[s]]$over > 0 || nrow(stmp) < samplesize) {
        if(design[[s]]$seltype != "Unequal") {
          if(nrow(stmp) < samplesize) {
            stmp$wgt <- n.desired * stmp$wgt / nrow(stmp)
          } else {
            stmp$wgt <- n.desired * stmp$wgt / samplesize
          }
        } else {
          if(nrow(stmp) < samplesize) {
            n.caty <- length(design[[s]]$caty.n)
            n.temp <- n.short / n.caty
            nc <- design[[s]]$caty.n - n.temp
          } else {
            nc <- design[[s]]$caty.n
          }
          for(i in names(n.desired)) {
            stmp$wgt[stmp$mdcaty == i] <- n.desired[i] *
              stmp$wgt[stmp$mdcaty == i] / nc[i]
          }
        }
      }
      
      # Add stratum sample to the output object
      
      if(first) {
        sites <- stmp
        levels(sites$stratum) <- strata.names
        first <- FALSE
      } else {
        sites <- rbind(sites, stmp)
      }
      SiteBegin <- SiteBegin + nrow(stmp)
      
      # End the loop for strata
      
    }
    
    # End the section for a finite population (discrete points)
    
  } else if(type.frame == "linear") {
    
    # Begin the section for a linear network
    
    first <- TRUE
    SiteBegin <- SiteBegin
    sf.object$length_mdm <- as.numeric(st_length(sf.object))
    
    # Begin the loop for strata
    
    for(s in strata.names) {
      
      cat(paste("\nStratum:", s, "\n"))
      
      # Create the sample frame
      
      temp <- sf.object[, stratum, drop = TRUE] == s
      if(sum(temp) == 0) {
        warning(paste("\nThe stratum column in the attributes data frame contains no values that match \nthe stratum named \"", s, "\" in the design list.\n", sep=""))
        next
      }
      
      sframe <- subset(sf.object, temp)
      if(design[[s]]$seltype == "Equal") {
        sframe$mdcaty <- "Equal"
      } else if(design[[s]]$seltype == "Unequal") {
        sframe$mdcaty <- factor(sframe[, mdcaty, drop = TRUE])
      } else if(design[[s]]$seltype == "Continuous") {
        sframe$mdcaty <- sframe[, mdcaty, drop = TRUE]
      } else {
        stop(paste("\nThe value provided for the type of random selection, \"", design[[s]]$seltype, "\", \nfor stratum \"", s, "\" is not valid.", sep=""))
      }
      
      # If seltype is not "Equal", ensure that mdcaty contains valid values
      
      if(design[[s]]$seltype == "Unequal") {
        if(any(is.na(sframe$mdcaty)))
          stop(paste("\nMissing values were detected among the unequal probability category values for \nstratum \"", s, "\".", sep=""))
      } else if(design[[s]]$seltype == "Continuous") {
        if(any(is.na(sframe$mdcaty)))
          stop(paste("\nMissing values were detected among the unequal probability category values for \nstratum \"", s, "\".", sep=""))
        if(!is.numeric(sframe$mdcaty))
          stop(paste("\nThe type of random selection for stratum \"", s, "\" is \"Continuous\", \nbut the unequal probability category values are not numeric.", sep=""))
        if(any(sframe$mdcaty < 0))
          stop(paste("\nNonpositive values were detected among the unequal probability category values \nfor stratum \"", s, "\".", sep=""))
      }
      
      # If seltype is "Unequal", ensure that caty.n is provided and that the names
      # in caty.n and the levels of mdcaty are equivalent
      
      if(design[[s]]$seltype == "Unequal") {
        if(is.null(design[[s]]$caty.n))
          stop(paste("The type of random selection was set to \"Unequal\", but caty.n was not \nprovided for stratum \"", s, "\".", sep=""))
        temp <- match(names(design[[s]]$caty.n),
                      levels(as.factor(sframe$mdcaty)), nomatch=0)
        if(any(temp == 0)) {
          temp.str <- vecprint(names(design[[s]]$caty.n)[temp == 0])
          stop(paste("\nThe following names in caty.n for stratum \"", s, "\" do not occur \namong the levels of the mdcaty variable in att.frame:\n", temp.str, sep=""))
        }
      }
      
      # Ensure that panel and caty.n contain valid values
      
      if(!is.numeric(design[[s]]$panel))
        stop(paste(" The design list must contain numeric values in the panel argument for \nstratum \"", s, "\".\n", sep=""))
      design[[s]]$panel <- round(design[[s]]$panel)
      design[[s]]$panel <- design[[s]]$panel[design[[s]]$panel > 0]
      if(length(design[[s]]$panel) == 0)
        stop(paste(" The design list does not not contain any valid values of the panel \nargument for stratum \"", s, "\".\n", sep=""))
      
      if(design[[s]]$seltype == "Unequal") {
        if(!is.numeric(design[[s]]$caty.n))
          stop(paste(" The design list must contain numeric values in the caty.n argument for \nstratum \"", s, "\".\n", sep=""))
        design[[s]]$caty.n <- round(design[[s]]$caty.n)
        design[[s]]$caty.n <- design[[s]]$caty.n[design[[s]]$caty.n > 0]
        if(length(design[[s]]$caty.n) == 0)
          stop(paste(" The design list does not not contain any valid values of the caty.n \nargument for stratum \"", s, "\".\n", sep=""))
      }
      
      # As necessary, remove rows from sframe that have values of mdcaty which are not
      # included among the names in caty.n
      
      if(design[[s]]$seltype == "Unequal") {
        temp <- sframe$mdcaty %in% names(design[[s]]$caty.n)
        if(any(!temp)) {
          sframe <- sframe[temp,]
        }
      }
      
      # Determine overall sample size for the stratum
      
      if(is.null(design[[s]]$over))
        design[[s]]$over <- 0
      if(design[[s]]$seltype != "Unequal") {
        samplesize <- sum(design[[s]]$panel)
        n.desired <- sum(samplesize, design[[s]]$over)
      } else {
        if(sum(design[[s]]$panel) != sum(design[[s]]$caty.n))
          stop("\nThe sum of panel sample sizes does not equal sum of caty.n sample sizes")
        samplesize <- sum(design[[s]]$caty.n)
        if(design[[s]]$over == 0) {
          n.desired <- design[[s]]$caty.n
        } else {
          over.n <- design[[s]]$over * design[[s]]$caty.n /
            sum(design[[s]]$caty.n)
          if(any(over.n != floor(over.n)))
            warning(paste("\nOversample size is not proportional to category sample sizes for stratum\n\"", s, "\".\n", sep=""))
          n.desired <- design[[s]]$caty.n + ceiling(over.n)
        }
      }
      
      # Calculate mdm - inclusion probabilities
      
      if(design[[s]]$seltype == "Equal")
        sframe$mdm <- mdmlin(sframe$len, sframe$mdcaty, c(Equal=n.desired))
      else if(design[[s]]$seltype == "Unequal")
        sframe$mdm <- mdmlin(sframe$len, sframe$mdcaty, n.desired)
      else
        sframe$mdm <- n.desired * sframe$mdcaty /
        sum(sframe$len * sframe$mdcaty)
      
      # Select the sample
      
      sf::st_agr(sframe) <- "constant"
      stmp <- grtslin(sframe, sum(n.desired), SiteBegin, shift.grid, startlev,
                      maxlev)
      
      # Add the stratum variable
      
      stmp$stratum <- as.factor(rep(s,nrow(stmp)))
      
      # Add panel and oversample structure
      
      stmp$panel <- rep("OverSamp",nrow(stmp))
      np <- c(0,cumsum(design[[s]]$panel))
      for(i in 1:length(design[[s]]$panel))
        stmp$panel[(np[i]+1):np[i+1]] <- names(design[[s]]$panel[i])
      
      # If an oversample is present, then adjust the weights
      
      if(design[[s]]$over > 0) {
        if(design[[s]]$seltype != "Unequal") {
          stmp$wgt <- n.desired * stmp$wgt / samplesize
        } else {
          nc <- design[[s]]$caty.n
          for(i in names(n.desired)) {
            stmp$wgt[stmp$mdcaty == i] <- n.desired[i] *
              stmp$wgt[stmp$mdcaty == i] / nc[i]
          }
        }
      }
      
      # Add stratum sample to the output data frame
      
      if(first) {
        sites <- stmp
        levels(sites$stratum) <- strata.names
        first <- FALSE
      } else {
        sites <- rbind(sites, stmp)
      }
      SiteBegin <- SiteBegin + nrow(stmp)
      
      # End the loop for strata
      
    }
    
    # End the section for a linear network
    
  } else if(type.frame == "area") {
    
    # Begin the section for a polygonal area
    
    first <- TRUE
    SiteBegin <- SiteBegin
    sf.object$area_mdm <- as.numeric(st_area(sf.object))
    
    # Begin the loop for strata
    
    for(s in strata.names) {
      
      cat(paste("\nStratum:", s, "\n"))
      
      # Create the sample frame
      
      temp <- sf.object[, stratum, drop = TRUE] == s
      if(sum(temp) == 0) {
        warning(paste("\nThe stratum column in the attributes data frame contains no values that match \nthe stratum named \"", s, "\" in the design list.\n", sep=""))
        next
      }
      
      sframe <- subset(sf.object, temp)
      if(design[[s]]$seltype == "Equal") {
        sframe$mdcaty <- "Equal"
      } else if(design[[s]]$seltype == "Unequal") {
        sframe$mdcaty <- factor(sframe[, mdcaty, drop = TRUE])
      } else if(design[[s]]$seltype == "Continuous") {
        sframe$mdcaty <- sframe[, mdcaty, drop = TRUE]
      } else {
        stop(paste("\nThe value provided for the type of random selection, \"", design[[s]]$seltype, "\", \nfor stratum \"", s, "\" is not valid.", sep=""))
      }
      
      # If seltype is not "Equal", ensure that mdcaty contains valid values
      
      if(design[[s]]$seltype == "Unequal") {
        if(any(is.na(sframe$mdcaty)))
          stop(paste("\nMissing values were detected among the unequal probability category values for \nstratum \"", s, "\".", sep=""))
      } else if(design[[s]]$seltype == "Continuous") {
        if(any(is.na(sframe$mdcaty)))
          stop(paste("\nMissing values were detected among the unequal probability category values for \nstratum \"", s, "\".", sep=""))
        if(!is.numeric(sframe$mdcaty))
          stop(paste("\nThe type of random selection for stratum \"", s, "\" is \"Continuous\", \nbut the unequal probability category values are not numeric.", sep=""))
        if(any(sframe$mdcaty < 0))
          stop(paste("\nNonpositive values were detected among the unequal probability category values \nfor stratum \"", s, "\".", sep=""))
      }
      
      # If seltype is "Unequal", ensure that caty.n is provided and that the names
      # in caty.n and the levels of mdcaty are equivalent
      
      if(design[[s]]$seltype == "Unequal") {
        if(is.null(design[[s]]$caty.n))
          stop(paste("The type of random selection was set to \"Unequal\", but caty.n was not \nprovided for stratum \"", s, "\".", sep=""))
        temp <- match(names(design[[s]]$caty.n),
                      levels(as.factor(sframe$mdcaty)), nomatch=0)
        if(any(temp == 0)) {
          temp.str <- vecprint(names(design[[s]]$caty.n)[temp == 0])
          stop(paste("\nThe following names in caty.n for stratum \"", s, "\" do not occur \namong the levels of the mdcaty variable in att.frame:\n", temp.str, sep=""))
        }
      }
      
      # Ensure that panel and caty.n contain valid values
      
      if(!is.numeric(design[[s]]$panel))
        stop(paste(" The design list must contain numeric values in the panel argument for \nstratum \"", s, "\".\n", sep=""))
      design[[s]]$panel <- round(design[[s]]$panel)
      design[[s]]$panel <- design[[s]]$panel[design[[s]]$panel > 0]
      if(length(design[[s]]$panel) == 0)
        stop(paste(" The design list does not not contain any valid values of the panel \nargument for stratum \"", s, "\".\n", sep=""))
      
      if(design[[s]]$seltype == "Unequal") {
        if(!is.numeric(design[[s]]$caty.n))
          stop(paste(" The design list must contain numeric values in the caty.n argument for \nstratum \"", s, "\".\n", sep=""))
        design[[s]]$caty.n <- round(design[[s]]$caty.n)
        design[[s]]$caty.n <- design[[s]]$caty.n[design[[s]]$caty.n > 0]
        if(length(design[[s]]$caty.n) == 0)
          stop(paste(" The design list does not not contain any valid values of the caty.n \nargument for stratum \"", s, "\".\n", sep=""))
      }
      
      # As necessary, remove rows from sframe that have values of mdcaty which are not
      # included among the names in caty.n
      
      if(design[[s]]$seltype == "Unequal") {
        temp <- sframe$mdcaty %in% names(design[[s]]$caty.n)
        if(any(!temp)) {
          sframe <- sframe[temp,]
        }
      }
      
      # Determine overall sample size for the stratum
      
      if(is.null(design[[s]]$over))
        design[[s]]$over <- 0
      if(design[[s]]$seltype != "Unequal") {
        samplesize <- sum(design[[s]]$panel)
        n.desired <- sum(samplesize, design[[s]]$over)
      } else {
        if(sum(design[[s]]$panel) != sum(design[[s]]$caty.n))
          stop("\nThe sum of panel sample sizes does not equal sum of caty.n sample sizes")
        samplesize <- sum(design[[s]]$caty.n)
        if(design[[s]]$over == 0) {
          n.desired <- design[[s]]$caty.n
        } else {
          over.n <- design[[s]]$over * design[[s]]$caty.n /
            sum(design[[s]]$caty.n)
          if(any(over.n != floor(over.n)))
            warning(paste("\nOversample size is not proportional to category sample sizes for stratum\n\"", s, "\".\n", sep=""))
          n.desired <- design[[s]]$caty.n + ceiling(over.n)
        }
      }
      
      # Calculate mdm - inclusion probabilities
      
      if(design[[s]]$seltype == "Equal")
        sframe$mdm <- mdmarea(sframe$area_mdm, sframe$mdcaty, c(Equal=n.desired))
      else if(design[[s]]$seltype == "Unequal")
        sframe$mdm <- mdmarea(sframe$area, sframe$mdcaty, n.desired)
      else
        sframe$mdm <- n.desired * sframe$mdcaty /
        sum(sframe$area * sframe$mdcaty)
      
      # Select the sample
      
      sf::st_agr(sframe) <- "constant"
      stmp <- grtsarea(sframe, sum(n.desired), SiteBegin, shift.grid,
                       startlev, maxlev, maxtry)
      
      # Determine whether the realized sample size is less than the desired size
      
      if(nrow(stmp) < sum(n.desired))
        warning(paste("\nThe size of the selected sample was less than the desired size for stratum \n\"", s, "\".\n", sep=""))
      
      # Add the stratum variable
      
      stmp$stratum <- as.factor(rep(s, nrow(stmp)))
      
      # Add panel and oversample structure
      
      stmp$panel <- as.character(rep("OverSamp",nrow(stmp)))
      n.panel <- length(design[[s]]$panel)
      if(nrow(stmp) < samplesize) {
        n.short <- samplesize - nrow(stmp)
        n.temp <- n.short / n.panel
        if(n.temp != floor(n.temp)) {
          n.temp <- c(ceiling(n.temp), rep(floor(n.temp), n.panel-1))
          i <- 1
          while(sum(n.temp) != n.short) {
            i <- i+1
            ntemp[i] <- n.temp[i] + 1
          }
        }
        np <- c(0, cumsum(design[[s]]$panel - n.temp))
      } else {
        np <- c(0, cumsum(design[[s]]$panel))
      }
      for(i in 1:n.panel)
        stmp$panel[(np[i]+1):np[i+1]] <- names(design[[s]]$panel[i])
      
      # If an oversample is present or the realized sample size is less than the
      # desired size, then adjust the weights
      
      if(design[[s]]$over > 0 || nrow(stmp) < samplesize) {
        if(design[[s]]$seltype != "Unequal") {
          if(nrow(stmp) < samplesize) {
            stmp$wgt <- n.desired * stmp$wgt / nrow(stmp)
          } else {
            stmp$wgt <- n.desired * stmp$wgt / samplesize
          }
        } else {
          if(nrow(stmp) < samplesize) {
            n.caty <- length(design[[s]]$caty.n)
            n.temp <- n.short / n.caty
            nc <- design[[s]]$caty.n - n.temp
          } else {
            nc <- design[[s]]$caty.n
          }
          for(i in names(n.desired)) {
            stmp$wgt[stmp$mdcaty == i] <- n.desired[i] *
              stmp$wgt[stmp$mdcaty == i] / nc[i]
          }
        }
      }
      
      # Add stratum sample to the output data frame
      
      if(first) {
        sites <- stmp
        levels(sites$stratum) <- strata.names
        first <- FALSE
      } else {
        sites <- rbind(sites, stmp)
      }
      SiteBegin <- SiteBegin + nrow(stmp)
      
      # End the loop for strata
      
    }
    
    # End the section for a polygonal area
    
  }
  
  # As necessary, terminate parallel processing
  
  if(type.frame != "finite") {
    stopCluster(cl)
  }
  
  # Add DesignID name to the numeric siteID value to create a new siteID
  
  sites$siteID <- as.character(gsub(" ","0", paste(DesignID,"-",
                                                   format(sites$siteID), sep="")))
  
  # Add Evaluation Status and Evaluation Reason variables to the output data frame
  
  sites$EvalStatus <- rep("NotEval", nrow(sites))
  sites$EvalReason <- rep(" ", nrow(sites))
  
  # Add attributes from sf.object that are not included in sites
  
  tm <- match(sites$id, sf.object$id)
  geom_name <- attr(sf.object, "sf_column")
  if(design[[s]]$seltype == "Equal") {
    td <- match(c(id, stratum, "length_mdm", "area_mdm", geom_name),
                names(sf.object), nomatch=0)
  } else {
    td <- match(c(id, stratum, mdcaty, "length_mdm", "area_mdm", geom_name),
                names(sf.object), nomatch=0)
  }
  temp <- names(sf.object)[-td]
  if(length(temp) > 0) {
    for(i in temp) {
      sites[, i] <- sf.object[tm, i, drop = TRUE]
    }
  }
  
  # Remove the id attribute from sites
  
  temp <- names(sites)
  temp <- temp[!(temp %in% c("id", geom_name))]
  sites <- subset(sites, select=temp)
  
  # Add row names to sites
  
  n <- nrow(sites)
  IDs <- as.character(1:n)
  row.names(sites) <- IDs
  
  # Assign attributes to sites
  
  ifelse(is.null(startlev),
         attr(sites, "startlev") <- "Not specified",
         attr(sites, "startlev") <- startlev)
  ifelse(is.null(maxlev),
         attr(sites, "maxlev") <- "Not specified",
         attr(sites, "maxlev") <- maxlev)
  attr(sites, "endlev") <- attributes(stmp)$nlev
  attr(sites, "shift.grid") <- shift.grid
  attr(sites, "do.sample") <- do.sample
  
  # If requested, create a shapefile containing the sample information
  
  if(shapefile == TRUE) {
    nc <- nchar(out.shape)
    if(substr(out.shape, nc-3, nc) != ".shp") {
      if(substr(out.shape, nc-3, nc-3) == ".") {
        out.shape <- paste(substr(out.shape, 1, nc-4), ".shp", sep="")
      } else {
        out.shape <- paste(out.shape, ".shp", sep="")
      }
    }
    if(out.shape %in% list.files()) {
      warning(paste("\nThe output shapefile named \"", out.shape, "\" already exists and was \noverwritten.\n", sep=""))
      st_write(sites, out.shape, quiet = TRUE, delete_dsn = TRUE)
    } else {
      st_write(sites, out.shape, quiet = TRUE)
    }
  }
  
  # Create an object of class SpatialDesign
  SpointsMat <- sf::st_coordinates(sites)
  rownames(SpointsMat) <- IDs
  dat <- sf::st_set_geometry(sites, NULL)
  sp_obj <- sp::SpatialPointsDataFrame(sp::SpatialPoints(SpointsMat),data=dat)
  rslt <- spsurvey::SpatialDesign(design = design, sp_obj = sp_obj)
  
  # Return the SpatialDesign object
  
  invisible(rslt)
}

pd_OR_noself_neighbor <- function(initial_prev, rate, time, n, d, neighbor_mat, initial_loc = "m"){
  y <- matrix(0, nrow = time, ncol = n)
  initial_farm <- ifelse(initial_loc == "m", sample(1:600, 1), 
                         ifelse(initial_loc == "l", sample(601:1200, 1), 
                                ifelse(initial_loc == "r", sample(1201:1740,1),
                                       ifelse(initial_loc == "ll", sample(1741:1920,1),
                                              ifelse(initial_loc == "lm", sample(1921:2160,1),
                                                     ifelse(initial_loc == "lr", sample(2161:2820,1),
                                                            ifelse(initial_loc == "ul", sample(2821:4440,1),
                                                                   ifelse(initial_loc == "um", sample(4441:5100,1),
                                                                          ifelse(initial_loc == "ur", sample(5101:6000,1),return(NULL))))))))))
  y[1, initial_farm] <- initial_prev
  odds_y <- matrix(0, nrow = time, ncol = n)
  odds_y[1,] <- y[1,]/(1-y[1,])
  
  for(t in 2:time){
    for(i in 1:n){
      odds_y[t, i] <- odds_y[t-1,i] + rate*sum(d[i,] * neighbor_mat[,i] * y[t-1,]) - rate * y[t-1,i]
      y[t,i] <- odds_y[t,i]/(1 + odds_y[t,i])
    }
  }
  prev <- apply(y, 1, FUN = function(x){sum(x >= 0.1)/n})
  return(list(y = y, prev = prev))
}





spatial_sampling_spread <- function(coord, time, d, neighbor_mat, initial_prev, rate, sample_size, cut_off_prev = 0.1, seed = 202117455, n_iter){
  N <- dim(coord)[1]
  p <- rep(sample_size/N, N)
  X = as.matrix(coord)
  cumm <- matrix(0, ncol = 7, nrow = time)
  res_method <- matrix(0, ncol = 6, nrow = time)
  s <- list()
  set.seed(seed)
  res <- pd_OR_noself_neighbor(initial_prev, rate, time, N, d, neighbor_mat, initial_loc = initial_loc)
  set.seed(NULL)
  Equaldsgn <- list(None = list(panel = c(PanelOne = sample_size), seltype = "Equal"))
  for(t in 1:n_iter){
    for(i in 1:time){
      GRTS_res <- grts(design = Equaldsgn, DesignID = "EQUAL", type.frame = "finite", src.frame = "att.frame", att.frame = coord, xcoord = "X", ycoord = "Y", shapefile = F)
      s[[1]] <- BalancedSampling::cube(p, X)
      s[[2]]<- BalancedSampling::lcube(p, X, cbind(p))
      s[[3]] <- BalancedSampling::lpm(p, X, h = sample_size)
      s[[4]] <- BalancedSampling::scps(p, X)
      tmp <- matrix(c(GRTS_res$xcoord, GRTS_res$ycoord), ncol = 2)
      s[[5]] <- apply(tmp, MARGIN = 1, FUN = function(x){which(x[1] == coord$X & x[2] == coord$Y)})
      s[[6]] <- sample(1:N, sample_size, replace = F)
      res_method[i,] <- c(any(rbinom(sample_size, 1, prob = ifelse(res$y[i, s[[1]]] >= cut_off_prev, 1, 0))),
                          any(rbinom(sample_size, 1, prob = ifelse(res$y[i, s[[2]]] >= cut_off_prev, 1, 0))),
                          any(rbinom(sample_size, 1, prob = ifelse(res$y[i, s[[3]]] >= cut_off_prev, 1, 0))),
                          any(rbinom(sample_size, 1, prob = ifelse(res$y[i, s[[4]]] >= cut_off_prev, 1, 0))),
                          any(rbinom(sample_size, 1, prob = ifelse(res$y[i, s[[5]]] >= cut_off_prev, 1, 0))),
                          any(rbinom(sample_size, 1, prob = ifelse(res$y[i, s[[6]]] >= cut_off_prev, 1, 0))))
    }
    cumm[,1:6] <- cumm[,1:6] + res_method
  }
  cumm <- cumm/n_iter
  cumm[,7] <- res$prev
  return(res)
}



spatial_balance <- function(coord, time, d, neighbor_mat, initial_prev, rate, sample_size,  initial_loc = "m", cut_off_prev = 0.1, seed = 202117455, n_iter){
  N <- dim(coord)[1]
  p <- rep(sample_size/N, N)
  X = as.matrix(coord)
  cumm <- matrix(0, ncol = 6, nrow = time)
  res_method <- matrix(0, ncol = 6, nrow = time)
  s <- list()
  set.seed(seed)
  #res <- pd_OR_noself_neighbor(initial_prev, rate, time, N, d, neighbor_mat, initial_loc = initial_loc)
  set.seed(NULL)
  Equaldsgn <- list(None = list(panel = c(PanelOne = sample_size), seltype = "Equal"))
  for(t in 1:n_iter){
    for(i in 1:time){
      GRTS_res <- grts(design = Equaldsgn, DesignID = "EQUAL", type.frame = "finite", src.frame = "att.frame", att.frame = coord, xcoord = "X", ycoord = "Y", shapefile = F)
      s[[1]] <- BalancedSampling::cube(p, X)
      s[[2]]<- BalancedSampling::lcube(p, X, cbind(p))
      s[[3]] <- BalancedSampling::lpm(p, X, h = sample_size)
      s[[4]] <- BalancedSampling::scps(p, X)
      tmp <- matrix(c(GRTS_res$xcoord, GRTS_res$ycoord), ncol = 2)
      s[[5]] <- apply(tmp, MARGIN = 1, FUN = function(x){which(x[1] == coord$X & x[2] == coord$Y)})
      s[[6]] <- sample(1:N, sample_size, replace = F)
      res_method[i,] <- c(BalancedSampling::sb(p,X,s[[1]]), BalancedSampling::sb(p,X,s[[2]]), BalancedSampling::sb(p,X,s[[3]]), BalancedSampling::sb(p,X,s[[4]]), BalancedSampling::sb(p,X,s[[5]]), BalancedSampling::sb(p,X,s[[6]]))
    }
    cumm[,1:6] <- cumm[,1:6] + res_method
  }
  cumm <- cumm/n_iter
  #cumm[,7] <- res$prev
  return(cumm)
}


rand.separated <- function(n,x0,x1,y0,y1,d){
  while(1){
    nums <- cbind(runif(n,x0,x1),runif(n,y0,y1))
    if(min(dist(nums)) >= d) return(round(nums, 1))
  }
}



spatial_sampling <- function(time, initial_prev, rate, sample_size,  initial_loc = "m", cut_off_prev = 0.1, seed = 202117455, n_iter){
  cumm <- matrix(0, ncol = 2, nrow = time)
  ul.nums <- rand.separated(1620, 0, 100, 130, 195, 0.1)
  um.nums <- rand.separated(660, 100, 200, 130, 195, 0.1)
  ur.nums <- rand.separated(900, 200, 300, 130, 195, 0.1)
  l.nums <- rand.separated(600, 0, 100, 65, 130, 0.1)
  m.nums <- rand.separated(600, 100, 200, 65, 130, 0.1)
  r.nums <- rand.separated(540, 200, 300, 65, 130, 0.1)
  ll.nums <- rand.separated(180, 0, 100, 0, 65, 0.1)
  lm.nums <- rand.separated(240, 100, 200, 0, 65, 0.1)
  lr.nums <- rand.separated(660, 200, 300, 0, 65, 0.1)
  
  coord <- rbind(m.nums, l.nums, r.nums, ll.nums, lm.nums, lr.nums, ul.nums, um.nums, ur.nums)
  sp_coord <- sp::SpatialPoints(coords = coord)
  d_geo <- as.matrix(dist(coord, diag = T))
  d <- exp(-d_geo/10)
  neighbor_mat <- (d_geo <= 5)
  coord <- data.frame(X = coord[,1], Y = coord[,2])
  N <- dim(coord)[1]
  p <- rep(sample_size/N, N)
  X = as.matrix(coord)
  res_method <- matrix(0, ncol = 1, nrow = time)
  s <- list()
  set.seed(seed)
  res <- pd_OR_noself_neighbor(initial_prev, rate, time, N, d, neighbor_mat, initial_loc = initial_loc)
  set.seed(NULL)
  for(t in 1:n_iter){
    for(i in 1:time){
      bas_res <- SDraw::bas.point(sp_coord, sample_size)
      s[[1]] <- as.numeric(as.character(bas_res$geometryID))
      res_method[i,] <- c(any(rbinom(sample_size, 1, prob = ifelse(res$y[i, s[[1]]] > cut_off_prev, 1, 0))))
    }
    cumm[,1] <- cumm[,1] + res_method
    cumm[,2] <- cumm[,2] + res$prev
  }
  cumm <- cumm/n_iter
  return(cumm)
}
