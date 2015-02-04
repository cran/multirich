#This file is part of the multirich package. See package information for license and details
#THIS CODE IS PROVIDED AS IS, WITHOUT ANY WARRANTY OF SUITABILITY FOR A PARTICULAR FUNCTION.

#' multirich (package)
#'
#' Package to calculate the number of unique trait combinations and scale the
#' number of unique trait combinations by the total number of trait combinations
#' possible as measures of multivariate richness.
#'
#' @details \tabular{ll}{ %%Note tabular{ll} is tabular{lowercase(LL)} not tabular{11}
#' Package: \tab multirich\cr
#' Type: \tab Package\cr
#' Version: \tab 2.0.2\cr
#' Date: \tab 2014-02-04\cr
#' License: \tab GPL-2 (or later)\cr
#' }
#' @author Alexander "Sasha" Keyel\cr
#' Maintainer: Sasha Keyel <skeyel@@gmail.com>
#' @references Keyel, A.C. and K. Wiegand. (in review) Validating the use of
#' Unique Trait Combinations for measuring multivariate functional richness.
#' @keywords package
#' @name multirich-package
#' @docType package
NULL
  
#' mvfd Main function to calculate functional richness
#'
#' Main function to calculate multivariate richness.
#' Goal is to mirror format of dbfd from Laliberte & Shipley 2011.
#'
#' Currently takes a record x trait matrix and an optional community matrix.
#'
#' @param in.mat A record x trait matrix.  Needs to be in matrix format, not as a dataframe
#' @param in.com A community x species matrix.  If only one community with all the species is considered, you can enter "none", and the code will auto-create the community needed for the script. 
#' @param unequal.abund A feature not currently in script, intended as an option to indicate whether abundance data should be incorporated.
#' @param resolution This input controls rounding of the data (categorization was acheived by rounding for simplicity).  0 indicates integers, 1 = 1 decimal place, 2 = 2 decimal places, -1 = 10's place., etc.)
#' @param st.range Option to control range standardization.  This was never properly scripted and should remain 0.  Any desired range transformations should be done prior to this script.
#' @param log.trans Option to control whether or not data are log transformed.  This has not been properly tested as I found it easier to log-transform the data manually in Excel.
#' @param col.mins If option is "use data" the function will get the minimum from the data.  Otherwise a vector of minimum values to use can be specified
#' @param col.maxs If option is "use data" the function will get the maximum values from the data.  Otherwise a vector of maximum values can be specified.
#' @param traitspace If set to "use data", the function will estimate the traitspace as the product of trait ranges.  Otherwise, the specified traitspace will be used for scaling (e.g., if you want to input a traitspace based on a convex hull)
#' @param calc.ovr An option to determine whether overlap is calculated.  This may be slow or buggy, so in some cases it may be easier to turn it off.
#' @param force.matrix an option to determine whether to try to force an input into matrix format
#' @author A.C. Keyel
#' @example
#' man-roxygen/mvfd_ex.r
#' @export mvfd
mvfd = function(in.mat,in.com = "none",unequal.abund = F,resolution = 0,st.range = 0,log.trans = 0,col.mins = "use data",col.maxs = "use data",traitspace = "use data",calc.ovr = 1, force.matrix = TRUE){

    #Check length of in.mat.  If it is length 0, return 0 with a warning.
    if (length(in.mat) == 0){
        warning("Species matrix is empty.  Script will return a value of 0.  If this is undesirable, make appropriate adjustments to your code.")
        return(list(utc = 0,sutc = NA,rmvo = 0,smvo = 0,cmvo = 0,meanmvo = 0,medmvo = 0,maxmvo = 0,minmvo = 0))    
        } 

    if (force.matrix == TRUE){
        in.mat = as.matrix(in.mat)
        }

    #Check to make sure in.mat is a suitably formatted matrix.
    matrix.check(in.mat)

    #Set up even abundances for one community if community vector is missing.
    if (length(in.com) == 1){ #This is to avoid getting warning messages when entering communities!
        if (in.com == "none"){
            in.com = matrix(rep(1,nrow(in.mat)),nrow = 1, dimnames = list(c("Community1"),rownames(in.mat)))        
            }
        }
    
    #Pre-process data
      # Warning: not all options are functional
    dpp = data.preprocess(in.mat,log.trans,st.range,col.mins,col.maxs)
    in.mat = dpp[[1]]
    col.mins = dpp[[2]]
    col.maxs = dpp[[3]]

    #Convert resolution into the format needed.
    if (length(resolution) == 1){
        cell.res = rep(resolution,(ncol(in.mat)*nrow(in.mat))) #Write resolution for every cell, for rounding .
        col.res = rep(resolution,ncol(in.mat)) #Assign a resolution value for every column, for calculation of traitspace
        }
        
    #Give an error message if resolution is not in a format currently supported by the code.
    if (length(resolution) != 1){
        col.res = resolution
        print("length(resolution != 1), and the code to handle that has not actually been scripted yet.")
        }

    #Convert data to categories for calculation of multivariate richness & functional redundancy
    cat.mat = df.categorize(in.mat,cell.res)

    #Set up column minimums & maximums
      #Done after categorization to make it so the traitspace includes the full categories.
    if (col.mins == "use data"){
        col.mins = apply(cat.mat,2,min) #Get minimum for each column
        }
    if (col.maxs == "use data"){
        col.maxs = apply(cat.mat,2,max) #Get max for each column
        }

    #Get size of traitspace for scaling purposes
    if (traitspace == "use data"){
        traitspace = get.traitspace(cat.mat,col.res,col.mins,col.maxs)
        }

    #Convert to functional units & calculate multivariate richness and multivariate overlap
    mvd.out = calc.mvd(cat.mat,in.com,traitspace,calc.ovr)
    
    return(mvd.out)
    }

#' Check that input to mvfd is in matrix format
#' 
#' @param in.mat An input to be tested for matrix formatting
matrix.check = function(in.mat){
    row.test = nrow(in.mat)
    col.test = ncol(in.mat)
    
    is.err = 0
    err.message = ""
    
    if (length(row.test) == 0 | length(col.test) == 0){
        is.err = 1
        this.err = 'Input is not in matrix format.  Note that for matrices with only one row, some R processes will convert them out of matrix format.\n'
        err.message = sprintf('%s%s',err.message, this.err)
        }
    
    if (typeof(in.mat) != "double"){
        is.err = 1
        this.err = sprintf('Input is of type %s. The input needs to be in matrix format, which is type double.\n',typeof(in.mat))
        err.message = sprintf('%s%s', err.message, this.err)
        }
    
    #Check that matrix has row and column names
    if (length(rownames(in.mat)) == 0){
        is.err = 1
        this.err = 'Input record x trait matrix must have row names\n'
        err.message = sprintf('%s%s', err.message, this.err)
        }
        
    if (length(colnames(in.mat)) == 0){
        is.err = 1
        this.err = 'Input record x trait matrix must have column names'
        err.message = sprintf('%s%s', err.message, this.err)
        }

    if (is.err == 1){ stop(err.message) }
    }

#' Preprocess data
#'
#' @param in.mat A record x trait matrix
#' @param log.trans Whether or not to do a log transformation. It is recommended to do all log-transformations outside of the multirich package.
#' @param st.range Must equal zero. This option may be added in future versions of the package
#' @param col.mins Minimum column values to use. "use data" will extract these from the dataset
#' @param col.maxs Maximum column values to use. "use data" will extract these from the dataset
data.preprocess = function(in.mat,log.trans = 0,st.range = 0,col.mins,col.maxs){

    ## Log-transform data using natural log if desired #NOTE: Currently log-transforms EVERYTHING!
    if (log.trans == 1 | log.trans == "e"){
        in.mat = log(in.mat)
        #Log-transform specified min/max values. If column minimums are based on the data, they will be assigned later
        if (col.mins != "use data"){ 
            col.mins = log(col.mins)
            }
        if (col.maxs != "use data"){
            col.maxs = log(col.maxs)
            }
        }

    ## Log-transform data using base 10
    if (log.trans == 10){
        in.mat = log10(in.mat)
        #Log10-transform specified min/max values.  If col mins and maxes use the default, they will be assigned later based on the logged data.
        if (col.mins != "use data"){
            col.mins = log10(col.mins)
            }
        if (col.maxs != "use data"){
            col.maxs = log10(col.maxs)
            }
        }

    ## Standardize data if desired
    
    #**# This part of the code is non-functional & any standardization must be done prior to running the main function.
    
    # #Standardize by range
    # if (st.range == 1){
    # 
    #     #**# This does not work:
    #     in.mat = apply(in.mat,2,range.standardize,col.mins,col.maxs)
    #     }
    #
    # #Standardize by mean and standard deviation
    # if (st.range == 2){
    #     #This is not yet scripted
    #     }
    
    return(list(in.mat,col.mins,col.maxs))
    }

## These functions are used by the data.preprocess function, but are currently non-functional
# ### Standardize data by range ###
# range.standardize = function(x,st.range,col.mins,col.maxs){
# 
#     col.range = col.maxs - col.mins #**# col.mins is a vector.  Need to get the right elements to match up.
#         
#     st.x = (x - col.mins) / col.range 
#     
#     return(st.x)
#     }
# 
# ### Standardize data by mean and standard deviation ###
# sd.standardize = function(inputs){
#     #**# Not yet scripted
#
#     return(st.x)
#     }


#' Calculate trait space (simple)
#'
#' Function to calculate the traitspace using minimum and maximum column values.
#' Creates a hypercubic traitspace.
#'
#' @param in.mat A record x trait matrix
#' @param col.res The number of decimal places to use for rounding (resolution) purposes
#' @param col.mins The minimum values to use for the trait space. "use data" extracts minimums from in.mat
#' @param col.maxs The maximum values to use for the trait space. "use data" extracts maximums from in.mat
#' 
get.traitspace = function(in.mat,col.res,col.mins,col.maxs){

    #Convert resolution from number of decimal places to powers of 10
    act.res = 10^(-col.res)

    #Compute number of categories for each column
    col.ranges = col.maxs - col.mins
    col.mod.ranges = col.ranges / act.res + 1 #Add 1 to include the lower endpoint

    #Get product of the column ranges to give entire trait space.
    traitspace = prod(col.mod.ranges) 
        
    return(traitspace)
    }


#'Categorize data
#'
#' Currently by rounding, which is a sub-optimal approach for many questions
#'
#' @param in.mat A record x trait matrix
#' @param cell.res The number of decimal places to round the data to.
df.categorize = function(in.mat, cell.res){
  
    # WARNING: THIS STEP WILL FAIL IF THE INPUT IS A DATAFRAME!  NEEDS TO BE IN MATRIX FORMAT.
    cat.mat = round(in.mat,cell.res) 

    return(cat.mat)
    }


#' Calculate multivariate richness and functional overlap
#'
#' Function to reduce to funcitonal units, and calculate multivariate richness
#' and functional overlap
#'
#' @param in.mat A record by trait matrix
#' @param in.com A record by community matrix
#' @param traitspace The total trait space for scaling purposes
#' @param calc.ovr An indicator for whether overlap metrics should be calculated
calc.mvd = function(in.mat,in.com,traitspace,calc.ovr){
    
    #Create output vectors and lists
    us.mvr.vec = c()  #Create a vector for unscaled multivariate richness    
    mvo.raw = list()     #create a vector for the raw overlap values for examining occupied niche space
    mvo.simple = c()  #Create a vector for simple overlap
    mvo.coverage = c() #Create a vector for coverage overlap
    mvo.mean = c()     #Create a vector for mean overlap
    mvo.median = c()    #Create a vector for median overlap
    mvo.max = c()       #Create a vector for max overlap
    mvo.min = c()      #Create a vector for min overlap
    
    #Calculate richness and overlap (optional) for each community by joining it to the record x trait matrix
    for (a.row in 1:nrow(in.com)){
        #Get row of community
        this.com = in.com[a.row, ]
        
        #Join with in.matrix
        n.cols = ncol(in.mat) + 1 #Original columns in in.mat plus the column we're adding
        r.nam = rownames(in.mat) #Row names from in.mat should == col names from in.com.
        c.nam = c(rownames(in.com)[a.row],colnames(in.mat))
        join.mat = matrix(c(this.com,in.mat), ncol = n.cols, dimnames = list(r.nam,c.nam))
        
        #Drop records where the species is not present
        com.mat = join.mat[join.mat[ ,1] != 0, ]
        
        #Need to drop the whole community column.  This step is embedded in a series of other steps required to maintain the matrix formatting because R converts one-row matrices into vectors.
        cm.r = rownames(com.mat)
        cm.c = colnames(com.mat)
        com.mat = matrix(com.mat, ncol = n.cols, dimnames = list(cm.r,cm.c)) #BS step to make R preserve matrix format when only 1 row remains

        #Drop the community column
          #**# NOTE in order to incorporate abundances, this portion of the code will need re-scripting
        com.mat = com.mat[ ,2:ncol(com.mat)]
        cm.c = cm.c[2:length(cm.c)] #drop the whole community column here as well.
        com.mat = matrix(com.mat, ncol = (n.cols - 1),dimnames = list(cm.r,cm.c)) #BS step to make R preserve matrix format when only 1 row remains.
        
        if (calc.ovr == 1){
            #Compute multivariate overlap
            dups = duplicated(com.mat) 
            this.mvo = calc.mvo(com.mat,dups)
              mvo.raw = append(mvo.raw,list(this.mvo[[1]]))
              mvo.simple = c(mvo.simple,this.mvo[[2]])
              mvo.coverage = c(mvo.coverage,this.mvo[[3]])
              mvo.mean = c(mvo.mean,this.mvo[[4]])
              mvo.median = c(mvo.median,this.mvo[[5]])
              mvo.max = c(mvo.max,this.mvo[[6]])
              mvo.min = c(mvo.min,this.mvo[[7]])
            }
        
        #Compute multivariate richness
        fu.mat = com.mat[!duplicated(com.mat), ]  #Create a matrix of just functional units
        fu.mat = matrix(fu.mat, ncol = (n.cols -1)) #Step to make R preserve matrix format when only 1 row remains.
        this.mvr = nrow(fu.mat)                     #Count the number of functional units
        us.mvr.vec = c(us.mvr.vec,this.mvr)        
        }

    sc.mvr.vec = us.mvr.vec / traitspace  #Compute scaled multivariate richness


    return(list(utc = us.mvr.vec,sutc = sc.mvr.vec,rmvo = mvo.raw,smvo = mvo.simple,cmvo = mvo.coverage,meanmvo = mvo.mean,medmvo = mvo.median,maxmvo = mvo.max,minmvo = mvo.min))    
    }


#' Partition multivariate richness Beta Component
#'
#' Function to calculate dissimilarity between a pair of communities, and
#' partition it into nestedness-related & turnover components  
#' Either the Sorensen index or the Jaccard index can be used to calculate dissimilarity
#'
#' @param in.mat A record x trait matrix.  Needs to be in matrix format, not as a dataframe
#' @param in.com A community x record matrix
#' @param index.rows A vector with 2 elements, that gives the row number for the
#' pair of communities of interest.
#' @param index.type specifies which index to use. Options are Sorensen (default) & Jaccard
#' @return \itemize{
#' \item aa Overlap between the two communities
#' \item dissimilarity Dissimilarity between the two communities
#' \item turnover This gives the turnover between the two communities.  To get the percent of dissimilarity due to turnover, divide by total dissimilarity.
#' \item Bnes Nestedness-related dissimilarity between the two communities.  To get the percent of dissimilarity due to this, divide by total dissimilarity.
#' }
#' @references Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity.  Global Ecology and Biogeography 19: 134-143\cr
#' Baselga, A. 2012. The relationship between species replacement, dissimilarity derived from nestedness, and nestedness.  Global Ecology and Biogeography 21: 1223-1232\cr
#' Villeger, S. Grenouillet, G., and Brosse, S. 2013.  Decomposing functional Beta-diversity reveals that low functional Beta-diversity is driven by low functional turnover in European fish assemblages.  Global Ecology and Biogeography 22: 671-681.
#' @author A.C. Keyel
#' @example man-roxygen/part_mvr_beta.r
#' @export part.mvr.beta
part.mvr.beta = function(in.mat,in.com,index.rows,index.type = "Sorensen"){
    
    #Get a,b,c
    #using same general approach as for mvr - could be optimized!
    #a = number of species present in both sites
    #b = number of species present in the second site but not the first
    #c = number of species present in the first site but not the second.

    #Get names for formatting step below 
    row.nam = rownames(in.mat)
    col.nam = colnames(in.mat)
    
    #Convert in.mat to functional units
    dup.index = duplicated(in.mat)
    fu.mat = in.mat[!dup.index, ]  #Create a matrix of just functional units

    #Get a matrix of just deleted rows & reformat to matrix format (needed if only 1 record)
    dup.mat = in.mat[dup.index, ]
    dup.r.nam = row.nam[dup.index]
    dup.mat = matrix(dup.mat, ncol = length(col.nam), dimnames = list(dup.r.nam,col.nam))    

    #Get a list of just deleted community records
    d.com = in.com[ ,dup.index]

    #Correct for R converting my matrix to a vector when there is only one row
    if (length(nrow(d.com)) == 0){
        d.com.names = names(d.com)
        d.com = matrix(d.com, nrow = length(d.com)) #Restore matrix format
        rownames(d.com) = d.com.names #Restore row names
        }

    dup.A = d.com[index.rows[1], ]
    dup.B = d.com[index.rows[2], ]

    #Remove duplicate names 
    row.nam = row.nam[!dup.index]

    #Reformat to fix matrix format when only one record  
    fu.mat = matrix(fu.mat, ncol = length(col.nam),dimnames = list(row.nam,col.nam)) #BS step to make R preserve matrix format when only 1 row remains.

    #Join the community data for the non-duplicate rows
    j.com = in.com[ ,!dup.index]
    
    #Get rows of communities to be compared
    com.A = j.com[index.rows[1], ]
    com.B = j.com[index.rows[2], ]
    
    n.cols = ncol(fu.mat) + 2 #Original columns in in.mat plus the column we're adding
    c.nam = c("A","B",col.nam)  
    join.mat = matrix(c(com.A,com.B,fu.mat), ncol = n.cols, dimnames = list(row.nam,c.nam))

    #If rows were dropped when converting to functional units, need to know which species are still represented
    if (nrow(dup.mat) > 0 ){

        #Go through deleted rows, and figure out which rows they correspond to.
        for (i in 1:nrow(dup.mat)){
            #Merge the deleted rows with the non-duplicate rows appropriately
            this.row = dup.mat[i, ]
            this.row = listtotext(this.row,",") #This is because R was only comparing the first elements of my vectors, and not the vectors themselves!
            this.dupa = dup.A[i]
            this.dupb = dup.B[i]
        
            #If both communities have no record from this row, go to the next row.
            if (this.dupa == 0 & this.dupb == 0){
                next 
                }
            
            for (j in 1:nrow(join.mat)){
                join.row = join.mat[j,3:ncol(join.mat) ] #exclude the first two elements from the check for equivalence.
                join.row = listtotext(join.row,",")
                join.A = join.mat[j,1]
                join.B = join.mat[j,2]
                
                #Warning: for vectors, R only compares the FIRST ELEMENTS!  STUPID R! I converted my vectors to strings above to make the comparison of all elements possible.
                if (this.row == join.row){
                                    
                    #Check & update com.A (only if it is not already 1)
                    if (join.A != 1){
                        if (this.dupa == 1){
                            join.mat[j,1] = 1 #If the community has a presence in the duplicate row, assign it to the merged row.
                            }
                        }
                    
                    #Check & update com.B (only if it is not already 1)
                    if (join.B != 1){
                        if (this.dupb == 1){
                            join.mat[j,2] = 1 #If the community has a presence in the duplicate row, update it.
                            }
                        }   
                    }
                }     
            }
        }
    

    #Drop records where neither species is present (#strictly speaking this isn't necesary, the coding below would not count these rows anyhow)
    com.mat = join.mat[join.mat[ ,1] != 0 | join.mat[ ,2] != 0, ]

    #need to count number of shared functional units, and number of functional units unique to each species.
    aa = 0     #a = number of species present in both sites
    bb = 0     #b = number of species present in the second site but not the first
    cc = 0     #c = number of species present in the first site but not the second.
    
    for (a.row in 1:nrow(com.mat)){
        this.row = com.mat[a.row, ]
        this.A = this.row[1]
        this.B = this.row[2]
        
        #If both are 1, add a shared species 
        if (this.A == 1 & this.B == 1){
            aa = aa + 1
            }
        
        if (this.A == 1 & this.B == 0){
            cc = cc + 1
            }
        
        if (this.B == 1 & this.A == 0){
            bb = bb + 1
            }
        
        }
            
    #Parameters to make the index work correctly.  Sorensen is the default and will occur unless Jaccard is specified.
    w = 2
    v = 1
    
    #Jaccard implementation
    if (index.type == "Jaccard"){
        w = 1
        v = 2      
        }
  
    #Calculate Dissimilarity   (#letters doubled to avoid conflict between c default r function c()
    dissimilarity.numerator = bb + cc
    dissimilarity.denominator = w * aa + bb + cc
    dissimilarity = dissimilarity.numerator / dissimilarity.denominator
      
    #Calculate turnover component
    turnover.numerator = v * min(bb,cc)
    turnover.denominator = aa + v * min(bb,cc)
    turnover = turnover.numerator / turnover.denominator
    
    #Calculate "nestedness-resultant" component
    Bnes.numerator.1 = max(bb,cc) - min (bb,cc)
    Bnes.denominator.1 = w * aa + bb + cc
    Bnes.numerator.2 = aa
    Bnes.denominator.2 = aa + v * min(bb,cc)
    Bnes = (Bnes.numerator.1 / Bnes.denominator.1) * (Bnes.numerator.2 / Bnes.denominator.2)

    return(list(aa,dissimilarity,turnover,Bnes))
    
    }

#' Calculate functional overlap
#'
#' @param in.mat a record x trait matrix
#' @param dups dups indicate which records in the record x trait matrix are duplicates 
calc.mvo = function(in.mat,dups){

    total.records = nrow(in.mat)

    #Simple redundancy
    simple.mvo = sum(dups) / total.records #Divide number of redundant rows over the total number of rows.  Gives an index of redundancy that varies between 0 and 1.  sum(dups) gives count of "TRUE" which is the number of redundant records


    br.nams = rownames(in.mat) #Step to get R to preserve matrix format when only one species is considered
    bc.nams = colnames(in.mat)
    b.nrow = nrow(in.mat)
    
    b.mat = in.mat[!dups, ] #Get matrix of non-duplicated values
    br.nams2 = br.nams[!dups]
    b.mat = matrix(b.mat, nrow = length(br.nams2), dimnames = list(br.nams2,bc.nams))
    
    t.mat = in.mat[dups, ]  #Get matrix of duplicated values
    
    #Convert the above matrices to text lists to allow equality testing.
    b.as.text = c()
    for (arow in 1:nrow(b.mat)){
        this.row = b.mat[arow, ]
        as.text = listtotext(this.row,",")
        b.as.text = c(b.as.text,as.text)
        }

    t.as.text = c("there are no duplicate values")
    #Need a check to prevent no duplicates from crashing the code!
    if (length(nrow(t.mat)) != 0){
        if(nrow(t.mat) != 0){
            t.as.text = c()
            for (arow in 1:nrow(t.mat)){
                this.row = t.mat[arow, ]
                as.text = listtotext(this.row,",")
                t.as.text = c(t.as.text,as.text)
                }
            }
        }
        
    #Go through matrix of non-duplicated values
    b.count = rep(0,length(b.as.text)) #Create a vector to count number of overlaps for the vector.
    names(b.count) = b.as.text
    for (i in 1:length(b.as.text)){
        item = b.as.text[i]
        
        #Go through matrix of duplicated values and count copies
        for (d in 1:length(t.as.text)){
            dupli = t.as.text[d]

            #If the duplicate matches the row, increment the duplicate count by one.
            if (dupli == item){
                b.count[i] = b.count[i] + 1
                }
            }
        }
    
    #Coverage overlap is length( b.count that is not 0) / length b.count
    coverage.mvo = length(b.count[b.count > 0]) / length(b.count)

    #mean, min and max can be done directly from b.count.
    min.mvo = min(b.count)
    max.mvo = max(b.count)
    mean.mvo = mean(b.count)
    med.mvo = median(b.count)
    
    mvo = list(b.count,simple.mvo,coverage.mvo,mean.mvo,med.mvo,max.mvo,min.mvo)
    return (mvo)
    }

#' List to text
#'
#' Take a list (or vector) and convert it to text separated by separator.
#' Does not work for nested lists
#' ORIGINALLY FROM MYR.R script
#'
#' @param inlist An unnested list
#' @param separator A separator to separate elements of the list.
listtotext = function(inlist,separator){

    #Check that there is no nesting in the list structure
    for (thing in inlist){
        if (length(thing) > 1){
            print("One or more list elements containied multiple elements.")
            print("This is not allowed, as it will lead to strange outputs.")
            return(NA)
            }
        }
        
    #Check that separator is input (otherwise it can let a typo through where just a one element list or a separator is given)
    sep.check = separator #sep.check is not used later.  It is here to make it so R will crash if it is not input.

    #Perform list to text function
    textout = sprintf("%s",inlist[1])
    listlen = length(inlist)
    if (listlen > 1){
        for (item in inlist[2:listlen]){
            textout = sprintf("%s%s%s",textout,separator,item)
            }
        }
    return(textout)
    }
