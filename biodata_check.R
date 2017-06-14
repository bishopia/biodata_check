#Title: biodata_check() function
#Author: Ian Bishop
#Date: 10 June 2017

biodata_check <- function(list_filename="", data_filename="", biodata_filename="", data_type=c("List", "Dataframe")) {

  #required libraries
  library(qdapTools)
  library(reshape2)
  
  #Import Biodata Complete Taxonomy  
  biodata_names <- read.csv(paste(biodata_filename, ".csv", sep=""), header=TRUE, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
  
  
  if (data_type=="List") {

    #Import user list of taxa
    user_list <- read.csv(paste(list_filename, ".csv", sep=""), header=TRUE, stringsAsFactors = FALSE)
  
    #Add user list to new conversion dataframe
    conversion_df <- user_list
    names(conversion_df)[1] <- "UserTaxa"
  
    #Compare user names to Biodata bench names, 
    conversion_df$CurrentBiodataName <- lookup_e(conversion_df[,1], biodata_names[,c("BenchTaxonName", "BiodataTaxonName")])
  
    #Add column for change notes
    conversion_df$ConversionNotes <- NA
    conversion_df[is.na(conversion_df$CurrentBiodataName)==FALSE,]$ConversionNotes <- "Name Updated to Current Biodata Name"
    conversion_df[is.na(conversion_df$CurrentBiodataName)==TRUE,]$ConversionNotes <- "Name Not Found In Biodata"
    conversion_df$ConversionNotes[as.character(conversion_df$CurrentBiodataName)==as.character(conversion_df$UserTaxa)] <- "Name Already Current"
  
    #Replace NAs in CurrentBiodataName column with original user name.
    conversion_df[is.na(conversion_df$CurrentBiodataName)==TRUE,]$CurrentBiodataName <- conversion_df[is.na(conversion_df$CurrentBiodataName)==TRUE,]$UserTaxa
  
    #confirm
    conversion_df <- data.frame(lapply(conversion_df, as.character), stringsAsFactors=FALSE)
  
    #export conversion_df to .csv
    # write.csv(conversion_df, "OUTPUT_converted_names.csv", row.names = FALSE)
    return(conversion_df)
  
  } else if(data_type=="Dataframe") {
  
    #Import user wide format tidy dataset  
    user_data <- read.csv(paste(data_filename, ".csv", sep=""), check.names = FALSE, stringsAsFactors=FALSE, header=TRUE, fileEncoding="UTF-8-BOM")
    
    #Add user list to new conversion dataframe
    working_wide <- user_data
  
    #convert to long format
    working_wide <- melt(working_wide, value.name="Count", variable.name="Taxon", id.vars="Sample")
  
    #add BiodataTaxon column of current Biodata names
    working_wide$BiodataTaxon <- lookup_e(working_wide$Taxon, biodata_names[,c("BenchTaxonName", "BiodataTaxonName")])
  
  
    
    #change sheet
    change_sheet <- working_wide[,c("Taxon", "BiodataTaxon")]
    change_sheet <- unique(change_sheet)
    #Add column for change notes
    change_sheet$ConversionNotes <- NA
    change_sheet[is.na(change_sheet$BiodataTaxon)==FALSE,]$ConversionNotes <- "Name Updated to Current Biodata Name"
    change_sheet[is.na(change_sheet$BiodataTaxon)==TRUE,]$ConversionNotes <- "Name Not Found In Biodata"
    change_sheet$ConversionNotes[as.character(change_sheet$BiodataTaxon)==as.character(change_sheet$Taxon)] <- "Name Already Current"
    #fill NAs in change_sheet$BiodataTaxon with Taxon
    change_sheet$BiodataTaxon[is.na(change_sheet$BiodataTaxon)] <- as.character(change_sheet$Taxon[is.na(change_sheet$BiodataTaxon)])
  
  
    
    #Fill in BiodataTaxon NAs with original taxon.
    working_wide$BiodataTaxon[is.na(working_wide$BiodataTaxon)] <- as.character(working_wide$Taxon[is.na(working_wide$BiodataTaxon)])
    
    #collected list of names that have been changed
    changed_names <- as.data.frame(unique(working_wide[working_wide$Taxon != working_wide$BiodataTaxon,]$Taxon))
    colnames(changed_names) <- "Changed Names"
    
    
    
    #Replace NAs in BiodataTaxon column with original Taxon.
    working_wide$Taxon <- working_wide$BiodataTaxon
    #Remove BiodataTaxon column
    working_wide$BiodataTaxon <- NULL
  
    #Collapse any species values that are now synonymized
    working_wide <- aggregate(Count ~ Sample + Taxon, FUN = sum, data=working_wide)
    
    #Convert from long format back to wide format
    working_wide <- dcast(working_wide, Sample ~ Taxon)
    
    #consider exporting this as XLSX with two datasheets
    #export conversion_df to .csv
    # write.csv(working_wide, "OUTPUT_renamed_data.csv", row.names = FALSE)
    # #export changed name list
    # write.csv(changed_names, "OUTPUT_name_changes.csv", row.names = FALSE)
    
    output<-list(UpdatedDataset = working_wide, ChangedNames = changed_names, NameTranslationSheet=change_sheet)
    return(output)
  
  } else {return(print("ERROR: Please choose either 'List' or 'Dataframe' for data_type"))}
}